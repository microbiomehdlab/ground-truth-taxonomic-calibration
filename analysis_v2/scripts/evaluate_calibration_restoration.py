#!/usr/bin/env python3
"""Compare uncorrected and calibrated disease models against unspiked baseline."""
from __future__ import annotations
import argparse, csv, hashlib, math
from collections import defaultdict
from pathlib import Path

FAMILY = ("cohort", "study", "analysis_population", "target_label", "assembly_arm",
          "profiler", "contrast")
KEY = FAMILY + ("dose_level", "feature")
REQUIRED = set(KEY) | {"effect", "q_value", "include", "spike_fraction_target"}


def num(x):
    z = float(x)
    if not math.isfinite(z): raise ValueError(f"non-finite number: {x}")
    return z


def render(x): return format(x, ".17g")


def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()


def read_models(path, contrast):
    result={}
    with path.open(newline="",encoding="utf-8") as f:
        r=csv.DictReader(f,delimiter="\t")
        if not REQUIRED<=set(r.fieldnames or []): raise ValueError("disease model columns incomplete")
        for row in r:
            if row["include"]!="1" or row["contrast"]!=contrast: continue
            k=tuple(row[x] for x in KEY)
            if k in result: raise ValueError(f"duplicate model row: {k}")
            result[k]=row
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--uncorrected-results",required=True,type=Path)
    p.add_argument("--corrected-results",required=True,type=Path)
    p.add_argument("--corrected-manifest",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    p.add_argument("--contrast",default="CRC_vs_Control")
    p.add_argument("--q-threshold",type=float,default=.05)
    a=p.parse_args()
    if a.outdir.exists() and any(a.outdir.iterdir()): raise SystemExit("[ERROR] OUTDIR must be new or empty")
    if not 0<a.q_threshold<1: raise SystemExit("[ERROR] invalid q threshold")
    try:
        u=read_models(a.uncorrected_results,a.contrast); c=read_models(a.corrected_results,a.contrast)
        if set(u)!=set(c): raise ValueError("uncorrected and corrected model row identities differ")
        targets={}
        with a.corrected_manifest.open(newline="",encoding="utf-8") as f:
            r=csv.DictReader(f,delimiter="\t")
            needed=set(FAMILY[:6])|{"target_feature"}
            if not needed<=set(r.fieldnames or []): raise ValueError("corrected manifest columns incomplete")
            for row in r:
                k=tuple(row[x] for x in FAMILY[:6])
                old=targets.setdefault(k,row["target_feature"])
                if old!=row["target_feature"]: raise ValueError("conflicting target feature")
        ledger=[]; stats=defaultdict(lambda:defaultdict(float))
        for k,ur in sorted(u.items()):
            family=k[:len(FAMILY)]; dose=k[len(FAMILY)]; feature=k[-1]
            if dose=="baseline":
                cr=c[k]
                if abs(num(ur["effect"])-num(cr["effect"]))>1e-10 or abs(num(ur["q_value"])-num(cr["q_value"]))>1e-10:
                    raise ValueError("calibration changed an unspiked baseline model")
                continue
            baseline_key=family+("baseline",feature)
            if baseline_key not in u: raise ValueError("positive row lacks baseline feature")
            br=u[baseline_key]; cr=c[k]
            be,ue,ce=map(num,(br["effect"],ur["effect"],cr["effect"]))
            bc,uc,cc=(num(br["q_value"])<=a.q_threshold,
                      num(ur["q_value"])<=a.q_threshold,num(cr["q_value"])<=a.q_threshold)
            target=targets.get(family[:6]); role="direct_target" if feature==target else "bystander"
            before,after=abs(ue-be),abs(ce-be)
            restoration=(1-after/before) if before else math.nan
            transition=("rescued_baseline" if bc and not uc and cc else
                        "harmed_baseline" if bc and uc and not cc else
                        "induced_removed" if not bc and uc and not cc else
                        "induced_remaining" if not bc and uc and cc else
                        "induced_new" if not bc and not uc and cc else
                        "baseline_retained" if bc and uc and cc else
                        "baseline_still_lost" if bc and not uc and not cc else "not_called")
            ledger.append(dict(zip(FAMILY,family),dose_level=dose,feature=feature,feature_role=role,
                spike_fraction_target=ur["spike_fraction_target"],baseline_effect=render(be),
                uncorrected_effect=render(ue),corrected_effect=render(ce),
                absolute_effect_error_before=render(before),absolute_effect_error_after=render(after),
                relative_effect_restoration="NA" if math.isnan(restoration) else render(restoration),
                baseline_called=int(bc),uncorrected_called=int(uc),corrected_called=int(cc),transition=transition))
            if role=="bystander":
                context=family+(dose,ur["spike_fraction_target"])
                s=stats[context];s["features"]+=1;s["baseline"]+=bc;s["uncorrected"]+=uc;s["corrected"]+=cc
                s[transition]+=1;s["before"]+=before;s["after"]+=after;s["improved"]+=after<before
        if not ledger: raise ValueError("no positive-dose comparisons")
        a.outdir.mkdir(parents=True)
        lp=a.outdir/"calibration_restoration_ledger.tsv"
        with lp.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=list(ledger[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(ledger)
        summary=[]
        for context,s in sorted(stats.items()):
            n=s["features"]; before=s["before"]; after=s["after"]
            summary.append(dict(zip(FAMILY+("dose_level","spike_fraction_target"),context),
                features_tested=int(n),baseline_biomarkers=int(s["baseline"]),
                uncorrected_biomarkers=int(s["uncorrected"]),corrected_biomarkers=int(s["corrected"]),
                baseline_biomarkers_rescued=int(s["rescued_baseline"]),baseline_biomarkers_harmed=int(s["harmed_baseline"]),
                induced_calls_removed=int(s["induced_removed"]),induced_calls_remaining=int(s["induced_remaining"]),
                induced_calls_created=int(s["induced_new"]),fraction_effect_errors_improved=render(s["improved"]/n),
                mean_absolute_effect_error_before=render(before/n),mean_absolute_effect_error_after=render(after/n),
                relative_mean_effect_restoration="NA" if before==0 else render(1-after/before)))
        sp=a.outdir/"calibration_restoration_summary.tsv"
        with sp.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(summary)
        note=a.outdir/"DEVELOPMENT_ONLY.txt";note.write_text("status=DEVELOPMENT_ONLY\ncausal_correction_claim=NO\n",encoding="utf-8")
        prov=a.outdir/"calibration_restoration.sha256"
        prov.write_text("".join(f"{digest(x)}  {x.resolve()}\n" for x in
            (a.uncorrected_results,a.corrected_results,a.corrected_manifest,lp,sp,note)),encoding="utf-8")
        (a.outdir/"SUCCESS").write_text(f"contexts={len(summary)}\nstatus=PASS\n",encoding="utf-8")
        print(f"[PASS] Calibration restoration evaluation: {a.outdir}")
    except (ValueError,OSError,KeyError) as e: raise SystemExit(f"[ERROR] {e}") from e

if __name__=="__main__": main()
