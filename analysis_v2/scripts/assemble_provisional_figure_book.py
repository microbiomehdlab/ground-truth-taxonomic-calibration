#!/usr/bin/env python3
"""Assemble a visibly provisional, checksum-tracked figure review book.

This is a presentation artifact, not a scientific analysis. It never promotes
an inventory entry to checked or infers a missing analysis from a filename.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import shutil
from pathlib import Path

FIELDS = ("figure", "panel", "title", "state", "asset", "seal", "note")
STATES = {"checked", "provisional", "placeholder"}
FORMATS = {".png", ".jpg", ".jpeg", ".svg", ".pdf"}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def inside(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    return True


def build(inventory: Path, root: Path, outdir: Path) -> None:
    if outdir.exists():
        raise ValueError(f"output directory exists: {outdir}")
    with inventory.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(FIELDS):
            raise ValueError("inventory columns must be exactly: " + ", ".join(FIELDS))
        rows = list(reader)
    if not rows:
        raise ValueError("empty figure inventory")
    seen = set()
    prepared = []
    for row in rows:
        identity = (row["figure"], row["panel"])
        if not all(identity) or identity in seen or row["state"] not in STATES:
            raise ValueError(f"invalid or duplicate figure panel: {identity}")
        seen.add(identity)
        asset_name = row["asset"].strip()
        seal_name = row["seal"].strip()
        if row["state"] == "placeholder" and (asset_name or seal_name):
            raise ValueError(f"placeholder has asset or seal: {identity}")
        if row["state"] == "checked" and (not asset_name or not seal_name):
            raise ValueError(f"checked panel lacks asset or seal: {identity}")
        asset = root / asset_name if asset_name else None
        seal = root / seal_name if seal_name else None
        if asset and (asset.suffix.lower() not in FORMATS or not inside(asset, root)):
            raise ValueError(f"invalid or external asset: {asset_name}")
        if seal and not inside(seal, root):
            raise ValueError(f"external seal: {seal_name}")
        available = bool(asset and asset.is_file() and asset.stat().st_size > 0)
        sealed = bool(seal and seal.is_file() and seal.stat().st_size > 0)
        if row["state"] == "checked" and not (available and sealed):
            raise ValueError(f"checked panel missing asset or seal: {identity}")
        prepared.append((row, asset, available, sealed))

    outdir.mkdir(parents=True)
    assets = outdir / "assets"
    assets.mkdir()
    manifest = []
    cards = []
    for index, (row, asset, available, sealed) in enumerate(prepared, start=1):
        state = row["state"]
        copied = ""
        checksum = ""
        if state != "placeholder" and available:
            copied = f"assets/{index:02d}_{asset.name}"
            shutil.copy2(asset, outdir / copied)
            checksum = digest(outdir / copied)
        status = state.upper() if copied else "MISSING PREVIEW" if state != "placeholder" else "PLACEHOLDER"
        manifest.append(dict(figure=row["figure"], panel=row["panel"], title=row["title"],
                             state=state, preview_status=status, source_asset=row["asset"],
                             source_seal=row["seal"], seal_present="1" if sealed else "0",
                             copied_asset=copied, copied_sha256=checksum, note=row["note"]))
        title = html.escape(f"{row['figure']} · {row['panel']} · {row['title']}")
        note = html.escape(row["note"])
        visual = (f'<img src="{html.escape(copied, quote=True)}" alt="{title}">' if copied and
                  Path(copied).suffix.lower() != ".pdf" else
                  f'<object data="{html.escape(copied, quote=True)}" type="application/pdf"></object>'
                  if copied else '<div class="empty">No vetted preview supplied</div>')
        cards.append(f'<section class="card {state}"><div class="badge">{html.escape(status)}</div>'
                     f'<h2>{title}</h2>{visual}<p>{note}</p></section>')
    with (outdir / "FIGURE_MANIFEST.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest)
    page = """<!doctype html><html lang="en"><meta charset="utf-8">
<title>Development-only figure book</title><style>
body{font:16px system-ui,sans-serif;max-width:1250px;margin:2rem auto;padding:0 1rem;color:#233}
h1{margin-bottom:.2rem}.warning{background:#fff1ca;padding:1rem;border-left:6px solid #a65}
.card{border:1px solid #ccd;border-radius:8px;padding:1rem;margin:1.5rem 0;break-inside:avoid}
.badge{font-weight:bold;color:#a44}.checked .badge{color:#176b56}
img,object{display:block;max-width:100%;width:100%;max-height:850px;object-fit:contain;background:#fafafa}
object{height:700px}.empty{padding:4rem;text-align:center;background:#f4f5f6;color:#667}
@media print{.card{page-break-inside:avoid}body{max-width:none}}
</style><h1>Provisional paper figure book</h1>
<p class="warning">DEVELOPMENT ONLY. Panels are not final manuscript evidence. A preview's presence
does not upgrade its scientific status. Read FIGURE_MANIFEST.tsv for source paths and checksums.</p>
""" + "\n".join(cards) + "</html>\n"
    (outdir / "index.html").write_text(page, encoding="utf-8")
    (outdir / "DEVELOPMENT_ONLY.txt").write_text("status\tDEVELOPMENT_ONLY\n", encoding="utf-8")
    (outdir / "source.sha256").write_text(f"{digest(inventory)}  {inventory.resolve()}\n", encoding="utf-8")
    (outdir / "SUCCESS").write_text(f"panels\t{len(rows)}\nstatus\tPASS\n", encoding="utf-8")
    print(f"[PASS] assembled {len(rows)} figure panels; {sum(bool(m['copied_asset']) for m in manifest)} previews")
    print(outdir / "index.html")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    try:
        build(args.inventory, args.project_root, args.outdir)
    except (OSError, ValueError) as exc:
        raise SystemExit(f"[ERROR] {exc}") from exc


if __name__ == "__main__":
    main()
