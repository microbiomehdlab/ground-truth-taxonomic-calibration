#!/usr/bin/env python3
import sys, argparse, gzip

def open_auto(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")

def load_taxdump(nodes_path, names_path):
    parent, rank, name = {}, {}, {}
    with open(nodes_path) as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 3: continue
            taxid, parent_taxid, r = parts[0], parts[1], parts[2]
            parent[taxid] = parent_taxid
            rank[taxid] = r
    with open(names_path) as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4: continue
            taxid, nm, _, cls = parts[0], parts[1], parts[2], parts[3]
            if cls == "scientific name":
                name[taxid] = nm
    return parent, rank, name

def ascend_to_rank(taxid, target_rank, parent, rank):
    cur = taxid
    seen = set()
    while cur and cur != "1" and cur not in seen:
        seen.add(cur)
        if rank.get(cur) == target_rank:
            return cur
        cur = parent.get(cur, "")
    return None

def lineage_path(taxid, parent, name):
    # Build textual path from root-ish to taxid
    chain = []
    cur = taxid
    seen = set()
    while cur and cur not in seen:
        seen.add(cur)
        nm = name.get(cur)
        chain.append(nm if nm else None)
        nxt = parent.get(cur)
        if not nxt or nxt == cur: break
        cur = nxt
    # chain currently from node up to root. Reverse and clean.
    chain = [c for c in chain if c][::-1]
    # Prefix taxon labels with their taxonomy source.
    return "NCBI; " + "; ".join(chain) + ";"

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kraken-out", required=True)
    ap.add_argument("--nodes", required=True)
    ap.add_argument("--names", required=True)
    ap.add_argument("--out-species", required=True)
    ap.add_argument("--out-genus", required=True)
    ap.add_argument("--out-all", required=True)
    return ap.parse_args()

def main():
    a = parse_args()
    parent, rank, name = load_taxdump(a.nodes, a.names)

    f_species = open(a.out_species, "w")
    f_genus   = open(a.out_genus, "w")
    f_all     = open(a.out_all, "w")

    with open_auto(a.kraken_out) as f:
        for line in f:
            # Kraken2 output format: <C/U>\t<read_id>\t<taxid>\t...
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3: continue
            status, rid, taxid = parts[0], parts[1], parts[2]

            if status == "U" or taxid == "0":
                # unclassified
                f_all.write(f"{rid}\t?\tNCBI;\n")
                continue

            # prefer species; fallback to genus; else unknown
            sid = ascend_to_rank(taxid, "species", parent, rank)
            gid = ascend_to_rank(taxid, "genus", parent, rank)

            wrote_any = False
            if gid:
                gp = lineage_path(gid, parent, name)
                f_genus.write(f"{rid}\tG\t{gp}\n")
                wrote_any = True
            if sid:
                sp = lineage_path(sid, parent, name)
                f_species.write(f"{rid}\tS\t{sp}\n")
                # and prefer species for ALL
                f_all.write(f"{rid}\tS\t{sp}\n")
            else:
                # No species; use genus in ALL if we had it
                if gid:
                    gp = lineage_path(gid, parent, name)
                    f_all.write(f"{rid}\tG\t{gp}\n")
                elif not wrote_any:
                    f_all.write(f"{rid}\t?\tNCBI;\n")

    for h in (f_species, f_genus, f_all):
        h.close()

if __name__ == "__main__":
    main()
