#!/usr/bin/env python3
"""Map sylph hits against IMG/VR sketches to host taxonomy.

Sylph emits one TSV row per (query_contig, IMG/VR_reference, ANI). We:
  1. Pick the best IMG/VR hit per query (max containment_index / ANI)
  2. Look up the host taxonomy of that IMG/VR sequence in the metadata TSV
  3. Extract the genus, emit:  vOTU_id, host_genus, host_score, method=sylph
"""
import argparse, csv, sys

def best_per_query(sylph_tsv):
    best = {}
    with open(sylph_tsv) as fh:
        header = None
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                continue
            row = dict(zip(header, f))
            q = row.get("Query_file") or row.get("Query") or row.get("query") or row.get("Genome_file")
            r = row.get("Genome_file") or row.get("Reference") or row.get("reference") or row.get("Contig_name")
            ani = row.get("Adjusted_ANI") or row.get("ANI") or row.get("Containment_ind") or "0"
            try:
                score = float(ani)
            except ValueError:
                score = 0.0
            if not q or not r: continue
            cur = best.get(q)
            if cur is None or score > cur[1]:
                best[q] = (r, score)
    return best

def load_imgvr_host(meta_tsv):
    out = {}
    with open(meta_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # Common IMG/VR column names — try a few
        host_col = None
        for cand in ("Host taxonomy", "Host_taxonomy_prediction",
                     "Host taxonomy prediction", "Host taxonomy (NCBI)"):
            if cand in header:
                host_col = header.index(cand); break
        if host_col is None:
            sys.exit(f"IMG/VR metadata has no host-taxonomy column. Header: {header[:5]}...")
        id_col = 0  # UViG / sequence id is column 0 in standard IMG/VR exports
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > host_col:
                out[f[id_col]] = f[host_col]
    return out

def extract_genus(host_tax):
    if not host_tax: return ""
    for tok in host_tax.split(";"):
        tok = tok.strip()
        if tok.startswith("g__") and len(tok) > 3:
            return tok[3:]
        if tok.lower().startswith("genus:"):
            return tok.split(":",1)[1].strip()
    # fallback: last non-empty token
    parts = [t.strip() for t in host_tax.split(";") if t.strip()]
    return parts[-1] if parts else ""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sylph", required=True)
    ap.add_argument("--imgvr_meta", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    best = best_per_query(args.sylph)
    host_lookup = load_imgvr_host(args.imgvr_meta)

    n_total = len(best); n_hit = 0
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["vOTU_id","host_genus","host_score","method"])
        for q, (r, score) in best.items():
            host_tax = host_lookup.get(r, "")
            genus = extract_genus(host_tax)
            if genus:
                n_hit += 1
                w.writerow([q, genus, f"{score:.4f}", "sylph"])
    print(f"[sylph_to_host] {n_hit}/{n_total} queries got a host genus call")

if __name__ == "__main__":
    main()
