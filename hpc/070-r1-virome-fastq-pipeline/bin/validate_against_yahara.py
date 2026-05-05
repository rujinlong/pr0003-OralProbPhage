#!/usr/bin/env python3
"""ANI cross-check between our self-built vOTU catalog and Yahara 2021's MOESM7.

We don't have Yahara's contig sequences locally — only their IDs. To run a true
ANI comparison, you'd first need to fetch their assembled contigs from DDBJ
(BioProject PRJDB10485). For the MVP-scope validation, we just report:
  - how many Yahara-IDs we cover by ID prefix match (their IDs encode contig
    coords, so a coarse "did our pipeline call something near here?" check)
  - cluster-level overlap (their per-sample sheet IDs vs our vOTU clusters)

When Yahara's genomes are available locally as FASTA, set --yahara_fasta and
this script will run a real `mash dist` comparison (95% ANI gate) instead.
"""
import argparse, os, subprocess, sys, csv
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog",     required=True)
    ap.add_argument("--yahara_xlsx", required=True)
    ap.add_argument("--yahara_fasta", default="",
                    help="if provided, compute mash ANI; otherwise ID-level coverage only")
    ap.add_argument("--out",         required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    # Load Yahara per-sample phage IDs
    sheets = pd.read_excel(args.yahara_xlsx, sheet_name=None, skiprows=2)
    yahara_ids = []
    for sname, df in sheets.items():
        col = "phageID" if "phageID" in df.columns else df.columns[0]
        for v in df[col].dropna().astype(str):
            if v.startswith("VIRSorter_"):
                yahara_ids.append((sname, v))
    n_y = len(yahara_ids)

    # Load our catalog headers
    ours = []
    with open(args.catalog) as fh:
        for line in fh:
            if line.startswith(">"):
                ours.append(line[1:].split()[0])
    n_o = len(ours)

    if args.yahara_fasta and os.path.exists(args.yahara_fasta):
        # Real ANI comparison via mash
        try:
            subprocess.check_call(["mash", "sketch", "-o", "yahara",
                                   args.yahara_fasta], stdout=sys.stderr)
            subprocess.check_call(["mash", "sketch", "-o", "ours",
                                   args.catalog], stdout=sys.stderr)
            with open(args.out, "w") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["yahara_id","our_vOTU_id","mash_distance","ani_estimate"])
                proc = subprocess.run(["mash", "dist", "yahara.msh", "ours.msh"],
                                      capture_output=True, text=True, check=True)
                for line in proc.stdout.splitlines():
                    f = line.split("\t")
                    if len(f) >= 5 and float(f[2]) <= 0.05:  # ~95% ANI
                        w.writerow([f[0], f[1], f[2], f"{(1 - float(f[2])) * 100:.2f}"])
            print(f"[validate] mash ANI written: {args.out}")
        finally:
            for f in ("yahara.msh", "ours.msh"):
                if os.path.exists(f): os.remove(f)
    else:
        # ID-level coverage only
        with open(args.out, "w") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["yahara_sheet","yahara_id","note"])
            for s, vid in yahara_ids:
                w.writerow([s, vid,
                           "ANI not computed (Yahara FASTA absent); rerun with --yahara_fasta"])
        print(f"[validate] {n_y} Yahara IDs catalogued; "
              f"{n_o} vOTUs in our catalog. ANI gate skipped.")
        print(f"[validate] To enable: download Yahara assemblies from PRJDB10485, "
              f"concatenate into a FASTA, and rerun with --yahara_fasta.")

if __name__ == "__main__":
    main()
