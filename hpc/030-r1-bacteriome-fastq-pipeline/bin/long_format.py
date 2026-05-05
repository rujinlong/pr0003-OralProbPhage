#!/usr/bin/env python3
"""Convert merge_metaphlan_tables.py wide TSV to long format suitable for DBI ingest."""
import sys
import pandas as pd

src, dst, bioproject = sys.argv[1], sys.argv[2], sys.argv[3]
df = pd.read_csv(src, sep="\t", skiprows=1)
df = df.rename(columns={df.columns[0]: "species_path"})
df["species"] = df["species_path"].str.split("|").str[-1]

long_df = df.melt(id_vars=["species_path", "species"],
                  var_name="run_id", value_name="rel_abund")
long_df = long_df[long_df["rel_abund"] > 0]
long_df["run_id"] = long_df["run_id"].str.replace(r"\.metaphlan$", "", regex=True)
long_df["bioproject"] = bioproject
long_df.to_csv(dst, sep="\t", index=False)
print(f"wrote {dst} ({len(long_df)} non-zero rows)")
