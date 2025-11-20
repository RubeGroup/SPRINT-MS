import argparse
import pandas as pd
import json
import numpy as np

parser = argparse.ArgumentParser(description="Extract columns from CSV using flat JSON config")

parser.add_argument("input_path", help="Path to the input CSV file")
parser.add_argument("--json", required=True, help="Path to JSON file with column lists")
parser.add_argument("--key", required=True, help="Key in JSON to extract (e.g., 'intensity', 'lfq')")
parser.add_argument("--output", required=False, help="Output file name (defaults to <key>.csv)")

args = parser.parse_args()

with open(args.json, 'r') as f:
    config = json.load(f)

columns_to_extract = config[args.key]

df_temp = pd.read_csv(args.input_path, sep='\t').replace(np.nan, 0.)

signal = df_temp[columns_to_extract]


output_path = args.output #or f"{args.key}.csv"
signal.to_csv(output_path, sep='\t')

