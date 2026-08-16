#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np

def main(input_file, output_file):
    records = []

    with open(input_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split()

            # 不足 5 列：无耦合信息，跳过
            if len(parts) < 5:
                continue

            group1, group2, transcript, direction, flag = parts[:5]

            try:
                flag = int(flag)
            except ValueError:
                continue

            records.append({
                "group1": group1,
                "group2": group2,
                "transcript": transcript,
                "flag": flag
            })

    if not records:
        raise ValueError("No valid records found.")

    df = pd.DataFrame(records)

    # === 在 group1 内，对 transcript 计算 final_score ===
    summary = (
        df.groupby(["group1", "transcript"])["flag"]
        .agg(
            P=lambda x: (x == 1).sum(),
            N=lambda x: (x == -1).sum()
        )
        .reset_index()
    )

    summary["T"] = summary["P"] + summary["N"]

    summary["direction_score"] = (summary["P"] - summary["N"]) / summary["T"]
    summary["final_score"] = summary["direction_score"] * np.log2(summary["T"] + 1)

    # === 转成热图矩阵：行=group1，列=transcript ===
    matrix = summary.pivot(
        index="group1",
        columns="transcript",
        values="final_score"
    ).fillna(0)

    matrix.to_csv(output_file, sep="\t")

    print(f"[OK] Final score matrix written to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate APA–expression coupling final score matrix (group1 × transcript)"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input Total.updown_elongeshort.list file"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output matrix file (TSV)"
    )

    args = parser.parse_args()
    main(args.input, args.output)

