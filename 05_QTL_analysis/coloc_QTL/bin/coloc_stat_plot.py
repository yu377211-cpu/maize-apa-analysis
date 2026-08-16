#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def main():
    parser = argparse.ArgumentParser(
        description="Plot Fig7A: 3aQTL vs 3aeQTL colocalization"
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input TSV/CSV file")
    parser.add_argument("-o", "--output_prefix", required=True,
                        help="Output file prefix (no extension)")
    parser.add_argument("--sep", default="\t",
                        help="File separator (default: tab)")
    parser.add_argument("--show_percent", action="store_true",
                        help="Show transcript overlap percentage")
    args = parser.parse_args()

    # 读取数据
    df = pd.read_csv(args.input, sep=args.sep)

    # 去掉列名可能的空格
    df.columns = df.columns.str.strip()

    # 自动匹配列名
    col_map = {
        "aQTL_transcript": ["aQTL_transcript", "3aQTL_transcript"],
        "aeQTL_transcript": ["aeQTL_transcript", "3aeQTL_transcript"]
    }

    def get_column(possible_names):
        for name in possible_names:
            if name in df.columns:
                return df[name]
        raise ValueError(f"None of these columns found: {possible_names}")

    tissues = df["Tissue"]
    aQTL = get_column(col_map["aQTL_transcript"])
    aeQTL = get_column(col_map["aeQTL_transcript"])
    tx_overlap = df["transcript_overlap"]
    var_overlap = df["variant_overlap"]

    x = np.arange(len(tissues))
    width = 0.25

    fig, ax1 = plt.subplots(figsize=(10, 6))

    # 主轴柱状图
    bars1 = ax1.bar(x - width, aQTL, width, label="3′aQTL transcripts")
    bars2 = ax1.bar(x, aeQTL, width, label="3′aeQTL transcripts")
    bars3 = ax1.bar(x + width, tx_overlap, width, label="Transcript overlap")

    ax1.set_ylabel("Number of transcripts")
    ax1.set_xticks(x)
    ax1.set_xticklabels(tissues, rotation=45, ha="right")

    # 标注 overlap 数值 + 百分比
    for i in range(len(tx_overlap)):
        value = tx_overlap.iloc[i]
        label = f"{value}"
        if args.show_percent:
            percent = (value / aQTL.iloc[i]) * 100 if aQTL.iloc[i] > 0 else 0
            label = f"{value}\n({percent:.1f}%)"
        ax1.text(x[i] + width,
                 value + max(aQTL)*0.01,
                 label,
                 ha="center",
                 va="bottom",
                 fontsize=8)

    # 副轴折线图
    ax2 = ax1.twinx()
    ax2.plot(x, var_overlap,
             marker='o',
             linestyle='-',
             linewidth=2,
             label="Shared QTL variants")
    ax2.set_ylabel("Number of shared QTL variants")

    # 合并图例
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2,
               loc="upper right",
               frameon=False)

    ax1.set_title("Limited colocalization between 3′aQTLs and 3′aeQTLs")

    fig.tight_layout()

    # 输出
    pdf_file = args.output_prefix + ".pdf"
    png_file = args.output_prefix + ".png"

    plt.savefig(pdf_file)
    plt.savefig(png_file, dpi=300)

    print(f"Saved: {pdf_file}")
    print(f"Saved: {png_file}")


if __name__ == "__main__":
    main()
