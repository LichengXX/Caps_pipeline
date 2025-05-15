#!/usr/bin/env bash
# split_by_chrom.sh
#
# 用法：./split_by_chrom.sh caps_candidates.tsv
# 将根据第二列 Chrom 将文件拆分为 Chrom1.tsv、Chrom2.tsv ... 等。
# 每个输出文件都包含原始文件的表头行。

set -euo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: $0 <input_tsv>"
  exit 1
fi

input="$1"
# 提取表头
header=$(head -n 1 "$input")

# 从第二行开始处理每一行
tail -n +2 "$input" | \
while IFS=$'\t' read -r marker chrom pos ref alt enz change seq; do
  out="${chrom}.caps.tsv"
  # 如果首次创建该染色体文件，则写入表头
  if [ ! -f "$out" ]; then
    printf "%s\n" "$header" > "$out"
  fi
  # 将当前行追加到对应文件
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$marker" "$chrom" "$pos" "$ref" "$alt" "$enz" "$change" "$seq" \
    >> "$out"
done

