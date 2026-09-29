#!/usr/bin/env bash
set -euo pipefail
workflow_dir="$(cd "$(dirname "$0")" && pwd)"
step=${1:?选择分析步骤：10_monocle2|11_monocle3|12_slingshot|13_tradeseq}
shift
case "$step" in
  10_monocle2|11_monocle3|12_slingshot|13_tradeseq) Rscript "$workflow_dir/$step.R" "$@" ;;
  *) printf '%s\n' '未知步骤' >&2; exit 2 ;;
esac
