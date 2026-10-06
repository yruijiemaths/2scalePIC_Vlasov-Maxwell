#!/usr/bin/env bash
set -euo pipefail

method="${PIC_METHOD:-TSI-Sim}"
grid="${PIC_GRID:-auto}"
build_dir="${BUILD_DIR:-build}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --method) method="$2"; shift 2 ;;
    --grid) grid="$2"; shift 2 ;;
    --build-dir) build_dir="$2"; shift 2 ;;
    --help)
      echo "Usage: $0 [--method NAME] [--grid auto|standard|half-grid] [--build-dir DIR]"
      exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

cmake -S . -B "$build_dir" -DPIC_METHOD="$method" -DPIC_GRID="$grid" -DCMAKE_BUILD_TYPE="${BUILD_TYPE:-Release}"
cmake --build "$build_dir" --parallel

echo "Built $method on $grid grid: $build_dir/bin/pic-vm.exe"

