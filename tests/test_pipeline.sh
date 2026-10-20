#!/usr/bin/env bash
set -euo pipefail
T=$(mktemp -d); python tests/make_toy.py "$T"
python src/run_ml.py --bracken "$T/bracken_out" --map "$T/sample_id_map.csv" --s3 "$T/s3.xlsx" --out "$T/res" --repeats 2 --perms 3
test -s "$T/res/summary.json" && test -s "$T/res/shap_taxonomy.csv" && echo "OK"
python src/extra_analyses.py --bracken "$T/bracken_out" --map "$T/sample_id_map.csv" --s3 "$T/s3.xlsx" --out "$T/res" --repeats 2
test -s "$T/res/extra.json" && test -s "$T/res/learning_curve.csv" && echo "OK extra"
