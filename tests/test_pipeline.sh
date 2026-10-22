#!/usr/bin/env bash
set -euo pipefail
T=$(mktemp -d); python tests/make_toy.py "$T"
python src/run_ml.py --bracken "$T/bracken_out" --map "$T/sample_id_map.csv" --s3 "$T/s3.xlsx" --out "$T/res" --repeats 2 --perms 3
test -s "$T/res/summary.json" && test -s "$T/res/shap_taxonomy.csv" && echo "OK"
python src/extra_analyses.py --bracken "$T/bracken_out" --map "$T/sample_id_map.csv" --s3 "$T/s3.xlsx" --out "$T/res" --repeats 2
test -s "$T/res/extra.json" && test -s "$T/res/learning_curve.csv" && echo "OK extra"
python src/make_figures.py --bracken "$T/bracken_out" --map "$T/sample_id_map.csv" --s3 "$T/s3.xlsx" --res "$T/res" --out "$T/fig"
test -s "$T/fig/roc.png" && echo "OK figures"
python src/sensitivity.py --bracken "$T/bracken_out" --map "$T/sample_id_map.csv" --s3 "$T/s3.xlsx" --res "$T/res" --repeats 2 && echo "OK sensitivity"
python src/compare_models.py --bracken "$T/bracken_out" --map "$T/sample_id_map.csv" --s3 "$T/s3.xlsx" --res "$T/res" --repeats 2 && echo "OK models"
python src/nested_ratio.py --bracken "$T/bracken_out" --map "$T/sample_id_map.csv" --s3 "$T/s3.xlsx" --res "$T/res" --repeats 2 --perms 3 && echo "OK nested_ratio"
python src/audit_evaluation.py --bracken "$T/bracken_out" --map "$T/sample_id_map.csv" --s3 "$T/s3.xlsx" --res "$T/res" --out "$T/res" --quick > /dev/null && test -s "$T/res/audit_evaluation.json" && echo "OK audit"
python src/freeze_model.py --bracken "$T/bracken_out" --map "$T/sample_id_map.csv" --s3 "$T/s3.xlsx" --out "$T/frozen.joblib" && python - "$T" <<'PY'
import sys, pandas as pd
t = sys.argv[1]; m = pd.read_csv(f"{t}/sample_id_map.csv")[["Run", "Group"]]; m.to_csv(f"{t}/labels.csv", index=False)
PY
python src/validate_external.py --model "$T/frozen.joblib" --bracken "$T/bracken_out" --labels "$T/labels.csv" --out "$T/ext" --perms 200 > /dev/null && test -s "$T/ext/external_validation.json" && echo "OK external"
