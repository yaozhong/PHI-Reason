#!/bin/bash
# Build then run the 6 ablation + 5 perturbation experiments, then analyze.
#
# Prerequisites:
#   1. Ollama running with the backbone loaded (see OLLAMA_URLS / OLLAMA_MODEL).
#   2. The FULL prompt set built at  ${PHI_EXP_ROOT}/${PHI_FULL_EXP}/data/inputs/*.txt
#      (produced by 01_profile_generation + 02_inference on the full evidence profile).
#   3. A gold mapping at ${PHI_GOLD_CSV}: phage,host or phage_id,host_species[,split]
#      (with or without header; default data/cherry_phage_host_pair.csv).
#
# Usage:
#   export PHI_EXP_ROOT=/path/to/experiments      # where E_reason_full/ lives
#   bash run_ablation_perturb.sh
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PHI_EXP_ROOT="${PHI_EXP_ROOT:-$(pwd)}"
export PHI_FULL_EXP="${PHI_FULL_EXP:-E_reason_full}"
PY="${PYTHON:-python3}"

echo "[build] deriving ablation/perturbation prompts from ${PHI_FULL_EXP}"
$PY "${SCRIPT_DIR}/build_ablation_perturb.py"

for e in E_abl_L0_base E_abl_L1_rbp E_abl_L2_rbp_blastn \
         E_abl_O1_no_rbp E_abl_O2_no_blastn E_abl_O3_no_25mer E_abl_O4_no_crispr \
         E_prt_P1_labelfree E_prt_P2_blastn_scram E_prt_P3_25mer_scram \
         E_prt_P4_crispr_scram E_prt_P5_all_scram; do
  echo "[run] $e"
  $PY "${SCRIPT_DIR}/run_reason.py" "$e"
done

echo "[analyze]"
$PY "${SCRIPT_DIR}/analyze_results.py"
