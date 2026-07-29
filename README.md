# PHI-Reason

Code for **PHI-Reason**: a frozen general-purpose LLM used as an evidence-integration
layer over named, perturbable structured biological evidence profiles for closed-set
candidate-conditioned phage–host species ranking.

This repository is the superset of
[`PHI-Reason-Species`](../PHI-Reason-Species): it keeps the same profile-building and
inference pipeline and adds the **evidence-perturbation platform** (`03_experiments/`) —
the ablation and label-perturbation experiments that measure how each named evidence
channel contributes to, is redundant with, or misleads the prediction under a fixed model.

---

## Repository structure

```
PHI-Reason/
├── config.sh                              # Path configuration — edit before first run
├── requirements.txt                       # Python dependencies
│
├── 00_annotation/                         # Upstream genome annotation wrappers
│   ├── run_prodigal.sh                    # Gene prediction (Prodigal)
│   ├── run_phrog.sh                       # PHROG functional annotation (DIAMOND)
│   ├── run_eggnog.sh                      # eggNOG-mapper annotation
│   ├── run_defensefinder.sh               # Defense-system detection (hosts)
│   └── run_blastn.sh                      # Whole-genome BLASTN (train DB → test queries)
│
├── 01_profile_generation/                 # Build phage and host text profiles
│   ├── phage_profiles/                    # base → +RBP (BLASTP) → +BLASTN context, id masking
│   ├── host_profiles/                     # base host profiles, taxonomy/knowledge, candidate list
│   └── blastp_context/                    # BLASTP phylogenetic context (numbered 01→05)
│
├── 02_inference/                          # LLM inference
│   ├── run_inference.py                   # Profile-assembly inference driver (async, resumable)
│   └── prompt/prompt_template.py          # System prompt and user message templates
│
├── 03_experiments/                        # Evidence-perturbation platform (this work)
│   ├── build_ablation_perturb.py          # Derive 6 ablation + 5 perturbation variants from FULL
│   ├── run_reason.py                      # Run one variant through the frozen backbone
│   ├── run_ablation_perturb.sh            # Build → run all 11 → analyze
│   └── analyze_results.py                 # Additive ladder / leave-one-out / perturbation / strata
│
└── data/                                  # Benchmark phage–host pair tables (no profile data)
    ├── cherry_phage_host_pair.csv
    ├── cherry_test_634.txt
    ├── vhdb_phage_host_pair.csv
    └── hic_phage_host_pair.csv
```

`00_annotation/`, `01_profile_generation/` and `02_inference/` are the profile-building
and inference pipeline (identical to `PHI-Reason-Species`). `03_experiments/` is the new
perturbation platform. Generated profile data are **not** included — build them with the
`00`/`01` pipeline (see `PHI-Reason-Species/README.md` for the full profile-building walkthrough).

---

## Evidence-perturbation platform (`03_experiments/`)

Operating on the pre-assembled FULL prompt files, `build_ablation_perturb.py` edits one
named evidence block at a time to produce 11 variants; the model, prompt and decoding are
unchanged.

**Ablations** (remove a channel)

| id | variant | content |
|----|---------|---------|
| L0 | `E_abl_L0_base`       | base only (no RBP / BLASTN / 25-mer / CRISPR) |
| L1 | `E_abl_L1_rbp`        | + RBP |
| L2 | `E_abl_L2_rbp_blastn` | + RBP + BLASTN |
| O4 | `E_abl_O4_no_crispr`  | FULL − CRISPR  (= additive L3, + 25-mer) |
| O1 | `E_abl_O1_no_rbp`     | FULL − RBP |
| O2 | `E_abl_O2_no_blastn`  | FULL − BLASTN |
| O3 | `E_abl_O3_no_25mer`   | FULL − 25-mer |

**Perturbations** (corrupt the host label inside a channel, keep everything else)

| id | variant | edit |
|----|---------|------|
| P1 | `E_prt_P1_labelfree`   | withhold host labels in BLASTN + 25-mer + CRISPR |
| P2 | `E_prt_P2_blastn_scram`| scramble BLASTN neighbour host |
| P3 | `E_prt_P3_25mer_scram` | scramble 25-mer host |
| P4 | `E_prt_P4_crispr_scram`| scramble CRISPR host |
| P5 | `E_prt_P5_all_scram`   | scramble all three |

Scrambling draws a replacement host uniformly from the in-prompt candidate catalogue
(fixed seed 42). All comparisons are apple-to-apple on the same evaluation set.

### Run

```bash
export PHI_EXP_ROOT=/path/to/experiments      # directory containing E_reason_full/
export PHI_GOLD_CSV=/path/to/gold.csv         # phage,host (species-level) per line
# Ollama must be serving the backbone (default endpoints :11435,:11436)
bash 03_experiments/run_ablation_perturb.sh
```

Key environment variables (all optional, sensible defaults): `PHI_EXP_ROOT`,
`PHI_FULL_EXP` (default `E_reason_full`), `PHI_GOLD_CSV`, `OLLAMA_URLS`, `OLLAMA_MODEL`,
`OLLAMA_SEED` (set for decode-variance / multi-seed runs).

---

## Requirements

```bash
pip install -r requirements.txt
```

System dependencies (Ollama, DIAMOND, BLAST+, Prodigal, eggNOG-mapper, PHROG, DefenseFinder)
match `PHI-Reason-Species`; see that repository's README for installation details.
