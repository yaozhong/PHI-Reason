#!/usr/bin/env python3
"""Ablation + perturbation analysis (all apple-to-apple on the same n=634).

Reports: the additive ladder (base -> +RBP -> +BLASTN -> +25mer -> +CRISPR),
leave-one-out marginals, label perturbations, and BLASTN-neighbour / CRISPR-coverage
strata. Species-level argmax top-1 throughout. Reads per-phage caches written by
run_reason.py (results/cache/*.json).

Environment overrides:
    PHI_EXP_ROOT   root holding experiment folders        (default: current dir)
    PHI_GOLD_CSV   gold mapping: phage,host or phage_id,host_species[,split] (with or
                   without header)  (default: data/cherry_phage_host_pair.csv)
    PHI_FULL_EXP   name of the FULL experiment            (default: E_reason_full)
"""
import json, glob, re, os, pathlib

ROOT = pathlib.Path(os.environ.get("PHI_EXP_ROOT", "."))
GOLD_CSV = pathlib.Path(os.environ.get("PHI_GOLD_CSV", "data/cherry_phage_host_pair.csv"))
FULLe = os.environ.get("PHI_FULL_EXP", "E_reason_full")


def load_gold(path):
    """Read phage->host, tolerating a header and an optional trailing split column."""
    g = {}
    for l in open(path):
        l = l.strip()
        if not l:
            continue
        p = l.split(",")
        if len(p) < 2 or p[0].lower() in ("phage", "phage_id") or p[1].lower() in ("host", "host_species"):
            continue
        g[p[0].strip()] = p[1].strip()
    return g


gold = load_gold(GOLD_CSV)
SRC = ROOT / FULLe / "data" / "inputs"
# Denominator = phages actually in the experiment (FULL inputs that have a gold label),
# NOT len(gold): a full gold CSV (e.g. Cherry 1,940 = train+test) far exceeds the 634 eval set.
N_TOTAL = sum(1 for f in SRC.glob("*.txt") if f.stem in gold) or len(gold)
BN = "## Phylogenetic Cluster Context (whole-genome BLASTN"
CR = "## CRISPR Spacer Matches"
CAND_RE = re.compile(r'^=== \d+ CANDIDATE HOSTS ===', re.M)   # candidate-count agnostic


def cand_pos(t, start=0):
    m = CAND_RE.search(t, start)
    return m.start() if m else len(t)


def block(t, h):
    i = t.find(h)
    if i < 0:
        return ""
    nx = re.search(r"\n## ", t[i + len(h):])
    j = i + len(h) + nx.start() + 1 if nx else cand_pos(t, i)
    return t[i:j]


# strata labels from the FULL prompts
grpA, hasCR = set(), set()
for f in SRC.glob("*.txt"):
    t = f.read_text()
    if re.search(r'infects ', block(t, BN)):
        grpA.add(f.stem)
    if CR in t and re.search(r'^\s*- ', block(t, CR), re.M):
        hasCR.add(f.stem)


def load(e):
    d = {}
    for f in glob.glob(str(ROOT / e / "results" / "cache" / "*.json")):
        try:
            j = json.load(open(f)); d[j.get("phage") or j.get("phage_id")] = j["scores"]
        except Exception:
            pass
    return d


def acc(d, subset=None):
    s1 = n = 0
    for a, sc in d.items():
        if a not in gold or not sc:
            continue
        if subset is not None and a not in subset:
            continue
        n += 1
        if sorted(sc, key=lambda h: (-sc[h], h))[0] == gold[a]:
            s1 += 1
    return n, (s1 / n if n else 0)


CACHE = {}
def G(e):
    if e not in CACHE:
        CACHE[e] = load(e)
    return CACHE[e]
def done(e):
    return acc(G(e))[0] >= N_TOTAL


FULL = acc(G(FULLe))[1]

print("=" * 60); print(f"[Table 1] Additive ladder  base->+RBP->+BLASTN->+25mer->+CRISPR (n={N_TOTAL})"); print("=" * 60)
prev = None
for t, e, c in [("L0", "E_abl_L0_base", "base only"), ("L1", "E_abl_L1_rbp", "+RBP"),
                ("L2", "E_abl_L2_rbp_blastn", "+BLASTN"), ("L3", "E_abl_O4_no_crispr", "+25mer"),
                ("L4", FULLe, "+CRISPR (=FULL)")]:
    n, s = acc(G(e))
    inc = f"{(s - prev) * 100:+.1f}" if prev is not None else "  -- "
    print(f"  {t} {c:16s} {('%.3f' % s if n >= N_TOTAL else f'{n}/{N_TOTAL} running'):>14s}   d{inc}")
    if n >= N_TOTAL:
        prev = s

print("\n" + "=" * 60); print(f"[Table 2] Leave-one-out marginal = FULL({FULL:.3f}) - drop-channel"); print("=" * 60)
for lab, e, ly in [("O2", "E_abl_O2_no_blastn", "BLASTN"), ("O1", "E_abl_O1_no_rbp", "RBP"),
                   ("O4", "E_abl_O4_no_crispr", "CRISPR"), ("O3", "E_abl_O3_no_25mer", "25mer")]:
    n, s = acc(G(e))
    print(f"  {lab} -{ly:7s} {('%.3f' % s + f'  marginal {(FULL - s) * 100:+.1f}pp') if n >= N_TOTAL else f'{n}/{N_TOTAL} running':32s}")

print("\n" + "=" * 60); print(f"[Table 3] Perturbation vs FULL={FULL:.3f} (n={N_TOTAL})"); print("=" * 60)
for lab, e, d in [("P1", "E_prt_P1_labelfree", "all labels withheld"), ("P2", "E_prt_P2_blastn_scram", "BLASTN scramble"),
                  ("P3", "E_prt_P3_25mer_scram", "25mer scramble"), ("P4", "E_prt_P4_crispr_scram", "CRISPR scramble"),
                  ("P5", "E_prt_P5_all_scram", "all scramble")]:
    n, s = acc(G(e))
    print(f"  {lab} {d:20s} {('%.3f' % s + f'  d{(s - FULL) * 100:+.1f}pp') if n >= N_TOTAL else f'{n}/{N_TOTAL} running':28s}")

print("\n" + "=" * 60); print(f"[Table 4] Strata  Group A (has BLASTN neighbour)={len(grpA)}  B={N_TOTAL - len(grpA)}  |  +/-CRISPR {len(hasCR)}/{N_TOTAL - len(hasCR)}"); print("=" * 60)
for lab, e in [("FULL", FULLe), ("-CRISPR", "E_abl_O4_no_crispr"), ("base", "E_abl_L0_base"), ("-BLASTN", "E_abl_O2_no_blastn")]:
    if not done(e):
        print(f"  {lab:14s} running"); continue
    d = G(e)
    nA, sA = acc(d, grpA); nB, sB = acc(d, set(gold) - grpA)
    nC, sC = acc(d, hasCR); nN, sN = acc(d, set(gold) - hasCR)
    print(f"  {lab:14s} A:{sA:.3f}(n{nA}) B:{sB:.3f}(n{nB}) | +CR:{sC:.3f}(n{nC}) -CR:{sN:.3f}(n{nN})")
