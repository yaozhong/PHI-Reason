#!/usr/bin/env python3
"""Derive the ablation (7) and perturbation (5) experiments (12 variants) from the FULL profile set.

Ablations remove one evidence channel at a time (additive ladder + leave-one-out);
perturbations corrupt the host label inside one channel (scramble / label-free) while
keeping every other field intact. The frozen model and prompt are unchanged; only the
named evidence blocks in the pre-assembled prompt files are edited.

Evidence blocks (matched by header in the assembled prompt):
    BLASTN  '## Phylogenetic Cluster Context (whole-genome BLASTN'   host: "infects X"
    25-mer  '## Genome-Level Shared 25-mers'                          host: "- Host_ (...)"
    CRISPR  '## CRISPR Spacer Matches'                                host: "- Host_ (...)"
    RBP     inline gene-line suffix ' <- RBP matches (BLASTP): ...'

Paths (override via environment):
    PHI_EXP_ROOT   directory holding the experiment folders (default: current dir)
    PHI_FULL_EXP   name of the FULL source experiment  (default: E_reason_full)

Each experiment writes prompt files to  <PHI_EXP_ROOT>/<name>/data/inputs/*.txt .
"""
import re, random, os, hashlib, pathlib

ROOT = pathlib.Path(os.environ.get("PHI_EXP_ROOT", "."))
FULL_EXP = os.environ.get("PHI_FULL_EXP", "E_reason_full")
SRC = ROOT / FULL_EXP / "data" / "inputs"          # FULL prompt files
CAND_RE = re.compile(r'^=== \d+ CANDIDATE HOSTS ===', re.M)   # candidate-count agnostic
SEED = 42


def cand_pos(t, start=0):
    """Start offset of the candidate-hosts section (any host count), or end of text."""
    m = CAND_RE.search(t, start)
    return m.start() if m else len(t)


def stem_seed(stem):
    """Deterministic per-file seed (hashlib; hash() is salted by PYTHONHASHSEED)."""
    return SEED + int(hashlib.md5(stem.encode()).hexdigest()[:8], 16)


import sys
_files = sorted(SRC.glob("*.txt"))
if not _files:
    sys.exit(f"no FULL prompt files in {SRC} (set PHI_EXP_ROOT / PHI_FULL_EXP)")
sample = _files[0].read_text()
if "=== USER MESSAGE ===" not in sample:
    sys.exit(f"FULL prompt missing '=== USER MESSAGE ===' marker: {_files[0]}")
CANDS = re.findall(r'^\s*\d+\.\s+(\S+_\S+)\s*(?:\||$)',
                   sample.split("=== USER MESSAGE ===")[1], re.M)

BN = "## Phylogenetic Cluster Context (whole-genome BLASTN"
KM = "## Genome-Level Shared 25-mers"
CR = "## CRISPR Spacer Matches"

# Warn loudly if an evidence channel is absent from ALL FULL prompts: its
# ablation/perturbation would otherwise be a silent no-op.
for _name, _marker in [("BLASTN", BN), ("25-mer", KM), ("CRISPR", CR), ("RBP", "← RBP matches")]:
    if not any(_marker in f.read_text() for f in _files):
        print(f"WARNING: evidence channel {_name} ('{_marker}') absent from all "
              f"{len(_files)} FULL prompts; its ablation/perturbation will be a no-op")


def strip_section(t, hdr):
    i = t.find(hdr)
    if i < 0:
        return t
    nxt = re.search(r"\n## ", t[i + len(hdr):])
    j = i + len(hdr) + nxt.start() + 1 if nxt else cand_pos(t, i)
    end = j
    while end < len(t) and t[end] == "\n":
        end += 1
    start = i
    while start > 0 and t[start - 1] == "\n":
        start -= 1
    return t[:start] + "\n\n" + t[end:]


def strip_rbp(t):
    return re.sub(r'\s*←\s*RBP matches \(BLASTP\):[^\n]*', '', t)


def block_span(t, hdr):
    i = t.find(hdr)
    if i < 0:
        return None
    nxt = re.search(r"\n## ", t[i + len(hdr):])
    j = i + len(hdr) + nxt.start() + 1 if nxt else cand_pos(t, i)
    return i, j


def transform_block(t, hdr, fn):
    sp = block_span(t, hdr)
    if not sp:
        return t
    i, j = sp
    return t[:i] + fn(t[i:j]) + t[j:]


def scram_blastn(blk, rng):   # "infects <to end of line>" (multi-word host names ok)
    return re.sub(r'(infects ).*',
                  lambda m: m.group(1) + rng.choice(CANDS).replace("_", " "), blk)


def scram_us(blk, rng):       # 25-mer / CRISPR: "- Host_ (" (underscored single token)
    return re.sub(r'^(\s*- )(\S+_\S+)(\s*\()',
                  lambda m: m.group(1) + rng.choice(CANDS) + m.group(3), blk, flags=re.M)


def lf_blastn(blk):
    return re.sub(r'(infects ).*', r'\1[host label withheld]', blk)


def lf_us(blk):
    return re.sub(r'^(\s*- )(\S+_\S+)(\s*\()', r'\1[host withheld]\3', blk, flags=re.M)


def build(name, fn):
    out = ROOT / name / "data" / "inputs"
    out.mkdir(parents=True, exist_ok=True)
    for f in out.glob("*.txt"):
        f.unlink()
    n = 0
    for f in SRC.glob("*.txt"):
        t = f.read_text()
        rng = random.Random(stem_seed(f.stem))
        (out / f.name).write_text(fn(t, rng))
        n += 1
    print(f"{name}: {n}")


# ===== Ablations (additive ladder L0-L2 + leave-one-out O1-O3) =====
build("E_abl_L0_base",       lambda t, r: strip_rbp(strip_section(strip_section(strip_section(t, BN), KM), CR)))
build("E_abl_L1_rbp",        lambda t, r: strip_section(strip_section(strip_section(t, BN), KM), CR))
build("E_abl_L2_rbp_blastn", lambda t, r: strip_section(strip_section(t, KM), CR))
build("E_abl_O1_no_rbp",     lambda t, r: strip_rbp(t))
build("E_abl_O2_no_blastn",  lambda t, r: strip_section(t, BN))
build("E_abl_O3_no_25mer",   lambda t, r: strip_section(t, KM))
build("E_abl_O4_no_crispr",  lambda t, r: strip_section(t, CR))   # = additive L3 (+25mer, FULL - CRISPR)

# ===== Perturbations (label-free P1 + per-channel scramble P2-P4 + all P5) =====
build("E_prt_P1_labelfree",    lambda t, r: transform_block(transform_block(transform_block(t, BN, lf_blastn), KM, lf_us), CR, lf_us))
build("E_prt_P2_blastn_scram", lambda t, r: transform_block(t, BN, lambda b: scram_blastn(b, r)))
build("E_prt_P3_25mer_scram",  lambda t, r: transform_block(t, KM, lambda b: scram_us(b, r)))
build("E_prt_P4_crispr_scram", lambda t, r: transform_block(t, CR, lambda b: scram_us(b, r)))
build("E_prt_P5_all_scram",    lambda t, r: transform_block(transform_block(transform_block(t, BN, lambda b: scram_blastn(b, r)), KM, lambda b: scram_us(b, r)), CR, lambda b: scram_us(b, r)))
print("done")
