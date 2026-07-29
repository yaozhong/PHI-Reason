#!/usr/bin/env python3
"""Derive the ablation (6) and perturbation (5) experiments from the FULL profile set.

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
import re, random, os, pathlib

ROOT = pathlib.Path(os.environ.get("PHI_EXP_ROOT", "."))
FULL_EXP = os.environ.get("PHI_FULL_EXP", "E_reason_full")
SRC = ROOT / FULL_EXP / "data" / "inputs"          # FULL prompt files
CANDHDR = "=== 223 CANDIDATE HOSTS ==="
SEED = 42

sample = open(next(iter(SRC.glob("*.txt")))).read()
CANDS = re.findall(r'^\s*\d+\.\s+(\S+_\S+)\s*(?:\||$)',
                   sample.split("=== USER MESSAGE ===")[1], re.M)

BN = "## Phylogenetic Cluster Context (whole-genome BLASTN"
KM = "## Genome-Level Shared 25-mers"
CR = "## CRISPR Spacer Matches"


def strip_section(t, hdr):
    i = t.find(hdr)
    if i < 0:
        return t
    nxt = re.search(r"\n## ", t[i + len(hdr):])
    j = i + len(hdr) + nxt.start() + 1 if nxt else t.find(CANDHDR, i)
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
    j = i + len(hdr) + nxt.start() + 1 if nxt else t.find(CANDHDR, i)
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
        rng = random.Random(SEED + hash(f.stem) % 100000)
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
