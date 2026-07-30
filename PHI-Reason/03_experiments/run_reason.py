#!/usr/bin/env python3
"""Run one ablation/perturbation experiment through the frozen backbone.

Reads pre-assembled prompt files from  <PHI_EXP_ROOT>/<exp>/data/inputs/*.txt ,
queries the model via the Ollama HTTP API, caches per-phage responses, and writes
species/genus ranking metrics to  <PHI_EXP_ROOT>/<exp>/results/ .

Usage:
    python run_reason.py <exp_name> [http://host:port]
    OLLAMA_SEED=43 python run_reason.py E_prt_P2_blastn_scram      # fix a decode seed

Environment overrides:
    PHI_EXP_ROOT   root holding experiment folders        (default: current dir)
    PHI_GOLD_CSV   gold mapping: phage,host or phage_id,host_species[,split] (with or
                   without header); species names underscore-joined  (default: data/cherry_phage_host_pair.csv)
    OLLAMA_URLS    comma-separated Ollama endpoints        (default: 11435,11436)
    OLLAMA_MODEL   model tag                                (default: qwen3-coder-next:q4_K_M)
    OLLAMA_SEED    decode seed (0 = backend default)        (default: 0)
"""
import sys, os, json, re, asyncio, aiohttp, pathlib
from json import JSONDecoder

ROOT = pathlib.Path(os.environ.get("PHI_EXP_ROOT", "."))
GOLD_CSV = pathlib.Path(os.environ.get("PHI_GOLD_CSV", "data/cherry_phage_host_pair.csv"))
MODEL = os.environ.get("OLLAMA_MODEL", "qwen3-coder-next:q4_K_M")
NUM_CTX, NUM_PREDICT = 40960, 4096
URLS = os.environ.get("OLLAMA_URLS", "http://127.0.0.1:11435,http://127.0.0.1:11436").split(",")
CONC_PER = 6
SEED = int(os.environ.get("OLLAMA_SEED", "0"))

exp = sys.argv[1]
if len(sys.argv) > 2 and sys.argv[2].startswith("http"):
    URLS = [sys.argv[2]]; CONC_PER = 8
EXPDIR = ROOT / exp
IND, RES = EXPDIR / "data" / "inputs", EXPDIR / "results"
CACHE = RES / "cache"
CACHE.mkdir(parents=True, exist_ok=True)

def load_gold(path):
    """Read phage->host, tolerating a header and an optional trailing split column."""
    g = {}
    for l in open(path):
        l = l.strip()
        if not l:
            continue
        p = l.split(",")
        if len(p) < 2 or p[0].lower() in ("phage", "phage_id") or p[1].lower() in ("host", "host_species"):
            continue  # header / malformed
        g[p[0].strip()] = p[1].strip()   # column 2 = host; ignore any split column
    return g


gold = load_gold(GOLD_CSV)
sample = next(IND.glob("*.txt")).read_text()
cand = sample.split("=== USER MESSAGE ===")[1]
ALL_HOSTS = list(dict.fromkeys(
    re.findall(r'^\s*\d+\.\s+(\S+_\S+)\s*(?:\[[^\]]+\])?\s*(?:\||$)', cand, re.M)))


def clean_host(x):
    return re.sub(r'\s*\[[^\]]*\]', '', str(x)).strip().replace(" ", "_").strip("_")


print(f"[{exp}] seed={SEED} candidates={len(ALL_HOSTS)} gold={len(gold)}", flush=True)


def parse_response(raw, gen_tokens):
    text = re.sub(r"<think>.*?</think>", "", raw.strip(), flags=re.DOTALL).strip()
    for t in ("<|endoftext|>", "<|im_start|>", "<|im_end|>"):
        if t in text:
            text = text.split(t, 1)[0].strip()
    scores = {h: 0.0 for h in ALL_HOSTS}
    cutoff = gen_tokens >= NUM_PREDICT - 50
    perr = False
    js = text.find("{"); payload = None
    if js != -1:
        try:
            payload, _ = JSONDecoder().raw_decode(text[js:])
        except Exception:
            perr = True
            for m in re.finditer(r'"host"\s*:\s*"([^"]+)"[^}]*?"score"\s*:\s*([\d.]+)', text[js:]):
                h = clean_host(m.group(1))
                try:
                    s = float(m.group(2))
                    if h in scores:
                        scores[h] = max(scores[h], min(1.0, max(0.0, s)))
                except Exception:
                    pass
    else:
        perr = True
    if isinstance(payload, dict):
        for row in payload.get("predictions", []):
            try:
                h = clean_host(row.get("host", "")); s = float(row.get("score", 0.0))
                if h in scores:
                    scores[h] = max(0.0, min(1.0, s))
            except Exception:
                pass
    reasoning = payload.get("reasoning", "") if isinstance(payload, dict) else ""
    return scores, cutoff, perr, reasoning


sem = {u: asyncio.Semaphore(CONC_PER) for u in URLS}


async def infer(session, acc):
    cf = CACHE / f"{acc}.json"
    if cf.exists():
        try:
            return acc, json.loads(cf.read_text())["scores"]
        except Exception:
            pass
    txt = (IND / f"{acc}.txt").read_text()
    system = txt.split("=== SYSTEM PROMPT ===", 1)[1].split("=== USER MESSAGE ===", 1)[0].strip()
    user = txt.split("=== USER MESSAGE ===", 1)[1].strip()
    url = URLS[hash(acc) % len(URLS)]
    payload = {"model": MODEL, "prompt": user, "system": system, "stream": False, "keep_alive": "120m",
               "options": {"temperature": 0.1, "seed": SEED, "num_predict": NUM_PREDICT, "num_ctx": NUM_CTX,
                           "num_gpu": 99, "stop": ["<|endoftext|>", "<|im_start|>", "<|im_end|>"]}}
    async with sem[url]:
        for attempt in range(3):
            try:
                async with session.post(f"{url}/api/generate", json=payload,
                                        timeout=aiohttp.ClientTimeout(total=1800)) as r:
                    d = await r.json()
                raw = d.get("response", ""); gt = d.get("eval_count", 0)
                scores, cutoff, perr, reasoning = parse_response(raw, gt)
                cf.write_text(json.dumps({"phage": acc, "scores": scores, "cutoff": cutoff,
                                          "gen_tokens": gt, "reasoning": reasoning}))
                return acc, scores
            except Exception as e:
                if attempt == 2:
                    print(f"  FAIL {acc}: {repr(e)[:80]}", flush=True)
                    return acc, {h: 0.0 for h in ALL_HOSTS}
                await asyncio.sleep(3)


async def main():
    accs = [f.stem for f in IND.glob("*.txt")]
    conn = aiohttp.TCPConnector(limit=0)
    async with aiohttp.ClientSession(connector=conn) as s:
        res = {}; done = 0
        for fut in asyncio.as_completed([infer(s, a) for a in accs]):
            a, sc = await fut; res[a] = sc; done += 1
            if done % 50 == 0:
                print(f"  {done}/{len(accs)}", flush=True)
    s1 = s5 = s10 = g1 = g5 = n = 0; mrr = 0.0; rows = []
    for a, sc in res.items():
        if a not in gold:
            continue
        t = gold[a]; tg = t.split("_")[0]; n += 1
        order = [h for h, _ in sorted(sc.items(), key=lambda kv: (-kv[1], kv[0]))]
        rank = order.index(t) + 1 if t in order else 999
        grank = min([i + 1 for i, h in enumerate(order) if h.split("_")[0] == tg] or [999])
        s1 += rank <= 1; s5 += rank <= 5; s10 += rank <= 10
        g1 += grank <= 1; g5 += grank <= 5; mrr += 1.0 / rank
        rows.append((a, t, order[0], rank, grank))
    m = {"exp": exp, "model": MODEL, "seed": SEED, "n": n, "num_ctx": NUM_CTX,
         "species_top1": round(s1 / n, 4), "species_top5": round(s5 / n, 4), "species_top10": round(s10 / n, 4),
         "genus_top1": round(g1 / n, 4), "genus_top5": round(g5 / n, 4), "species_mrr": round(mrr / n, 4)}
    json.dump(m, open(RES / "metrics.json", "w"), indent=1)
    with open(RES / "top_predictions.tsv", "w") as f:
        f.write("phage\ttrue\tpred\tsp_rank\tgenus_rank\n")
        for a, t, p, r, gr in sorted(rows):
            f.write(f"{a}\t{t}\t{p}\t{r}\t{gr}\n")
    print(f"[{exp}] sp@1={m['species_top1']} genus@1={m['genus_top1']} (n={n}, seed={SEED})", flush=True)


asyncio.run(main())
