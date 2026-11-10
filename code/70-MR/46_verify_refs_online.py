# -*- coding: utf-8 -*-
# 46_verify_refs_online.py -- verify all 38 references against authoritative
# sources (NCBI eutils esummary for PMIDs, CrossRef for DOIs). Read-only.
# Compares: journal abbreviation, year, volume, issue, pages, title, first
# author surname. Prints PASS/DRIFT per field. Never edits files.
import io, re, json, time, urllib.request, sys

REF = r"D:/GEO-limma,肌少症1/70-MR/paper/References_draft.md"
s = io.open(REF, encoding="utf-8").read()

refs = {}
for m in re.finditer(r"^(\d+)\.\s+(.+)$", s, re.M):
    n = int(m.group(1))
    if 1 <= n <= 39 and n not in refs:
        refs[n] = m.group(2)

assert len(refs) == 39, "expected 39 refs, got %d" % len(refs)

def grab(txt, pat):
    m = re.search(pat, txt)
    return m.group(1).strip() if m else None

pmids, dois = {}, {}
for n, t in refs.items():
    pmids[n] = grab(t, r"PMID\s+(\d+)")
    dois[n] = grab(t, r"DOI\s+(\S+?)\s*(?:\*|$)")

# ---- first-pass parse of our own fields from the citation string ----
def cite_fields(t):
    # journal: token before ". YEAR;" pattern -> "J abbrev. YEAR;VOL(ISS):PAGES"
    m = re.search(r"([A-Z][^.]*?)\.\s*(\d{4})", t)
    jr = m.group(1).strip() if m else None
    yr = m.group(2) if m else None
    vi = re.search(r";(\d+)(?:\(([^)]+)\))?[:.]", t)
    vol = vi.group(1) if vi else None
    iss = vi.group(2) if vi else None
    pg = re.search(r"(?:\(\d+\))?[:.](\d+(?:[–-]\d+)?|e\d+)\.?\s", t + " ")
    # more robust: look for :pages or .pages after vol(issue)
    pg = None
    m2 = re.search(r"\)\s*[:.]?\s*([e]?\d+(?:[–-]\d+)?)\s*\.", t)
    if m2: pg = m2.group(1)
    else:
        m3 = re.search(r";(\d+):([e]?\d+[–-]?\d*)\.", t)
        if m3: pg = m3.group(2)
    au = re.match(r"[^.]+\.\s", t)
    first_author = re.match(r"^([A-ZÀ-ž' -]+?),", t)
    fa = first_author.group(1) if first_author else None
    # title: text between first author block ". " and ". JOURNAL"
    tm = re.match(r"^[^.]+\.\s+(.+?)\.\s+[A-Z][A-Za-z ]+?\.\s*\d{4}", t)
    title = tm.group(1) if tm else None
    return dict(journal=jr, year=yr, vol=vol, issue=iss, pages=pg,
                first_author=fa, title=title)

def norm_title(x):
    x = re.sub(r"\*|et al\.?$", "", x or "")
    return re.sub(r"\s+", " ", x).strip().lower()

def http_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": "ref-check/1.0"})
    return json.load(urllib.request.urlopen(req, timeout=45))

# ---- eutils esummary for all PMIDs (one batch call) ----
ids = [pmids[n] for n in sorted(refs) if pmids[n]]
print("PMIDs present:", len(ids), "| DOIs present:", len([d for d in dois.values() if d]))
es = http_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=" + ",".join(ids) + "&retmode=json")["result"]

issues = []
for n in sorted(refs):
    t = refs[n]
    ours = cite_fields(t)
    pm = pmids[n]
    if not pm:
        issues.append((n, "NO-PMID (verify via DOI/CrossRef)"))
        continue
    r = es[pm]
    # check first author surname
    a0 = r["authors"][0]["name"]  # "Surname Initials"
    fa_ours = (ours["first_author"] or "").replace(" ", "")
    fa_pub = a0.split(" ")[0].replace(" ", "")
    fa_ok = fa_pub.lower().startswith(re.sub(r"[^A-Za-zÀ-ž]", "", fa_ours)[:4].lower())
    # title match (our title = text after authors up to journal)
    title_ours = norm_title(ours["title"])
    title_pub = norm_title(r["title"])
    title_ok = title_ours[:40] in title_pub or title_pub[:40] in title_ours
    # journal
    jr_pub = r.get("source") or ""
    jr_ok = jr_pub.lower().startswith((ours["journal"] or "")[:6].lower())
    # year
    yr_pub = (r.get("pubdate") or "")[:4]
    yr_ok = yr_pub == ours["year"]
    # vol/issue/pages
    vol_ok = str(r.get("volume") or "") == str(ours["vol"] or "")
    iss_ok = (str(r.get("issue") or "") == str(ours["issue"] or "")) if ours["issue"] else True
    pg_pub = (r.get("pages") or "").replace("-", "–")
    pg_ours = (ours["pages"] or "").replace("-", "–")
    pg_ok = pg_ours[:5] in pg_pub or pg_pub[:5] in pg_ours if pg_ours else True
    flags = []
    if not fa_ok: flags.append("AUTHOR1 ours=%s pub=%s" % (fa_ours, fa_pub))
    if not title_ok: flags.append("TITLE ours=%s... pub=%s..." % (title_ours[:45], title_pub[:45]))
    if not jr_ok: flags.append("JOURNAL ours=%s pub=%s" % (ours["journal"], jr_pub))
    if not yr_ok: flags.append("YEAR ours=%s pub=%s" % (ours["year"], yr_pub))
    if not vol_ok: flags.append("VOL ours=%s pub=%s" % (ours["vol"], r.get("volume")))
    if not iss_ok: flags.append("ISSUE ours=%s pub=%s" % (ours["issue"], r.get("issue")))
    if not pg_ok: flags.append("PAGES ours=%s pub=%s" % (pg_ours, pg_pub))
    if flags:
        issues.append((n, "; ".join(flags)))
    else:
        print("PASS ref %2d  (%s)" % (n, jr_pub))

print("\n==== ISSUES ====")
for n, msg in issues:
    print("ref %2d: %s" % (n, msg))
