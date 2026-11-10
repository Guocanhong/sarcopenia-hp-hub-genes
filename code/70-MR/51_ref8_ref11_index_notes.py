# -*- coding: utf-8 -*-
"""51_ref8_ref11_index_notes.py v2 — real-character anchors."""
import io

p = r"D:/GEO-limma,肌少症1/70-MR/paper/References_draft.md"
s = io.open(p, encoding="utf-8").read()

pairs = [
    ("verified via CrossRef and the publisher page, 2026-11-07; volume/article number "
     "not yet assigned on the publisher citation line \u2014 complete at proof stage)*",
     "verified via CrossRef and the publisher page, 2026-11-07; volume/article number not "
     "yet assigned on the publisher citation line \u2014 complete at proof stage; "
     "PubMed/Europe PMC searches on 2026-11-08 return no record \u2014 PubMed indexing of "
     "this recently published article is still pending, so no PMID exists yet)*"),
    ("page; preprint \u2014 not yet peer-reviewed published)*",
     "page and required by the API terms of use; bioRxiv preprint \u2014 not peer-review "
     "published and therefore not indexed in PubMed, by design; DOI verified via CrossRef: "
     "posted 2020-08-10, 14 authors)*"),
]

n = 0
for old, new in pairs:
    if new in s and old not in s:
        print("SKIP (applied):", old[:55])
        continue
    assert old in s, "NOT FOUND: " + old[:70]
    assert s.count(old) == 1, "NOT UNIQUE: " + old[:70]
    s = s.replace(old, new, 1)
    n += 1

io.open(p, "w", encoding="utf-8", newline="").write(s)
print("patched:", n)
