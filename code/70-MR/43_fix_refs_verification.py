# -*- coding: utf-8 -*-
# 43_fix_refs_verification.py -- apply citation metadata verified against
# NCBI eutils (PMID 37086739, 15962564, 27625185) and the Nature/CrossRef
# records for DOI 10.1038/s41467-026-78442-8. Idempotent; verified 2026-11-07.
import io, os

P = r"D:/GEO-limma,肌少症1/70-MR/paper/References_draft.md"

def sub(s, old, new, tag):
    if new in s:
        print("SKIP (already applied):", tag); return s
    assert old in s, "anchor missing: " + tag
    assert s.count(old) == 1, "anchor not unique: " + tag
    print("OK:", tag)
    return s.replace(old, new, 1)

s = io.open(P, encoding="utf-8").read()

# --- ref 3: real record PMID 37086739 (authors 2-3, title, DOI all corrected;
#     vol/issue/pages 8(6):553-564 confirmed)
s = sub(s,
    "3. Li Y, Yin J, Fu W, et al. Global prevalence and trends of *Helicobacter pylori* "
    "infection among asymptomatic individuals: a systematic review and meta-analysis. "
    "Lancet Gastroenterol Hepatol. 2023;8(6):553\u2013564. PMID 37086739. "
    "DOI 10.1016/S2468-1253(23)00050-1 *(journal/order of authors 2\u20133 to confirm at "
    "typesetting from the article PDF)*",
    "3. Li Y, Choi H, Leung K, Jiang F, Graham DY, Leung WK. Global prevalence of "
    "*Helicobacter pylori* infection between 1980 and 2022: a systematic review and "
    "meta-analysis. Lancet Gastroenterol Hepatol. 2023;8(6):553\u2013564. PMID 37086739. "
    "DOI 10.1016/S2468-1253(23)00070-5 *(metadata verified against PubMed/eutils, "
    "2026-11-07)*",
    "ref 3")

# --- ref 8: title/journal/authors confirmed (CrossRef + publisher page);
#     article number not yet assigned by the publisher as of 2026-11-07
s = sub(s,
    "Nat Commun. 2026;17 (published 2026-10-15). DOI 10.1038/s41467-026-78442-8 "
    "*(article number/page range to fill at typesetting)*",
    "Nat Commun. 2026;17 (published 2026-10-15). DOI 10.1038/s41467-026-78442-8 "
    "*(title/journal/authors verified via CrossRef and the publisher page, 2026-11-07; "
    "volume/article number not yet assigned on the publisher citation line \u2014 "
    "complete at proof stage)*",
    "ref 8")

# --- ref 21: PMID corrected 27504818 -> 27625185 (eutils title match)
s = sub(s,
    "Genet Epidemiol. 2016;40(7):597\u2013608. PMID 27504818 "
    "*(PMID to confirm at typesetting; DOI verified)*. DOI 10.1002/gepi.21998",
    "Genet Epidemiol. 2016;40(7):597\u2013608. PMID 27625185. "
    "DOI 10.1002/gepi.21998 *(PMID corrected and verified via PubMed/eutils, "
    "2026-11-07: 27504818 belonged to an unrelated spectroscopy paper)*",
    "ref 21")

# --- ref 24: vol/issue/pages corrected 133(4):617-624 -> 133(3):569-572 (eutils)
s = sub(s,
    "Epidemiol Infect. 2005;133(4):617\u2013624. PMID 15962564. "
    "DOI 10.1017/S095026880400319X *(mixed-evidence counterpoint for circulating "
    "cytokines; volume/pages to confirm at typesetting \u2014 PubMed shows Epidemiol "
    "Infect, not the earlier-search \"Epidemiol Infect 2005\"; "
    "DOI 10.1017/s095026880400319x)*",
    "Epidemiol Infect. 2005;133(3):569\u2013572. PMID 15962564. "
    "DOI 10.1017/S095026880400319X *(volume/pages corrected and verified via "
    "PubMed/eutils, 2026-11-07; mixed-evidence counterpoint for circulating cytokines)*",
    "ref 24")

io.open(P, "w", encoding="utf-8", newline="").write(s)

# readback
s = io.open(P, encoding="utf-8").read()
for must in ["Li Y, Choi H, Leung K, Jiang F, Graham DY, Leung WK",
             "S2468-1253(23)00070-5", "PMID 27625185", "2005;133(3):569\u2013572",
             "not yet assigned on the publisher citation line"]:
    assert must in s, "readback fail: " + must
for gone in ["Yin J, Fu W", "00050-1", "PMID 27504818", "133(4):617"]:
    assert gone not in s, "old text remains: " + gone
print("readback OK")
