# -*- coding: utf-8 -*-
# 44_american_spelling.py -- Medicine (Baltimore) follows the AMA Manual of
# Style, which prescribes American spelling. Convert British -> American in
# the manuscript body files and the build script. Reference titles are left
# untouched (verified: zero British spellings in References_draft.md).
# "Age Ageing" (journal name) is protected: only lowercase "ageing" converts.
# Idempotent.
import io, glob, re

FILES = [
    r"D:/GEO-limma,肌少症1/70-MR/paper/Abstract_draft.md",
    r"D:/GEO-limma,肌少症1/70-MR/paper/Introduction_draft.md",
    r"D:/GEO-limma,肌少症1/70-MR/paper/Methods_draft.md",
    r"D:/GEO-limma,肌少症1/70-MR/paper/Results_draft.md",
    r"D:/GEO-limma,肌少症1/70-MR/paper/Discussion_draft.md",
    r"D:/GEO-limma,肌少症1/70-MR/paper/Tables_draft.md",
    r"D:/GEO-limma,肌少症1/70-MR/paper/Supplementary_Methods_S1_draft.md",
    r"D:/GEO-limma,肌少症1/70-MR/code/19_build_submission_md.py",
]

# (pattern with word boundary, replacement; {0,1} preserves case via lambda)
MAP = [
    (r"titre(s?)", r"titer\1"),
    (r"harmonised", r"harmonized"),
    (r"harmonisation", r"harmonization"),
    (r"standardised", r"standardized"),
    (r"analysed", r"analyzed"),
    (r"residualised", r"residualized"),
    (r"residualisation", r"residualization"),
    (r"characterised", r"characterized"),
    (r"characterisation", r"characterization"),
    (r"organised", r"organized"),
    (r"summarised", r"summarized"),
    (r"recognised", r"recognized"),
    (r"emphasise", r"emphasize"),
    (r"maximise", r"maximize"),
    (r"hospitalisation", r"hospitalization"),
    (r"signalling", r"signaling"),
    (r"programmes?", lambda m: "programs" if m.group(0) == "programmes" else "program"),
    (r"ageing", r"aging"),          # lowercase only; "Age Ageing" protected
    (r"centre(s?)", r"center\1"),
    (r"fibres", r"fibers"),
]

total = 0
for path in FILES:
    s = io.open(path, encoding="utf-8").read()
    n0 = s
    for pat, rep in MAP:
        # NOTE: repl must be a string (re.sub processes backrefs) or a callable
        # that already returns the FINAL literal (callables never get backref
        # expansion -- the Nov-7 '\1' pollution bug came from wrapping a
        # backref string in a lambda).
        s = re.sub(r"\b" + pat + r"\b", rep, s)
    # case-preserving pass for capitalized forms (Harmonisation etc.)
    for pat, rep in MAP:
        cap_pat = pat[0].upper() + pat[1:]
        cap_rep = rep if callable(rep) else rep[0].upper() + rep[1:]
        if isinstance(rep, str):
            s = re.sub(r"\b" + cap_pat + r"\b", cap_rep, s)
    if s != n0:
        io.open(path, "w", encoding="utf-8", newline="").write(s)
        diff = sum(1 for a, b in [(n0, s)][0:1])
        print("CONVERTED:", path.split("/")[-1] if "/" in path else path.split("\\")[-1])
        total += 1
    else:
        print("no change :", path.split("/")[-1] if "/" in path else path.split("\\")[-1])

# restore protected journal name (capitalised pass may have touched it)
for path in FILES:
    t = io.open(path, encoding='utf-8').read()
    if 'Age Aging' in t:
        io.open(path, 'w', encoding='utf-8', newline='').write(t.replace('Age Aging', 'Age Ageing'))
        print('restored Age Ageing in', path)

# readback: no British spellings remain in body (Age Ageing allowed)
BAD = re.compile(r"\b(titre|harmonised|harmonisation|standardised|analysed|residualised|"
                 r"residualisation|characterised|characterisation|organised|summarised|"
                 r"recognised|emphasise|maximise|hospitalisation|signalling|programme|"
                 r"programmes|ageing|centre|centres|fibres)\b")
for path in FILES:
    s = io.open(path, encoding="utf-8").read()
    hits = BAD.findall(s)
    assert not hits, "leftover: %s in %s" % (hits, path)
print("readback OK; converted files:", total)
