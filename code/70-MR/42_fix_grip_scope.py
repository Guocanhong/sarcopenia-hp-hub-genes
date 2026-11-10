# -*- coding: utf-8 -*-
# 42_fix_grip_scope.py -- disclose that the antigen-spectrum line used only the
# right-hand grip GWAS (left-hand grip not tested in this line). Idempotent.
import io, os

PAPER = r"D:/GEO-limma,肌少症1/70-MR/paper"
OLD = ("For the antigen-spectrum line (6 exposures \u00d7 4 outcomes = 24 tests), "
       "Bonferroni correction was applied (\u03b1 = 0.05/24 = 0.0021)")
NEW = ("For the antigen-spectrum line, the six antigen exposures were tested against "
       "four of the five outcome GWASs described above \u2014 the right-hand grip GWAS "
       "served as the bilateral-strength outcome and the left-hand grip GWAS was not "
       "included in this line \u2014 giving 24 tests; Bonferroni correction was applied "
       "(\u03b1 = 0.05/24 = 0.0021)")

def patch(path, old, new):
    s = io.open(path, encoding="utf-8").read()
    if new in s:
        print("SKIP (already applied):", os.path.basename(path))
        return
    assert old in s, "anchor not found in " + path
    assert s.count(old) == 1, "anchor not unique in " + path
    io.open(path, "w", encoding="utf-8", newline="").write(s.replace(old, new, 1))
    print("OK:", os.path.basename(path))

patch(os.path.join(PAPER, "Methods_draft.md"), OLD, NEW)

# readback
s = io.open(os.path.join(PAPER, "Methods_draft.md"), encoding="utf-8").read()
assert "the left-hand grip GWAS was not included in this line" in s
assert "6 exposures \u00d7 4 outcomes = 24 tests), Bonferroni" not in s
print("readback OK")
