# -*- coding: utf-8 -*-
# 28_steiger_writeback.py -- Write Steiger backfill results into the manuscript:
# (1) fix 27_steiger_backfill.R relative-path bug; (2) 19_build_submission_md.py:
#     add walking steiger files to the Table 5 dict + update Table 5 note;
# (3) Results 3.1/3.5 + drafting note; (4) Tables_draft.md Table S2: fill antigen
#     Steiger TRUE, add 2 reverse walking rows, update note.
import io, re, glob

def rw(p, t): io.open(p, 'w', encoding='utf-8').write(t)

def patch(path, pairs, label):
    t = io.open(path, encoding='utf-8').read()
    for old, new in pairs:
        if t.count(old) != 1:
            print('%s :: count=%d for %r' % (label, t.count(old), old[:80]))
            return False
        t = t.replace(old, new)
    rw(path, t)
    print(label, 'OK (%d patches)' % len(pairs))
    return True

ok = True

# ---------- 1) fix 27_steiger_backfill.R ----------
ok &= patch(r'D:/GEO-limma,\u808c\u5c11\u75c71/70-MR/code/27_steiger_backfill.R'.encode().decode('unicode_escape')
            if False else 'D:/GEO-limma,\u808c\u5c11\u75c71/70-MR/code/27_steiger_backfill.R',
    [( '  fwrite(as.data.frame(s), out_name)',
       '  fwrite(as.data.frame(s), file.path(OUT, out_name))')],
    '27_steiger_backfill.R')

# ---------- 2) 19_build_submission_md.py ----------
B = 'D:/GEO-limma,\u808c\u5c11\u75c71/70-MR/code/19_build_submission_md.py'
ok &= patch(B, [
    ('for fn in ["line3_steiger_grip_right_to_igg_pos.csv", "line3_steiger_grip_right_to_selfreport.csv",\n'
     '           "line3_steiger_alm_to_igg_pos.csv", "line3_steiger_alm_to_selfreport.csv"]:',
     'for fn in ["line3_steiger_grip_right_to_igg_pos.csv", "line3_steiger_grip_right_to_selfreport.csv",\n'
     '           "line3_steiger_alm_to_igg_pos.csv", "line3_steiger_alm_to_selfreport.csv",\n'
     '           "line3_steiger_walking_to_igg_pos.csv", "line3_steiger_walking_to_selfreport.csv"]:'),
    ('The grip \\u2192 seropositivity signal failed the Steiger directionality check (FALSE); the walking pace \\u2192 seropositivity signal survives Bonferroni correction across the six reverse tests (\\u03b1 = 0.0083) but is unsupported by the weighted-median estimator (P = 0.054), absent for self-reported infection, and lacks a computable directionality test for the walking-pace pairs (Results 3.5).',
     'The grip \\u2192 seropositivity signal failed the Steiger directionality check (FALSE), and so did the walking pace \\u2192 seropositivity signal (Steiger P = 0.348, backfilled 2026-11-07 from the frozen harmonised data); that signal also survives Bonferroni correction across the six reverse tests (\\u03b1 = 0.0083) but is unsupported by the weighted-median estimator (P = 0.054) and absent for self-reported infection (Results 3.5).'),
], '19_build_submission_md.py')

# ---------- 3) Results_draft.md ----------
R = 'D:/GEO-limma,\u808c\u5c11\u75c71/70-MR/paper/Results_draft.md'
ok &= patch(R, [
    # 3.1 sentence: antigen line now also confirmed
    ('and Steiger directionality tests confirmed the hypothesised exposure-to-outcome orientation for the forward infection-status and gene-expression analyses (Table S2).',
     'and Steiger directionality tests confirmed the hypothesised exposure-to-outcome orientation for all forward-line analyses (infection status, antigen spectrum and gene expression; Table S2).'),
    # 3.5 walking pair sentence: add Steiger FALSE
    ('Genetic liability to faster usual walking pace was associated with lower odds of seropositivity (\\u03b2 = \\u22120.300, SE = 0.104, P = 0.0041, 56 instruments; weighted median \\u03b2 = \\u22120.317, P = 0.054, directionally concordant), with no association for self-reported infection (\\u03b2 = \\u22120.003, P = 0.24).',
     'Genetic liability to faster usual walking pace was associated with lower odds of seropositivity (\\u03b2 = \\u22120.300, SE = 0.104, P = 0.0041, 56 instruments; weighted median \\u03b2 = \\u22120.317, P = 0.054, directionally concordant), but the Steiger directionality test did not support the hypothesised direction for this pair (correct direction = FALSE, Steiger P = 0.35). There was no association of walking pace with self-reported infection (\\u03b2 = \\u22120.003, P = 0.24; Steiger direction = TRUE).'),
    # 3.5 summary paragraph
    ('That signal is, however, unsupported by the weighted-median estimator (P = 0.054), absent for self-reported infection (P = 0.24), and not corroborated by a directionality test, which is unavailable for the walking-pace pairs given the small seropositivity outcome sample (Table 5).',
     'That signal is, however, unsupported by the weighted-median estimator (P = 0.054), absent for self-reported infection (P = 0.24), and fails the Steiger directionality check (correct direction = FALSE, Steiger P = 0.35; Table 5).'),
    # drafting note status
    ('3.5 now states this explicitly with WM P = 0.054 and no-Steiger caveat (Table 5 walking rows have no Steiger entry; recomputation pending).',
     '3.5 states this explicitly with WM P = 0.054; Steiger backfilled 2026-11-07 (walking \\u2192 seropositivity FALSE, P = 0.348; walking \\u2192 self-report TRUE) and Table 5 walking rows now carry Steiger entries.'),
], 'Results_draft.md')

# ---------- 4) Tables_draft.md Table S2 ----------
T = 'D:/GEO-limma,\u808c\u5c11\u75c71/70-MR/paper/Tables_draft.md'
t = io.open(T, encoding='utf-8').read()
# fill Steiger TRUE for the 24 antigen rows (5th column: Egger -, Q value, Steiger -)
n = 0
for m in re.finditer(r'\| Antigen \| (\w+) -> (\w+) \| - \| ([\w.e-]+) \| - \| \|', t):
    row_old = m.group(0)
    row_new = '| Antigen | %s -> %s | - | %s | TRUE | Steiger backfilled 2026-11-07 |' % (m.group(1), m.group(2), m.group(3))
    t = t.replace(row_old, row_new); n += 1
print('antigen rows filled:', n)
ok &= (n == 24)
# add reverse walking rows before the MR-PRESSO row
anchor = '| Forward (MR-PRESSO) | Seropositivity -> ALM'
new_rows = ('| Reverse | Usual walking pace -> Anti-helicobacter pylori IgG seropositivity | - | - | FALSE | Steiger P = 0.348, backfilled 2026-11-07 |\n'
            '| Reverse | Usual walking pace -> Non-cancer illness code, self-reported: helicobacter pylori | - | - | TRUE | Steiger backfilled 2026-11-07 |\n')
assert t.count(anchor) == 1
t = t.replace(anchor, new_rows + anchor)
# update the S2 note
old_note = 'Note: \\u2013\\u2013 entries are sensitivity statistics that were not computed or not estimable for the pair (Steiger for the antigen-level and reverse walking-pace pairs is pending recomputation; Egger intercepts require at least three instruments).'
old_note = old_note.replace('\\u2013', '\u2013')
new_note = 'Note: \u2013\u2013 entries are sensitivity statistics that were not computed or not estimable for the pair (Egger intercepts require at least three instruments). Steiger directionality for the antigen-level and reverse walking-pace pairs was backfilled on 2026-11-07 by re-running directionality_test on the frozen harmonised datasets (code/27_steiger_backfill.R; per-pair outputs line1c_steiger_*.csv and line3_steiger_walking_*.csv); all 24 antigen pairs and the walking \u2192 self-report pair support the hypothesised direction, while walking \u2192 seropositivity does not (Steiger P = 0.348).'
assert t.count(old_note) == 1
t = t.replace(old_note, new_note)
rw(T, t)
print('Tables_draft.md OK')

print('ALL OK' if ok else 'PROBLEMS FOUND')
