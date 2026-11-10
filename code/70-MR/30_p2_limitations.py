# -*- coding: utf-8 -*-
# 30_p2_limitations.py -- P2 final item: expand Discussion 4.5 limitations with
# (a) BIA-derived ALM outcome limitation (new "Ninth"), (b) candidate-selection
# degrees of freedom (expanded "Finally"), and (c) fix the 4.6 "study., so" glitch.
# Reuses existing ref 16 (Pei 2020) only -- no new entries, no renumbering.
import io

P = 'D:/GEO-limma,\u808c\u5c11\u75c71/70-MR/paper/Discussion_draft.md'
t = io.open(P, encoding='utf-8').read()

# ---------- (a) insert BIA limitation before "Finally" ----------
anchor_a = ('Eighth, self-reported infection could not be analysed as a forward exposure '
            '(fewer than three independent variants), leaving a minor scope gap that is '
            'transparently reported. Finally,')
assert t.count(anchor_a) == 1, 'anchor A: %d' % t.count(anchor_a)

bia = (
 "Ninth, the ALM outcome GWAS derives lean mass from bioelectrical impedance analysis "
 "(BIA) rather than dual-energy X-ray absorptiometry (Pei et al., 2020) {{CITE 16: Pei "
 "2020 (PMID 33097823) \u2014 BIA-derived ALM GWAS, UK Biobank N = 450,243}}. BIA estimates "
 "lean mass from whole-body electrical conductivity, is sensitive to hydration status "
 "and fluid distribution, and is less accurate than DXA in older adults and at higher "
 "adiposity; because chronic infection and systemic inflammation can in principle "
 "perturb fluid balance, a component of the estimated effect could conceivably reflect "
 "inflammation-associated hydration change rather than true lean-mass loss. The "
 "estimate should therefore be read as applying to ALM as measured by BIA; the "
 "compensating advantage is the sample size that BIA enables (N = 450,243, far "
 "exceeding current DXA-based lean-mass GWAS), and replication against DXA-derived "
 "lean-mass GWAS as these resources scale is a priority. Finally,")
t = t.replace(anchor_a, bia)
assert t.count('is a priority. Finally, the hub-gene candidates') == 1, 'joint failed'

# ---------- (b) expand the "Finally" sentence ----------
old_finally = ('Finally, the hub-gene candidates derive from our integrative transcriptomic '
               'screening analysis (Methods 2.1), whose full pipeline is documented in '
               'Supplementary Methods S1; other inflammatory genes outside this candidate '
               'set were not examined.')
assert t.count(old_finally) == 1, 'anchor B: %d' % t.count(old_finally)
new_finally = (
 "Finally, the hub-gene candidates derive from our integrative transcriptomic screening "
 "analysis (Methods 2.1), whose full pipeline is documented in Supplementary Methods S1. "
 "This nomination layer involves researcher degrees of freedom that are not covered by "
 "the MR multiple-testing correction: DEGs were called at nominal P < 0.05 without a "
 "fold-change threshold, WGCNA modules depend on soft-thresholding and cut-height "
 "parameters, and the four-gene set emerged from requiring concordance between two of "
 "three machine-learning feature rankings and then intersecting the two disease lines. "
 "The gene-expression Bonferroni family therefore protects only within the nominated "
 "set; inflammatory genes outside it were not tested, and the nomination step itself "
 "was not pre-registered or independently replicated. The screening analysis did include "
 "external-cohort validation of the nominated genes' discriminatory performance and "
 "single-cell and protein-level localisation (Supplementary Methods S1), and NFKBIA "
 "survived its pre-specified Bonferroni family while the other three genes were honestly "
 "null. Nevertheless, independent replication of the nomination step \u2014 for example, "
 "reproducing the infection-response and muscle-depletion hub genes in an external "
 "transcriptomic cohort \u2014 would materially strengthen the gene-expression line, and "
 "the results at this layer should be read as hypothesis-generating rather than as an "
 "agnostic discovery screen.")
t = t.replace(old_finally, new_finally)

# ---------- (c) fix 4.6 punctuation glitch ----------
old_c = 'is the contribution of the present study., so extending this association map to muscle mass is a natural next step rather than an outlier claim.'
assert t.count(old_c) == 1, 'anchor C: %d' % t.count(old_c)
t = t.replace(old_c, 'is the contribution of the present study; extending this association map to muscle mass is a natural next step rather than an outlier claim.')

# ---------- drafting notes ----------
note_anchor = '- **Word count:** \u2248 2,130 words'
assert t.count(note_anchor) == 1
t = t.replace(note_anchor,
 '  13. BIA + nomination-freedom limitations \u2192 **ref 16 only** (Pei, already cited in Methods; inserted as Ninth + expanded Finally on 2026-11-07; 4.6 "study., so" glitch fixed).\n'
 '- **Word count:** \u2248 2,450 words')

io.open(P, 'w', encoding='utf-8').write(t)

# verify
t2 = io.open(P, encoding='utf-8').read()
for probe in ['Ninth, the ALM outcome GWAS derives lean mass from bioelectrical impedance',
              'far exceeding current DXA-based lean-mass GWAS',
              'researcher degrees of freedom that are not covered by',
              'hypothesis-generating rather than as an agnostic discovery screen',
              'the contribution of the present study; extending this association map',
              'ref 16 only']:
    print('OK ' if probe in t2 else 'MISSING', probe[:60])
