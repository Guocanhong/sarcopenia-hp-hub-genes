# -*- coding: utf-8 -*-
# 29_magnitude_gap_paragraph.py -- Insert the observational-vs-MR magnitude
# discussion paragraph at the end of Discussion 4.2 paragraph 1. Reuses existing
# refs 7 (Davey Smith & Ebrahim 2003) and 33 (Burgess et al. 2016) only -- no new
# entries, no renumbering. Also appends a drafting-notes line.
import io

P = 'D:/GEO-limma,\u808c\u5c11\u75c71/70-MR/paper/Discussion_draft.md'
t = io.open(P, encoding='utf-8').read()

anchor = ('the present genetic design, together with the null reverse-MR line, '
          'addresses both concerns directly.')
assert t.count(anchor) == 1, 'anchor not unique: %d' % t.count(anchor)

para = (

 "\n\nThe magnitude of our estimates is far smaller than the observational "
 "association, and this apparent discrepancy is expected rather than "
 "contradictory, because the two quantities are not commensurable. The NHANES "
 "III estimate is a cross-sectional odds ratio for a binary syndromic "
 "definition \u2014 confirmed sarcopenia combines low mass with weak grip and slow "
 "gait \u2014 contrasting seropositive with seronegative individuals, whereas our "
 "estimates are continuous per-unit effects on the mass component alone: the SD "
 "change in appendicular lean mass per unit increase in genetically predicted "
 "liability to seropositivity. Three further factors make a smaller genetic "
 "estimate the anticipated result. First, the observational contrast remains "
 "exposed to residual confounding by socioeconomic status, smoking and "
 "comorbidity, and to reverse causation; eliminating these biases is precisely "
 "what the genetic design is for (Davey Smith and Ebrahim, 2003) {{CITE 7: "
 "Davey Smith G, Ebrahim S. Int J Epidemiol 2003;32(1):1\u201322, PMID 12689998 \u2014 "
 "MR mitigates confounding and reverse causation}}. Second, a single "
 "serological measurement cannot distinguish current from past infection and "
 "carries within-person measurement variability, which both blurs and can "
 "inflate observed cross-sectional contrasts. Third, two-sample IVW estimates "
 "from instruments explaining a modest fraction of exposure variance are "
 "attenuated toward the null, so our \u03b2 of \u22120.061 SD per log-odds should be "
 "read as a conservative lower bound (Burgess et al., 2016) {{CITE 33: Burgess "
 "S, Davies NM, Thompson SG. Genet Epidemiol 2016;40(7):597\u2013608, DOI "
 "10.1002/gepi.21998 \u2014 weak-instrument attenuation and overlap-bias bounds}}; "
 "on this scale, the difference in genetic liability between seropositive and "
 "seronegative individuals corresponds to a shift in appendicular lean mass "
 "that is a modest fraction of an SD, not the twofold odds of the binary "
 "observational contrast. What triangulation requires is consistency in "
 "direction across designs with different bias structures, not equality of "
 "magnitude \u2014 and on that criterion the observational and genetic lines agree."
)

t = t.replace(anchor, anchor + para)

# drafting-notes bookkeeping
note_anchor = '- **Word count:** \u2248 1,900 words'
assert t.count(note_anchor) == 1
t = t.replace(note_anchor,
 '  12. Observational-vs-MR magnitude gap \u2192 **refs 7, 33 only** (no new entries; paragraph appended to 4.2 after the Wu & Chen sentence on 2026-11-07; arguable points deliberately conservative: no quantitative conversion attempted, "conservative lower bound" from weak-instrument attenuation per Burgess 2016).\n'
 '- **Word count:** \u2248 2,130 words after the 2026-11-07 magnitude-gap paragraph (Medicine has no hard Discussion limit)')

io.open(P, 'w', encoding='utf-8').write(t)

# verify
t2 = io.open(P, encoding='utf-8').read()
for probe in ['far smaller than the observational association',
              'conservative lower bound (Burgess et al., 2016)',
              '(Davey Smith and Ebrahim, 2003)',
              'not equality of magnitude',
              'refs 7, 33 only']:
    print('OK ' if probe in t2 else 'MISSING', probe[:60])
print('word count:', len(t2.split()))
