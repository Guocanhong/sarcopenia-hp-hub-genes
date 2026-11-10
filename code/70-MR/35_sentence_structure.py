# -*- coding: utf-8 -*-
"""35_sentence_structure.py — sentence-structure pass (2026-11-07).

Targets: mechanical parallelism (Results 3.5 "Genetic liability to..." x3),
overloaded long sentences split (3.3/3.4/4.3/Abstract/titre para), duplicate
stats removed from 3.5 summary paragraph, repeated templates varied (4.1, 4.6),
filler removed (Intro "upstream, modifiable", "MR-PRESSO-free", "model here").
No citation tokens, no frozen numbers touched. Idempotent: skips patches whose
old text is absent but new text already present.
"""
import io, os

P = r'D:/GEO-limma,肌少症1/70-MR/paper/'
log = []

def patch(fn, old, new, tag):
    p = os.path.join(P, fn)
    t = io.open(p, encoding='utf-8').read()
    if old in t:
        assert t.count(old) == 1, 'AMBIGUOUS %s [%s]' % (fn, tag)
        t = t.replace(old, new)
        io.open(p, 'w', encoding='utf-8').write(t)
        log.append('%s [%s] APPLIED' % (fn, tag))
    elif new in t:
        log.append('%s [%s] already applied' % (fn, tag))
    else:
        raise AssertionError('NOT FOUND %s [%s] :: %s...' % (fn, tag, old[:80]))

# ---------- Results 3.5: break the threefold "Genetic liability to..." run ----------
R35_OLD = ('Genetic liability to higher right-hand grip strength was nominally associated '
 'with lower odds of seropositivity (IVW random-effects log-OR = −0.124, SE = 0.048, '
 'P = 0.010, 169 instruments; weighted median log-OR = −0.145, P = 0.052), but the Steiger '
 'directionality test could not confirm the hypothesised direction for this pair (snpR² '
 'exposure 0.018 vs outcome 0.019; correct direction = FALSE, Steiger P = 0.74), which is '
 'expected given the small outcome sample (N = 8,735) and the low variance explained. '
 'Genetic liability to ALM showed no association with seropositivity (β = −0.018, P = 0.22, '
 '659 instruments; Steiger direction = TRUE) and neither grip strength nor ALM was '
 'associated with self-reported infection (P = 0.56 and P = 0.15 respectively; 100 and 333 '
 'instruments; both Steiger direction = TRUE). Genetic liability to faster usual walking '
 'pace was associated with lower odds of seropositivity (β = −0.300, SE = 0.104, '
 'P = 0.0041, 56 instruments; weighted median log-OR = −0.317, P = 0.054, of the same '
 'sign), but the Steiger directionality test did not support the hypothesised direction '
 'for this pair (correct direction = FALSE, Steiger P = 0.35). There was no association of '
 'walking pace with self-reported infection (β = −0.003, P = 0.24; Steiger direction = TRUE).')
R35_NEW = ('Genetically predicted higher right-hand grip strength was nominally associated '
 'with lower odds of seropositivity (IVW random-effects log-OR = −0.124, SE = 0.048, '
 'P = 0.010, 169 instruments; weighted median log-OR = −0.145, P = 0.052), but the Steiger '
 'directionality test could not confirm the hypothesised direction for this pair (snpR² '
 'exposure 0.018 vs outcome 0.019; correct direction = FALSE, Steiger P = 0.74), as '
 'expected given the small outcome sample (N = 8,735) and the low variance explained. '
 'ALM liability, in turn, was unassociated with seropositivity (β = −0.018, P = 0.22, '
 '659 instruments; Steiger direction = TRUE), and neither grip strength nor ALM was '
 'associated with self-reported infection (P = 0.56 and P = 0.15, respectively; 100 and '
 '333 instruments; both Steiger direction = TRUE). Genetically predicted faster usual '
 'walking pace was likewise associated with lower odds of seropositivity (β = −0.300, '
 'SE = 0.104, P = 0.0041, 56 instruments; weighted median log-OR = −0.317, P = 0.054, of '
 'the same sign), but the Steiger test again did not support the hypothesised direction '
 '(correct direction = FALSE, Steiger P = 0.35). Walking pace showed no association with '
 'self-reported infection (β = −0.003, P = 0.24; Steiger direction = TRUE).')

# ---------- Results 3.5 summary paragraph: drop stats already given above ----------
R35B_OLD = ('Of the two nominal signals, grip strength → seropositivity (P = 0.010) failed '
 'the Steiger directionality check (correct direction = FALSE), whereas walking pace → '
 'seropositivity (P = 0.0041) does survive Bonferroni correction across the six reverse '
 'tests (α = 0.0083). That signal is, however, unsupported by the weighted-median '
 'estimator (P = 0.054), absent for self-reported infection (P = 0.24), and fails the '
 'Steiger directionality check (correct direction = FALSE, Steiger P = 0.35; Table 5).')
R35B_NEW = ('Of the two nominal signals, grip strength → seropositivity failed the Steiger '
 'directionality check, whereas walking pace → seropositivity does survive Bonferroni '
 'correction across the six reverse tests (α = 0.0083) but lacks support from the '
 'weighted-median estimator, the self-reported outcome and the directionality test '
 '(Table 5).')

patch('Results_draft.md', R35_OLD, R35_NEW, 'R35-parallel')
patch('Results_draft.md', R35B_OLD, R35B_NEW, 'R35-dup-stats')

# ---------- Results 3.4: split 83-word heterogeneity sentence ----------
patch('Results_draft.md',
 'which widens the random-effects confidence interval; the IVW random-effects estimate is therefore',
 'which widens the random-effects confidence interval. The IVW random-effects estimate is therefore',
 'R34-split-het')

# ---------- Results 3.4: split C1QA/CTNNBIP1 compound sentence ----------
patch('Results_draft.md',
 'or low grip strength (IVW random-effects log-OR = −0.022, P = 0.12), and CTNNBIP1 (3 harmonised instruments) was null',
 'or low grip strength (IVW random-effects log-OR = −0.022, P = 0.12). CTNNBIP1 (3 harmonised instruments) was null',
 'R34-split-genes')

# ---------- Results 3.3: split GroEL overdispersion sentence; drop filler "here" ----------
patch('Results_draft.md',
 'the multiplicative random-effects model here shrinks rather than inflates the standard '
 'error relative to the fixed-effect model (random-effects SE 0.0039 versus fixed-effects '
 'SE 0.0072), the same overdispersion artefact documented for CTNNBIP1 in 3.4.',
 'the multiplicative random-effects model shrinks rather than inflates the standard error '
 'relative to the fixed-effect model (random-effects SE 0.0039 versus fixed-effects SE '
 '0.0072). This is the same overdispersion artefact documented for CTNNBIP1 in 3.4.',
 'R33-split-overdisp')

# ---------- Results 3.3: split GroEL decision-rule sentence ----------
patch('Results_draft.md',
 'the GroEL→ALM estimate is reported as suggestive — directionally concordant with the '
 'seropositivity signal, but not presented as a Bonferroni-confirmed finding.',
 'the GroEL→ALM estimate is reported as suggestive; it is directionally concordant with '
 'the seropositivity signal but not Bonferroni-confirmed.',
 'R33-split-suggestive')

# ---------- Results 3.2: split titre-age-structure sentence ----------
patch('Results_draft.md',
 'rather than by the magnitude of the antibody response, although the age-structure '
 'difference between the two exposure GWAS (the titre GWAS measured in ALSPAC children '
 'aged 7 years, seropositivity in UK Biobank adults; Methods 2.2.1) prevents a clean '
 'dissociation of the two explanations.',
 'rather than by the magnitude of the antibody response. The age-structure difference '
 'between the two exposure GWAS (the titre GWAS measured in ALSPAC children aged 7 years, '
 'seropositivity in UK Biobank adults; Methods 2.2.1), however, prevents a clean '
 'dissociation of the two explanations.',
 'R32-split-titre')

# ---------- Abstract: split the relative-chain sentence ----------
patch('Abstract_draft.md',
 'anti-GroEL antibody levels showed a nominally significant random-effects association '
 'with lower ALM (β = −0.021 SD, P = 3.1 × 10⁻⁸), which was not supported by the '
 'fixed-effect or weighted-median estimators and is therefore reported as suggestive '
 'under the uniform decision rule.',
 'anti-GroEL antibody levels showed a nominally significant random-effects association '
 'with lower ALM (β = −0.021 SD, P = 3.1 × 10⁻⁸). The fixed-effect and weighted-median '
 'estimators did not support this signal, which is therefore reported as suggestive under '
 'the uniform decision rule.',
 'ABS-split-groel')

# ---------- Discussion 4.1: vary the "was causally associated with lower ALM" template ----------
patch('Discussion_draft.md',
 'the master inhibitor and feedback target of the NF-κB pathway, was causally associated with lower ALM (β = −0.044 SD units, P = 4.9 × 10⁻⁴; weighted median P = 2.3 × 10⁻⁸)',
 'the master inhibitor and feedback target of the NF-κB pathway, lowered ALM (β = −0.044 SD units, P = 4.9 × 10⁻⁴; weighted median P = 2.3 × 10⁻⁸)',
 'D41-vary-template')

# ---------- Discussion 4.3: split GTEx tail sentence ----------
patch('Discussion_draft.md',
 'so a muscle-resolved cis-instrument set of usable strength is not currently extractable '
 'and the tissue-resolved extension remains future work conditional on larger muscle eQTL '
 'resources.',
 'so a muscle-resolved cis-instrument set of usable strength is not currently extractable. '
 'The tissue-resolved extension therefore remains conditional on larger muscle eQTL '
 'resources.',
 'D43-split-gtex')

# ---------- Discussion 4.6: "causal effect ... causal support" repetition ----------
patch('Discussion_draft.md',
 'with causal support for an NF-κB–centric inflammatory programme (NFKBIA)',
 'with convergent evidence for an NF-κB–centric inflammatory programme (NFKBIA)',
 'D46-causal-repeat')

# ---------- Introduction: filler "upstream, modifiable" repetition ----------
patch('Introduction_draft.md',
 'No MR study has yet traced sarcopenia back to such an upstream, modifiable exposure.',
 'No MR study has yet traced sarcopenia back to such an exposure.',
 'INT-trim-filler')

# ---------- Introduction: fix "MR-PRESSO-free" coined compound ----------
patch('Introduction_draft.md',
 'with sensitivity analyses (MR-Egger, weighted median, weighted mode, Cochran\'s Q, '
 'MR-PRESSO-free leave-one-out and Steiger directionality) applied throughout.',
 'with a standard sensitivity battery (MR-Egger, weighted median, weighted mode, '
 'Cochran\'s Q, leave-one-out and Steiger directionality) applied throughout.',
 'INT-mrpresso-free')

print('\n'.join(log))
print('%d patches processed' % len(log))
