# -*- coding: utf-8 -*-
"""32_fix_edit_review.py — 执行 Academic_edit_review_2026-11-07.md 的 E1–E31 修复。
顺序 patch + 每条断言 count==1（E26 为文件内 2 处）。不改动冻结数值；不增删引用。
E10/E16 按报告结论保留不改；E25 由 E3 覆盖（randomisation 短语整体删除）。
"""
import io

P = 'D:/GEO-limma,肌少症1/70-MR/paper/'
log = []

def patch(fn, pairs):
    p = P + fn
    t = io.open(p, encoding='utf-8').read()
    for old, new, tag in pairs:
        if old in t:
            c = t.count(old)
            assert c == 1, 'FAIL %s [%s] count=%d :: %s...' % (fn, tag, c, old[:70])
            t = t.replace(old, new)
            log.append('%s [%s] OK' % (fn, tag))
        else:
            # already applied in a previous partial run; verify new state exists
            assert (new in t) or (new == '' and old not in t), \
                'FAIL %s [%s] neither old nor new present :: %s...' % (fn, tag, old[:70])
            log.append('%s [%s] ALREADY-APPLIED (skipped)' % (fn, tag))
    io.open(p, 'w', encoding='utf-8').write(t)

# ---------- Abstract ----------
patch('Abstract_draft.md', [
    ('bilateral grip strength, usual walking pace and low grip strength',
     'grip strength, usual walking pace and low grip strength', 'E26'),
    ('with concordant alternative estimators and clean sensitivity profiles',
     'with concordant alternative estimators and supportive sensitivity analyses', 'E4'),
])

# ---------- Introduction ----------
patch('Introduction_draft.md', [
    ('Sarcopenia — the progressive loss of skeletal muscle mass and strength — affects',
     'Sarcopenia, the progressive loss of skeletal muscle mass and strength, affects', 'E13a'),
    ('This gap underscores the need to identify modifiable upstream risk factors for '
     'sarcopenia whose causal relevance can be established with robust evidence.',
     'Identifying modifiable upstream factors is therefore a priority.', 'E2'),
    ('was associated with a more than twofold risk of confirmed sarcopenia '
     '(OR = 2.05, 95% CI 1.70–2.48)',
     'was associated with approximately twofold odds of confirmed sarcopenia '
     '(OR 2.05, 95% CI 1.70–2.48)', 'E30'),
    ('and — intriguingly — antibodies', 'and, intriguingly, antibodies', 'E13b'),
    ('bilateral grip strength, usual walking pace and low grip strength',
     'grip strength, usual walking pace and low grip strength', 'E26'),
    ('formalised this concern and demonstrated that reproducible causal inference '
     'requires targeted hypothesis-driven designs with systematic validation across '
     'complementary data sources',
     'formalised this concern and argued that reproducible causal inference benefits '
     'from targeted, hypothesis-driven designs with systematic validation across '
     'complementary data sources', 'E24'),
])

# ---------- Methods ----------
patch('Methods_draft.md', [
    ('with genetic variants as instrumental variables (IVs) serving as natural '
     'randomisation.', 'with genetic variants as instrumental variables (IVs).', 'E3+E25'),
    ('All participating cohorts had received appropriate ethical approval and written '
     'informed consent; no new human data were collected for the present study. ',
     '', 'E5'),
])

# ---------- Results ----------
patch('Results_draft.md', [
    # E21
    ('Together, the divergence between the seropositivity and IgG-titre instruments '
     'indicates that the causal signal is driven by infection status (seropositivity) '
     'rather than by the magnitude of the antibody response — noting additionally that '
     'the titre GWAS was measured in ALSPAC children aged 7 years whereas seropositivity '
     'was measured in UK Biobank adults (Methods 2.2.1), so the two exposures also '
     'differ in age structure.',
     'Together, the divergence between the seropositivity and IgG-titre instruments is '
     'consistent with a signal driven by infection status (seropositivity) rather than '
     'by the magnitude of the antibody response, although the age-structure difference '
     'between the two exposure GWAS (the titre GWAS measured in ALSPAC children aged 7 '
     'years, seropositivity in UK Biobank adults; Methods 2.2.1) prevents a clean '
     'dissociation of the two explanations.', 'E21'),
    # E29
    ('Bonferroni threshold α = 0.05/24 = 2.1 × 10⁻³',
     'Bonferroni threshold α = 0.05/24 = 0.0021', 'E29'),
    # E31 (3.3 first occurrence)
    ('the anti-GroEL (chaperonin-60/HSP60) immune response',
     'the anti-GroEL (*H. pylori* chaperonin-60, the functional homologue of human '
     'HSP60) immune response', 'E31'),
    # E8
    ('The MR-Egger estimate was likewise directionally concordant but imprecise',
     'The MR-Egger estimate agreed in direction but was imprecise', 'E8'),
    ('provides suggestive, directionally concordant evidence',
     'provides suggestive evidence', 'E8'),
    # E13f
    ('and — applying the same decision rule as for the OMP antigen signal — the estimate',
     'and, applying the same decision rule as for the OMP antigen signal, the estimate',
     'E13f'),
    # E19
    ('showed no association with grip strength (β = 0.007, P = 0.086)',
     'showed no association with right grip strength (β = 0.007, P = 0.086)', 'E19'),
    # E28
    ('the random-effects IVW estimate is therefore the appropriately conservative '
     'primary estimate',
     'the IVW random-effects estimate is therefore the appropriately conservative '
     'primary estimate', 'E28'),
    # E8
    ('The MR-Egger estimate was directionally concordant but non-significant',
     'The MR-Egger estimate agreed in direction but was non-significant', 'E8'),
    # E15
    ('Thus, among the four candidate hub genes, only NFKBIA — the inhibitor of the '
     'NF-κB pathway — exhibited a causal signal, and it pointed to the same outcome '
     '(ALM) as the infection exposures; the single nominal CTNNBIP1–low grip estimate '
     'did not survive fixed-effects or weighted-median checks and is reported as '
     'suggestive only (3.4).',
     'Thus, among the four candidate hub genes, only NFKBIA exhibited a causal signal, '
     'directed at the same outcome (ALM) as the infection exposures.', 'E15'),
    # E27
    ('weighted median β = −0.145, P = 0.052',
     'weighted median log-OR = −0.145, P = 0.052', 'E27'),
    # E6
    ('Overall, the reverse line found no consistent causal effect of muscle phenotypes '
     'on *H. pylori* infection across the two infection definitions. Of the two nominal '
     'signals,',
     'Of the two nominal signals,', 'E6'),
    # E14
    ('Genetic liability to ALM was not associated with seropositivity',
     'Genetic liability to ALM showed no association with seropositivity', 'E14'),
    # E27 + E8
    ('weighted median β = −0.317, P = 0.054, directionally concordant',
     'weighted median log-OR = −0.317, P = 0.054, of the same sign', 'E27+E8'),
    # E11 + E8
    ('reported as suggestive under the uniform decision rule and directionally '
     'concordant with the seropositivity signal',
     'reported as suggestive and pointing in the same direction as the seropositivity '
     'signal', 'E11+E8'),
    # E13g
    ('converge on a single outcome — appendicular lean mass:',
     'converge on a single outcome, appendicular lean mass:', 'E13g'),
])

# ---------- Discussion ----------
patch('Discussion_draft.md', [
    # E1 + E17
    ('In this two-sample Mendelian randomization study, we leveraged genetic variation '
     'as a natural experiment to test the causal relationship between *Helicobacter '
     'pylori* infection and sarcopenia-related traits. Three independent lines of '
     'evidence converged on a single outcome — appendicular lean mass (ALM). First,',
     'Two independent lines of evidence converged on a single outcome, appendicular '
     'lean mass (ALM), and a third, suggestive line pointed in the same direction. '
     'First,', 'E1+E17'),
    # E13c
    ('NFKBIA — the master inhibitor and feedback target of the NF-κB pathway — was '
     'causally associated',
     'NFKBIA, the master inhibitor and feedback target of the NF-κB pathway, was '
     'causally associated', 'E13c'),
    # E31 (4.1 second occurrence)
    ('the anti-GroEL (chaperonin-60/HSP60) response',
     'the anti-GroEL response', 'E31'),
    # E18 + E8
    ('No comparable effects were observed for grip strength, usual walking pace, '
     'antibody titre per se, the other five antigens, or the other three candidate hub '
     'genes, indicating that the robust signals are outcome-specific (ALM) and '
     'pathway-specific (NFKBIA), with the antigen-level GroEL estimate providing '
     'suggestive directional concordance.',
     'No comparable estimates were observed for grip strength, usual walking pace, '
     'antibody titre per se, the other five antigens, or the other three candidate hub '
     'genes; these nulls are not proof of absence (Section 4.4), but the robust signals '
     'were confined to ALM and NFKBIA, with the antigen-level GroEL estimate providing '
     'suggestive directional support.', 'E18'),
    # E7
    (' — but such associations are vulnerable to confounding by socioeconomic status, '
     'comorbidity and reverse causation (muscle-wasting states altering infection '
     'susceptibility); the present genetic design, together with the null reverse-MR '
     'line, addresses both concerns directly.',
     '.', 'E7'),
    # E9
    ('concluded that direct causal evidence for this link "has not been fully '
     'explored" (Liu and Peng, 2025)',
     'concluded that direct causal evidence for this link remained lacking (Liu and '
     'Peng, 2025)', 'E9'),
    # E13d
    ('binary syndromic definition — confirmed sarcopenia combines low mass with weak '
     'grip and slow gait — contrasting',
     'binary syndromic definition (confirmed sarcopenia combines low mass with weak '
     'grip and slow gait) contrasting', 'E13d'),
    # E31 (4.2 direction-of-homology fix)
    ('GroEL is the *H. pylori* homologue of chaperonin-60 (HSP60) and shares '
     'substantial sequence homology with human mitochondrial HSP60, making it',
     'GroEL is the chaperonin-60 of *H. pylori* and shares substantial sequence '
     'homology with human mitochondrial HSP60, making it', 'E31b'),
    # E22
    ('indicating that the systemic acute-phase component of the inflammatory response, '
     'like the antigen-specific component dissected here, carries a deleterious muscle '
     'signal',
     'indicating that the systemic acute-phase component of the inflammatory response '
     'carries a deleterious muscle signal, in parallel with the suggestive '
     'antigen-specific estimate dissected here', 'E22'),
    # E13e
    ('The direction of our estimate — greater predicted NFKBIA expression, lower ALM — '
     'is therefore consistent',
     'The direction of our estimate (greater predicted NFKBIA expression, lower ALM) '
     'is therefore consistent', 'E13e'),
    # E28
    ('the random-effects IVW estimate is the appropriately conservative primary '
     'estimate',
     'the IVW random-effects estimate is the appropriately conservative primary '
     'estimate', 'E28'),
    # E20
    ('A related consideration is statistical power: with 11 instruments for '
     'seropositivity and modest per-instrument variance explained, our design had power '
     'to detect effects of the observed magnitude on ALM, but effects on grip strength '
     'would need to be of comparable size to be detectable, and genuinely weaker '
     'functional effects could have been missed.',
     'A related consideration is statistical power: with 11 instruments for '
     'seropositivity and modest per-instrument variance explained, the precision of the '
     'primary estimate (SE = 0.017 SD units) implies 80% power at α = 0.05 for effects '
     'of approximately 0.05 SD or larger; the observed ALM effect (β = −0.061) was '
     'therefore within the detectable range, whereas genuinely smaller effects on grip '
     'strength or walking pace could have been missed.', 'E20'),
    # E12
    ('Ninth, the ALM outcome GWAS derives lean mass',
     'Eighth, the ALM outcome GWAS derives lean mass', 'E12'),
    # E13h
    ('distinct bias structures — agreement across lines with different assumptions is '
     'stronger evidence than any single estimate;',
     'distinct bias structures, where agreement across lines with different assumptions '
     'is stronger evidence than any single estimate;', 'E13h'),
    # E8 (4.6, keep Abstract's identical phrase)
    ('suggestive, directionally concordant evidence implicating the anti-GroEL immune '
     'response',
     'suggestive evidence implicating the anti-GroEL immune response', 'E8'),
    # E23
    ('and NFKBIA-centred NF-κB signalling is a candidate druggable node',
     'and NFKBIA-centred NF-κB signalling offers a biologically tractable target for '
     'pharmacological modulation', 'E23'),
])

print('\n'.join(log))
print('ALL %d PATCHES APPLIED' % len(log))
