# -*- coding: utf-8 -*-
"""Remove empty rhetoric per 2026-11-07 user directive:
template background intros, unsupported subjective evaluations,
repeated viewpoints, no-information summary sentences, over-decorated
adjectives. No citations touched, no frozen numbers touched.
Idempotent: patches already applied are skipped."""
import io

P = 'D:/GEO-limma,肌少症1/70-MR/paper/'
log = []

def patch(fn, old, new, tag):
    p = P + fn
    t = io.open(p, encoding='utf-8').read()
    if new in t and old not in t:
        log.append('%s [%s] SKIP (already applied)' % (fn, tag))
        return
    c = t.count(old)
    assert c == 1, 'FAIL %s [%s] count=%d :: %s...' % (fn, tag, c, old[:80])
    t = t.replace(old, new)
    io.open(p, 'w', encoding='utf-8').write(t)
    log.append('%s [%s] OK' % (fn, tag))

# ============ Introduction ============

# I1: template opener + no-info closer -> pose the research question
patch('Introduction_draft.md',
 'Despite its clinical importance, the mechanisms driving age-related muscle loss remain incompletely elucidated: chronic low-grade inflammation, mitochondrial dysfunction and oxidative stress are repeatedly implicated, yet approved pharmacological treatments are lacking {{CITE: Jang SY, Choi KM. Impact of adipose tissue and lipids on skeletal muscle in sarcopenia. J Cachexia Sarcopenia Muscle 2025;16(4):e70000, PMID 40641114}}. Identifying modifiable upstream factors is therefore a priority.',
 'Chronic low-grade inflammation, mitochondrial dysfunction and oxidative stress are repeatedly implicated in age-related muscle loss, yet approved pharmacological treatments are lacking {{CITE: Jang SY, Choi KM. Impact of adipose tissue and lipids on skeletal muscle in sarcopenia. J Cachexia Sarcopenia Muscle 2025;16(4):e70000, PMID 40641114}}; whether a modifiable upstream driver of these processes contributes to muscle loss is unknown.',
 'I1')

# I2: template transition sentence deleted (next sentence carries the evidence)
patch('Introduction_draft.md',
 'Although traditionally regarded as a gastric pathogen, chronic *H. pylori* infection is increasingly recognised for its extra-gastric manifestations. A growing observational literature',
 'A growing observational literature',
 'I2')

# I3: over-decorated adverb
patch('Introduction_draft.md',
 'and, intriguingly, antibodies against bacterial heat-shock proteins',
 'and antibodies against bacterial heat-shock proteins',
 'I3')

# I4a: vague concern -> concrete finding attributed to the two cited studies
patch('Introduction_draft.md',
 'Yet the recent proliferation of "broad-panel" MR screens — testing hundreds of circulating proteins or metabolites against dozens of outcomes — has been accompanied by growing concern that many nominally significant signals do not survive rigorous scrutiny. Two 2026 methodological studies in *Nature Communications* formalised this concern and argued that reproducible causal inference benefits from targeted, hypothesis-driven designs with systematic validation across complementary data sources (Zhang et al., 2026; Machnik et al., 2026)',
 'Yet broad-panel MR screens — testing hundreds of circulating proteins or metabolites against dozens of outcomes — often produce nominally significant signals that fail rigorous scrutiny: two 2026 methodological studies in *Nature Communications* showed this directly and argued for targeted, hypothesis-driven designs with systematic validation across complementary data sources (Zhang et al., 2026; Machnik et al., 2026)',
 'I4a')

# I4b: unsupported paradigm claim -> grounded recommendation
patch('Introduction_draft.md',
 '}}. Against this background, single-exposure MR studies anchored in specific biological hypotheses, triangulated across independent evidence lines, represent the more reliable paradigm.',
 '}}. Single-exposure designs anchored in specific biological hypotheses and triangulated across independent evidence lines follow this recommendation.',
 'I4b')

# ============ Results ============

# R1: emphasis adverb with no content
patch('Results_draft.md',
 'Notably, this exposure exhibited pronounced between-instrument heterogeneity',
 'This exposure exhibited pronounced between-instrument heterogeneity',
 'R1')

# R2: bold topic sentence immediately repeated verbatim by the body sentence
patch('Results_draft.md',
 '**Seropositivity is causally associated with lower appendicular lean mass.** Genetic liability to anti-*H. pylori* IgG seropositivity was causally associated with lower appendicular lean mass (IVW random-effects β = −0.061 SD units per log-odds of seropositivity',
 '**Seropositivity is causally associated with lower appendicular lean mass.** Genetic liability to seropositivity lowered ALM (IVW random-effects β = −0.061 SD units per log-odds of seropositivity',
 'R2')

# R3: same verbatim repetition in 3.4
patch('Results_draft.md',
 'Greater genetically predicted expression of NFKBIA (NF-κB inhibitor alpha) in whole blood was causally associated with lower ALM (IVW random-effects β',
 'Greater genetically predicted whole-blood expression of NFKBIA (NF-κB inhibitor alpha) lowered ALM (IVW random-effects β',
 'R3')

# R4: "Taken together" summary repeating the two preceding paragraphs; keep
#     the chaperonin-60 definition and the provisional caveat in one sentence
patch('Results_draft.md',
 'Taken together, the antigen-spectrum analysis provides suggestive evidence that the forward-line infection signal may be carried in part by the anti-GroEL (*H. pylori* chaperonin-60, the functional homologue of human HSP60) immune response: the general seropositivity signal and the GroEL-specific estimate converge on ALM in sign and outcome, while antibody responses to the other five antigens do not reproduce it. With only four harmonised GroEL instruments and no estimator corroboration beyond the random-effects model, this antigen-level localisation remains provisional (Discussion 4.3).',
 'The GroEL estimate (*H. pylori* chaperonin-60, the functional homologue of human HSP60) rests on four harmonised instruments and no estimator corroboration beyond the random-effects model, so antigen-level localisation of the infection signal remains provisional (Discussion 4.3).',
 'R4')

# R5: no-information section-closing sentence (3.6 carries the convergence claim)
patch('Results_draft.md',
 'DTL could not be tested (fewer than three instruments even at P < 1 × 10⁻⁵; Methods 2.3).\n\nThus, among the four candidate hub genes, only NFKBIA exhibited a causal signal, directed at the same outcome (ALM) as the infection exposures.',
 'DTL could not be tested (fewer than three instruments even at P < 1 × 10⁻⁵; Methods 2.3).',
 'R5')

# ============ Discussion ============

# D1: tail of 4.1 repeats the First/Second/Third enumeration above it
patch('Discussion_draft.md',
 'these nulls are not proof of absence (Section 4.4), but the robust signals were confined to ALM and NFKBIA, with the antigen-level GroEL estimate providing suggestive directional support.',
 'these nulls are not proof of absence (Section 4.4).',
 'D1')

# D2: unsupported self-evaluation ("complement, and move beyond") + decoration
patch('Discussion_draft.md',
 'Our findings complement, and move beyond, a rapidly expanding MR literature on sarcopenia-related traits that now spans lifestyle factors and vitamin D (Zhou et al., 2026)',
 'Recent MR studies of sarcopenia-related traits cover lifestyle factors and vitamin D (Zhou et al., 2026)',
 'D2')

# D3: generalisation from a single example ("broad panels yield mixed and often
#     counterintuitive signals") -> the concrete counterintuitive result itself
patch('Discussion_draft.md',
 'screened 91 circulating inflammatory factors and identified a protective effect of TNF-β on muscle phenotypes, illustrating that broad panels of downstream cytokines yield mixed and often counterintuitive signals (Mei et al., 2026)',
 'screened 91 circulating inflammatory factors and reported a protective effect of TNF-β, a counterintuitive result for a classical pro-inflammatory cytokine (Mei et al., 2026)',
 'D3')

# D4: layer comparison restated nearly verbatim from Introduction para 4;
#     keep the layer contrast, drop the repeated eradicable clause
patch('Discussion_draft.md',
 'whereas the present study traces the axis to its upstream driver (infection status), dissects the antigen-specific immune response that carries the signal, and instruments the pathway\'s own feedback node (NFKBIA) — and, unlike circulating cytokines, the upstream exposure we identify is detectable and eradicable. To our knowledge',
 'whereas the present study tests the upstream driver (infection status), the antigen-specific immune response carrying the signal, and the pathway\'s feedback node (NFKBIA). To our knowledge',
 'D4')

# D5: subjective metaphor + unsupported assertion before the cited biology
patch('Discussion_draft.md',
 'The antigen-level result also provides a biologically sharper picture than the general serology alone. Among six *H. pylori* antigens, only the GroEL-specific immune response produced a nominally significant random-effects signal, directionally concordant with the general serology. This specificity is not random with respect to biology: GroEL is',
 'Among six *H. pylori* antigens, only the GroEL-specific immune response produced a nominally significant random-effects signal, directionally concordant with the general serology. GroEL is',
 'D5')

# D6: self-praise adjective ("coherent"); keep the epistemic marker
patch('Discussion_draft.md',
 'A coherent, though partly hypothesised, mechanistic chain integrates the convergent signals.',
 'The convergent signals admit a partly hypothesised mechanistic chain.',
 'D6')

# D7: emphasis adverb
patch('Discussion_draft.md',
 'the principal cytoplasmic inhibitor of NF-κB and, critically, an NF-κB target gene itself',
 'the principal cytoplasmic inhibitor of NF-κB and an NF-κB target gene itself',
 'D7')

# D7b: rhetorical framing before cited evidence
patch('Discussion_draft.md',
 'These inflammatory routes are not merely historical constructs: muscle-resident NF-κB signalling arising from multiple cell populations is now shown to actively orchestrate distinct catabolic states',
 'Muscle-resident NF-κB signalling arising from multiple cell populations actively orchestrates distinct catabolic states',
 'D7b')

# D9: paragraph repeating Results 3.4 (conservative RE + concordant WM);
#     keep the one new interpretive point, cross-reference the rest
patch('Discussion_draft.md',
 'The substantial between-instrument heterogeneity for NFKBIA→ALM (Cochran\'s Q = 27.34, df = 4, P = 1.7 × 10⁻⁵) is itself informative: the five cis-eQTLs likely capture more than one regulatory mechanism (e.g. direct expression effects and effects on transcript isoform usage), and the IVW random-effects estimate is the appropriately conservative primary estimate. That the weighted-median estimate (P = 2.3 × 10⁻⁸) sits between the random- and fixed-effect IVW estimates and agrees in direction supports that the causal signal is not an artefact of a single heterogeneous variant.',
 'The substantial between-instrument heterogeneity for NFKBIA→ALM (Cochran\'s Q = 27.34, df = 4, P = 1.7 × 10⁻⁵) likely reflects the five cis-eQTLs capturing more than one regulatory mechanism (e.g. direct expression effects and effects on transcript isoform usage); the conservative random-effects primary estimate and the concordant weighted-median estimate address this heterogeneity (Results 3.4).',
 'D9')

# D10: meta-commentary opener ("deserves comment") duplicating the section heading
patch('Discussion_draft.md',
 'The specificity of the causal effect for ALM — rather than grip strength or walking pace — deserves comment. ALM is a direct "quantity" measure',
 'ALM is a direct "quantity" measure',
 'D10')

# D11: filler list intro under an existing "Strengths." heading
patch('Discussion_draft.md',
 '**Strengths.** This study has several design features that guard against the known failure modes of MR: (i)',
 '**Strengths.** The design guards against known MR failure modes: (i)',
 'D11')

# D12: Nat Comms argument restated from Introduction para 3
patch('Discussion_draft.md',
 'This architecture operationalises the design principles highlighted by two recent methodological studies in *Nature Communications*, which showed that broad-panel MR screens generate many signals that fail rigorous validation and argued for targeted, hypothesis-driven designs with systematic cross-checking across complementary evidence sources (Zhang et al., 2026; Machnik et al., 2026)',
 'These features follow the design principles of two recent methodological studies in *Nature Communications*, which showed that broad-panel MR screens generate many signals that fail rigorous validation (Zhang et al., 2026; Machnik et al., 2026)',
 'D12')

# D13: vague premise opener ("epidemiologically grounded and increasingly causal")
patch('Discussion_draft.md',
 'The premise that *H. pylori* acts as a systemic-risk infection is epidemiologically grounded and increasingly causal: meta-analyses associate seropositivity with atherosclerotic cardiovascular disease (Wang et al., 2020)',
 'Beyond the gut, meta-analyses associate *H. pylori* seropositivity with atherosclerotic cardiovascular disease (Wang et al., 2020)',
 'D13')

# D14: self-congratulatory clause already covered by 4.2 first-MR claim
patch('Discussion_draft.md',
 ', both published in this journal; extending the same design to muscle phenotypes is the contribution of the present study. Future work should',
 ', both published in this journal. Future work should',
 'D14')

print('\n'.join(log))
print('TOTAL: %d patches' % len(log))
