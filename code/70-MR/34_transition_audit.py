# -*- coding: utf-8 -*-
"""Paragraph-transition and connective audit (2026-11-07 user directive):
delete empty connectives, fix logic gaps, make paragraphs respond to
predecessors. Enumeration connectives (First/Second/.../Finally in genuine
lists) and real however/by-contrast contrasts are retained.
Idempotent."""
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

# T1 (Intro, paragraph gap): para 1 ends "whether a modifiable upstream
# driver ... is unknown"; para 2 used to open with bare prevalence.
# Bridge: nominate H. pylori as a candidate for that driver (hypothesis
# framing, not asserted causality), so para 2 answers para 1's question.
patch('Introduction_draft.md',
 '*Helicobacter pylori* is one of the most prevalent chronic infections globally, colonising roughly 40–50% of the world\'s population with marked geographic and ancestral variation {{CITE:',
 'One candidate for such a driver is *Helicobacter pylori*, one of the most prevalent chronic infections globally: it colonises roughly 40–50% of the world\'s population with marked geographic and ancestral variation {{CITE:',
 'T1')

# T2 (Results 3.2, empty summary connective): "Together," adds nothing;
# the sentence is itself the synthesis of the two preceding facts.
patch('Results_draft.md',
 'Together, the divergence between the seropositivity and IgG-titre instruments',
 'The divergence between the seropositivity and IgG-titre instruments',
 'T2')

# T3 (Methods 2.4, duplicated clause): "Directional pleiotropy was further
# evaluated by the MR-Egger intercept test" repeats the parenthetical
# "(intercept detects directional pleiotropy)" two sentences earlier.
patch('Methods_draft.md',
 'Cochran\'s Q statistic assessed between-SNP heterogeneity. Directional pleiotropy was further evaluated by the MR-Egger intercept test, and the correct causal direction was verified with Steiger directionality tests',
 'Cochran\'s Q statistic assessed between-SNP heterogeneity. The correct causal direction was verified with Steiger directionality tests',
 'T3')

# T4 (Discussion 4.2, paragraph gap): the magnitude paragraph ends on
# "the observational and genetic lines agree"; the antigen paragraph then
# opened abruptly with "Among six antigens...". Factual bridge naming what
# the antigen-spectrum line does.
patch('Discussion_draft.md',
 'Among six *H. pylori* antigens, only the GroEL-specific immune response produced a nominally significant random-effects signal, directionally concordant with the general serology. GroEL is',
 'The antigen-spectrum analysis localises the infection signal further: among six *H. pylori* antigens, only the GroEL-specific immune response produced a nominally significant random-effects signal, directionally concordant with the general serology. GroEL is',
 'T4')

# T5 (Discussion 4.3, empty transitional sentence): "The NFKBIA result slots
# directly into this pathway." adds nothing the next sentence doesn't carry;
# the paragraph now opens directly with NFKBIA/IκBα, which responds to the
# preceding paragraph's NF-κB/MuRF1 ending by content.
patch('Discussion_draft.md',
 'The NFKBIA result slots directly into this pathway. NFKBIA encodes IκBα',
 'NFKBIA encodes IκBα',
 'T5')

print('\n'.join(log))
print('TOTAL: %d patches' % len(log))
