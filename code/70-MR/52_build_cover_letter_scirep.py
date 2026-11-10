# -*- coding: utf-8 -*-
"""Build Cover Letter docx for Scientific Reports — MR manuscript.
Adapted from 50_build_cover_letter.py (Medicine version)."""
import docx
from docx.shared import Pt, Cm
from docx.enum.text import WD_LINE_SPACING

doc = docx.Document()

sec = doc.sections[0]
sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
sec.top_margin = sec.bottom_margin = Cm(2.5)
sec.left_margin = sec.right_margin = Cm(2.5)

st = doc.styles['Normal']
st.font.name = 'Times New Roman'
st.font.size = Pt(12)
st.element.rPr.rFonts.set(docx.oxml.ns.qn('w:eastAsia'), 'Times New Roman')
st.paragraph_format.space_after = Pt(10)
st.paragraph_format.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE


def add_runs(p, text, bold=False, italic=False, size=None):
    parts = text.split('*')
    for j, seg in enumerate(parts):
        if not seg:
            continue
        r = p.add_run(seg)
        r.bold = bold
        r.italic = italic or (j % 2 == 1)
        if size:
            r.font.size = Pt(size)


def para(text, bold=False, italic=False, align=None, size=None, space_after=10):
    p = doc.add_paragraph()
    if align is not None:
        p.alignment = align
    p.paragraph_format.space_after = Pt(space_after)
    add_runs(p, text, bold=bold, italic=italic, size=size)
    return p


para('November 9, 2026', space_after=16)

para('Editor-in-Chief')
para('*Scientific Reports* (Springer Nature)', space_after=16)

para('Subject: Submission of Original Research \u2014 "Causal effects of '
     '*Helicobacter pylori* infection on sarcopenia-related traits: a two-sample '
     'Mendelian randomization study implicating NF-\u03baB signaling in appendicular '
     'lean mass"', bold=True, space_after=16)

para('Dear Editor,', space_after=12)

body = [
    'We are pleased to submit the above manuscript for consideration as Original '
    'Research in *Scientific Reports*. The study is reported in accordance with the '
    'STROBE-MR guideline, and a completed checklist is uploaded with the submission.',

    'Chronic *Helicobacter pylori* infection has been repeatedly linked '
    'epidemiologically to lower grip strength, reduced appendicular lean mass, and '
    'sarcopenia in older adults, but these observations are cross-sectional and '
    'vulnerable to reverse causation and residual confounding \u2014 precisely the '
    'biases that a genetic design is built to mitigate. Whether the '
    'infection-muscle association has a causal component, and through which '
    'inflammatory pathway it might act, has not been tested genetically. We '
    'therefore conducted a two-sample Mendelian randomization (MR) program '
    'triangulating four independent lines of evidence: (i) seropositivity '
    '(N = 8,735) and circulating anti-*H. pylori* IgG titer (N = 4,683) against '
    'five muscle-outcome GWASs; (ii) six antigen-specific antibody levels (CagA, '
    'catalase, GroEL, OMP, UreA, VacA) against four of these outcomes under a '
    'Bonferroni grid; (iii) four blood hub genes (C1QA, NFKBIA, DTL, CTNNBIP1) '
    'using cis-eQTL instruments; and (iv) reverse-direction analyses with power '
    'quantification. Sensitivity estimators, heterogeneity diagnostics, Steiger '
    'directionality, and leave-one-out analyses were applied uniformly, and every '
    'signal was classified by a decision rule fixed in advance: only signals '
    'corroborated by the fixed-effect and weighted-median estimators are reported '
    'as supported, whereas random-effects-only signals are reported as suggestive.',

    'Genetically predicted seropositivity lowered appendicular lean mass '
    '(\u03b2 = \u22120.061 SD units per log-odds of seropositivity; P = 3.7 \u00d7 10\u207b\u2074; '
    '11 instruments), with concordant fixed-effect, weighted-median, and MR-Egger '
    'estimates. The remaining lines were individually suggestive but directionally '
    'convergent: the six-antigen line nominated OMP, and the cis-eQTL line linked '
    'higher predicted NFKBIA expression \u2014 an NF-\u03baB-related transcript \u2014 '
    'to lower lean mass: evidence that the pathway participates in the phenotype '
    'rather than of one specific direction of NF-\u03baB signaling. We report these '
    'signals with their decision-rule classifications stated explicitly rather '
    'than asserted as confirmed causal effects, and we discuss why seropositivity, '
    'an eradicable upstream exposure, is the actionable target of this axis.',

    'We believe the manuscript suits *Scientific Reports*, which publishes '
    'technically sound research judged on methodological rigor and completeness of '
    'reporting rather than perceived novelty. All analyses used publicly available '
    'GWAS summary statistics retrieved through the IEU OpenGWAS database; no new '
    'human or animal experimentation was performed. Instrument lists, harmonized '
    'effect estimates, and the full analysis code are deposited at '
    'https://github.com/Guocanhong/sarcopenia-hp-hub-genes (folder `70-MR`) and '
    'will be archived with a Zenodo DOI upon acceptance. In line with journal '
    'policy, we disclose that Tencent WorkBuddy (a GLM large language model, '
    'Tencent) was used to assist with language editing and analysis-code '
    'development; the manuscript contains a corresponding disclosure, and the '
    'authors reviewed and take full responsibility for all content.',

    'This manuscript is original, has not been published previously, and is not '
    'under consideration elsewhere. All authors have approved the submission and '
    'declare no conflicts of interest. This work was supported by the Natural '
    'Science Foundation of Fujian Province of China (Grant No. 2024J011631); the '
    'funder had no role in study design, data collection and analysis, decision to '
    'publish, or preparation of the manuscript.',

    'Thank you for your consideration.',
]
for b in body:
    para(b, space_after=12)

para('Sincerely,', space_after=24)

para('Guocan Hong, MD', bold=True, space_after=2)
para('Corresponding author, on behalf of all authors', space_after=2)
para('Department of Health Management, Affiliated Mindong Hospital,', space_after=2)
para('Fujian Medical University', space_after=2)
para('89 Heshan Road, Fu\u2019an 355000, Fujian Province, China', space_after=2)
para('Email: 95493638@qq.com', space_after=16)

out = r'D:/MRout/70-MR/submission/Cover_Letter_SciRep.docx'
doc.save(out)
print('saved:', out)
