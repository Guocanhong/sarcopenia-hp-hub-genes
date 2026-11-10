# -*- coding: utf-8 -*-
"""48_editor_compliance.py — Medicine 编辑部审查补丁（幂等，带断言）

E1: 补 Acknowledgements / Author contributions / Funding / Conflicts of interest 四节
E2: ALSPAC 官方致谢措辞（2024-2026 多刊标准文本，实锤核实）
E3: Methods 2.6 加 STROBE-MR 报告规范声明（checklist 单独建补充文件）
A2: Limitations 补抗原滴度相关性条目（第十）
A3: Limitations 补性别分层/非线性不可行条目（第十一）
B1: Limitations 补 traits≠综合征定义条目（第九）
"""
import io
import os

PAPER = r"D:/GEO-limma,肌少症1/70-MR/paper"
DISC = os.path.join(PAPER, "Discussion_draft.md")
METH = os.path.join(PAPER, "Methods_draft.md")

LIMIT_NEW = (
    "Ninth, the outcomes analyzed are individual sarcopenia-related traits "
    "(grip strength, appendicular lean mass, walking pace, low grip strength) rather "
    "than a consensus-defined sarcopenia diagnosis, which requires a combined judgment "
    "of muscle mass and function; the estimates therefore characterize trait-level "
    "causal effects and should not be extrapolated to clinical sarcopenia "
    "classification. Tenth, antigen-specific antibody titers are correlated within "
    "individuals, so the six antigen lines are not mutually independent; antigen-level "
    "attribution (for example, OMP versus GroEL) cannot be adjudicated without "
    "multivariable MR, which is infeasible at the available instrument counts, and the "
    "antigen-spectrum findings should be read as evidence that the anti-*H. pylori* "
    "humoral response as a whole matters rather than as localization to a single "
    "antigen. Eleventh, the number of instruments per exposure (3\u201315 variants) "
    "was too small to support sex-stratified or exposure\u2013outcome nonlinearity "
    "analyses, which therefore could not be performed. "
)

DECL = """# Declarations

## Acknowledgements

We are extremely grateful to all the families who took part in the ALSPAC study, the midwives for their help in recruiting them, and the whole ALSPAC team, which includes interviewers, computer and laboratory technicians, clerical workers, research scientists, volunteers, managers, receptionists and nurses. The UK Medical Research Council and Wellcome (Grant ref: 217065/Z/19/Z) and the University of Bristol provide core support for ALSPAC; a comprehensive list of grants funding is available on the ALSPAC website (http://www.bristol.ac.uk/alspac/external/documents/grant-acknowledgements.pdf). The present analysis used only publicly deposited summary statistics; we thank the participants of ALSPAC, the UK Biobank and all contributing GWAS consortia, and the IEU OpenGWAS team. [Writing assistance disclosure to be completed before submission.]

## Author contributions

[To be completed before submission. Suggested CRediT assignment: Conceptualization, Methodology, Formal analysis, Writing \u2014 original draft: [initials]; Data curation, Software: [initials]; Supervision, Writing \u2014 review and editing: [initials].]

## Funding

[Grant/funding information to be completed before submission. If none, insert: This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.]

## Conflicts of interest

[To be confirmed before submission. Suggested standard text: The authors have no conflicts of interest to declare.]
"""


def patch(path, old, new, tag):
    s = io.open(path, encoding="utf-8").read()
    if new in s:
        print("SKIP (already applied):", tag)
        return
    assert s.count(old) == 1, "anchor not unique/missing in %s :: %s" % (path, tag)
    s = s.replace(old, new, 1)
    io.open(path, "w", encoding="utf-8", newline="").write(s)
    print("OK:", tag)


# --- B1/A2/A3: Limitations 第九/十/十一条（插在 Finally 收束句之前） ---
patch(
    DISC,
    "Finally, the hub-gene candidates derive from our integrative transcriptomic screening analysis",
    LIMIT_NEW + "Finally, the hub-gene candidates derive from our integrative transcriptomic screening analysis",
    "Limitations Ninth/Tenth/Eleventh",
)

# --- E3: Methods 2.6 STROBE-MR 声明 ---
patch(
    METH,
    "no additional ethical approval was required for the present secondary analysis of aggregate data.",
    "no additional ethical approval was required for the present secondary analysis of aggregate data. "
    "Reporting complies with the STROBE-MR guideline; a completed checklist is provided in the supplementary material.",
    "STROBE-MR statement in 2.6",
)

# --- E1+E2: Declarations 四节（插在 4.6 结尾与起草注记分隔线之间） ---
patch(
    DISC,
    "longitudinal lean-mass measurements.\n\n---\n",
    "longitudinal lean-mass measurements.\n\n---\n\n" + DECL + "\n---\n",
    "Declarations block (Acknowledgements/Contributions/Funding/COI)",
)

print("All patches processed.")
