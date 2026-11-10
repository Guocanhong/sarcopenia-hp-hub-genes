# -*- coding: utf-8 -*-
"""49_authors_funding.py — 按上篇（论文_Medicine_submission.docx, 60-Review/Medicine_Submission_Package）
的作者/单位/通讯/基金/贡献分工/COI/AI 声明，填入本篇 MR 稿的标题页与 Declarations 占位。
来源（实锤）:
  - 上篇标题页: manuscript_for_review.txt (60-Review/Medicine_Submission_Package)
  - 上篇 Declarations: 24-Declarations/Declarations_英文初稿_v1.md
幂等 + 断言。"""
import io

FILES_TITLE = r"code/19_build_submission_md.py"
FILES_DECL = r"paper/Discussion_draft.md"

AUTHORS = ("**Authors:** Zhenzhen Cai^1^, Shaojian Chen^2^, Jiangchao Qian^1^, "
           "Wen Lin^2^, Guocan Hong^2,\\*^")
AFFIL = ("**Affiliations:** ^1^Department of Pain Management, Affiliated Mindong Hospital, "
         "Fujian Medical University, 89 Heshan Road, Fu'an 355000, Fujian Province, China; "
         "^2^Department of Health Management, Affiliated Mindong Hospital, "
         "Fujian Medical University, 89 Heshan Road, Fu'an 355000, Fujian Province, China")
CORR = ("**Correspondence:** Guocan Hong, Department of Health Management, "
        "Affiliated Mindong Hospital, Fujian Medical University, 89 Heshan Road, "
        "Fu'an 355000, Fujian Province, China. Email: 95493638@qq.com")

CONTRIB = ("Zhenzhen Cai: conceptualization, methodology, formal analysis, data curation, "
           "visualization, and writing — original draft. Shaojian Chen: data curation, "
           "investigation, and validation. Jiangchao Qian: data curation, software, and "
           "visualization. Wen Lin: investigation and data curation. Guocan Hong: "
           "conceptualization, supervision, project administration, funding acquisition, and "
           "writing — review and editing. All authors contributed to data interpretation, "
           "read and approved the final version of the manuscript, and agree to be accountable "
           "for all aspects of the work.")

FUNDING = ("This work was supported by the Natural Science Foundation of Fujian Province of "
           "China (grant number 2024J011631, \\u201cSystematic study on the correlation between "
           "*Helicobacter pylori* infection and sarcopenia and observation of its clinical "
           "application\\u201d). The funder had no role in study design, data collection and "
           "analysis, decision to publish, or preparation of the manuscript.")

COI = "The authors declare no competing interests."

AI_DISC = ("During the preparation of this work (2026), the authors used Tencent WorkBuddy "
           "(a GLM large language model, Tencent) to assist with language editing, manuscript "
           "drafting and revision, and analysis code development. No generative AI tool was "
           "used to collect, generate, or analyze the experimental data, and no AI-generated "
           "images are included in this manuscript. The study was designed, all analyses were "
           "performed and validated, and all AI-assisted output was reviewed and edited by the "
           "authors, who take full responsibility for the accuracy and integrity of the "
           "published work. The AI tool was not used to generate scientific results or images "
           "and is not an author of this manuscript.")


def patch(path, pairs):
    s = io.open(path, encoding="utf-8").read()
    n = 0
    for old, new in pairs:
        if new in s and old not in s:
            print("SKIP (already applied):", path, "::", new[:55])
            continue
        assert old in s, path + " :: NOT FOUND :: " + old[:70]
        assert s.count(old) == 1, path + " :: not unique :: " + old[:70]
        s = s.replace(old, new, 1)
        n += 1
        print("OK:", new[:60])
    io.open(path, "w", encoding="utf-8", newline="").write(s)
    return n


# ---- 1) title page in build script ----
t_pairs = [
    ("**Authors:** [Author list and affiliations to be finalised before submission — corresponding author: Guo Canhong]", AUTHORS),
    ("**Affiliations:** [Department / Institution / City, China — to be finalised]", AFFIL),
    ("**Correspondence:** [Name, Department, Address, Email — to be finalised]", CORR),
]
n1 = patch(FILES_TITLE, t_pairs)

# ---- 2) Declarations in Discussion draft ----
d_pairs = [
    # Author contributions placeholder
    ("[To be completed before submission. Suggested CRediT assignment: Conceptualization, "
     "Methodology, Formal analysis, Writing — original draft: [initials]; Data curation, "
     "Software: [initials]; Supervision, Writing — review and editing: [initials].]",
     CONTRIB),
    # Funding placeholder
    ("[Grant/funding information to be completed before submission. If none, insert: This "
     "research did not receive any specific grant from funding agencies in the public, "
     "commercial, or not-for-profit sectors.]",
     FUNDING),
    # COI placeholder
    ("[To be confirmed before submission. Suggested standard text: The authors have no "
     "conflicts of interest to declare.]",
     COI),
    # Writing assistance disclosure in Acknowledgements
    ("[Writing assistance disclosure to be completed before submission.]",
     AI_DISC),
]
n2 = patch(FILES_DECL, d_pairs)

print("\ntotal patched:", n1 + n2)

# readback
s = io.open(FILES_TITLE, encoding="utf-8").read()
assert "Guocan Hong^2" in s and "95493638" in s
d = io.open(FILES_DECL, encoding="utf-8").read()
assert "2024J011631" in d and "no competing interests" in d
assert "Zhenzhen Cai: conceptualization" in d
assert "[To be completed before submission" not in d
assert "[Grant/funding information" not in d
assert "[Writing assistance disclosure" not in d
print("readback OK")
