# 1-parse_hpa.py  从 HPA 子页面提取 IHC 图 URL + 抗体 + 强度
import re

RAW = "raw"
genes = ["ENSG00000173372-C1QA", "ENSG00000100906-NFKBIA", "ENSG00000178585-CTNNBIP1"]
tissues = ["skeletal_muscle", "stomach", "spleen"]
tlabel = {"skeletal_muscle": "Skeletal muscle", "stomach": "Stomach", "spleen": "Spleen"}

for g in genes:
    for t in tissues:
        html = open(f"{RAW}/{g}_{t}.html", encoding="utf-8", errors="ignore").read()
        # 抗体标记位置（HPAxxxxx / CABxxxxx 作为区块头）
        ab_marks = [(m.start(), m.group(0)) for m in re.finditer(r"\b(?:HPA|CAB)\d{6}\b", html) ]
        # 全尺寸图链接（href 里的非 medium .jpg）
        hrefs = [(m.start(), m.group(1)) for m in re.finditer(r'href="(//images\.proteinatlas\.org/[^"]+\.jpg)"', html)
                 if "medium" not in m.group(1)]
        seen, rows = set(), []
        for pos, u in hrefs:
            if u in seen: continue
            seen.add(u)
            ab = next((a for p, a in reversed(ab_marks) if p < pos), "?")
            rows.append((ab, u))
        print(f"\n=== {g} | {tlabel[t]} : {len(rows)} imgs ===")
        for ab, u in rows:
            print(f"  [{ab}] https:{u}")
