<<<ANSWER>>>
[立场] 修正
核心结论：支持 P0 修复数据真实性与身份信息、P1 抽取共享模块、P2 建立核对协议与分阶段重跑的分层方案；但反对把 `revision_20261116/` 纳入本轮修改，反对把 `01_mr_estimators.py` 自检块强制移到 `tests/`，反对在未检索的情况下扩大 README 声明的修改范围。

[理由]
1. 硬编码结果数字：各方一致认定 `code/07_figures/02_figure1_design.py` L169–171 的 PTS 三元组必须从 `results/tables/concordance_by_cohort.tsv` 读入中间变量，标签由 rho 格式化生成；必须显式规定舍入规则（0.0976 → +0.10 / 0.202 → +0.20 / 0.313 → +0.31），保证图面已发表数字不变。`02_select_instruments.py` docstring 中的 “23.28%” 必须删除或改为指向结果表。
2. 身份信息：所有正式包脚本文件头的作者姓名、邮箱必须清除；`07_figures/01_main_figures_1_2_3.R` 中的期刊相关措辞必须改为中性描述。`06_EC_IR/paper_ms/skeleton_ebiomedicine.tex` 属文稿且不在本轮合并/重跑范围，本轮不改文件名。
3. 随机数：
   - bootstrap / 置换 / 模块抽样（7 个脚本）保留，种子必须来自 `config.SEED` 或 `--seed`；
   - `10_positive_control.py` 作为阳性对照 / 功效模拟保留，但输出必须隔离并明确标注模拟性质，不进主结果表；
   - `01_mr_estimators.py` 的 `if __name__ == "__main__"` 自检块（seed=1、TRUE=0.30、n=400）是单元自检，不写入 results；可接受方案为 (a) 保留 `__main__` 但改用 `config.SEED`，或 (b) 移到 `tests/` 作为独立测试。两种方案不改变正式结果。
   - `revision_20261116/proxy_uncertainty.py` 的 `default_rng(20261116)`：该目录不在主 agent 核实的正式可复现包范围内，本轮不改其种子，仅确认其输出不进入正式结果路径。
4. 去重策略：反对为减少文件数而合并 FLOW 阶段脚本或把 45 个脚本合成巨文件；应抽取真重复函数到共享模块。建议：
   - `code/common/io_utils.py`：`resolve_input`、`remote_size`
   - `code/common/locus_utils.py`：`_find_member`、`genes_near_loci`
   - `code/common/coloc_math.py`：`_logsum`
   - `code/common/module_sampling.py`：`zmat`、`build_background`、`draw_module`
   其中 `zmat`/`build_background`/`draw_module` 涉及零模型与得分的统计契约，抽取前必须先 diff + 回归测试；若实现不等价，则保留两份并改名，禁止强行统一。`_parse_args` 同名不合并。
5. README 声明：原句 “No analysis script contains a hard-coded path, accession or threshold. The only constant fixed in source is the random seed.” 中，“唯一固定常量是随机种子”已被 PTS 硬编码和 `01_mr_estimators.py` 的 seed=1 证伪，必须修正。但 “无硬编码 path/accession/threshold” 部分是否成立尚未核验，不应在未检索全部脚本的情况下直接改写。最终措辞应在主 agent 完成检索后确定。
6. 核对与重跑：必须建立文稿-表格对照表；已确认 concordance_by_cohort.tsv、incremental_value.tsv、table1_coloc_summary.tex 三组数字与 README 四舍五入一致。README 中 ΔAUC 的 95% CI（−0.014 to +0.031）来源不在 incremental_value.tsv 该行，标为“未核验”，必须由主 agent 定位到具体 bootstrap 输出文件。重跑分阶段：阶段 0 不依赖大文件的自检；阶段 1 依赖本地数据，缺文件即停并报告缺失清单；禁止用模拟数据替代。

[异议]
- 对 GLM / DeepSeek / Nemotron 提议把 `revision_20261116/proxy_uncertainty.py` 的随机种子改为 `config.SEED`：异议。该目录不在本轮正式可复现包范围（正式包为 `code/{00_setup,...,07_figures}` + config），修改其种子属于范围蔓延。正确做法是：明确排除在正式结果路径外，不动代码。
- 对 Mistral 提议把 `01_mr_estimators.py` 的 `__main__` 自检块强制移到 `tests/test_mr_estimators.py`：部分异议。移到 tests/ 是一种可接受方案，但不是唯一方案；保留 `__main__` 自检并改用 `config.SEED` 同样合规。应把两种方案都写入计划，由主 agent 根据仓库结构选定，不应硬性规定必须新建 tests 目录。
- 对各方普遍把 README 声明整句改写的倾向：异议。在未全面检索前不能确认 “无硬编码 path/accession/threshold” 也违例；本轮只能确认 “唯一固定常量是种子” 这部分不成立。README 修改范围应限定为先修正已证伪部分，其余措辞留待检索后定稿。
- 对 Nemotron 关于 `10_positive_control.py` “不写结果表”的表述：该脚本当前输出路径未核验，不能默认它已混入正式表。应增加核验任务：确认其输出是否已在 `results/tables/` 之外；若已在正确位置，则无需移动，只需加标注。
- 对 DeepSeek 把 `revision_20261116/` 中文注释翻译列为 P1：异议。同前，该目录不在本轮范围，不应修改其注释。

[修正]
无修改（本轮为首次发言，无上一轮立场）。

[收敛] NO
各方在核心方向（硬编码数字、身份信息、共享模块而非合并脚本、分阶段重跑）上高度一致，但仍有三个实质分歧未解决：(1) `revision_20261116/` 是否应改种子（我方主张不改）；(2) `01_mr_estimators.py` 自检块是否必须移到 `tests/`（我方主张两种方案均可）；(3) README 声明修改范围是否应限定为已证伪部分（我方主张先检索再定全文措辞）。需下一轮就这三点达成一致后，方可 `[收敛] YES`。
<<<END>>>