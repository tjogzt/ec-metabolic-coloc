<<<ANSWER>>>
[结论] ec-metabolic-coloc 正式包应保留阶段化结构，通过抽取共享模块、移除硬编码结果数字与身份信息、统一随机数种子、建立表-文稿核对协议并分阶段重跑验证，方可进入正式复现；本轮仅制定计划，不修改仓库，亦未达成 `[收敛] YES`。

[理由]
1. 范围与文稿界定：仅 `ec-metabolic-coloc/code/{00_setup,02_colocalisation,03_burden,04_singlecell,05_mendelian_randomisation,06_cross_cancer,07_figures}` + config 在合并/重跑范围；`06_EC_IR`、`02_literature`、`03_proposal`、`05_comorb` 及 `paper_ms` 不在范围；`skeleton_ebiomedicine.tex` 属文稿文件名含期刊名，本轮不改。
2. 随机数裁决：统计程序（bootstrap/置换/模块抽样）保留，种子须来自 `config.SEED` 或 `--seed`；`proxy_uncertainty.py` 的 `np.random.default_rng(20261116)` 必须改为读取 `config.SEED`；`10_positive_control.py` 的模拟效应量保留为阳性对照，不写结果表；`01_mr_estimators.py` 的 `if __name__=="__main__"` 自检（seed=1，TRUE=0.30，n=400）保留为单元自检，不进入正式结果路径；R 侧 `set.seed(SEED)` 保留。
3. 硬编码结果数字必须移出：`02_figure1_design.py` L169-171 的 PTS 列表（74.0/0.313/+0.31/adipose 等）须从 `results/tables/concordance_by_cohort.tsv` 读入，禁止嵌入作图代码；`02_select_instruments.py` L59 附近 docstring 中的 “23.28%” 是注释性结果数字，应删除或改为描述性文字；README 中“无分析脚本含硬编码路径/accession/threshold”的声明与事实不符，应修正。
4. 代码质量与合并策略：反对把 45 个脚本合成少数巨文件，反对为减少文件数而合并 FLOW 阶段（会破坏 `docs/FLOW.md` 运行顺序）；应抽取真重复函数到共享模块：`_find_member`/`genes_near_loci` → `loci_utils.py`；`_logsum` → `coloc_math.py`；`resolve_input`/`remote_size` → `io_utils.py`；`zmat`/`build_background`/`draw_module` → `cross_cancer_utils.py`；`_parse_args` 两个 singlecell 脚本同名不代表应合并，逻辑一致时可抽至 `cli.py`，否则保留各自。
5. 身份信息清除：所有正式脚本文件头的 `Author : Tao Zhu (zhutao@tjh.tjmu.edu.cn)` 及类似邮箱/姓名必须清除；`07_figures/01_main_figures_1_2_3.R` 中 “the journal requirement for composite figures” 等期刊相关表述应改为中性描述；注释统一为英文，添加文件头与必要行间注释。
6. 数据核对协议：必须建立 tex/README 句子与 tsv 文件的逐项映射，由主 agent 批准后执行：`concordance_by_cohort.tsv` ↔ GSE13070/GSE18732/GSE15773 的 rho/p_rho/responsive concordance 及 README 四舍五入表述；`incremental_value.tsv` ↔ delta_auc_IR/p_IR/n/events 及 README ΔAUC/p/145 tumours/45 events，并标出 95% CI (−0.014 to +0.031) 的来源文件（可能在 incremental_boot 或另一张表），不得假定正确；`table1_coloc_summary.tex` ↔ 503 tests、PP.H4>0.8 共 3（liver 1、muscle 2）、Max PP.H4=0.962，README liver PP.H4=0.931 与表一致，0.962 为 muscle 最大值，README 未误指为 liver。本轮不声称已完成全文数字核对。

[风险与不确定]
- 原始数据不在 git 中，`data/README.md` 仅列 accession，完整 colocalisation / TCGA / 单细胞重跑可能因缺本地文件无法执行，重跑计划的可行性不确定。
- `incremental_value.tsv` 中未包含 95% CI，CI 来源文件未定位，若主 agent 未找到对应表，README 中的 CI 表述无法核对。
- 重复函数是否为“真重复”尚需逐函数比对签名与行为，合并前需单元测试确保输出一致，跨文件同名不一定应合并。
- `proxy_uncertainty.py` 种子硬写 20261116 是否影响已发表结果，需确认其是否写入 results 表；若已写入，改种子可能改变数字。
- `02_figure1_design.py` 从 tsv 读入后，若表格列名或数值精度与作图坐标不完全匹配，示意图可能出现微小偏差，需确认是否允许。
- 身份信息清除若遗漏注释或字符串，可能违反脱敏要求；R 脚本中的 `set.seed` 变量来源需确认与 Python `config.SEED` 一致。
- 未在上文出现的数字一律标“未核验”，不能编造文献或数据库条目；多重检验、批次/组成混杂、独立实验单位等问题在现有代码审查中未逐项展开，属证据缺口。

[建议]
1. P0（阻塞项，用户批准前不改代码）：提交修改计划清单（修改项、理由、优先级、预期效果、是否改变已发表数字）供用户批准；执行硬编码数字移出（`02_figure1_design.py` 从 `results/tables/concordance_by_cohort.tsv` 读取 rho/cohort；`02_select_instruments.py` 删除 docstring 中的 “23.28%”；修正 README 中与事实不符的“无硬编码”声明）；统一随机数种子（`proxy_uncertainty.py` 改用 `config.SEED`，确认其他脚本种子来自 `config.SEED` 或 `--seed`）；清除身份信息（所有正式脚本文件头删除作者/邮箱；`07_figures/01_main_figures_1_2_3.R` 删除期刊表述；`skeleton_ebiomedicine.tex` 属文稿，本轮不动）。
2. P1（代码结构）：新建 `code/utils/`，抽取 `_find_member`/`genes_near_loci` → `loci_utils.py`；`_logsum` → `coloc_math.py`；`resolve_input`/`remote_size` → `io_utils.py`；`zmat`/`build_background`/`draw_module` → `cross_cancer_utils.py`；`_parse_args` 若逻辑一致可抽至 `cli.py`，否则保留各自单cell脚本。明确不合并 45 脚本为巨文件，不合并 FLOW 阶段，不删除带种子的 bootstrap/置换/模块抽样，不新增模拟数据路径。反对的做法：为了减少文件数而合并 FLOW 阶段或把 45 个脚本合成少数巨文件；删除所有带种子的 bootstrap/置换/模块抽样；把示意图数字改成从错误的表读取，或继续嵌在作图代码里；新增模拟数据路径或为“跑通”而替换真实数据；在未核验前声称已完成全文数字核对。
3. P1（数据核对）：建立 tex/README 与 tsv 的逐项映射表（`concordance_by_cohort.tsv` ↔ GSE13070/GSE18732/GSE15773 rho/p_rho/responsive concordance 及 README 四舍五入；`incremental_value.tsv` ↔ delta_auc_IR/p_IR/n/events 及 README ΔAUC/p/145 tumours/45 events，并标注 95% CI −0.014 to +0.031 的来源文件待查，不得假定正确；`table1_coloc_summary.tex` ↔ 503 tests、PP.H4>0.8 共 3（liver 1、muscle 2）、Max PP.H4=0.962，README liver PP.H4=0.931 与表一致，0.962 为 muscle 最大值）。由主 agent 批准后执行交叉核对并输出不一致对照表；本轮不声称已完成全文核对。
4. P2（重跑与注释）：统一注释为英文，添加文件头与必要行间注释；保留 `10_positive_control.py` 作为阳性对照/功效模拟（明确不进入正式结果表），保留 `01_mr_estimators.py` 的 `if __name__=="__main__"` 自检（seed=1，TRUE=0.30，n=400，不写 results 表）；重跑分三阶段——阶段一无大文件依赖的自检；阶段二依赖本地数据，缺文件停止并报告；阶段三完整 colocalisation/TCGA/单细胞重跑，缺原始数据时停止，禁止为“跑通”而引入模拟数据或替换真实数据。
<<<END>>>