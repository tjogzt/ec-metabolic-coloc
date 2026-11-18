<<<ANSWER>>>
[结论]
本轮应形成一份“P0 修复数据真实性与身份信息、P1 抽取共享模块并补随机数 provenance、P2 建立文稿–表格核对协议”的修改计划；反对为减少文件数而合并 FLOW 阶段文件、反对删除带种子的 bootstrap/置换、反对用模拟数据补缺失原始数据；本轮只定计划，不改仓库、不重跑、不收敛。

[理由]
1. 硬编码违规点已定位且修复后不改变已发表数字：`07_figures/02_figure1_design.py` L169–171 的 PTS 数值（0.313/0.202/0.098）与 `results/tables/concordance_by_cohort.tsv` 一致，因此把它们改为从该表读入、标签由 rho 计算生成（round 至 2 位得到 “+0.31/+0.20/+0.10”），图面不变，只改变数据来源；同文件 15.0/16.6 等版式坐标保留。`05_mendelian_randomisation/02_select_instruments.py` docstring 里的 “23.28%” 属注释内结果数字，应删掉或改为指向 results 表，不能留在源码里。
2. 随机数三类裁决：① bootstrap/置换/模块抽样（`02_null_calibration.py`、`05_build_programme.py`、`07_incremental_value.py`、`12_permutation_null.py`、`15_locus_screen.py`、`03_score_across_cancers.py`、`04_normal_tissue_null.py`）全部保留，前提是种子统一来自 `config.SEED` 或 `--seed`，并把所用 seed 写入输出表或 run manifest 做 provenance；② `10_positive_control.py` 的功效模拟保留（已有、不新增），但输出必须隔离到独立目录（如 `results/positive_control/`）并在文稿中标注为模拟，不得混入正式结果表；③ `01_mr_estimators.py` 的 `__main__` 自检（`default_rng(1)`、TRUE=0.30、n=400）应移到 `tests/` 作为单元测试，seed 不再写死为 1，真值作为测试 fixture 保留在测试文件内，不属于分析路径。
3. 重复函数应“抽共享模块”而不是“合并脚本”：`_logsum` → 建议 `code/common/numerics.py`；`_find_member`/`genes_near_loci` → `code/common/locus_utils.py`；`resolve_input`/`remote_size` → `code/common/io_utils.py`；`zmat`/`build_background`/`draw_module` → `code/common/module_sampling.py`（注意 `04_normal_tissue_null.py` 是 null 路径，抽共享后其随机性仍须走 seed 并写 manifest）。`_parse_args` 同名不合并。任何抽取前必须先 diff 两份实现，确认逻辑一致；若有边界差异，以不改变现有输出为准，不能盲合。
4. README 声明必须改写：原句 “No analysis script contains a hard-coded path, accession or threshold. The only constant fixed in source is the random seed.” 被三类事实反驳——figure 脚本嵌入结果数字、自检脚本写死 seed=1、docstring 含结果数字。应改为准确声明，例如只声称“分析输入均来自 data/ 与 results/ 中间表；随机种子统一来自 config.SEED；模拟仅限标注明确的 positive control 与 tests/”。
5. 核对协议应设计为“句子→表”矩阵：README/tex 中每个数字声明映射到具体 tsv 的文件+列+行+四舍五入规则+状态。已知三组一致（四舍五入口径）：concordance_by_cohort.tsv 的 rho=0.0976/0.202/0.313、p_rho=0.00412、61.62%（n_responsive=99）↔ README 的 +0.10/+0.20/+0.31 与 61.6%；incremental_value.tsv 的 delta_auc_IR=0.000389、p_IR=0.15、n=145、events=45 ↔ README 的 ΔAUC=+0.0004、p=0.15、145、45；table1_coloc_summary.tex 的 503 tests / 3 个 PP.H4>0.8（liver 1、muscle 2）/ 最大 0.962 ↔ README 的 503、3、liver 0.93（表内 liver max 0.931）。README 的 95% CI −0.014 to +0.031 来源不在 incremental_value.tsv 本行，标“未核验”，必须由主 agent 定位到具体的 bootstrap 输出表后才能认定。
6. 重跑必须分阶段并设停止规则：Stage 0 不依赖大文件——import 全部模块、config 加载、tests/ 单测、以及只依赖 `results/tables/` 的 figure 脚本；Stage 1 依赖本地原始数据（colocalisation/TCGA/单细胞/MR），先核对 `data/README.md` 的 accession 对应文件是否存在，缺失则停止并输出缺失清单，禁止用模拟数据替代；Stage 2 执行核对协议并输出对照表。每阶段写 seed 与软件版本到 run manifest。

[风险与不确定]
- 共享模块抽取的最大风险是两份同名实现存在隐性差异（边界条件、容差），盲合并会悄悄改变结果；缓解办法是先 diff，并在 Stage 0/1 重跑中对改动前后输出做逐表比对（完全一致或在数值容差内），不一致即回退。
- README 的 ΔAUC 95% CI（−0.014 to +0.031）来源不明，属于证据缺口，未核验；若最终定位不到支撑表，应从 README 删除或改为“未核验”，不能假定正确。
- 原始数据不在 git 内，完整重跑可能因缺文件而无法完成；这是可复现性风险，只能报告缺失清单，不能靠模拟数据“跑通”。
- 用户“减少文件数”与“防腐化”存在张力：新增 3–4 个共享模块会让文件数不降反略增，但能消除跨文件重复；我的判断是不追求把 45 个脚本硬压成少数巨文件，净目标是“重复代码归零、阶段文件不动”。
- 重跑验证时还应检查输出是否含多重检验校正 p、效应量 CI、批次/队列分层与独立实验单位（patient/library/cell）口径，避免伪重复混入。

[建议]
P0（必须先做，预期不改变已发表数字）：
- P0-1 `code/07_figures/02_figure1_design.py`：PTS 改为读取 `results/tables/concordance_by_cohort.tsv`，按 cohort（adipose/diabetes/clamp）映射 rho，标签由 rho 四舍五入生成；改文件、因硬编码结果数字违规、效果是图随表自动更新、不改变数字。
- P0-2 全部正式包 .py/.R 文件头：删除 `Author: Tao Zhu (zhutao@tjh.tjmu.edu.cn)` 及一切作者/邮箱/期刊措辞；R 文件中 “journal requirement for composite figures” 改为中性表述；改所有 code/ 下文件、因身份信息清除要求、效果是源码去身份化、不改变数字。
- P0-3 `05_mendelian_randomisation/01_mr_estimators.py`：把 `__main__` 自检移到 `tests/test_mr_estimators.py`，seed 改由 `config.SEED` 或测试参数给出，不再写死 1；改该文件+新增测试文件、因写死 seed 违反统一种子规则、效果是自检可复现且不混入分析路径、不改变结果表。
- P0-4 `05_mendelian_randomisation/02_select_instruments.py`：删除 docstring 中 “23.28%…” 结果数字，改为指向对应 results 表；以及 `ec-metabolic-coloc/README.md` 声明句改写；改这两处、因 README 声明与事实不符、效果是声明可核验、不改变数字。

P1（计划批准后实施）：
- P1-1 抽取共享模块：`code/common/numerics.py`（`_logsum`）、`code/common/locus_utils.py`（`_find_member`、`genes_near_loci`）、`code/common/io_utils.py`（`resolve_input`、`remote_size`）、`code/common/module_sampling.py`（`zmat`、`build_background`、`draw_module`）；先 diff 验证一致再迁移，各脚本改为 import；改涉及的 6 个文件+新增 4 个模块、因跨文件重复是腐化温床、效果是单一实现源、预期不改变输出（需重跑比对确认）。
- P1-2 随机数 provenance：所有 bootstrap/置换/模块抽样脚本把实际使用的 seed 写入输出表新增列或 `results/run_manifest.tsv`；`10_positive_control.py` 输出移到 `results/positive_control/` 并在 README/文稿标注“simulated”；改 B 节列出的 7 个统计脚本+1 个模拟脚本、因可复现性要求、效果是任何结果可追溯到 seed、不改变数值。
- P1-3 核对协议：建立核对矩阵（README/tex 句子 → tsv 文件+列+行 → 四舍五入规则 → 一致/不一致/未核验），由主 agent 执行；首批必查项为 concordance、incremental_value、table1_coloc_summary 三组，并把 ΔAUC 的 95% CI 来源定位任务列为必完成项（未核验直至定位）。

P2（后续）：
- P2-1 全面统一注释为英文（含 `revision_20261116/` 三份脚本的中文注释，若它们被并入正式包；否则保持在修订目录隔离）；清掉残留身份信息。
- P2-2 文稿文件名含刊名（`06_EC_IR/paper_ms/skeleton_ebiomedicine.tex`）单独表态：本轮不动，因其不在代码审查范围且改名会破坏引用路径，留待文稿阶段由主 agent 决定。
- P2-3 重跑三阶段执行（Stage 0 自检+figures → Stage 1 依赖数据阶段带缺失即停 → Stage 2 核对输出对照表），每阶段写 manifest；禁止新增模拟数据路径。

明确反对的做法：① 反对为减少文件数而合并 `docs/FLOW.md` 的阶段文件或把 45 个脚本合成巨文件；② 反对删除/削弱带种子的 bootstrap、置换与模块抽样；③ 反对把示意图数字改从 concordance 以外的表读取（必须来自 `concordance_by_cohort.tsv`）；④ 反对用模拟/合成数据替代缺失的原始数据来“跑通”；⑤ 反对本轮立刻修改仓库或宣布 `[收敛] YES`——辩论需多轮、全体一致后才允许收敛。
<<<END>>>