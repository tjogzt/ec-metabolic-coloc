# 任务

# 任务：对现有分析代码做系统性审查，并收敛成一份可执行的修改方案

审查对象是仓库 `/Users/taozhu/my researches/gongbing` 中的正式分析包
`ec-metabolic-coloc/`（约 45 个 .py/.R，约 14520 行）。
同仓库的 `06_EC_IR/`、`02_literature/`、`03_proposal/`、`05_comorb/` 是早期探索、
文稿草稿和第三方包（如 scTenifoldKnk），**不在本次合并/重跑范围内**。
文稿数字来源：`06_EC_IR/paper_ms/` 与 `01_manuscript.pdf` 所对应的表格，
以及 `ec-metabolic-coloc/README.md`、`docs/FLOW.md`、`results/tables/`。

用户硬性要求（辩论必须逐条回应，不能用“建议性”语气把硬要求软化掉）：

1. 代码质量：是否精简、是否有冗余/重复模块、是否应合并以减少文件数、防止腐化。
2. 数据真实性：禁止非种子随机数、硬编码结果数字、模拟/虚假数据。
   数据必须先读入中间变量再传递，禁止把结果数字直接嵌进分析代码。
3. 运行验证：审查后要完整重跑全部代码并修到能跑通。本轮只制定计划，不执行重跑。
4. 数据核对：把代码输出整理成对照表，与文稿逐项交叉核对，标出不一致。
5. 注释：统一英文；文件头注释 + 必要行间注释；清除代码中的期刊名、作者姓名、邮箱等身份信息。
6. 辩论必须多轮，直到实际到场的全体成员对修改方案一致，才允许 `[收敛] YES`。
7. 主 agent 形成修改计划（修改项、理由、优先级、预期效果）交用户批准后才能改代码。
   **本轮禁止提议立刻改仓库。**

## 已由主 agent 在本地核实的事实（可直接引用，不要再发明文件或数字）

### A. 范围
- 正式可复现包：`ec-metabolic-coloc/code/{00_setup,02_colocalisation,03_burden,04_singlecell,05_mendelian_randomisation,06_cross_cancer,07_figures}` + `config/config.py` + `config/config.R`。
- `revision_20261116/` 三份修订脚本含中文注释，且 `proxy_uncertainty.py` 使用 `np.random.default_rng(20261116)`，种子未走 `config.SEED`。
- README 声称：“No analysis script contains a hard-coded path, accession or threshold. The only constant fixed in source is the random seed.” 这句与下面的发现不完全相符，辩论时要判断这句话是否应改。

### B. 随机数（全部带种子，没有发现裸 `np.random.rand` / `random.random`）
用途分三类，请分别裁决“保留 / 移出正式结果路径 / 必须改写”：
- 统计程序（bootstrap、置换、模块抽样）：`04_singlecell/02_null_calibration.py`、`03_burden/05_build_programme.py`、`03_burden/07_incremental_value.py`、`02_colocalisation/12_permutation_null.py`、`02_colocalisation/15_locus_screen.py`、`06_cross_cancer/03_score_across_cancers.py`、`06_cross_cancer/04_normal_tissue_null.py`。种子来自 `config.SEED` 或 `--seed`。
- 模拟功效 / 阳性对照的效应量：`02_colocalisation/10_positive_control.py` 明确写 “only simulate effect sizes”，位点、等位基因、标准误来自真实数据。默认 seed=`config.SEED`。
- 合成自检，且种子写死为 1：`05_mendelian_randomisation/01_mr_estimators.py` 的 `if __name__ == "__main__"`，`rng = np.random.default_rng(1)`，真值 `TRUE = 0.30`，`n = 400`。这是单元自检，不写 results 表。
- R：`config/config.R` 与 `07_figures/01_main_figures_1_2_3.R` 有 `set.seed(SEED)`。

### C. 硬编码结果数字（已定位）
`code/07_figures/02_figure1_design.py` 约 L169–171 把文稿数字直接画进示意图：
```
PTS = [(74.0, 0.313, "+0.31", ... "adipose"),
       (83.0, 0.202, "+0.20", ... "diabetes"),
       (92.0, 0.098, "+0.10", ... "clamp")]
```
这些数字与 `results/tables/concordance_by_cohort.tsv` 的 rho 一致（0.313、0.202、0.098），但是**嵌在作图代码里**，不是从该表读入。
同文件其余 15.0 / 16.6 等是版式坐标，不是科学结果。

`05_mendelian_randomisation/02_select_instruments.py` L59 附近 docstring 写了 “23.28% of its genome-wide significant rows”，这是注释里的结果数字，不是计算输入。

### D. 重复函数（跨文件同名，需判断是真重复还是巧合）
- `_find_member`, `genes_near_loci`：`02_colocalisation/06_define_risk_loci.py` 与 `13_gtex_v11_reader.py`
- `_logsum`：`05_coloc_abf_engine.py` 与 `08_prior_sensitivity.py`
- `resolve_input`, `remote_size`：`10_positive_control.py` 与 `15_locus_screen.py`
- `zmat`, `build_background`, `draw_module`：`06_cross_cancer/03_score_across_cancers.py` 与 `04_normal_tissue_null.py`
- `_parse_args`：两个 singlecell 脚本（同名不代表应合并）

用户希望减少文件数。正式包已经按阶段编号，合并会破坏 `docs/FLOW.md` 的运行顺序。请明确：哪些函数应抽到共享模块、哪些文件**不应**合并，以及共享模块的建议文件名。不要提议把 45 个脚本合成少数几个巨文件。

### E. 身份信息（必须清除的范围请你们定）
几乎每个正式脚本文件头都有：
`Author   : Tao Zhu (zhutao@tjh.tjmu.edu.cn)`
`07_figures/01_main_figures_1_2_3.R` 还有 “the journal requirement for composite figures”（无具体刊名）。
`06_EC_IR/paper_ms/skeleton_ebiomedicine.tex` 文件名含期刊名，但那是文稿不是代码；本轮是否动文稿文件名，请单独表态。

### F. 已核对的文稿关键数字（来自 results 表，不是模型记忆）
`results/tables/concordance_by_cohort.tsv`：
- GSE13070 clamp：rho=0.0976，p_rho=0.00412，responsive concordance=61.62%（n_responsive=99）
- GSE18732 T2D：rho=0.202
- GSE15773 adipose：rho=0.313
README 写成 +0.10 / +0.20 / +0.31 与 61.6%，与表一致（四舍五入）。

`results/tables/incremental_value.tsv`：
delta_auc_IR=0.000389，p_IR=0.15，n=145，events=45。
README：ΔAUC=+0.0004，p=0.15，145 tumours，45 events。一致（四舍五入）。
README 还写 95% CI −0.014 to +0.031；该 CI **不在** incremental_value.tsv 这一行里，可能在 incremental_boot / 另一张表。核对计划必须把 CI 的来源文件标出来，不能假定它正确。

`paper_ms/tables/table1_coloc_summary.tex`：五组织检验合计 503，PP.H4>0.8 共 3（liver 1，muscle 2），最大 PP.H4=0.962。
README：503 tests，3 exceed PP.H4>0.8，liver PP.H4=0.93。0.93 与表中 liver Max PP.H4=0.931 一致；0.962 是 muscle 的最大值，README 未把它说成 liver。

本轮**不要**声称已经完成全文数字核对。请设计核对协议：哪些 tex/README 句子对哪张 tsv，由主 agent 在你批准后执行。

### G. 重跑约束（计划必须尊重）
- 原始数据不在 git 里，`data/README.md` 只列 accession。完整 colocalisation / TCGA / 单细胞重跑可能缺本地文件。
- 计划必须分成：不依赖大文件就能跑的自检；依赖本地数据、缺文件就停止并报告的阶段；禁止为了“跑通”而改用模拟数据。
- 用户规则：分析必须用真实数据；需要模拟数据时必须先征得用户同意。功效模拟（10_positive_control）已存在，是否保留要你们裁决，但不能新增模拟数据路径。

## 你们必须交付的修改方案结构

在 [建议] 里给出优先级 P0/P1/P2 的修改项，每项含：改哪些文件、为什么、预期效果、是否会改变已发表数字。
明确写出你**反对**的做法（例如：为了减少文件数合并 FLOW 阶段；删除所有带种子的 bootstrap；把示意图数字改成从错误的表读取）。

文献和未在上文出现的数字一律标“未核验”，不要编造。




# 你的任务

你是 **Kimi（moonshotai/kimi-k2.7-code）**，在 BioHex 六合团队中负责：代码实现 / 工程复现。
请**独立**给出你的最佳分析；不要假设其他成员会同意你，也不要等待他人。

严格按以下结构输出（不要添加其它章节）：

<<<ANSWER>>>
[结论] 一句话结论
[理由] 你的核心论据（可列 2–6 点）
[风险与不确定] 主要风险、混杂、证据缺口
[建议] 可执行下一步（可列 1–4 点）
<<<END>>>
