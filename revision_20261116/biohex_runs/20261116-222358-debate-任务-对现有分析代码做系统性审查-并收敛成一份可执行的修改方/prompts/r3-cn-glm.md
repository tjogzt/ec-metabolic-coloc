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




# 上一轮各方立场

### cn-deepseek DeepSeek (deepseek/deepseek-v4-pro)

<<<ANSWER>>>
[立场] 修正：P0 修复硬编码数字与身份信息已无争议；种子、共享模块抽取、舍入规则三处核心分歧必须按“先核验再裁决”的门槛统一，本轮仍不能收敛。

[理由]
1. **revision_20261116/proxy_uncertainty.py 种子**：GLM 的两步门槛（先证明输出未进入文稿/结果表，再统一到 config.SEED）是唯一同时满足“禁止非种子随机数”与“不改变已发表数字”的方案。Nemotron 的“必须改”和 Grok 的“禁止改”都是一刀切，前者可能踩红线，后者可能留下永久例外。Kimi 主张“不改”但未设核验步骤，若该脚本输出已混入正式结果，则“不改”同样会维持一个未走 config.SEED 的种子，与 README 声明冲突。因此，正确裁决是：由主 agent 核验该脚本输出是否进入 results/tables/ 或文稿数字；若未进入，改种子为 config.SEED；若已进入，保留原种子并标注为遗留例外，不得改动。该目录在 ec-metabolic-coloc/ 内，属于审查范围，不能以“不在正式包”为由排除。

2. **MR 自检 seed=1**：该自检块不写结果表，seed=1 是已播种的随机数，不违反“禁止非种子随机数”。将其改为 config.SEED 会破坏测试预言的确定性，且无任何收益；移到 tests/ 增加文件并可能改变导入路径，同样无必要。保留 seed=1 作为单元测试 fixture，仅在 README 声明中承认存在此例外。Kimi 的“两种方案均可”模糊了风险，Mistral 的“移到 tests/ 并改种子”是过度工程，应予反对。

3. **zmat/build_background/draw_module 抽取**：这些函数是跨癌得分与零模型之间的统计契约，即使文本相同，调用上下文（输入矩阵维度、抽样单位）的微小差异也可能导致数值偏移。Grok 的“默认不共享，回归一致才抽”是唯一安全策略。Kimi、Mistral、Nemotron 的“diff 后直接抽”缺少数值验证门槛，可能引入难以察觉的 p 值变化。我们修正为：先逐行 diff，若等价则对现有 results/tables/ 做回归测试，输出完全一致才抽取到 code/common/module_sampling.py；不等价则保留两份并改名以暴露差异。

4. **共享模块组织**：单一 shared_utils.py 将位点逻辑、数值计算、IO、模块抽样混入一个文件，是用新耦合换旧重复。采纳 GLM 的分文件方案：code/common/{locus_utils, numerics, io_utils, module_sampling}.py。Kimi 的命名建议（io_utils, locus_utils, coloc_math, module_sampling）与此一致，无实质分歧。

5. **舍入冲突**：图中 PTS 的 0.098 来自表内 0.0976 的三位舍入，标签 +0.10 来自两位舍入。读表改造时必须锁定规则：图坐标使用三位精度（与现有 0.098 一致），标签使用两位精度并带符号（与现有 +0.31/+0.20/+0.10 一致）。改造后必须用黄金对照复现现有输出，失败即停。Grok 和 GLM 已指出此坑，其他三方遗漏。

6. **skeleton_ebiomedicine.tex 改名**：该文件属 06_EC_IR/ 文稿，明确不在本轮范围。我们撤回此前将其列为 P2 的主张，本轮不动。

7. **README 声明**：原句“唯一常数是随机种子”已被 PTS 硬编码和 MR 自检 seed=1 证伪，必须修正。但“无硬编码 path/accession/threshold”部分尚未核验，不能提前宣称合规。修正顺序：先由主 agent 对全仓做 grep 审计，确认除已定位项外无其他硬编码，再撰写替换句。Kimi 的“先修正已证伪部分”与 GLM 的“后置改写”本质一致，可合并为：先临时注释掉被证伪的句子，待审计完成后替换为准确声明。

[异议]
- **对 Nemotron**：(1) 将 proxy_uncertainty.py 改种子列为无条件 P0，未设核验门槛，可能改变已发表数字，违反硬要求。(2) 主张 _parse_args 抽取，但 CLI 解析器是脚本契约，即使文本相似也不应共享，否则两个 singlecell 脚本的参数集互相牵制。(3) 未处理 0.0976→0.098 的舍入冲突。(4) README 的 liver 值写为 0.931，实际表内为 0.931，但 README 写 0.93，这是四舍五入一致，不是错误。
- **对 Grok**：(1) 对 proxy 种子“禁止改”的表述过于绝对，应改为“核验通过前禁止改”，核验通过后统一到 config.SEED 是应做的，否则种子政策留一个永久例外。(2) 其余无异议，整体最严谨。
- **对 Kimi**：(1) 主张 revision 不改种子，但未区分“核验后决定”，若该脚本输出已进入正式结果，则不改会维持不一致；应接受两步门槛。(2) 主张 MR 自检两种方案均可，但我们认为保留 seed=1 是唯一正确方案，改或移都是无收益变更，应明确反对另两种方案。(3) 图中未写死 61.6% 是事实，但这不是 Kimi 的核心论点，不影响其立场。
- **对 Mistral**：(1) 主张 MR 自检移到 tests/ 并改种子，我们反对，理由如上。(2) 对 proxy 种子的“先核验再裁决”与我们一致，无异议。(3) 对 zmat 的“默认不共享”与我们一致，无异议。
- **对 GLM**：(1) 整体无重大错误，两步门槛、分文件共享模块、舍入黄金对照、README 后置改写均采纳。(2) 唯一商榷：GLM 主张 MR 自检“不移也不改”，我们完全同意；但其表述为“保留为 fixture”，未明确反对 Kimi/Mistral 的改种子方案，我们在此明确反对任何改动。

[修正]
相对上一轮（GLM 转述的 DeepSeek 立场），我们做出以下修正：
1. **proxy_uncertainty.py 种子**：从“P2 改为 config.SEED”修正为“先核验是否进入文稿，再决定改或不改”，优先级 P1（核验）→ P2（修改，若允许）。
2. **共享模块**：从“单一 code/shared_utils.py”修正为“按主题分文件 code/common/{locus_utils, numerics, io_utils, module_sampling}.py”，并增加“zmat 类函数默认不共享，回归一致才抽”的门槛。
3. **skeleton_ebiomedicine.tex**：撤回 P2 重命名主张，本轮不动。
4. **舍入冲突**：新增“锁定舍入规则 + 黄金对照”的硬约束，读表改造必须复现现有标签和坐标。
5. **MR 自检**：明确反对任何改动（改种子或移 tests/），保留 seed=1。
6. **README 声明**：从“直接改写整句”修正为“先审计全仓，再替换声明”。

[收敛] NO
未解决的核心分歧：
1. **proxy_uncertainty.py 种子**：Nemotron 仍持“必须改”（未设门槛），Grok 持“禁止改”（未区分核验后），需全体接受“先核验再裁决”的两步门槛。
2. **MR 自检 seed=1**：Kimi 持“两种方案均可”，Mistral 持“移到 tests/ 并改种子”，需统一到“保留 seed=1，不移不改”。
3. **zmat 类函数抽取门槛**：Kimi、Nemotron 可能仍持“diff 后直接抽”，需接受“回归一致才抽，默认不共享”。
4. **README 声明修改范围**：Kimi 主张先修正已证伪部分，GLM 主张后置改写，实质一致但表述需统一为“先临时注释被证伪句，审计后替换全文”。
以上任一条未获全体明确同意前，不得 `[收敛] YES`。
<<<END>>>

### cn-kimi Kimi (moonshotai/kimi-k2.7-code)

<<<ANSWER>>>
[立场] 修正：接受“默认不共享、等价才抽取”的严格门槛；MR 合成自检保留 seed=1 作为测试 fixture 并隔离；revision_20261116/ 默认不改动；图读表必须带舍入黄金对照；共享模块按主题拆分；本轮仍不收敛。

[理由]
1. 硬要求的优先级应是“不改变已发表数字”高于“种子统一”。因此，任何会改变随机流、函数输出或数值的改动，必须先证明其输出未进入 results/tables/ 或文稿，并通过回归测试。
2. `code/07_figures/02_figure1_design.py` 中的科学数字（rho 与标签）必须从 `results/tables/concordance_by_cohort.tsv` 读入中间变量；版式坐标（74.0/83.0/92.0、15.0/16.6 等）不是科学结果，可保留。读表后必须复现现有图上的 0.098 与标签 “+0.31/+0.20/+0.10”；若不能复现，停止修改。舍入规则建议锁定为：y 坐标显示三位小数（`{rho:.3f}`），标签显示两位带符号（`{rho:+.2f}`）。
3. 身份信息清除范围：删除所有正式 `.py/.R` 文件头的作者姓名、邮箱；R 中 “journal requirement” 改为中性英文；不引入新的团体名。`06_EC_IR/paper_ms/skeleton_ebiomedicine.tex` 不在正式包范围，本轮不动。
4. 共享模块抽取须设双重门槛：先逐行 diff，再对现有 `results/tables/` 做数值回归（容差 0）。建议按主题拆分为 `code/common/{locus_utils.py, numerics.py, io_utils.py, module_sampling.py}`，禁止单一 `shared_utils.py` 混装异质函数。`zmat`/`build_background`/`draw_module` 默认保留两份，不等价则改名以暴露差异；`_logsum`、`_find_member`/`genes_near_loci`、`resolve_input`/`remote_size` 同样处理。`_parse_args` 不抽取。
5. 随机数三类裁决：统计推断程序（bootstrap/置换/模块抽样）保留；`10_positive_control.py` 保留但输出隔离、标注 simulated、不扩展、不填缺失数据；`01_mr_estimators.py` 的 `__main__` 合成自检是单元测试 fixture，保留 `seed=1`、`TRUE=0.30`、`n=400`，仅移到 `tests/` 或显式标注 synthetic self-test，不写 results，不改种子；`revision_20261116/` 不在正式包，隔离并审计 provenance，默认不改种子，若用户后续纳入再按两步门槛处理。
6. README 的 “only constant fixed in source is the random seed” 声明已被证伪，必须重写；但重写须在全仓 path/accession/threshold 硬编码审计完成后进行，不能提前写另一句绝对声明。
7. 重跑计划分三阶段：A. 不依赖大文件的自检/单元测试；B. 依赖本地数据的阶段，先跑输入清单预检，缺文件即停并报告；C. 禁止用模拟数据补洞。完整重跑受原始数据可用性约束，计划必须明确这一点。
8. 核对协议：建立机器可读对照表（claim_id / 文本位置 / 声称值 / 舍入规则 / 源 tsv 文件列 / 状态）。`README` 中 95% CI（−0.014 to +0.031）标 `source_missing`，定位来源文件前禁止重跑 `03_burden/07_incremental_value.py` 等 bootstrap 脚本，以免覆盖现有结果。`config.py` 与 `config.R` 的 SEED 一致性未核验，须加审计。

[异议]
- 对 DeepSeek：单一 `code/shared_utils.py` 把位点逻辑、数值、IO、抽样混装，是用新耦合换旧重复；未处理 0.0976 → 0.098 → +0.10 的舍入冲突；提议用 “BioHex Consortium” 替换作者行是引入新身份信息；将 `skeleton_ebiomedicine.tex` 改名超出本轮范围；对 `revision_20261116/` 改种子未设“是否进入文稿”的隔离门槛。
- 对 GLM：基本赞同其“两步门槛”和主题分文件方案；但“revision 证明隔离后改种子”若写入本轮计划，应列为 P2 且仅在用户明确扩大范围后执行，因为该目录不在正式可复现包内。
- 对 Mistral：主张将 MR 自检 seed=1 改为 `config.SEED` 并移入 `tests/`；我们认为 seed=1 是测试 fixture 的预言常量，改种子无收益且破坏测试确定性，移动入口可接受但不得改种子。
- 对 Grok：赞同其等价性测试、舍入黄金对照、README 后置改写等严格立场；但“分阶段重跑是请求用户改写硬要求 3”过于绝对——在原始数据不在 git 的前提下，唯一诚实的执行方案是预检输入+缺数据停报，这仍是执行而非削弱“禁止模拟数据”的硬要求。
- 对 Nemotron：将 `revision_20261116/proxy_uncertainty.py` 改种子列为 P0 不可接受，因为该脚本范围与 provenance 均未核验，改种子可能改变已发表数字；未处理图内舍入冲突；对 README liver 0.93 与表 0.931 的质疑不成立，二者是四舍五入关系。

[修正]
- 撤回“MR 自检 `seed=1` 应改为 `config.SEED`”的上一轮立场，改为保留 seed=1 作为 fixture，仅隔离/标注。
- 撤回“同名函数 diff 后即可抽取共享”的立场，改为“默认不共享，diff + 数值回归一致才抽取，否则改名保留两份”。
- 明确接受“图读表必须带黄金对照、先锁定舍入规则”。
- 明确 `revision_20261116/` 不在本轮正式修改范围，默认隔离审计，不改动种子。
- 将共享模块组织从单一 `_shared/` 改为按主题 `code/common/`。

[收敛] NO
P0 卫生项（硬编码结果数字读表、身份信息清除、README 声明修正）已趋一致，但 revision 种子处置、MR 自检种子、zmat 类函数共享门槛三处，其他成员上一轮仍有“必改/禁改/默认不共享”等对立表述；我方本轮修正后尚未获得他们确认，故不能判定实际到场全体成员已达成一致。
<<<END>>>

### intl-mistral Mistral Large 4 (mistralai/mistral-large-4-0)

<<<ANSWER>>>
[立场] 修正 —— 赞同五方已形成的共同框架（不合并 FLOW 阶段、保留全部带种子的统计推断、PTS 改为读表、清除身份信息、分阶段重跑、本轮不改仓库）；我在四个未决点上给出统一裁决规则：MR 自检种子保留为 fixture、revision 目录本轮零改动、共享模块抽取须过“文本 diff + 输出回归 + RNG 调用序列”三重门槛、README 改写后置于全量审计。以上未获全体接受前不收敛。

[理由]
1. **PTS 硬编码（P0）**：`02_figure1_design.py` L169–171 的科学数字必须经中间变量来自 `results/tables/concordance_by_cohort.tsv`，且连接键必须显式——表行是 GSE accession（GSE15773/GSE18732/GSE13070），图上标签是 adipose/diabetes/clamp，二者映射关系要写在代码或小配置里，读表头校验、连接失败即停，禁止换表。**舍入锁定**（此点多方遗漏，只有 xai 指出）：表值 0.0976 在现图内画的是 0.098（三位舍入，作图坐标）、标签是 "+0.10"（两位带符号）。改造后必须“作图坐标用 `f"{rho:.3f}"`、标签用 `f"{rho:+.2f}"`”，并以黄金对照断言生成标签严格等于 "+0.31"/"+0.20"/"+0.10"、坐标严格等于 0.313/0.202/0.098；任一不符则停止、不发布，不得强行“从表读取”后放行一个变了的图。版式坐标 74.0/83.0/92.0 与 15.0/16.6 是排版常量不是科学结果，可保留。
2. **随机数三类裁决**：(a) 7 个统计脚本（bootstrap/置换/模块抽样）全部保留，种子来自 `config.SEED` 或 `--seed`——共识；(b) `10_positive_control.py` 保留为阳性对照/功效模拟，但其输出路径是否已在 `results/tables/` 之外属**未核验**，须先核验再决定是否迁移；若迁移，输出必须隔离并标注 simulated；不新增任何模拟路径；(c) `01_mr_estimators.py` 的 `__main__` 自检：`default_rng(1)` 是**有种子**的，不违反“禁止非种子随机数”；`TRUE=0.30`、`n=400` 是测试 fixture 而非分析常量。裁决：种子保持 1 不改（改成 `config.SEED` 会改变抽样，若断言绑定该次抽样则直接破坏自检的可判定性），把自检块移到 `tests/test_mr_estimators.py` 并在 README 注明“该测试使用固定 fixture 种子，不属于分析管线”。`revision_20261116/` 明确不在本轮范围：零代码改动，仅在计划中注明“遗留修订目录，未纳入正式 FLOW”；`proxy_uncertainty.py` 的种子统一留待用户未来决定是否归档、且先证明其输出未进入文稿（两步门槛），本轮不得以任何优先级改动它。
3. **去重而不合并**：禁止合并 FLOW 阶段脚本、禁止把 45 个脚本合成巨文件。真重复函数抽到按主题划分的共享模块：`code/common/locus_utils.py`（`_find_member`, `genes_near_loci`）、`code/common/numerics.py`（`_logsum`）、`code/common/io_utils.py`（`resolve_input`, `remote_size`）、`code/common/module_sampling.py`（`zmat`, `build_background`, `draw_module`）。抽取门槛三重：文本 diff 相同 + 受影响脚本对现有 `results/tables/` 的输出逐值回归一致（容差 0）+ 随机数调用序列一致（同种子下调用次数/顺序不同即不同输出）。`zmat`/`build_background`/`draw_module` 是零模型与得分的统计契约函数，**默认不共享**；回归不通过就保留两份并改名以暴露差异。`_parse_args` 不抽取。诚实表述硬要求 1：抽取会净增约 4 个文件，文件数的“精简”来自消除重复而非合并，净增必须由用户批准，不能写成“已满足精简”。
4. **身份信息（P0）**：所有正式脚本文件头的 `Author: Tao Zhu (zhutao@tjh.tjmu.edu.cn)` **删除**而非替换为团体名（替换仍是身份信息）；R 文件中 “journal requirement for composite figures” 改为中性英文描述；`06_EC_IR/paper_ms/skeleton_ebiomedicine.tex` 是文稿、不在本轮范围，文件名不动——单独向用户表态。
5. **README 声明**：原句后半句 “The only constant fixed in source is the random seed” 已被 PTS 硬编码和自检 fixture 种子证伪；前半句 “无硬编码 path/accession/threshold” **未核验**。因此 README 改写列为 P1，必须在主 agent 对 `code/**` 完成 path/accession/threshold/科学数字的全量 grep 审计并出具报告之后进行，新措辞须可验证（如“分析路径内所有随机性经 config.SEED 播种；分析代码中不嵌入结果数值”），最终措辞待审计后定稿。另须审计 `config.py` 的 SEED 与 `config.R` 的 SEED 是否一致（**未核验**）。
6. **23.28% 注释**：`02_select_instruments.py` docstring 中的结果数字只删除，不改指向未定位的表。
7. **核对协议**：建立“文稿/README 句子 → 具体 tsv/tex 单元格”对照表，读单元格原文而非记忆值。已一致：concordance rho 0.313/0.202/0.0976 与 README +0.31/+0.20/+0.10（四舍五入一致）；incremental_value.tsv 的 delta_auc=0.000389、p=0.15、n=145、events=45 与 README +0.0004/p=0.15 一致；table1_coloc_summary.tex 的 503 tests / 3 个 PP.H4>0.8 / 最大 0.962（muscle）、liver 0.931 与 README 0.93 一致。**未核验**：README 中 ΔAUC 的 95% CI −0.014 to +0.031 不在 incremental_value.tsv 该行，来源文件待定位——在定位并验证之前，**禁止重跑** `03_burden/07_incremental_value.py` 及任何 bootstrap 脚本，以免覆盖现有结果表。
8. **分阶段重跑计划**：阶段 A——不依赖大文件的自检（import、单元测试、config 一致性）；阶段 B——依赖本地数据的阶段，缺文件即停止并输出缺失清单，禁止用模拟数据补洞；阶段 C——禁止新增模拟数据路径。重跑检查单纳入独立实验单位、伪重复、批次/组成混杂、多重检验、效应量与 CI——作为审查项，不预设为已证实缺陷。

[异议]
- **对 DeepSeek**：抽取后输出逐值比对容差 0 正确，但须补充 RNG 调用序列一致性门槛；把 revision 目录的种子修改列入计划项与“本轮范围不含 revision_20261116”自相矛盾；若曾主张用团体名替换作者行，那是引入新身份信息，应删除而非替换。
- **对 GLM**：两步门槛裁决、主题分文件、舍入黄金对照均采纳；但 skeleton_ebiomedicine.tex 改名列入 P2 越出本轮范围，反对；MR 自检“seed=1 是 fixture、不改种子”的论证正确，采纳。
- **对 Kimi**：同意 README 修改范围限于已证伪部分、审计后定稿；但“自检改用 config.SEED”作为可选方案应放弃——fixture 种子改动可能破坏绑定该抽样的断言，计划里只保留“移到 tests/ 且保持 seed=1”；同意 revision 目录不动代码。
- **对 xai (Grok)**：整体最严谨，采纳其黄金对照舍入锁、CI 不预设文件名、config 跨语言种子未核验、硬要求 1 不能宣称满足、版式坐标不标未核验。一点补充：仅说“从图读表”不足以保证图面不变，必须明确“作图坐标三位舍入、标签两位舍入”的双规则。
- **对 Nemotron**：把 `proxy_uncertainty.py` 改种子列为无条件 P0 逻辑不成立——该目录不在本轮范围，且未证明其输出未进文稿，改种子可能改变已发表数字；反对其 `_parse_args` 抽取到 cli.py。

[修正] 无修改（本轮为我方首次完整表态）。相对五方共同共识，我方新增并坚持的裁决是：MR 自检 seed=1 保留为 fixture 不改种子；PTS 改造附“三位坐标/两位标签”双舍入黄金对照；共享模块抽取增加 RNG 调用序列一致性门槛，`zmat` 三函数默认不共享；硬要求 1 须诚实表述为“净增约 4 个共享文件，待用户批准”。

[收敛] NO。仍有实质分歧未获全体接受：(1) `proxy_uncertainty.py` 种子——DeepSeek（P2 改）、Nemotron（P0 改）vs 多数（本轮不改、两步门槛）；(2) MR 自检种子——Kimi 仍保留“改用 config.SEED”选项 vs 我方与 GLM 的“保留 seed=1”；(3) `zmat`/`build_background`/`draw_module` 抽取门槛表述（diff 即抽 vs 回归+RNG 序列一致才抽、默认不抽）；(4) README 改写时机与范围。任一条未全体同意前不得 `[收敛] YES`。
<<<END>>>

### intl-xai Grok 4.7 (x-ai/grok-4.7)

<<<ANSWER>>>
[立场] 修正。卫生修复可以进入待批准计划，但种子一律不改写、统计函数默认不共享、重跑不得覆盖现有表；本轮不改仓库，核心分歧未消除，不收敛。

[理由]
待批准计划只做不改变已发表数字的卫生与审计。执行门槛：用户批准前不改仓库；任何数值回归对不上、或缺本地数据因而无法回归，都是失败并停止，不是放行。

P0（预期不改变已发表数字；黄金对照失败则该项不落地）
- `code/07_figures/02_figure1_design.py`：只把科学字段从 `results/tables/concordance_by_cohort.tsv` 读入中间变量。连接键用已核实对应关系（adipose/GSE15773、T2D/GSE18732、clamp/GSE13070），禁止按行序猜测，禁止改从其他表读取。74.0/83.0/92.0 与 15.0/16.6 保持版式常量，不从结果表派生。锁定格式：坐标 `f"{rho:.3f}"` 必须得到 0.313、0.202、0.098；标签 `f"{rho:+.2f}"` 必须得到 +0.31、+0.20、+0.10。禁止把原始 0.0976 直接画进图。任一字符串或坐标不一致即停。
- 正式路径 `code/{00_setup..07_figures}`、`config/config.py`、`config/config.R`、`07_figures/01_main_figures_1_2_3.R`：删除作者姓名与邮箱；期刊措辞改为中性英文。不改统计逻辑。
- `05_mendelian_randomisation/02_select_instruments.py`：删除 docstring 中的 “23.28%”，改为指向结果表，不写入新数字。
- 重跑协议写入计划而非本轮执行：不覆盖 `results/tables/`；输出到独立目录再 diff；diff 不一致只报告，不回写、不改统计代码去凑文稿。缺原始数据即停并列出缺失，禁止用模拟数据补洞。在定位到 ΔAUC 95% CI（−0.014 to +0.031）的真实来源文件之前，禁止重跑 `03_burden/07_incremental_value.py` 及任何 bootstrap/置换脚本。该 CI 现状为 source_missing，不得假定正确。

P1（仍不改已发表数字）
- 先对正式 `code/**` 检索 path、accession、threshold 和嵌在源码里的科学结果数字，再整句替换 README 的绝对声明。禁止先写“已合规”。替换方向：正式结果路径的种子来自 `config.SEED` 或 `--seed`；科学结果从中间表读入；自检预言和遗留修订种子不是该声明的例外补丁，而是声明必须承认的非结果常量。声明应改，因为现句已被 PTS 硬编码和非 `config.SEED` 种子证伪。
- 同名函数先逐行 diff。文本相同才允许抽到 `code/common/io_utils.py`（`resolve_input`、`remote_size`）、`code/common/locus_utils.py`（`_find_member`、`genes_near_loci`）、`code/common/coloc_math.py`（`_logsum`）。`zmat`、`build_background`、`draw_module` 默认不共享；仅当文本相同且受影响现表逐值一致（容差 0）才可新建 `code/common/module_sampling.py`。缺数据无法回归则不抽取。不等价则保留两份并改名，禁止强行统一。`_parse_args` 不合并。阶段脚本 00–07 不合并，不合成巨文件。精简度量是重复行，不是文件数。
- 核对协议由主 agent 在批准后执行，本轮不声称全文已核对。机器表列：claim_id、文本位置、声称值、舍入规则、源文件与列、状态（match / rounding_match / mismatch / source_missing / 未核验）。可预填的仅限已核实项：concordance 三 rho 与 61.62%/61.6%；incremental 的 0.000389/+0.0004、p=0.15、n=145、events=45；coloc 503、PP.H4>0.8 共 3、liver Max 0.931 对 README 0.93、muscle 最大 0.962。其余标未核验。
- `10_positive_control.py`：保留为已有、仅模拟效应量的非结果脚本；不扩展、不填缺失、不进文稿来源表。输出路径未核验，先核验；若已混入正式表，只隔离并标注 simulated，不重新生成。保留须在批准时单独征得同意。不新增模拟路径。
- `config.py` 与 `config.R` 的 SEED 是否相同未核验。只核对；不一致则报告并停，禁止自动统一。
- `revision_20261116/`：在审查范围内，不在正式 FLOW 重跑范围内。必须核验三份脚本输出是否进入 `results/tables/` 或文稿数字。无论结果如何，本计划都不把 `default_rng(20261116)` 改成 `config.SEED`。20261116 已是种子，不违反“禁止非种子随机数”；改种子会破坏该脚本原输出的可复现性。若已进入文稿，改种子永久禁止。若未进入，仍不改种子，只排除出正式结果路径并记录遗留。

P2
- 修订脚本中文注释可另列注释级英译，且 diff 不得触碰随机数那一行；未获准碰该目录则不动。
- `01_mr_estimators.py` 的 `__main__` 保留 `seed=1`、`TRUE=0.30`、`n=400`，加英文横幅标明合成自检、不写 results、不进入 FLOW。不强制迁入 `tests/`（仓库是否已有测试入口未核验）。禁止把 `TRUE=0.30` 改从结果表读取。
- `06_EC_IR/paper_ms/skeleton_ebiomedicine.tex` 本轮不改名。

反对：合并 FLOW 阶段或把约 45 个脚本合成巨文件；删除带种子的 bootstrap/置换；从错误的表读示意图数字；用原始 0.0976 替换图上 0.098；把 revision 种子改写列为 P0 或“隔离后必须统一”；把 MR 自检种子改成 `config.SEED`；为跑通而新增模拟数据；重跑时覆盖现表；SEED 不一致时自动选边；未过双重门槛就抽取 `zmat` 一类函数；用改代码去挽救 README 原句。

风险：读表改造若格式未锁定会改图；共享统计函数若调用次数或全局状态不同会改 p 值；重跑覆盖会毁掉文稿对照基线。三者都用失败即停控制，而不是用“预期一致”放过。

[异议]
- 对 DeepSeek：容差 0 和缺 CI 来源前禁止重跑增量脚本，采纳。强制迁入 `tests/test_mr_estimators.py` 不采纳：测试入口是否存在未核验，隔离不等于新建目录。缺数据时的回归必须算失败，不能只写“有差异才回退”。
- 对 GLM：两步门槛的核验步采纳；“证明隔离后必须改成 `config.SEED`”不采纳。这是溯源破坏，不是用户规则。用户禁止的是非种子随机数，不是遗留种子与 `config.SEED` 不同。主题分文件和默认不共享 `zmat`，采纳。MR 种子维持 1，采纳；强制迁入 `tests/` 不采纳。
- 对 Kimi：不把 revision 种子改掉，方向对；但“目录完全不出现在计划里”漏掉了用户要求的裁决。自检“两种方案交主 agent 选”不采纳，计划必须单选。README 不能只改半句后继续留着可能仍为假的绝对声明，应整句后置到检索之后。阳性对照输出路径未核验，采纳。最新文本仍未明确禁止把原始 rho 当作图坐标。
- 对 Mistral：图内黄金对照、默认不共享 `zmat`、README 后置，采纳。“自检改 `config.SEED` 零风险”不成立：断言内容未核验，`config.SEED` 的具体值也未核验。修订脚本“未进文稿就改种子”与 GLM 同类错误，不采纳。
- 对 Nemotron：先核验 revision 输出是否进文稿，这点对；核验不是改种子的授权。把重复函数抽取写成已共识、而不写等价门槛，对 `zmat` 过宽。把“计划中包含一次未来的文件移动”说成超出“本轮不重跑”，混淆了计划与执行。README 冲突应改声明，不应改自检预言去凑原句。0.0976/0.098 已由多方覆盖，该条异议已过时。

[修正]
相对上一轮被引用的“未证明隔离前禁止改种子”：改为无论隔离结果如何都不改 20261116，核验只决定如何标注遗留，不授权改种子。补充容差 0、缺数据即不抽取、精确格式串、禁止画原始 0.0976、重跑不覆盖现表、SEED 不一致不自动统一、自检不强制迁目录。

[收敛] NO
未决：revision 种子在确认未进文稿后是否仍可改（我方：不可）；MR 的 `seed=1` 是否改为 `config.SEED`（我方：不改）；自检是否必须迁入 `tests/`（我方：不必须）。这三项未获全体明确同意前不得收敛，也不得把未决项写成已经批准的改法。
<<<END>>>

### intl-nvidia Nemotron 3.5 Super (nvidia/nemotron-3.5-super-vl)

<<<ANSWER>>
[立场] 修正 —— 撤回先前对 `proxy_uncertainty.py` 无条件 P0 改种子、`_parse_args` 抽取及 README liver 值表述；采纳“先核验隔离、再统一种子”的两步门槛，默认不动 revision 代码；共享模块默认不抽，须 diff + 零容差数值回归；MR 自检保留 seed=1 作为 fixture。

[理由]
1. `proxy_uncertainty.py` 位于 `revision_20261116/`，不在正式可复现包范围；未证明其输出未进入文稿/结果表前改种子即可能改变已发表数字，违反硬要求。本轮只能将其标记为遗留、隔离、不改代码，待主 agent 核验后由用户裁决是否统一到 `config.SEED`。
2. `_parse_args` 是脚本契约，同名不代表应合并；抽取会引入循环导入和参数耦合，无收益，应保持两份独立解析器。
3. fact G 明确：README 写 liver PP.H4=0.93，表中 liver Max PP.H4=0.931，四舍五入一致；0.962 是 muscle 最大值。若我先前表述为 README=0.931 是错误，应更正。
4. 共享模块抽取必须以“相同输入→相同输出（容差为 0）”的回归测试为前提，`zmat`/`build_background`/`draw_module` 涉及零模型与得分，默认保留两份。
5. 图内 PTS 数字必须从 `results/tables/concordance_by_cohort.tsv` 经中间变量读入，并锁定舍入规则：坐标三位（0.098）、标签两位带符号（+0.10），黄金对照失败则停止。

[异议]
- DeepSeek：(1) 将 `proxy_uncertainty.py` 改种子列为 P2，未设“先证明输出未进文稿”的门槛；(2) 提议单一 `shared_utils.py` 装下 8 个异质函数是反模式，应按主题分文件；(3) 提议重命名 `skeleton_ebiomedicine.tex` 越出本轮范围；(4) 未处理 0.0976（表）与 0.098（代码）的舍入冲突，未说明种子来源。
- GLM：(1) 两步门槛原则正确，但未明确写出“本轮不修改 `revision_20261116/` 任何代码”；(2) 对 MR 自检 seed=1 的表述偏向保留，但未明确拒绝将其改为 `config.SEED` 作为 P0；(3) 未明确要求 `10_positive_control.py` 输出须先核验为不在 `results/tables/` 才能标记为隔离。
- Kimi：(1) 对 MR 自检给出“两种方案均可”但未在计划中二选一，留有模糊；(2) 未明确说明 `revision_20261116/` 本轮不改代码，仅说“不在范围”；(3) 未处理 0.0976/0.098 单元格原文差异。
- Mistral：(1) 倾向将 MR 自检种子改为 `config.SEED`，但该块不写结果、改种子无收益且破坏测试预言的可判定性，应保留 seed=1；(2) “先核验再裁决” revision 种子是可以接受的折中，但应补充“本轮无论核验结果如何均不修改其代码”；(3) 未要求按单元格原文核对 0.0976 与 0.098，未禁止假定 CI 来源文件名为 `incremental_boot.tsv`。
- Grok：(1) 将“图中写死 61.6%”归咎于 Kimi，但 Kimi 立场中未提此项，属误读；(2) 说 Kimi 要求作图脚本额外读入 `n_responsive`，Kimi 亦未提出，属误读；(3) 正确指出共同遗漏（单元格原文、CI 文件名假设、config.SEED 跨语言一致性），但未给出具体解决动作。

[修正] 相对上一轮：
1. 撤回 `proxy_uncertainty.py` 无条件 P0 改种子；改为：本轮不改 `revision_20261116/` 任何代码，仅标记为遗留并隔离，待主 agent 核验其输出是否进入正式结果后，由用户裁决是否统一到 `config.SEED`。
2. 撤回 `_parse_args` 抽取至 `cli.py`；同意保持两份独立解析器，不合并。
3. 更正 README liver 值表述：README 为 0.93，表为 0.931，四舍五入一致；不得说成 README=0.931。
4. 采纳共享模块默认不抽、须 diff + 零容差回归通过才抽，`zmat`/`build_background`/`draw_module` 默认保留两份。
5. 采纳 MR 自检保留 `seed=1`，可移至 `tests/` 但不改种子，不作为 P0。
6. 采纳 README 修改后置于全量硬编码审计后，仅先修正已证伪的“唯一常量是随机种子”子句。
7. 采纳图内舍入锁定：坐标 `f"{rho:.3f}"`、标签 `f"{rho:+.2f}"`，黄金对照失败即停。

[收敛] NO。虽大部分核心分歧已趋于一致（revision 种子、README 时机、舍入锁定、共享模块门槛、`_parse_args` 不抽），但 MR 自检种子（保留 seed=1 vs 改 `config.SEED`）仍未获得全体明确同意；且共享模块的目录命名与文件拆分方案（单文件 vs 按主题分文件）尚未定稿。未见所有实际到场成员在本轮对整合方案做出无异议确认，故不得写 `[收敛] YES`。
<<<END>>>

# 你的任务

你是 **GLM（z-ai/glm-5.3-prime）**，负责：方法架构 / 中文综合。
请逐条阅读其他方的立场：

1. 明确指出他方论证中的错误、遗漏、不严谨或不可核验之处；没有则写“无异议”。
2. 明确你是否修改自己的立场，以及修改了什么。
3. 不要重复任务背景，不要客套。
4. 只有在你认为**所有实际到场的成员**已就核心分歧达成一致时，才写 `[收敛] YES`；
   有实质未解决分歧时写 `[收敛] NO`。

严格按以下结构输出：

<<<ANSWER>>>
[立场] 支持 / 反对 / 修正（选一个，后接一句话结论）
[理由] 你的核心论据
[异议] 对其他方的具体质疑；无则写“无异议”
[修正] 你相对上一轮修改了什么；无则写“无修改”
[收敛] YES / NO
<<<END>>>
