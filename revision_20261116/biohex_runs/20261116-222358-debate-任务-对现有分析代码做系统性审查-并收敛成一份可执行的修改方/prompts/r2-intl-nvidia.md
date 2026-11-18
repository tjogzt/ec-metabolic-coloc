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
[立场] 修正。正式包的修改方案必须比现有提议更严格：共享模块抽取须经等价性测试且默认不合并统计相关函数；MR 自检保留原样但隔离出分析路径；revision 脚本隔离不动；README 声明须在全面检索后重写；舍入政策须锁定并验证；本轮不收敛，待上述分歧解决后再议。

[理由]
1. **共享模块抽取的风险被低估**：`zmat`、`build_background`、`draw_module` 等函数直接参与零模型构建与得分计算，任何实现差异（浮点运算顺序、边界处理、随机数调用次数）都会在相同种子下产生不同输出，从而改变已发表的 p 值或置信区间。glm、kimi、mistral、nvidia 均提议抽取这些函数，但未要求先通过“相同输入→相同输出”的回归测试。xai 正确指出必须等价测试通过才许动，否则保持两份。我采纳这一立场，并进一步要求：即使 diff 显示文本相同，也必须在抽取后对受影响脚本的输出表做逐值比对（容差为 0），任何差异即回退。这比“先 diff 再抽取”更安全。
2. **MR 自检种子的处理**：glm 和 kimi 提议将 `seed=1` 改为 `config.SEED`，但这会改变自检的预言值（尽管不写结果表），且该自检是门控测试，改种子可能掩盖与原始开发环境的浮点差异。mistral 提议移到 `tests/` 并改种子，同样有风险。xai 提议保留 `seed=1` 与 `TRUE=0.30` 作为预言，仅隔离出分析路径，这是最保守且不改变任何行为的选择。我采纳：将该自检移到 `tests/test_mr_estimators.py`，保持原种子与真值，并在 README 中注明该测试使用固定种子 1，不属于分析管线。
3. **revision 脚本的归属**：nvidia 提议将 `proxy_uncertainty.py` 的种子改为 `config.SEED`，但该脚本可能已被用于生成某些中间结果（未核验）。xai 正确指出在未证明其未进入文稿前禁止改种子。我采纳：`revision_20261116/` 整体隔离，不纳入正式 FLOW，不修改任何代码，仅在计划中标注为“遗留修订，未纳入正式结果路径，待用户裁决是否归档或删除”。
4. **README 声明的修改时机**：glm、kimi、mistral 均提议在 P0 修改 README，但未要求先完成对 path/accession/threshold 的全面检索。xai 强调必须等检索完成后再写，否则可能再次出现声明与事实不符。我采纳：README 修改降为 P1，且必须在主 agent 完成全仓库硬编码检索（包括路径、accession、阈值、版式常量中的科学数字）并出具报告后，再拟定新声明。
5. **舍入政策冲突**：所有他方均未充分讨论 `02_figure1_design.py` 中 0.098（三位）与标签 `+0.10`（两位）的舍入不一致问题。表值为 0.0976，硬编码 0.098 是四舍五入到三位，标签是两位。若简单从表读取并格式化为两位，图上数字将从 0.098 变为 0.10，改变已发表示意图。xai 指出这是冲突。我要求：在 P0 修改该脚本时，必须显式定义舍入规则——图上数字保留三位小数（`f"{rho:.3f}"`），标签保留两位带符号（`f"{rho:+.2f}"`），并验证生成值与硬编码完全一致（0.098 与 +0.10）。若不一致，停止修改并报告，不得强行“从表读取”。
6. **核对协议中 CI 来源**：所有方均指出 95% CI 来源未核验，但未明确禁止在定位前重跑增量价值脚本。xai 正确警告重跑可能改变 CI。我要求在核对协议中明确：在 CI 来源文件定位并验证之前，禁止重跑 `03_burden/07_incremental_value.py` 或任何 bootstrap 相关脚本，以免覆盖现有结果表。
7. **反对做法的一致性**：各方在反对合并阶段脚本、删除 bootstrap、新增模拟数据上一致，无异议。

[异议]
- **对 glm**：未考虑舍入冲突和共享模块的等价性测试要求，对 README 修改时机过于乐观。
- **对 kimi**：同上，且提议将 `proxy_uncertainty.py` 种子改为 `config.SEED` 有风险。
- **对 mistral**：提议将 MR 自检移到 `tests/` 并改种子，可能改变预言值；未强调等价性测试。
- **对 nvidia**：提议抽取 `_parse_args` 到 `cli.py` 不必要，且可能引入循环导入；提议改 revision 种子有风险。
- **对 xai**：无实质异议，但其“黄金对照失败则停”的表述需明确为“生成值与硬编码不一致则停，不得发布”，我已在立场中细化。

[修正] 无修改（本轮首次发言）。

[收敛] NO。核心分歧未解决：共享模块抽取的等价性测试标准、MR 自检种子是否保留、revision 脚本是否改种子、README 修改时机、舍入政策锁定。需各方就上述五点达成一致后方可收敛。
<<<END>>>

### cn-glm GLM (z-ai/glm-5.3-prime)

<<<ANSWER>>>
[立场] 修正 —— 五方在“不合并 FLOW 阶段、保留全部带种子的统计推断、硬编码结果数字必须改为读表、身份信息必须清除、本轮不改仓库”上已实质一致；但三处直接冲突（revision 种子、MR 自检种子、zmat 抽取门槛）尚未被全体接受统一裁决规则，收敛条件不满足。

[理由]
1. **种子裁决的正确标准是“是否进入已发表数字”，不是“是否统一到 config.SEED”**。用户硬规则同时要求“禁止非种子随机数”与“不得改变已发表数字”。`revision_20261116/proxy_uncertainty.py` 的 `default_rng(20261116)`：在未证明其输出未进入文稿/结果表之前，改种子本身就是一次结果变更，必须禁止；证明隔离之后，统一到 `config.SEED` 才是安全操作。这是一个两步门槛，不是“改”或“不改”的单选题。`01_mr_estimators.py` 的 `seed=1`、`TRUE=0.30`、`n=400` 是单元自检的测试预言（fixture）：用户规则禁止的是**非种子**随机数，seed=1 是种子，不违规；把它移到 `tests/` 入口是合理的，但把种子改成 `config.SEED` 没有必要，反而破坏“自检在固定预言下通过”的可判定性。测试 fixture 里的常数不属于 README 所说的“分析源码中固定常数”。
2. **共享模块抽取必须设双重门槛**：文本 diff 相同 + 对现有 `results/tables/` 回归输出完全一致。`zmat`/`build_background`/`draw_module` 是跨癌得分脚本与零模型脚本之间的契约函数，两份实现即使文本相同，调用上下文（输入矩阵、抽样单位）也可能不同；默认不共享，仅当回归一致才允许抽到共享模块，否则保留两份并改名以暴露差异。单一 `shared_utils.py` 把位点逻辑、数值、IO、模块抽样混进一个文件，是用新耦合换旧重复，应改为按主题分文件。
3. **舍入冲突必须显式处理**：clamp rho=0.0976 在 README 是 +0.10（两位）、在图内是 0.098（三位）。把 PTS 改为读表时，必须先锁定标签生成规则（“+0.XX” 两位舍入），并用黄金对照复现现有字符串 “+0.31/+0.20/+0.10”；复现失败即停止，不得改图。否则“读表改造”本身就会改图上数字——这是多数方（DeepSeek、Kimi、Mistral、Nemotron）都漏掉的坑，只有 Grok 指出。
4. **README 声明句的替换必须后置于全量审计**：先对 `code/**` 做 path/accession/threshold 的 grep 审计，确认除已定位项外无其他硬编码，再写替换句。先写“已合规”后补审计，顺序颠倒。
5. 核对协议、CI 标记 source_missing、分阶段重跑（A 自检 / B 缺文件即停 / C 禁止模拟补洞）、`10_positive_control.py` 保留但输出隔离并标注 simulated、`skeleton_ebiomedicine.tex` 本轮不动——这五点五方已一致，直接进入计划。

[异议]
- **对 DeepSeek**：(a) 把 `proxy_uncertainty.py` 改种子列为 P2，未设“先证明输出未进文稿”的门槛，可能改变已发表数字，违反硬要求；且 `revision_20261116/` 明确不在本轮合并/重跑范围，不应出现在正式修改项里。(b) 单一 `code/shared_utils.py` 装下 8 个异质函数是反模式。(c) P2 重命名 `skeleton_ebiomedicine.tex` 超出本轮范围（任务明说 06_EC_IR 不在范围）。(d) 未处理 0.0976→+0.10 vs 0.098 的舍入冲突，未提种子 provenance。
- **对 Kimi**：(a) P1-2 要求把 MR 自检 seed=1 改为 `config.SEED`，混淆了“分析种子政策”与“测试 fixture 常数”；seed=1 不违反用户规则，改它是无收益的变更。(b) `_shared/` 下划线前缀目录名不利于作为包导入，建议 `code/common/`（命名小事，但既然要定就定一个可 import 的）。(c) 同样未处理舍入冲突。
- **对 Mistral**：(a) P0-3 “seed 不再写死为 1”与 Kimi 同错。(b) 完全未回应 `revision_20261116/` 的种子 20261116 问题——这是任务给定事实，遗漏。(c) 其余（manifest、positive_control 隔离、三阶段重跑、回归比对）是本轮最完整的工程方案，我采纳。
- **对 Grok**：整体最严谨，采纳其门槛化裁决、舍入黄金对照、README 后置改写、zmat 默认不共享。两点商榷：(a) “revision 禁止改种子”应表述为“门槛未过前禁止”，门槛过了（证明隔离）之后统一到 config.SEED 是应做的，不是永久禁止——否则种子政策留一个永久例外。(b) `_parse_args` 完全禁止合并我同意，但无需在计划里反复强调，一次说清即可。
- **对 Nemotron**：(a) 最严重的错误：把“`proxy_uncertainty.py` 改用 `config.SEED`”列为 **P0**，与 Grok 的“禁止改种子”直接冲突，且未证明该脚本输出未进文稿就改种子，可能直接改变已发表数字——P0 是“必须先做”的级别，把一个可能改结果的操作放在 P0 违反本轮“不改变已发表数字”的底线。(b) `_parse_args` “逻辑一致可抽至 cli.py”：CLI 解析器是脚本契约，即使逻辑相似也不值得抽，抽了反而让两个 singlecell 脚本的参数集互相牵制。(c) 未处理舍入冲突、无种子 manifest。

[修正]
本轮为首次表态，无上一轮立场。相对五方现有共识文本，我施加的修正是：(1) 种子问题从“改/不改”改为两步门槛裁决（先证明隔离，再统一）；(2) MR 自检种子维持 seed=1 作为 fixture，仅移入口径为 tests/，不改种子；(3) 共享模块按主题分文件（建议 `code/common/{locus_utils, numerics, io_utils, module_sampling}.py`），抽取以 diff+回归一致为双重门槛，`zmat`/`build_background`/`draw_module` 默认不共享；(4) 读表改造必须附标签黄金对照，锁定舍入规则；(5) README 替换句后置于全量硬编码审计。

[收敛] NO
未解决的核心分歧（需下一轮逐条表决）：
1. `proxy_uncertainty.py` 种子：Nemotron（P0 改）/ DeepSeek（P2 改）vs Grok（禁止改）——待全体接受“两步门槛”裁决。
2. MR 自检 seed=1：Kimi、Mistral（改为 config.SEED）vs Grok（保留为预言）——待表决。
3. `zmat`/`build_background`/`draw_module` 抽取：Kimi、Mistral、Nemotron、DeepSeek（diff 后抽）vs Grok（回归一致才抽，默认不抽）——待表决。
4. `skeleton_ebiomedicine.tex` 改名：DeepSeek 单方主张，其余四方反对——待 DeepSeek 撤回或说明理由。
5. 共享模块文件组织（单文件 vs 主题分文件）与目录命名——待统一。
以上任一条未获全体明确同意前，不得 `[收敛] YES`，也不得向用户提交含未决项的修改计划。
<<<END>>>

### cn-kimi Kimi (moonshotai/kimi-k2.7-code)

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

### intl-mistral Mistral Large 4 (mistralai/mistral-large-4-0)

<<<ANSWER>>
[立场] 修正 — P0 卫生修复项（硬编码结果数字移出、身份信息清除、README 声明修正）四方已实质一致；但在 revision_20261116 种子处置、zmat 类函数是否默认共享、图内舍入规则三处仍有实质分歧，本轮未达收敛。

[理由]
1. **P0 三项无争议，直接进入计划**：`02_figure1_design.py` 的 PTS 三元组必须改为从 `results/tables/concordance_by_cohort.tsv` 读入中间变量再绘图（反对从任何其他表读取）；全部正式脚本文件头删除作者姓名/邮箱；README 的 "only constant fixed in source is the random seed" 声明必须修正——采纳 Grok 的补充：替换句须等全仓 path/accession/threshold 检索完成后再落笔，不得提前宣称合规。
2. **随机数三类裁决可统一**：统计程序（bootstrap/置换/模块抽样，种子走 `config.SEED` 或 `--seed`）保留；`10_positive_control.py` 保留但隔离——输出不得进入文稿来源表、不得扩展、不得填缺失数据，不新增任何模拟路径；`01_mr_estimators.py` 的 `__main__` 自检移到测试入口，显式标注 synthetic self-test，不写 results，种子统一走 `config.SEED`（该块不产出结果表，改种子零风险）。
3. **revision_20261116 是最大分歧点**：Nemotron 主张"必须改种子"，Grok 主张"禁止改种子、先核验是否进入文稿"。我的统一方案：先由主 agent 核验三份修订脚本的输出是否进入 `results/tables/` 或文稿数字——若未进入，种子改走 `config.SEED`、中文注释译英；若已进入，保留原种子并标注遗留，**不得改**（改种子会改变已发表数字，踩红线）。一刀切的"必改"或"必不改"都不可接受。
4. **共享模块必须设门槛**：采纳 Grok 的做法——所有同名函数先逐行 diff + 对现有 results 表做数值回归，通过才抽取；`zmat`/`build_background`/`draw_module` 默认不共享，不等价则改名保留两份（否则可能改变 p 值）；`_parse_args` 不合并。采纳 GLM 的度量说明：向用户明确"精简"的度量是重复代码行数而非文件数，阶段编号 00–07 绝不合并。
5. **图读表须带黄金对照**：采纳 Grok 的补充——读表后的加载器必须复现现有 "+0.31/+0.20/+0.10" 标签和 0.098 三位显示；0.0976 → 图坐标 0.098（三位）vs 标签 +0.10（两位）是舍入政策，须先锁定舍入规则，对照失败则停止、不得发布改动。
6. **证据缺口列入核对协议**：README 的 95% CI（−0.014 to +0.031）标 `source_missing`，检索 `results/tables/` 定位真实来源文件，找不到保持未核验；`config.py` 与 `config.R` 的 SEED 一致性加断言（GLM 提出，未核验）；重抽样独立实验单位（patient/library/cell）做审计登记，本轮不改统计逻辑（Grok/Nemotron 均指出，属证据缺口）。

[异议]
- 对 Nemotron：`proxy_uncertainty.py` "必须改为 `config.SEED`" 有风险——若该脚本输出已进入已发表表，改种子即改变已发表数字；应先核验再裁决，不能"必须改"。另，"阶段三完整重跑"易被误读为本轮执行，任务明确本轮只定计划、不重跑。
- 对 DeepSeek：建议重命名 `skeleton_ebiomedicine.tex` 越出本轮范围（`06_EC_IR/` 明确不在合并/重跑范围，且要求单独表态），应登记为文稿侧待办而非 P2 代码项；其将 proxy 种子改为 `config.SEED` 同样须先核验是否进入文稿。
- 对 GLM 与 Kimi：整体无异议；两者的核对协议（GLM 层 3 对照表、Kimi 的 `docs/reconciliation_plan.tsv`）应合并为一个机器可读对照表，列含 claim_id / 文本位置 / 声称值 / 舍入规则 / 源 tsv 文件列 / 状态（match / rounding_match / mismatch / source_missing / 未核验）。
- 对 Grok：整体赞同其门槛式、失败即停的风格；唯一小分歧是 MR 自检种子——Grok 主张保留 `seed=1` 不改写，我倾向改为 `config.SEED`，因该块不写 results、改动零风险，且仓库自述要求"唯一常数是种子"的统一性；两方案皆可，但须二选一写进计划，不能模糊。

[修正] 首轮立场（无上一轮发言）。基于审阅，我将三处立场收紧：① revision_20261116 种子从"直接改 config.SEED"调整为"先核验是否进入文稿再裁决"；② `zmat`/`build_background`/`draw_module` 从"可抽取共享"调整为"默认不共享，数值回归通过才抽，否则改名保留两份"；③ 为图读表增加"锁定舍入规则 + 黄金对照失败即停"的硬约束。

[收敛] NO — 未决实质分歧：① revision_20261116 种子处置（Nemotron"必改" vs Grok"禁改" vs 我的"先核验再裁决"）；② zmat 类函数是否默认共享；③ 图内舍入规则（两位标签 vs 三位坐标）未锁定；④ `config.py`/`config.R` 的 SEED 一致性未核验。P0 清单（硬编码数字、身份信息、README 声明）本身已达成一致。
<<<END>>

### intl-xai Grok 4.7 (x-ai/grok-4.7)

<<<ANSWER>>>
[立场] 修正：方向可保留，但五份方案都不能原样批准；若干“必须改”会改数字、改错事实或把硬要求悄悄改写。
[理由] 可执行裁决只保留这些：作图科学数字必须经中间变量来自 `results/tables/concordance_by_cohort.tsv`，且先读表头、写明 cohort 连接键，并用断言锁住 adipose/diabetes/clamp 与舍入后的 +0.31/+0.20/+0.10，连接失败就停，禁止换表；版式坐标（74/83/92 与 15.0/16.6）不是科学结果，可留在代码里。作者姓名和邮箱删除，不替换成新的团体名；R 里无刊名的 “journal requirement” 改为中性英文即可。`23.28%` 只删除，不指向未定位的表。FLOW 阶段脚本不合并；同名函数先 diff，等价才抽到 `code/common/{numerics,locus_utils,io_utils,sampling}.py`，不等价就改名保留；`_parse_args` 不抽。带种子的 bootstrap/置换/模块抽样保留。不新增模拟路径，不把缺失原始数据换成模拟。本轮不改仓库。硬要求 1（减文件数）的正式裁决是拒绝合并，抽取会使文件数增加约 4，必须由用户批准这一否定，不能写成“已满足精简”。硬要求 2 对既有 `10_positive_control.py` 和 MR `__main__` 合成自检只能作为待批准例外，不能宣称已经合规。硬要求 3 与“原始数据不在 git、缺文件即停”冲突，分阶段停报是在请求用户改写该硬要求，不是执行该硬要求。
[异议] DeepSeek：把同名函数写成“显然复制粘贴”未 diff；用 “BioHex Consortium” 替换作者行是新身份信息；“必须立即修改”违反本轮不改仓库；把 `skeleton_ebiomedicine.tex` 列入修改项越出正式包；README 改写把未被证伪的“无硬编码 path/accession/threshold”一并推翻，并写成另一句未证实的绝对声明；未处理 clamp 的 0.0976 与代码 0.098 不能互相替代。GLM：把 `default_rng(1)` 改成 `config.SEED` 列为 P0 无依据——这是已播种、不写 results 的自检，且未证明没有绑死 seed=1 的断言；“不纳入 revision 却仍改其种子”自相矛盾；其余关于 SEED 跨语言未核验、不能为跑绿而改写硬要求 3，成立。Kimi：代码里是 0.098 不是 0.0976；图中并未显示写死 61.6%，此条是编造；要求图脚本额外读入 concordance 与 n_responsive 会改变图，超出已定位缺陷；15.0/16.6 已被核实为版式坐标，不能再标成未核验；`Created:` 不在清除范围内。Mistral：在未证明阳性对照写入主结果表之前，强制迁到 `results/positive_control/` 可能打断 FLOW；给已发表 tsv 加 seed 列会改表结构，provenance 只能做旁路 manifest；把自检搬进新的 `tests/` 增加文件且同样可能改行为；“图面不变”在 y 使用全精度 rho 时不成立。Nemotron：`proxy_uncertainty.py` 的种子改写被自己标为可能改变已发表数字，仍列入无条件 P0，逻辑不成立；`_parse_args` 不应抽取；README 的 liver 值是 0.93 不是 0.931。各方共同遗漏：事实 C 写表内 rho 为 0.098，事实 F 写 0.0976，核对必须读单元格原文；CI −0.014 to +0.031 不得预设文件名 `incremental_boot.tsv`；`config.py` 与 `config.R` 的 SEED 是否相同未核验；重跑检查单应覆盖独立实验单位、伪重复、批次/组成混杂、多重检验、效应量与 CI，但这是审查项，不是已证实缺陷。
[修正] 无修改（首轮，无既有立场）。
[收敛] NO
<<<END>>>

# 你的任务

你是 **Nemotron 3.5 Super（nvidia/nemotron-3.5-super-vl）**，负责：数据分析 / 长上下文。
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
