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

### cn-glm GLM (z-ai/glm-5.3-prime)

<<<ANSWER>>>
[结论] 该包整体质量较高（种子齐全、无裸随机数、结果数字基本来自 results 表），但存在四类必须修复的问题：作图代码内嵌结果数字、README 关于“唯一常量是种子”的声明与事实不符、跨文件重复函数未抽共享模块、文件头身份信息未清除；修改方案应以“不改变任何已发表数字”为红线，分 P0（数字溯源与身份信息）、P1（去重与注释统一）、P2（重跑验证与核对协议）执行，本轮只出计划不改仓库。

[理由]
1. **随机数三类裁决（回应硬要求 2）**：
   - 统计程序类（7 个脚本，种子走 `config.SEED` 或 `--seed`）：**保留**。bootstrap/置换/模块抽样是方法本身，种子化保证可复现，不属于“非种子随机数”违规。
   - `10_positive_control.py`（只模拟效应量，位点/SE 来自真实数据）：**保留但加护栏**——在文件头与输出表元数据中显式标注“simulated effect sizes, positive control, not primary result”，并确保其输出不进入主结果表路径。这是阳性对照，不是伪造数据。
   - `05_mr/01_mr_estimators.py` 的 `rng = default_rng(1)` + `TRUE=0.30`：**改写为 P1**——把种子改为 `config.SEED`（或 `--seed` 参数），并明确该块是 `__main__` 单元自检、不写 results。种子写死为 1 虽不进结果，但违反“唯一常量是种子”的仓库自述，且 `revision_20261116/proxy_uncertainty.py` 的 `default_rng(20261116)` 同理应统一到 `config.SEED`（若该脚本要纳入正式包；若维持“修订脚本不在合并范围”，则标注为遗留、不进 FLOW）。
2. **硬编码结果数字（回应硬要求 2 的另一半）**：`07_figures/02_figure1_design.py` L169–171 的 PTS 三元组是**必须改的 P0**。正确修法：从 `results/tables/concordance_by_cohort.tsv` 读入 rho，经中间变量（如 `df_fig = read_tsv(...)` → `pts = df_fig[["x_pos","rho","label"]]` 或按 cohort 键映射）传入绘图函数；标签 "+0.31" 等由 `f"+{rho:.2f}"` 生成。**反对**把示意图改成从错误的表（如 incremental_value.tsv）读取。`02_select_instruments.py` docstring 里的 "23.28%" 属注释性结果引用，P1 改为“见 results/tables/ 对应表”或删除具体数字，避免注释与数据漂移。
3. **README 声明与事实的冲突**：README 称"No analysis script contains a hard-coded path, accession or threshold. The only constant fixed in source is the random seed."——这与 (a) `02_figure1_design.py` 内嵌 0.313/0.202/0.098、(b) `01_mr_estimators.py` 的 `default_rng(1)` 与 `TRUE=0.30`、(c) `revision_20261116/proxy_uncertainty.py` 的写死种子，均不完全相符。**P0 修改 README**：要么改声明为准确表述（如"analysis constants are sourced from config and results tables; figure annotations are read from results tables"），要么先修代码再保留原声明。二者取其一，不能声明与代码长期背离。
4. **去重与文件合并（回应硬要求 1）**：用户要减文件数，但**反对**合并 FLOW 阶段脚本——阶段编号是运行顺序的机器可读契约，合并会破坏 `docs/FLOW.md` 并增加单文件腐化风险。正确做法是**抽共享模块而非减脚本数**：
   - 新建 `code/common/`（或 `code/_shared/`）：`coloc_math.py`（收 `_logsum`）、`loci_utils.py`（收 `_find_member`、`genes_near_loci`）、`io_utils.py`（收 `resolve_input`、`remote_size`）、`module_sampling.py`（收 `zmat`、`build_background`、`draw_module`）。各原脚本改为 import，函数体删除。
   - `_parse_args` 同名是 argparse 惯例，**不合并**。
   - 判断依据：`_logsum`/`resolve_input`/`zmat` 等是逐字符或语义级重复（同一算法多处拷贝 = 腐化温床）；`_find_member` 需主 agent diff 确认两处实现是否一致，若语义不同则各自保留并改名（如 `_find_member_loci` vs `_find_member_gtex`），**不强行统一**。
   - 净效果：文件数不降反可能持平，但重复代码行数下降、单一事实来源建立——这才是“精简”的正确度量，需向用户明确说明“减文件数”与“减重复”在此处应取后者。
5. **身份信息（回应硬要求 5）**：所有 `.py` 文件头的 `Author : Tao Zhu (zhutao@tjh.tjmu.edu.cn)` 是 P0 清除项（含邮箱）。`01_main_figures_1_2_3.R` 的 "the journal requirement for composite figures" 无具体刊名，P1 改为中性表述 "composite figure layout convention"。`06_EC_IR/paper_ms/skeleton_ebiomedicine.tex` 文件名含期刊名，但属文稿且不在本轮合并/重跑范围——**本轮不动**，仅在计划中登记为“文稿侧待办，需用户单独批准”，避免范围蔓延。
6. **重跑与核对协议（回应硬要求 3、4、6）**：分三层——
   - **层 1（无大文件自检，P2 执行）**：`01_mr_estimators.py` 的 `__main__` 自检、`02_figure1_design.py` 从 tsv 重读后重绘、所有脚本的 `--help`/import 冒烟测试、config 一致性检查（config.py 与 config.R 的 SEED 是否同值——这一点上文未给出，需主 agent 核验）。
   - **层 2（依赖本地数据）**：按 FLOW 顺序逐阶段跑，缺文件即停并输出缺失 accession 清单（对照 `data/README.md`），**禁止**以模拟数据补洞。
   - **层 3（文稿核对协议）**：建对照表，列为 {文稿句子/表格单元格 → 源 tsv 文件与行列 → 一致/四舍五入一致/不一致}。已知的必查项：README 的 95% CI (−0.014, +0.031) 必须定位到具体表（候选：incremental_boot 类输出，**未核验**，不得假定）；table1 的 503/3/0.962/0.931 与 coloc 汇总表逐格核对；concordance 三 rho 与 61.62% 的 README 表述。核对由主 agent 在方案批准后执行，本轮不声称已完成。

[风险与不确定]
- **数字不变性风险**：抽共享模块若两处 `_logsum`/`_find_member` 实现有细微差异（浮点容差、边界处理），统一后可能改变下游小数位。缓解：抽取前先 diff 两实现；抽取后对受影响表做逐值比对（P2 层 1 的比对脚本），任何数值变动即回滚并改为“保留两份、改名区分”。
- **`revision_20261116/` 归属未定**：三份修订脚本是否纳入正式包影响 README 声明的修改范围与 FLOW 更新，需用户裁决；我倾向不纳入（早期修订痕迹），仅修其种子写法以防未来被引用。
- **CI 来源缺口**：README 的 ΔAUC 置信区间来源文件未定位，若最终找不到生成该 CI 的表，则 README 该句需标注“来源待补”或删除——这是潜在的真实不一致点，核对协议必须覆盖。
- **config.py 与 config.R 的 SEED 一致性未核验**：若两者不同值，Python 与 R 随机路径不可比，需在层 1 自检中加入断言。
- **重跑可行性**：原始数据不在 git，层 2 大概率部分阶段无法执行；计划必须预设“缺数据 → 报告而非造数据”的终止路径，且用户硬要求 3 的“修到能跑通”应解释为“能跑通的部分跑通、缺数据的部分明确报告缺失”，这一点需在计划中向用户显式确认，避免被理解为“无论如何都要跑绿”。

[建议]
1. **P0（改前必须用户批准，预期不改变任何已发表数字）**：
   - `07_figures/02_figure1_design.py`：PTS 三元组改为从 `results/tables/concordance_by_cohort.tsv` 读入中间变量再绘图；标签字符串由 rho 格式化生成。
   - 全部正式脚本文件头：删除 Author/邮箱行，替换为英文用途说明 + 依赖说明。
   - `README.md`：修正"No hard-coded constant"声明使其与代码事实一致（在上述两项完成后）。
   - `05_mr/01_mr_estimators.py`：`default_rng(1)` → `config.SEED`，docstring 标注单元自检性质。
2. **P1（结构性，预期不改变数字，但需比对验证）**：
   - 新建 `code/common/`（`coloc_math.py`、`loci_utils.py`、`io_utils.py`、`module_sampling.py`），迁移 D 项重复函数；迁移前 diff 确认实现一致，迁移后跑数值比对。
   - 注释统一英文：`revision_20261116/` 中文注释（若保留该目录则翻译；若不纳入则不动）；`02_select_instruments.py` docstring 的 "23.28%" 改为指向结果表。
   - `01_main_figures_1_2_3.R` 去期刊措辞；`proxy_uncertainty.py` 种子统一到 `config.SEED`。
3. **P2（验证与核对，本轮只定协议不执行）**：
   - 层 1 自检清单（含 config SEED 一致性断言、figure 重绘比对、共享模块数值比对脚本）。
   - 层 2 数据依赖重跑计划：按 FLOW 阶段列出所需输入文件与 accession，缺即停即报。
   - 层 3 文稿-代码对照表模板：{tex/README 位置， 源 tsv， 核对状态}，重点覆盖 ΔAUC 的 CI 来源、table1 五组织汇总、concordance 三 rho。
4. **明确反对的做法**（写入计划供辩论）：
   - 反对合并 FLOW 阶段脚本以凑“减文件数”指标；
   - 反对删除或去种子化任何 bootstrap/置换程序；
   - 反对把 figure1 示意数字改从非 concordance 表读取；
   - 反对为“跑通”引入任何新模拟数据路径（`10_positive_control` 保留但加显式标注护栏）；
   - 反对本轮直接改仓库——所有 P0/P1/P2 均为待批准计划项。
   辩论收敛条件：全体到场成员对上述 P0–P2 清单与反对清单无异议后，方可 `[收敛] YES`；任何成员对共享模块归属或 README 措辞有实质分歧，则该单项降级为“待用户裁决”而非强行收敛。
<<<END>>>

### cn-kimi Kimi (moonshotai/kimi-k2.7-code)

<<<ANSWER>>>
[结论]  
正式包 `ec-metabolic-coloc/` 需要在**不破坏 FLOW 阶段顺序**的前提下，把跨文件重复函数抽到共享模块、清除代码中的硬编码结果数字与身份信息，并建立“结果数字 → 源表”的可追溯核对协议；当前尚未达到 `[收敛] YES`，需由主 agent 就下列修改项逐条组织辩论。

[理由]  
1. **硬编码结果数字是 P0 级违规**：`code/07_figures/02_figure1_design.py` 把 `concordance_by_cohort.tsv` 的 rho（0.313、0.202、0.0976）直接写在作图代码里，违反“数据必须先读入中间变量再传递”；`05_mendelian_randomisation/02_select_instruments.py` 的 docstring 也嵌入了 “23.28%” 这一结果数字。两者都必须改掉。  
2. **共享模块比合并文件更合理**：`_find_member`、`genes_near_loci`、`_logsum`、`resolve_input`、`zmat` 等函数在多个阶段重复出现；把它们抽到 `code/_shared/` 下的主题模块，可在不破坏 `docs/FLOW.md` 运行顺序的前提下减少冗余，也便于单元测试。  
3. **随机数整体可保留，但种子来源要统一**：统计程序里的 bootstrap / 置换 / 模块抽样均来自 `config.SEED` 或 CLI `--seed`，属于合法、可复现的随机过程，应保留；`02_colocalisation/10_positive_control.py` 的效应量模拟是既有功效/阳性对照，可保留但必须在文档中明确标注为 simulation；`05_mendelian_randomisation/01_mr_estimators.py` 的单元自检把种子写死为 1，应改写为默认 `config.SEED`，但不影响已发表结果表。  
4. **README 的自我声明需要修正**：README 称 “No analysis script contains a hard-coded path, accession or threshold. The only constant fixed in source is the random seed.” 这一说法与已发现的硬编码结果数字、版式常量等不完全相符，应改为更准确的限定语句。  
5. **身份信息清除范围应限于代码与运行文档**：所有正式 `.py`/`.R` 文件头的 `Author: Tao Zhu (...)` 和邮箱必须删除；`07_figures/01_main_figures_1_2_3.R` 中 “the journal requirement for composite figures” 这类指向投稿语境的注释也应删除；`06_EC_IR/paper_ms/skeleton_ebiomedicine.tex` 是文稿文件而非分析代码/运行文档，本轮不建议改动其文件名。  
6. **重跑必须分阶段且禁止用模拟数据补洞**：原始数据不在 git 中，计划应区分为（a）不依赖大文件的静态/单元自检、（b）依赖本地数据的阶段并在缺文件时停止报告、（c）禁止为跑通而改用模拟数据；`10_positive_control.py` 作为已有模拟路径可保留，但不得新增模拟数据分支。  

[风险与不确定]  
- **硬编码数字的精确映射未核验**：已知 `02_figure1_design.py` 的三个 rho 与 `concordance_by_cohort.tsv` 一致，但其余 15.0/16.6 等版式坐标是否也隐含科学结果、以及该文件是否还有其他硬编码数字，需要主 agent 通读确认。  
- **CI 来源未确认**：README 中 incremental value 的 95% CI（−0.014 to +0.031）不在 `incremental_value.tsv` 的同一行，是否来自 `incremental_boot.tsv` 或另一张表尚未核验；核对协议必须把它标为待确认项。  
- **重复函数语义是否完全一致**：`_logsum`、`resolve_input` 等同名函数跨文件出现，需要主 agent 逐行 diff 判断是真重复还是因历史原因有细微差异，避免抽共享模块时改变行为。  
- **revision 脚本是否纳入正式路径**：`revision_20261116/proxy_uncertainty.py` 使用 `np.random.default_rng(20261116)` 且含中文注释；如果该脚本属于正式重跑路径，则种子应走 `config.SEED` 且注释应改为英文；若仅为一次性修订稿，可降级为 P2。  
- **“改变已发表数字”风险**：把硬编码数字改为读表理论上**不应**改变数字，但若读表时列名、四舍五入或行筛选与当前硬编码不完全一致，示意图可能产生微小差异；需用版本控制快照比对。  

[建议]  
**P0（必须先做，否则影响数据真实性与可复现性）**  
- **P0-1 消除硬编码结果数字**  
  - 改 `code/07_figures/02_figure1_design.py`：从 `results/tables/concordance_by_cohort.tsv` 读入 clamp/T2D/adipose 的 rho、concordance、n_responsive，再绘制；禁止把 0.313/0.202/0.098/61.6% 等科学数字直接写在代码里。  
  - 改 `code/05_mendelian_randomisation/02_select_instruments.py`：把 docstring 中的 “23.28%” 替换为指向生成该统计量的脚本/结果表的说明，或删除。  
  - 预期效果：满足用户“数据必须先读入中间变量再传递”的硬要求；不改变已发表数字（读同一源表）。  
- **P0-2 清除代码中的身份信息**  
  - 删除所有正式 `.py`/`.R` 文件头的 `Author: Tao Zhu`、`zhutao@tjh.tjmu.edu.cn`、`Created:` 等身份字段；保留功能性的 `Purpose/Input/Output/Usage` 英文注释。  
  - 删除 `07_figures/01_main_figures_1_2_3.R` 中 “the journal requirement for composite figures” 这类投稿语境注释。  
  - 预期效果：满足匿名化与英文注释要求；不改变任何结果。  
- **P0-3 建立结果数字核对协议（由主 agent 执行）**  
  - 创建 `docs/reconciliation_plan.tsv`，列：`manuscript_location`（README 句子 / tex 表 / 图注）、`claimed_value`、`source_table`、`source_column`、`source_row_key`、`status（一致/不一致/未核验）`、`notes`。  
  - 必须包含的待确认项：README 中 incremental value 的 95% CI 来源（`incremental_value.tsv` 还是 `incremental_boot.tsv` 或其他）。  
  - 先核对 `concordance_by_cohort.tsv`、`incremental_value.tsv`、`table1_coloc_summary.tex` 与 README 的对应关系，再扩展到所有 `results/tables/*.tsv`。  
  - 预期效果：为后续“代码输出 ↔ 文稿”逐项交叉核对提供机器可读的追踪表；不改变代码。  

**P1（应在 P0 完成后尽快做，提升工程健壮性）**  
- **P1-1 抽取共享模块以减少重复**  
  - 新建 `code/_shared/loci_utils.py`：放入 `_find_member`、`genes_near_loci`。  
  - 新建 `code/_shared/stats_utils.py`：放入 `_logsum` 等数值辅助函数。  
  - 新建 `code/_shared/io_utils.py`：放入 `resolve_input`、`remote_size`。  
  - 新建 `code/_shared/module_utils.py`：放入 `zmat`、`build_background`、`draw_module`。  
  - 原文件改为 `from _shared.xxx_utils import ...`，保持 `docs/FLOW.md` 的阶段编号与调用顺序不变。  
  - 预期效果：减少重复代码、降低腐化风险；若 diff 确认语义一致，则不改变结果。  
- **P1-2 统一随机种子来源并文档化**  
  - 改 `05_mendelian_randomisation/01_mr_estimators.py` 的 `if __name__ == "__main__"` 自检块：把 `rng = np.random.default_rng(1)` 改为默认 `np.random.default_rng(config.SEED)` 或接受 `--seed` 参数；保留 `TRUE=0.30`、`n=400` 作为自检参数，但加注释说明这是 synthetic self-test，不写 results。  
  - 若 `revision_20261116/proxy_uncertainty.py` 属于正式路径，同样把种子改走 `config.SEED`。  
  - 在 README 或 `docs/FLOW.md` 增加一节，列出所有使用随机数的脚本、用途、种子来源。  
  - 预期效果：消除“硬编码种子 1”与 20261116 的例外；不改变正式结果表（自检不写表）。  
- **P1-3 修正 README 关于硬编码常量的声明**  
  - 把 README 中 “No analysis script contains a hard-coded path, accession or threshold. The only constant fixed in source is the random seed.” 改为准确版本，例如保留前半句（路径/accession/threshold 不硬编码），后半句改为 “The only fixed scientific constants in analysis scripts are the random seed and CLI defaults; all published numbers are read from results tables or generated by the pipeline.”  
  - 预期效果：与代码实际状态一致，避免误导。  
- **P1-4 英文化与注释规范化**  
  - 把 `revision_20261116/` 三份脚本的中文注释翻译为英文（若保留在仓库中）。  
  - 统一文件头格式：`"""Purpose: ... | Inputs: ... | Outputs: ... | Usage: ..."""`，禁止出现作者、邮箱、期刊名。  
  - 预期效果：满足注释统一要求；不改变结果。  

**P2（可做可后做，低风险优化）**  
- **P2-1 静态检查与轻量级重跑计划**  
  - 制定 `docs/reproduce_checklist.md`：Phase A 为 `python -m py_compile`、`R CMD check`、自包含 `__main__` 单元自检；Phase B 为检查 `data/README.md` 所列原始文件是否存在，缺文件则停止并报告；Phase C 为完整重跑；Phase D 为与 `results/tables/` 和文稿比对。  
  - 预期效果：让后续重跑有明确入口，避免为跑通而造假数据。  
- **P2-2 文稿文件名**  
  - `06_EC_IR/paper_ms/skeleton_ebiomedicine.tex` 含期刊名，但它属于文稿而非分析代码/运行文档；本轮不建议重命名，若用户坚持可在文稿提交前单独处理。  
  - 预期效果：避免本轮范围蔓延。  

**我明确反对的做法**  
- 把 `00_setup` 到 `07_figures` 的 45 个脚本合并成少数几个巨文件：这会破坏 `docs/FLOW.md` 的运行顺序、增加调试难度、模糊阶段责任。  
- 删除所有带种子的 bootstrap / 置换 / 模块抽样：这些是合法统计推断所必需，删除会削弱结论的统计基础。  
- 为了“跑通”而新增模拟数据分支或把缺失的原始数据替换为合成数据：违反用户“分析必须用真实数据”的硬规则。  
- 把 `02_figure1_design.py` 的硬编码数字改为从错误的表读取：必须先由核对协议确认 `concordance_by_cohort.tsv` 就是唯一正确源表，再修改代码。  
- 本轮直接修改仓库：必须先由主 agent 把上述计划整理为修改计划，经用户批准后方可执行。  
<<<END>>>

### intl-mistral Mistral Large 4 (mistralai/mistral-large-4-0)

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

### intl-xai Grok 4.7 (x-ai/grok-4.7)

<<<ANSWER>>>
[结论] 正式包只允许做“不改统计量的完整性修复”：结果数字必须从表读入、非结果随机路径必须移出正式结果、缺数据必须失败停止；禁止为减文件数合并阶段脚本，禁止本轮改仓库。
[理由]
1. 减文件数不是防腐化。腐化来自两套科学逻辑和嵌进代码的结果数字，不来自 `docs/FLOW.md` 的阶段编号。同名函数先当巧合：`_find_member`/`genes_near_loci`、`_logsum`、`resolve_input`/`remote_size` 只有在 diff 与数值等价测试通过后才可抽到共享模块；`zmat`/`build_background`/`draw_module` 是得分与零模型的契约，未证明等价前合并会改 p 值，禁止动。`_parse_args` 禁止合并。阶段脚本 00–07、图脚本与分析脚本、MR 自检与正式 MR，全部不得合并成巨文件。
2. 随机数裁决：bootstrap/置换/模块抽样（种子来自 `config.SEED` 或 `--seed`）保留在正式路径，这是带种子的推断，不是虚假数据；删除会毁掉已报告的不确定性。`10_positive_control.py` 移出正式结果路径：脚本可留，输出不得进入文稿表，不得扩展、不得用来填缺失数据。`01_mr_estimators.py` 的 `seed=1`、`TRUE=0.30`、`n=400` 是门控自检，移到测试入口，禁止写入 `results/`。`revision_20261116/` 的 `default_rng(20261116)` 移出正式结果路径；未证明它未进入文稿之前，禁止改种子。
3. `02_figure1_design.py` 的 PTS（0.313/0.202/0.098 与 “+0.31/+0.20/+0.10”）是嵌死的结果，必须经中间变量从 `concordance_by_cohort.tsv` 读入。15.0/16.6 等版式坐标不是科学结果，不得数据化。删除字面量之前，加载器必须复现现有标签字符串；复现失败则停，不得改图。`02_select_instruments.py` 注释中的 “23.28%” 必须删除，不得重算后写回新数字。
4. README“唯一写死的常数是随机种子、无硬编码 path/accession/threshold”已被 PTS 与非 `config.SEED` 种子证伪。该句必须改；path/accession/threshold 是否还有别的硬编码仍是未核验，不得靠记忆补命中。文稿数字本轮只建立对照协议，不得声称已完成全文核对。
5. 身份信息：正式包每个脚本头的作者姓名与邮箱，以及 `01_main_figures_1_2_3.R` 中指向期刊要求的句子，必须清除并改为英文文件头（目的、输入、输出、种子来源、不做什么）。方法名、队列 accession、GTEx 等不是身份，不得误删。`06_EC_IR/paper_ms/skeleton_ebiomedicine.tex` 不在本次代码合并/重跑范围，本轮禁止改文件名。
6. “重跑到能跑通”不得解释为退出码 0。缺原始数据时必须停并报告；禁止新模拟路径。功效模拟已存在，裁决是隔离而非新增，也不是静默删除后用别的模拟顶上。任何会改检验统计量、零模型、种子、重抽样单位或阈值的编辑，都是结果变更，超出本轮卫生修复，须另报用户。
[风险与不确定]
- 95% CI −0.014 to +0.031 不在 `incremental_value.tsv` 该行；来源文件未核验。在定位前重跑 `07_incremental_value.py` 或改 bootstrap 种子，可能改已写进 README 的区间。
- 0.0976 显示为 +0.10（两位）与图中 0.098（三位）是舍入政策冲突，不是已证实的计算错误。未锁定舍入规则就“从表读取”，会改图上数字。
- 无法确认：两套同名函数是否文本/语义相同；`config.py` 与 `config.R` 的 SEED 是否相等；`10_positive_control.py` 是否覆盖实证表；`revision_20261116` 是否喂给文稿；除已定位项外是否还有硬编码 path/accession/threshold。
- 单细胞/跨癌的重抽样单位若是 cell 而非 patient/library，p 值可能偏乐观。本轮不得借去重“修正”单位，否则发表数字改变。跨平台 BLAS/库版本可造成种子程序的浮点差；不得为追绿而改代码。
- 领域共识：缺数据不得用模拟填平。我的推断：PTS 的 y 值是舍入后的 rho 被当作示意图坐标。无法确认的条目不得写成已核对。
[建议]
1. P0（用户批准前不得改仓库；批准后先做这些，且不得改变已发表统计量）。`code/07_figures/02_figure1_design.py`：PTS 科学数字改为读 `results/tables/concordance_by_cohort.tsv` 到中间变量再传入作图；标签规则必须复现 “+0.31”“+0.20”“+0.10”；版式坐标不动；黄金对照失败则停。预期：图与表单源；已发表数字不变，仅当对照失败时不得发布。正式包 `code/**` 与 `config/**` 文件头，以及 `07_figures/01_main_figures_1_2_3.R`：删除作者姓名、邮箱和期刊要求措辞，改为英文头注释。预期：无数字变化。`revision_20261116/`：排除在正式结果路径外，禁止把种子改成 `config.SEED`。`10_positive_control.py`：输出不得进入 `results/tables` 的文稿来源表。`01_mr_estimators.py` 的 `__main__` 自检：不得写 results。README 绝对句：标为不合规；替换句必须等 path/accession/threshold 检索完成后再写，不得提前声称已合规。预期：文档主张改变，表内数字不变。
2. P1（等价测试通过才许动；失败则保持原文件）。仅在 diff 证明相同后抽取：`resolve_input`/`remote_size` → `code/lib/io_inputs.py`；`_logsum` → `code/lib/logsumexp.py`（须含边界数值测试）；`_find_member`/`genes_near_loci` → `code/lib/locus_index.py`。`zmat`/`build_background`/`draw_module` 默认不共享；仅当得分脚本与零模型脚本对现有表的回归完全一致时，才允许 `code/lib/module_null.py`，否则改名以暴露差异并保持两份。禁止合并 FLOW 阶段脚本。MR 自检移出分析入口到测试入口，保留 `seed=1` 与真值 0.30 作为预言，不改写。删除 “23.28%” 注释，不替换新百分比。增加 SEED 一致性检查：`config.py` 与 `config.R` 不一致则失败，不得静默统一。预期：行为不变；若回归表不一致，停止，不得为通过而改统计代码。
3. P2 与明确反对。P2：`revision_20261116/` 若仍留在仓库，只做英文注释与身份清除，不重跑、不改种子；文稿文件名本轮不动，另案由用户决定。反对：为减少文件数合并 00–07；删除全部带种子 bootstrap/置换；把示意图改成从错误表或从 README 字符串读取；把 15.0/16.6 版式坐标改成数据驱动；新增任何模拟数据路径；用 `10_positive_control` 填 colocalisation/TCGA/单细胞缺文件；在未定位 CI 来源时重跑增量价值；改 PP.H4 阈值、先验或重抽样单位；把 45 个脚本收成少数巨文件；本轮直接改仓库；把 `skeleton_ebiomedicine.tex` 改名并入代码 PR。
4. 批准后由主 agent 执行的核对与重跑，本轮不跑。对照表列：claim_id、文本位置、声称值、舍入规则、truth_file、truth_column、truth_value、status（match / rounding_match / mismatch / source_missing / 未核验）、是否会改变已发表数字。必填行只使用已核实事实：clamp rho 0.0976、p_rho 0.00412、61.62%（n_responsive=99）对 README +0.10 与 61.6%，并单列图内 0.098 对 +0.10 的舍入；T2D rho 0.202 对 +0.20；adipose rho 0.313 对 +0.31；delta_auc_IR 0.000389、p_IR 0.15、n=145、events=45 对 README +0.0004、p=0.15、145、45；CI −0.014 to +0.031 标 source_missing，检索 `results/tables/` 后记录真实文件名，找不到则保持未核验，不得假定正确；coloc 503、PP.H4>0.8 共 3（liver 1，muscle 2）、liver 最大 0.931 对 README 0.93、muscle 最大 0.962 不得记成 liver。tex 必须回链到生成它的 tsv，禁止 tex 对 tex。其余 README/tex 句子标未核验，由检索枚举，不得凭记忆加数字。重跑分层：A 不依赖大文件（导入、SEED 一致性、MR 自检、对已存在 tsv 的图加载器黄金对照），缺 tsv 即停；B 依赖本地数据的 coloc/burden/单细胞/MR/跨癌，缺文件非零退出并报告路径，禁止替代输入；C 禁止项为任何新模拟或用功效模拟充当 B。成功标准是 A 通过且 B 要么在预先声明的容差内复现现表，要么缺文件停止；数值不符是发现，不是改到变绿的许可。
<<<END>>>

### intl-nvidia Nemotron 3.5 Super (nvidia/nemotron-3.5-super-vl)

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

# 你的任务

你是 **DeepSeek（deepseek/deepseek-v4-pro）**，负责：Leader / 深度推理 / 统计推断。
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
