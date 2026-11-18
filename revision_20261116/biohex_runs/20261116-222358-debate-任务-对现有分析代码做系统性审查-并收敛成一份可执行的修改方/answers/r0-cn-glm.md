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