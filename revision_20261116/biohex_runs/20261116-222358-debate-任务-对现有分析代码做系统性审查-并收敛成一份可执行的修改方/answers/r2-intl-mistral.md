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