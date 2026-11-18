<<<ANSWER>>>
[立场] 修正 — 撤回我上一轮两项有风险的 P0 提议（`default_rng(1)` 改 `config.SEED`、revision 脚本种子统一），采纳“等价性门槛 + 黄金对照 + 先核验后裁决”框架；剩余分歧收窄到 MR 自检的位置与 revision 失败分支两处，尚未收敛。

[理由]
1. **MR 自检的最终裁决（我方修正后的立场）**：`05_mendelian_randomisation/01_mr_estimators.py` 的 `__main__` 块**原地保留**，`seed=1`、`TRUE=0.30`、`n=400` 一律不改，仅在块前加英文注释：“synthetic self-test with fixed seed 1; not part of the analysis pipeline; writes no results tables”。理由：(a) 该块不写 results，改种子收益为零；(b) 改种子可能掩盖跨环境浮点差异、破坏可能绑定 seed=1 的期望值断言（xai 的论证成立，我上一轮低估了）；(c) 移入新建 `tests/` 既增加文件数（与用户减文件意愿相反），又可能因导入路径变化改变行为。此方案满足 mistral 的“二选一”（选保留 seed=1 这一支）、kimi 的两可、nvidia 的必要性质疑、xai 的全部要求；唯一残留反对是 deepseek 的 tests/ 迁移——但 deepseek 的核心诉求是“保种子 + 与分析路径隔离”，原地 + `__main__` 门控 + 显式注释已实现隔离，请 deepseek 确认可接受等价实现。
2. **revision_20261116 统一裁决**：整体排除在正式包外，本轮**不改任何代码**（含种子与中文注释）；登记核验任务——确认三份脚本输出是否进入 `results/tables/` 或文稿数字；后续处置（归档/删除/种子统一）一律交用户批准。此方案吸收 mistral 的“先核验再裁决”，但删去“核验通过即自动改种子”分支：该目录不在正式包内，任何修改都属范围蔓延，须用户点头。nvidia 的“必须改”与自身“可能改变已发表数字”的警告自相矛盾（xai 已指出），应放弃。
3. **图脚本修复协议（P0）**：`02_figure1_design.py` 的 PTS 科学数字改为从 `results/tables/concordance_by_cohort.tsv` 读入中间变量；连接键为 cohort 标签（adipose/diabetes/clamp 映射必须断言）；舍入规则锁定——图注/坐标三位小数（0.098），标签带符号两位（+0.10）；加载后与现硬编码值做黄金对照，不一致即停、不得发布。版式坐标（74/83/92、15.0/16.6）保留。**前置任务**：先读表单元格原文——事实 C 写代码为 0.098、事实 F 写表值为 0.0976，二者是否同一精度未定，写加载器前必须确认（xai 的提醒成立）。`23.28%` docstring 数字只删除，不指向未定位的表。
4. **共享模块（P1）**：仅当逐行 diff 完全一致 + 受影响脚本输出表零容差回归通过后才抽取；`zmat`/`build_background`/`draw_module` 默认保留两份，不等价则改名；`_parse_args` 永不合并（nvidia 的 `cli.py` 提议有循环导入风险，放弃）。目标文件 `code/common/{io_utils, locus_utils, numerics, sampling}.py`。必须向用户如实说明：文件数净增约 4，“精简”的度量是消除的重复行数而非文件数，阶段编号 00–07 不合并——这是对硬要求 1 字面表述的否定，须用户批准，不得包装成“已满足”。
5. **README 声明（P1，非 P0）**：先由主 agent 全仓检索硬编码 path/accession/threshold/科学数字并出具报告，再重写；本轮只确认“唯一固定常量是随机种子”子句已被 PTS 硬编码证伪，其余子句未证伪不得一并推翻，也不得写新的绝对声明。
6. **重跑与核对**：阶段 0 自检（不依赖大文件）；阶段 1 依赖本地数据，缺文件即停并报缺失清单；禁止模拟数据补洞。采纳 deepseek 的硬约束：CI −0.014~+0.031 的来源文件定位并验证之前，**禁止重跑** `03_burden/07_incremental_value.py` 及任何 bootstrap 脚本，防止覆盖现有结果表。核对协议合并为一张机器可读对照表（claim_id / 文本位置 / 声称值 / 舍入规则 / 源 tsv 文件与列 / 状态 match·rounding_match·mismatch·source_missing·未核验）；`config.py` 与 `config.R` 的 SEED 一致性列为核验项（未核验，不预设相同）。
7. **身份信息（P0）**：作者姓名/邮箱行**直接删除**，不替换为任何团体署名（"BioHex Consortium" 同样是身份信息，deepseek 此议撤回）；R 脚本 "journal requirement" 改中性英文；`skeleton_ebiomedicine.tex` 属文稿、不在本轮范围，登记为文稿侧待办并单独向用户表态，不列入代码修改项。
8. **`10_positive_control.py`**：保留（既有功效模拟，用户规则允许经同意的模拟，此为存量）；先核验其输出实际位置，未证明混入主结果表之前**不迁移目录**（mistral 的预迁移有打断 FLOW 风险）；加模拟性质标注；不新增任何模拟数据路径。

[异议]
- **对 deepseek**：(a) "BioHex Consortium" 替换作者行是引入新身份信息，应纯删除；(b) tests/ 迁移增加文件且可能改导入行为，原地+注释已达成其隔离目标；(c) 将 `skeleton_ebiomedicine.tex` 重命名列入修改项越出正式包范围；(d) README 降 P1 正确，但其“全面检索后重写”须明确排除未被证伪子句的一并推翻。其“CI 定位前禁重跑”与零容差回归两项我方采纳。
- **对 kimi**：(a) 称图中写死 61.6% —— 已核实事实中 PTS 仅含 rho 三元组，61.62% 只出现在表与 README，此说法无依据（xai 已指出，属编造）；(b) 称代码含 0.0976 —— 代码是 0.098；(c) 要求图脚本额外读入 n_responsive 会改变图面、超出已定位缺陷，属范围蔓延；(d) 15.0/16.6 已核实为版式坐标，不应再标未核验。其余（revision 不动、README 限范围、两可方案）我方赞同。
- **对 mistral**：(a) MR 自检改 `config.SEED` 撤回支持——零收益、非零风险，且 mistral 自己已声明 Grok 的保留方案“皆可”，应直接选定保留 seed=1；(b) 预先将 positive_control 迁至 `results/positive_control/` 应改为先核验输出位置；(c) 给已发表 tsv 加 seed 列会改表结构，provenance 只能做旁路 manifest（xai 的反对成立）；(d) revision“核验通过即改种子”分支应改为“核验后交用户裁决”。
- **对 nvidia**：(a) `proxy_uncertainty.py` “必须改种子”与自认“可能改变已发表数字”逻辑矛盾；(b) `_parse_args` 抽取撤回；(c) 对 `TRUE=0.30` 的追问合理，裁决为：保留为自检预言值并注释说明，不属硬编码结果数字（不进任何结果路径）。
- **对 xai**：整体最严谨，但一处错误——其称“分阶段停报是在请求用户改写硬要求 3”。任务约束 G 由用户本人明文规定“缺文件就停止并报告、禁止为跑通改用模拟数据”，分阶段方案正是**执行**该约束而非改写它；xai 将用户自己的约束误判为对硬要求的软化。其余（0.098/0.0976 须读单元格原文、CI 不得预设文件名、SEED 跨语言未核验、文件数净增须用户批准）我方全部采纳。

[修正]
1. 撤回上一轮将 `default_rng(1)`→`config.SEED` 列为 P0 的提议，改为原地保留 seed=1 + 英文注释隔离。
2. 撤回对 `revision_20261116/proxy_uncertainty.py` 种子统一的提议（与“排除该目录”自相矛盾，xai 指出成立），改为整体不动 + 核验登记 + 用户裁决。
3. README 修改从 P0 降为 P1，且限定为“检索完成后仅重写已证伪子句”。
4. 新增采纳：图读表的黄金对照与舍入规则锁定（deepseek/mistral/xai）、CI 定位前禁重跑（deepseek）、positive_control 先核验输出位置再谈迁移（xai）、对照表机器可读化并合并 GLM 层 3 与 kimi 的 reconciliation_plan（mistral）。
5. 明确向用户声明：文件数将净增约 4，精简度量是重复行数；此为对硬要求 1 字面表述的否定，须批准。

[收敛] NO
剩余实质分歧仅两处：(1) MR 自检块位置——deepseek 主张迁 `tests/`，xai 与我方主张原地保留（种子=1 已无分歧，mistral 声明两案皆可）；(2) revision_20261116 核验后的失败分支——mistral 主张“未进文稿即改种子”，其余四方主张“任何修改交用户批准”。请 deepseek 确认原地+注释是否满足其隔离诉求、mistral 确认放弃自动改种子分支后，下一轮可收敛。
<<<END>>>