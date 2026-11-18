<<<ANSWER>>>
[立场] 修正 —— 采纳其他五方对我先前立场的合理批评，撤回 tests/ 迁移、团体署名、单一共享文件、revision 种子修改等提议；当前唯一未决的核心分歧是 Kimi 对 MR 自检种子的模糊立场，需其明确选择保留 seed=1 后方可收敛。

[理由]
1. **MR 自检种子**：`01_mr_estimators.py` 的 `__main__` 块是合成自检，不写结果表，`seed=1` 是测试 fixture 的预言常量。改种子零收益，且可能破坏绑定该次抽样的断言，使自检失去可判定性。原地保留 + 英文注释隔离已满足“与分析路径隔离”的要求，无需迁入 `tests/`（仓库是否已有测试入口未核验，迁入可能改变导入行为）。Kimi 仍持“两种方案均可”，未明确选择保留 seed=1，这是当前唯一阻碍收敛的实质分歧。
2. **revision_20261116/**：本轮不改任何代码，仅核验其输出是否进入 `results/tables/` 或文稿数字，后续处置（归档/删除/种子统一）一律交用户批准。我撤回先前将其种子修改列入计划项的提议，因为该目录不在正式包范围，任何修改都属范围蔓延。
3. **共享模块**：撤回单一 `shared_utils.py` 的提议，改为按主题分文件：`code/common/locus_utils.py`（`_find_member`, `genes_near_loci`）、`code/common/numerics.py`（`_logsum`）、`code/common/io_utils.py`（`resolve_input`, `remote_size`）、`code/common/module_sampling.py`（`zmat`, `build_background`, `draw_module`，默认不共享，须过 diff + 零容差数值回归 + RNG 调用序列一致性三重门槛才抽取）。`_parse_args` 不合并。
4. **身份信息**：撤回用 “BioHex Consortium” 替换作者行的提议，改为纯删除作者姓名与邮箱行；R 脚本 “journal requirement” 改为中性英文；`skeleton_ebiomedicine.tex` 本轮不动。
5. **图读表**：`02_figure1_design.py` 的科学数字从 `results/tables/concordance_by_cohort.tsv` 读入，连接键用已核实的 GSE accession 映射，舍入锁定为坐标 `f"{rho:.3f}"`（得 0.313/0.202/0.098）、标签 `f"{rho:+.2f}"`（得 +0.31/+0.20/+0.10），黄金对照失败即停。版式坐标保留。
6. **README**：后置于全量硬编码审计后重写，仅修正已证伪的“唯一常量是随机种子”子句，不推翻未证伪部分，不写新的绝对声明。
7. **重跑与核对**：分阶段，缺数据即停，禁止模拟补洞；CI 来源定位前禁止重跑 `03_burden/07_incremental_value.py` 等 bootstrap 脚本；核对协议为机器可读对照表，已核验项可预填，其余标未核验。
8. **阳性对照**：`10_positive_control.py` 保留，先核验输出位置，若混入正式表则隔离标注，不新增模拟路径。

[异议]
- **对 Kimi**：MR 自检种子“两种方案均可”的表述必须收敛为“保留 seed=1”。该块不写结果，改种子零收益且可能破坏绑定该抽样的断言；若保留“改用 config.SEED”选项，则计划存在歧义，无法通过黄金对照验证。请 Kimi 明确选择保留 seed=1。其余立场（revision 不动、README 限范围、共享模块门槛等）无异议。
- **对 Grok**：Grok 对 Kimi 的“图中写死 61.6%”和“要求读入 n_responsive”的异议，在 Kimi 最新立场中已不存在，可能为误读，但无伤大雅。Grok 其余论证严谨，采纳。
- **对 Mistral**：Mistral 指出我先前“把 revision 目录的种子修改列入计划项”与范围矛盾，已修正；其三重门槛、默认不共享 `zmat`、双舍入黄金对照等均采纳。
- **对 Nemotron**：Nemotron 指出我先前提议单一 `shared_utils.py` 是反模式，已修正为按主题分文件；其“先核验隔离、再统一种子”的两步门槛原则采纳，但本轮不改 revision 代码。
- **对 GLM**：无异议，其两步门槛、主题分文件、原地保留 MR 自检、CI 定位前禁重跑等均采纳。

[修正]
1. 撤回强制迁入 `tests/` 的提议，改为原地保留 `seed=1` + 英文注释隔离。
2. 撤回用团体名替换作者行的提议，改为纯删除。
3. 撤回单一 `shared_utils.py`，改为按主题 `code/common/{locus_utils, numerics, io_utils, module_sampling}.py`。
4. 撤回对 `revision_20261116/` 任何种子修改的提议，改为本轮零改动，仅核验隔离，用户裁决。
5. 撤回重命名 `skeleton_ebiomedicine.tex` 的提议。
6. 明确 MR 自检 `seed=1` 必须保留，不提供“改用 config.SEED”选项。

[收敛] NO
原因：Kimi 仍未明确放弃“MR 自检改用 config.SEED”的选项，该分歧未消除。待 Kimi 确认保留 seed=1 后，所有核心分歧将消除，可收敛。
<<<END>>>