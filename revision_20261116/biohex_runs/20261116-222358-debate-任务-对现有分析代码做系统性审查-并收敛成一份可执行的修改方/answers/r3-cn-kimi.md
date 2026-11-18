<<<ANSWER>>>
[立场] 修正（首次完整表态）— 支持“先核验再修改”框架；MR 自检原地保留 seed=1；共享模块按主题拆分但设置严格等价门槛；本轮不改 `revision_20261116/` 任何代码；README 须在全仓硬编码审计后整句重写。

[理由]
1. **`revision_20261116/proxy_uncertainty.py` 种子**：该目录在正式可复现包 `code/{00..07}` 之外，本轮修改范围排除。`np.random.default_rng(20261116)` 是带种子随机数，不违反“禁止非种子随机数”。主 agent 应先核验其输出是否进入 `results/tables/` 或文稿数字：若已进入，该种子永久保留；若未进入，未来是否统一到 `config.SEED` 必须由用户批准，不能在本计划里预设“核验通过即自动改”。
2. **MR 自检 `seed=1`**：`05_mendelian_randomisation/01_mr_estimators.py` 的 `if __name__ == "__main__"` 块是单元自检 fixture，`TRUE=0.30`、`n=400`、`seed=1` 均不写 results 表。改为 `config.SEED` 会改变抽样实现并可能破坏断言；移到新建 `tests/` 会增加文件数、可能改变导入路径且仓库是否已有测试入口未核验。最优方案是原地保留并加英文注释，明确其为“synthetic self-test, fixed seed, not in analysis pipeline”。
3. **共享模块拆分**：按主题拆分为 `code/common/{io_utils.py, locus_utils.py, numerics.py, module_sampling.py}`。抽取必须同时满足：逐行 diff 等价、零容差数值回归通过、同种子下 RNG 调用序列一致。`zmat`/`build_background`/`draw_module` 涉及跨癌零模型与得分的统计契约，默认保留两份；仅当三重门槛全过才抽入 `module_sampling.py`，否则保留并改名以暴露差异。`_parse_args` 是各脚本 CLI 契约，不抽取。阶段脚本 `00–07` 不合并。
4. **图脚本 PTS 硬编码**：`02_figure1_design.py` 中的科学数字改为从 `results/tables/concordance_by_cohort.tsv` 经中间变量读入，连接键显式映射（adipose/GSE15773、diabetes/GSE18732、clamp/GSE13070），失败即停。锁定舍入：图坐标 `f"{rho:.3f}"` 必须复现 0.313/0.202/0.098，标签 `f"{rho:+.2f}"` 必须复现 +0.31/+0.20/+0.10；任一不符停止。版式坐标 74.0/83.0/92.0、15.0/16.6 保留。
5. **README 声明**：原句“唯一固定常量是随机种子”已被 PTS 硬编码和 MR seed=1 证伪；前半句“无硬编码 path/accession/threshold”尚未核验。必须先由主 agent 对 `code/**` 全量 grep 审计，再整句替换为可验证声明，禁止提前写“已合规”。
6. **身份信息**：所有正式脚本文件头的 `Author: Tao Zhu (...)` 直接删除，不替换为任何团体署名；R 脚本中 “journal requirement” 改为中性英文。`06_EC_IR/paper_ms/skeleton_ebiomedicine.tex` 是文稿，不在本轮代码修改范围。
7. **核对协议**：建立机器可读对照表（claim_id / 文本位置 / 声称值 / 舍入规则 / 源文件与列 / 状态）。已核实项可预填：concordance rho、incremental ΔAUC/p/n/events、coloc 503/PP.H4>0.8/liver max 0.931。README 中 ΔAUC 95% CI −0.014~+0.031 来源文件未定位前标 `source_missing`，并禁止重跑 `03_burden/07_incremental_value.py` 及任何 bootstrap/置换脚本。
8. **重跑计划**：阶段 A 为不依赖大文件的自检（import、config 一致性、图脚本读表黄金对照）；阶段 B 为依赖本地数据的阶段，缺文件即停并输出缺失清单；全程禁止用模拟数据补洞、禁止覆盖现有 `results/tables/`。

[异议]
- **对 DeepSeek**：(1) 将 `proxy_uncertainty.py` 改种子列为 P2 并设“核验后决定”仍隐含未来自动修改，我反对，应明确“用户批准”门槛；(2) 用 “BioHex Consortium” 替换作者行是引入新身份信息，应纯删除；(3) 单一 `shared_utils.py` 是反模式，按主题分文件可避免异质函数耦合；(4) 其余（CI 定位前禁重跑、零容差回归、默认不共享 zmat）采纳。
- **对 GLM**：整体最一致。一点补充：共享模块实现时须检查循环导入，建议各 common 模块不依赖 `config.py` 以外的业务脚本；MR 自检“原地+注释”方案与我相同。
- **对 Mistral**：(1) 将 MR 自检移到 `tests/test_mr_estimators.py` 会增加文件数并可能改变导入路径，我主张明确原地保留，不应作为可选方案保留；(2) “三位坐标/两位标签”双舍入规则、RNG 调用序列门槛、README 后置审计均采纳。
- **对 Grok**：(1) “无论隔离结果如何都不改 20261116”过于绝对；若主 agent 证明其输出未进入任何结果/文稿路径，且用户明确批准，未来统一种子是可接受的（但本轮绝对不改）；(2) 不强制迁 tests/、不覆盖现表、SEED 不一致不自动统一均采纳。
- **对 Nemotron**：(1) 已撤回无条件改种子和 `_parse_args` 抽取，与我一致；(2) “MR 自检可移至 tests/”作为选项仍造成模糊，我主张在计划中单选“原地保留”；(3) 0.93/0.931 四舍五入一致、共享模块默认不抽、README 审计后重写均采纳。

[修正] 无（首次完整表态，无上一轮立场可修正）。

[收敛] NO
未解决核心分歧：
1. **MR 自检位置**：Mistral 主张迁至 `tests/`，DeepSeek/GLM/Grok/Nemotron 与我主张原地保留（或作为可选项），需明确单选。
2. **`revision_20261116` 种子未来处置**：Grok 主张“永不改”，DeepSeek 隐含“核验后自动改”，我/GLM/Mistral/Nemotron 主张“核验后用户批准”，表述未完全统一。
3. **共享模块文件名**：DeepSeek/GLM/Mistral 用 `numerics.py` 或 `coloc_math.py`、`sampling.py` 或 `module_sampling.py`，需定稿。
以上任一条未获全体明确同意前，不得 `[收敛] YES`。
<<<END>>>