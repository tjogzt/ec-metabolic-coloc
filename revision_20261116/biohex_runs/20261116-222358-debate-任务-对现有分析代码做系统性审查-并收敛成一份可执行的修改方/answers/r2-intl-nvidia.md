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