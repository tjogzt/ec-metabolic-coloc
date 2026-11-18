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