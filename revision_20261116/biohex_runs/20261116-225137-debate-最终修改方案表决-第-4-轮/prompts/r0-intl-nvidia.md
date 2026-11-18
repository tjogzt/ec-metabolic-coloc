# 任务

# 最终修改方案表决（第 4 轮）

前 3 轮已充分交锋。以下是对全部争议的**定稿文本**，主 agent 已就两个程序分歧作出裁决。
请只回答：你是否同意这份定稿方案；如不同意，指出哪一条必须改、改成什么。
不要重开已解决的科学/范围问题，不要新增未核验事实。

## 裁决 1：MR 自检 seed=1 原地保留
`05_mendelian_randomisation/01_mr_estimators.py` 的 `if __name__ == "__main__"` 是合成自检。
保留 `default_rng(1)`、`TRUE=0.30`、`n=400`，**不迁移 tests/**，加英文横幅注明
“synthetic self-test; fixed fixture seed; not in the analysis pipeline; writes no results”。
理由：改为 config.SEED 会改变抽样并可能使断言失效；仓库是否有测试入口未核验；
新增 tests/ 与“减少文件数”目标冲突；横幅已实现隔离。这是规则问题，不再两选一。

## 裁决 2：revision_20261116 本轮零代码改动
`proxy_uncertainty.py` 的 `default_rng(20261116)` 是**有种子**随机数，不违反硬要求。
主 agent 已核实其输出写到 `revision_20261116/results/proxy_uncertainty.tsv`，
不写正式 `results/tables/`。本轮：不改代码、不改种子、不翻译注释，仅登记为
“遗留修订目录，不在正式 FLOW”。未来是否统一到 config.SEED：由用户在看到核验报告后决定；
计划中不预置“核验通过即自动改”。Grok “永久禁止”的表述不采纳（用户是最终裁决者），
但在本轮计划里等价于不改。

## 定稿方案

### P0 — 卫生修复（预期不改变已发表数字；黄金对照失败即停，不回写）
1. `07_figures/02_figure1_design.py` L169–171 的科学数字改为经中间变量从
   `results/tables/concordance_by_cohort.tsv` 读入。显式连接键：
   adipose↔GSE15773、diabetes↔GSE18732、clamp↔GSE13070；禁止按行序推断。
   舍入双规则：坐标 `f"{rho:.3f}"` 必须等于 0.313/0.202/0.098；标签 `f"{rho:+.2f}"`
   必须等于 +0.31/+0.20/+0.10。版式常量 74.0/83.0/92.0、15.0/16.6 保留。
   不得把原始 0.0976 画进图。
2. 所有正式 `.py/.R` 文件头**删除** `Author : Tao Zhu (zhutao@tjh.tjmu.edu.cn)`
   （删除，不替换成任何团体名）；R 中 “journal requirement” 改中性英文。
   范围：`code/00_setup`–`07_figures`、`config/config.py`、`config/config.R`。
3. `05_mendelian_randomisation/02_select_instruments.py` docstring 中的 “23.28%” 删除，
   不写入替代数字。
4. 阳性对照/功效模拟：`02_colocalisation/10_positive_control.py` 保留。主 agent 已核实
   `results/tables/` 下**没有** `coloc_power.tsv` 或 `positive_control_t2d.tsv`，故未混入正式表。
   执行时仍须在输出侧标注 simulated；不新增模拟路径；不用于文稿数字。
5. 增量分析冻结：主 agent 已定位 README 的 95% CI（−0.014, +0.031）对应
   `results/tables/incremental_boot.tsv` 的 `boot_lo=-0.01378`、`boot_hi=0.03086`，
   但该行 `delta_auc=-0.00156`，与 README 的 +0.0004 不是同一次 CV 实现；
   `skeleton_v3.tex` L744–745 已写明“点估计在交叉验证下变号”，README 没写。
   定稿：README 必须补上交叉验证变号说明，禁止把 +0.0004 与该 CI 说成同一次估计。
   在核对协议完成前，不重跑 `03_burden/07_incremental_value.py` 及任何 bootstrap/置换脚本。

### P1 — 结构与审计（不改数字）
6. 注释规范：正式代码统一英文；文件头说明用途/输入/输出，行间只注释非显然逻辑。
   已有中文注释的只有 `revision_20261116/` 三份脚本，按裁决 2 本轮不动。
7. 去重而不合并：新建 `code/common/`，按主题分文件：
   `locus_utils.py`（`_find_member`, `genes_near_loci`）、
   `numerics.py`（`_logsum`）、`io_utils.py`（`resolve_input`, `remote_size`）、
   `module_sampling.py`（`zmat`, `build_background`, `draw_module`，仅在通过门槛时）。
   抽取三重门槛：逐行文本等价 + 受影响脚本对现有表逐值回归（容差 0）+ 同种子 RNG
   调用序列一致。`zmat` 三函数默认不共享；不过门槛则保留两份并改名以暴露差异。
   `_parse_args` 不抽取。`00`–`07` 阶段脚本不合并、不合成巨文件。
   诚实说明：净增约 4 个文件；精简体现在消除重复行，不是文件数。
8. README：先对 `code/**` 全量 grep（path/accession/threshold/科学数字），再整句替换
   “The only constant fixed in source is the random seed”。新措辞须可验证，例如
   “分析路径内所有随机性由 config.SEED 或 --seed 控制；分析代码不内嵌结果数值”。
   主 agent 已核实 `config.py` 与 `config.R` 的 SEED 均为 49（一致，仅报告，不自动改）。
9. 核对协议：机器可读表（claim_id / 文本位置 / 声称值 / 舍入规则 / 源文件与列 /
   状态 match|rounding_match|mismatch|source_missing|unverified）。
   已核实可预填：concordance rho 0.313/0.202/0.0976 对 +0.31/+0.20/+0.10；
   conc_resp 61.62% 对 61.6%；incremental +0.000389 对 +0.0004、p=0.15、n=145、
   events=45；coloc 503、PP.H4>0.8 共 3、liver 0.931 对 0.93、muscle 0.962。
   **主 agent 新发现**：`06_EC_IR/paper_ms/skeleton_ebiomedicine.tex` L74/L293/L409
   写 51.9%（447/862, p=0.15）且 L74 写 p=0.29，与表（51.74%、p=0.1616）及
   `skeleton_v3.tex`（51.7%、446/862、p=0.16）不一致；标 mismatch，禁止把旧稿当投稿版。

### P2 — 重跑台账
10. 阶段 A（不依赖大文件）：import 检查、config 跨语言一致性、图脚本黄金对照。
    阶段 B（依赖本地数据）：缺文件即停并输出缺失清单，禁止模拟补洞。
    阶段 C：全部重跑输出到独立目录，与现有 `results/tables/` 逐值 diff；不一致只报告，
    不覆盖、不回写、不为了对上文稿改统计代码。
11. `revision_20261116/` 与 `06_EC_IR/` 探索脚本不纳入正式重跑；`06_EC_IR/paper_ms/`
    只做核对，不改代码。

## 表决要求
- 同意：`[立场] 支持` + `[收敛] YES`。
- 不同意：写明必须改的条目与改法；`[收敛] NO`。
- 不要因为“还有未执行的审计”而投 NO——审计是批准后的执行步骤，不是方案缺陷。
- 不要新增未在上文出现的文件、数字或文献。




# 你的任务

你是 **Nemotron 3.5 Super（nvidia/nemotron-3.5-super-vl）**，在 BioHex 六合团队中负责：数据分析 / 长上下文。
请**独立**给出你的最佳分析；不要假设其他成员会同意你，也不要等待他人。

严格按以下结构输出（不要添加其它章节）：

<<<ANSWER>>>
[结论] 一句话结论
[理由] 你的核心论据（可列 2–6 点）
[风险与不确定] 主要风险、混杂、证据缺口
[建议] 可执行下一步（可列 1–4 点）
<<<END>>>
