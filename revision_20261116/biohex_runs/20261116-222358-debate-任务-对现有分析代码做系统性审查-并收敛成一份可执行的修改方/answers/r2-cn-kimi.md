<<<ANSWER>>>
[立场] 修正：接受“默认不共享、等价才抽取”的严格门槛；MR 合成自检保留 seed=1 作为测试 fixture 并隔离；revision_20261116/ 默认不改动；图读表必须带舍入黄金对照；共享模块按主题拆分；本轮仍不收敛。

[理由]
1. 硬要求的优先级应是“不改变已发表数字”高于“种子统一”。因此，任何会改变随机流、函数输出或数值的改动，必须先证明其输出未进入 results/tables/ 或文稿，并通过回归测试。
2. `code/07_figures/02_figure1_design.py` 中的科学数字（rho 与标签）必须从 `results/tables/concordance_by_cohort.tsv` 读入中间变量；版式坐标（74.0/83.0/92.0、15.0/16.6 等）不是科学结果，可保留。读表后必须复现现有图上的 0.098 与标签 “+0.31/+0.20/+0.10”；若不能复现，停止修改。舍入规则建议锁定为：y 坐标显示三位小数（`{rho:.3f}`），标签显示两位带符号（`{rho:+.2f}`）。
3. 身份信息清除范围：删除所有正式 `.py/.R` 文件头的作者姓名、邮箱；R 中 “journal requirement” 改为中性英文；不引入新的团体名。`06_EC_IR/paper_ms/skeleton_ebiomedicine.tex` 不在正式包范围，本轮不动。
4. 共享模块抽取须设双重门槛：先逐行 diff，再对现有 `results/tables/` 做数值回归（容差 0）。建议按主题拆分为 `code/common/{locus_utils.py, numerics.py, io_utils.py, module_sampling.py}`，禁止单一 `shared_utils.py` 混装异质函数。`zmat`/`build_background`/`draw_module` 默认保留两份，不等价则改名以暴露差异；`_logsum`、`_find_member`/`genes_near_loci`、`resolve_input`/`remote_size` 同样处理。`_parse_args` 不抽取。
5. 随机数三类裁决：统计推断程序（bootstrap/置换/模块抽样）保留；`10_positive_control.py` 保留但输出隔离、标注 simulated、不扩展、不填缺失数据；`01_mr_estimators.py` 的 `__main__` 合成自检是单元测试 fixture，保留 `seed=1`、`TRUE=0.30`、`n=400`，仅移到 `tests/` 或显式标注 synthetic self-test，不写 results，不改种子；`revision_20261116/` 不在正式包，隔离并审计 provenance，默认不改种子，若用户后续纳入再按两步门槛处理。
6. README 的 “only constant fixed in source is the random seed” 声明已被证伪，必须重写；但重写须在全仓 path/accession/threshold 硬编码审计完成后进行，不能提前写另一句绝对声明。
7. 重跑计划分三阶段：A. 不依赖大文件的自检/单元测试；B. 依赖本地数据的阶段，先跑输入清单预检，缺文件即停并报告；C. 禁止用模拟数据补洞。完整重跑受原始数据可用性约束，计划必须明确这一点。
8. 核对协议：建立机器可读对照表（claim_id / 文本位置 / 声称值 / 舍入规则 / 源 tsv 文件列 / 状态）。`README` 中 95% CI（−0.014 to +0.031）标 `source_missing`，定位来源文件前禁止重跑 `03_burden/07_incremental_value.py` 等 bootstrap 脚本，以免覆盖现有结果。`config.py` 与 `config.R` 的 SEED 一致性未核验，须加审计。

[异议]
- 对 DeepSeek：单一 `code/shared_utils.py` 把位点逻辑、数值、IO、抽样混装，是用新耦合换旧重复；未处理 0.0976 → 0.098 → +0.10 的舍入冲突；提议用 “BioHex Consortium” 替换作者行是引入新身份信息；将 `skeleton_ebiomedicine.tex` 改名超出本轮范围；对 `revision_20261116/` 改种子未设“是否进入文稿”的隔离门槛。
- 对 GLM：基本赞同其“两步门槛”和主题分文件方案；但“revision 证明隔离后改种子”若写入本轮计划，应列为 P2 且仅在用户明确扩大范围后执行，因为该目录不在正式可复现包内。
- 对 Mistral：主张将 MR 自检 seed=1 改为 `config.SEED` 并移入 `tests/`；我们认为 seed=1 是测试 fixture 的预言常量，改种子无收益且破坏测试确定性，移动入口可接受但不得改种子。
- 对 Grok：赞同其等价性测试、舍入黄金对照、README 后置改写等严格立场；但“分阶段重跑是请求用户改写硬要求 3”过于绝对——在原始数据不在 git 的前提下，唯一诚实的执行方案是预检输入+缺数据停报，这仍是执行而非削弱“禁止模拟数据”的硬要求。
- 对 Nemotron：将 `revision_20261116/proxy_uncertainty.py` 改种子列为 P0 不可接受，因为该脚本范围与 provenance 均未核验，改种子可能改变已发表数字；未处理图内舍入冲突；对 README liver 0.93 与表 0.931 的质疑不成立，二者是四舍五入关系。

[修正]
- 撤回“MR 自检 `seed=1` 应改为 `config.SEED`”的上一轮立场，改为保留 seed=1 作为 fixture，仅隔离/标注。
- 撤回“同名函数 diff 后即可抽取共享”的立场，改为“默认不共享，diff + 数值回归一致才抽取，否则改名保留两份”。
- 明确接受“图读表必须带黄金对照、先锁定舍入规则”。
- 明确 `revision_20261116/` 不在本轮正式修改范围，默认隔离审计，不改动种子。
- 将共享模块组织从单一 `_shared/` 改为按主题 `code/common/`。

[收敛] NO
P0 卫生项（硬编码结果数字读表、身份信息清除、README 声明修正）已趋一致，但 revision 种子处置、MR 自检种子、zmat 类函数共享门槛三处，其他成员上一轮仍有“必改/禁改/默认不共享”等对立表述；我方本轮修正后尚未获得他们确认，故不能判定实际到场全体成员已达成一致。
<<<END>>>