<<<ANSWER>>>
[结论] 当前代码包在可复现性、数据真实性、代码整洁性上存在多处与用户硬性要求及 README 声明不符的缺陷，必须进行系统性修改。核心问题包括：硬编码结果数字、身份信息泄露、重复函数未抽取、README 声明不准确。随机数使用整体合理，但需明确分类并文档化。需设计严格的数据核对与分阶段重跑协议，且修改不得改变已发表数字。

[理由]
1. **硬编码结果数字直接违反数据真实性要求**：`07_figures/02_figure1_design.py` 将三个队列的 rho 值（0.313, 0.202, 0.098）硬编码在作图代码中，而非从 `results/tables/concordance_by_cohort.tsv` 读入。即使当前数值与表一致，未来任何重跑或数据更新都会导致代码与结果脱节，这是不可接受的腐化风险。`05_mendelian_randomisation/02_select_instruments.py` 的 docstring 中写死 “23.28%” 同样属于结果数字嵌入源码，必须清除。
2. **身份信息必须清除**：几乎所有正式脚本文件头包含作者姓名和邮箱，违反匿名化要求。`07_figures/01_main_figures_1_2_3.R` 中的 “journal requirement” 虽无具体刊名，但为彻底起见应改为通用描述。`06_EC_IR/paper_ms/skeleton_ebiomedicine.tex` 文件名含期刊名，虽非代码，但建议重命名以完全去标识化。
3. **重复函数增加维护成本且易导致不一致**：跨文件的同名函数（如 `_find_member`, `_logsum`, `zmat` 等）显然是复制粘贴，应抽取到共享工具模块（如 `code/shared_utils.py`），而非合并阶段文件。合并会破坏 `docs/FLOW.md` 定义的运行顺序，必须反对。
4. **README 声明与事实不符**：README 声称 “No analysis script contains a hard-coded path, accession or threshold. The only constant fixed in source is the random seed.” 但 `02_figure1_design.py` 硬编码了结果数字，`02_select_instruments.py` 的 docstring 含结果数字，`revision_20261116/` 中种子未走 `config.SEED`。该声明必须修正为准确描述。
5. **随机数使用整体合理，但需明确分类**：统计程序（bootstrap、置换）使用 `config.SEED` 或 `--seed`，完全可复现，应保留。功效模拟（`10_positive_control.py`）基于真实数据结构模拟效应量，属于方法学验证，应保留但需在文档中明确其模拟性质，且不纳入主要结果。合成自检（`01_mr_estimators.py` 的 `__main__`）写死种子为 1，仅用于单元测试，不输出结果，保留无碍。R 中的 `set.seed(SEED)` 正常。
6. **数据核对与重跑计划缺失**：当前无系统性核对协议，增量价值 CI 来源不明（`incremental_value.tsv` 中无 CI 列），重跑约束未定义。必须设计分阶段重跑计划，优先运行不依赖大数据的自检，再运行依赖本地数据的脚本，缺文件时明确报告，严禁用模拟数据替代。

[风险与不确定]
- 修改 `02_figure1_design.py` 从 tsv 读取数据时，需确保读取逻辑与表格结构完全匹配，否则可能引入绘图错误。当前硬编码值与表格一致，修改后不应改变已发表数字，但需验证。
- 提取共享函数时，需逐一核对同名函数是否真正等价（实现细节可能因上下文微调），若存在差异则不能强行合并，需保留各自版本并重命名以避免混淆。
- 清除身份信息可能涉及修改大量文件头，需确保不误删必要的版权或许可证声明。
- 重跑计划依赖本地数据文件，若 `data/` 目录中缺少某些 accession 的原始数据，部分阶段将无法执行，需提前列出所需文件清单并评估可获取性。
- `incremental_value.tsv` 中缺少 CI，README 中的 CI 可能来自 bootstrap 结果文件（如 `incremental_boot.tsv` 或类似），需定位确切来源，否则无法完成核对。

[建议]
**P0（必须立即修改，影响结果真实性与合规）**
- **修改 `code/07_figures/02_figure1_design.py`**：移除 L169–171 硬编码的 PTS 列表，改为从 `results/tables/concordance_by_cohort.tsv` 读取 rho 和 cohort 信息。预期效果：消除硬编码结果数字，确保图形与数据表动态绑定。不会改变已发表数字（当前值一致）。
- **清除所有正式脚本文件头中的身份信息**：删除 `Author: Tao Zhu (zhutao@tjh.tjmu.edu.cn)` 行，可替换为项目名或 “BioHex Consortium”。涉及约 45 个文件。预期效果：满足匿名化要求。
- **修正 `README.md` 中的不实声明**：将 “No analysis script contains a hard-coded path, accession or threshold” 改为准确描述，例如 “Analysis scripts read all data from external files; no result numbers are hard-coded in source (with the exception of layout coordinates in figures).” 预期效果：文档与代码实际状态一致。
- **移除 `code/05_mendelian_randomisation/02_select_instruments.py` docstring 中的结果数字**：将 “23.28% of its genome-wide significant rows” 改为引用表格或删除具体百分比。预期效果：避免注释中的过时结果数字。

**P1（应尽快修改，提高可维护性与可复现性）**
- **提取重复函数到共享模块**：创建 `code/shared_utils.py`，将已验证等价的 `_find_member`, `genes_near_loci`, `_logsum`, `resolve_input`, `remote_size`, `zmat`, `build_background`, `draw_module` 移入，并更新所有引用文件的 import。注意：`_parse_args` 虽同名但可能因参数不同而不应合并。预期效果：减少代码重复，降低维护成本。不会改变结果。
- **统一注释为英文，清除中文注释**：`revision_20261116/` 中的三个脚本含中文注释，应翻译为英文。预期效果：符合国际化规范。
- **设计并执行数据核对协议**：列出所有 README 和 tex 文稿中的关键数字（如 rho, PP.H4, ΔAUC, CI 等），映射到对应的 `results/tables/` 文件及行列，由主 agent 逐项核对。特别需定位增量价值 CI 的来源文件（可能为 `incremental_boot.tsv` 或 `03_burden/07_incremental_value.py` 的输出）。预期效果：确保文稿数字与代码输出完全一致。
- **制定分阶段重跑计划**：
  - 阶段 1（自检）：运行所有含 `if __name__ == "__main__"` 的单元测试（如 `01_mr_estimators.py`），验证基础函数。
  - 阶段 2（小数据验证）：运行不依赖大型原始数据的脚本（如 `02_colocalisation/12_permutation_null.py` 若仅需汇总统计，或 `03_burden/05_build_programme.py` 若输入为已处理中间文件）。
  - 阶段 3（全数据重跑）：在确认本地数据文件齐全后，按 `docs/FLOW.md` 顺序完整重跑。缺文件时生成缺失报告，严禁用模拟数据替代。预期效果：验证代码可运行性，发现潜在依赖问题。

**P2（可选改进）**
- **重命名 `06_EC_IR/paper_ms/skeleton_ebiomedicine.tex`**：移除期刊名，改为通用名如 `skeleton_main.tex`。预期效果：完全去期刊标识。
- **审查 `revision_20261116/` 脚本的随机种子**：将 `proxy_uncertainty.py` 中的 `np.random.default_rng(20261116)` 改为使用 `config.SEED` 或通过命令行参数传入，以统一可复现性标准。

**明确反对的做法**
- 反对为了减少文件数而合并不同 FLOW 阶段的脚本，这会破坏模块化和运行顺序。
- 反对删除任何用于统计推断的 bootstrap/置换随机数，它们是标准方法，且种子固定可复现。
- 反对将示意图数字改为从错误的表读取，必须严格从对应的结果表读入。
- 反对在重跑时因数据缺失而改用模拟数据，必须如实报告缺失并停止。
<<<END>>>