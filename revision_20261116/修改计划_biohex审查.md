# BioHex 六方代码审查 — 修改计划（待批准）

日期：2026-11-16
审查对象：`ec-metabolic-coloc/` 正式分析包（45 个 .py/.R，约 14,520 行）
文稿：`06_EC_IR/paper_ms/`、`01_manuscript.pdf`、`supp_a1a4.pdf`
状态：**未修改任何分析代码**；仅执行了语法/导入验证。

## 一、团队运行记录

- 模式：BioHex debate（6 方，OpenRouter）。
- 第一轮辩论：3 轮交叉质证，6/6 每轮到场；分歧集中在 2 个程序问题。
- 裁决后第二轮表决：2 轮，6/6 全部 `[收敛] YES`。
- 缺席：无。未用三方冒充四方。
- 档案：`revision_20261116/biohex_runs/`（transcript、raw、leader-synthesis）。
- 两处主 agent 裁决：
  1. MR 自检 `seed=1` 原地保留，加英文横幅，不迁 `tests/`。
  2. `revision_20261116/` 本轮零代码改动，处置权交用户。

## 二、修改计划（按优先级）

| # | 修改项 | 文件 | 理由 | 优先级 | 预期效果 | 风险控制 |
|---|---|---|---|---|---|---|
| 1 | 图 1 示意图的科学数字改为从结果表读入 | `code/07_figures/02_figure1_design.py` L169–171 | 当前把 0.313/0.202/0.098 与 `+0.31/+0.20/+0.10` 直接写在绘图代码里，违反“数据先入中间变量” | P0 | 图与 `concordance_by_cohort.tsv` 强绑定，改表即改图，不再手抄 | 显式键映射 adipose↔GSE15773、diabetes↔GSE18732、clamp↔GSE13070；坐标 `f"{rho:.3f}"`、标签 `f"{rho:+.2f}"`；黄金对照不符即停 |
| 2 | 给图脚本加 `main()` / `if __name__ == "__main__"` 保护 | `02_figure1_design.py` | 当前模块级代码在 import 时直接作图（本次审查已触发一次重写，见第五节） | P0 | 导入不再产生副作用，便于测试与复用 | 不改任何绘图数值 |
| 3 | 删除作者姓名与邮箱 | `code/00_setup`–`07_figures`、`config/config.py`、`config/config.R` 共 42 个源文件 | 身份信息不应进入代码 | P0 | 匿名化；符合投稿/公开要求 | 纯删除，不替换为团体名；不碰统计逻辑 |
| 4 | 删除注释中的结果数字 “23.28%” | `code/05_mendelian_randomisation/02_select_instruments.py` L59 | 注释硬编码结果，且不可溯源 | P0 | 注释不再携带会过期的结果 | 只删除，不填替代数字 |
| 5 | 增量分析：README 补“交叉验证点估计变号”说明，并拆分点估计与 CI 两条 claim | `README.md` | README 把 +0.0004 与 CI −0.014/+0.031 并列，但该 CI 来自 `incremental_boot.tsv` 中 delta=−0.00156 的另一实现；`skeleton_v3.tex` 已写明变号，README 没写 | P0 | 避免把两次 CV 实现说成同一次估计 | 核对完成前冻结重跑 `07_incremental_value.py` 与所有 bootstrap/置换脚本 |
| 6 | 抽取真正重复的函数到 `code/common/` | 新建 `locus_utils.py`、`numerics.py`、`io_utils.py`、`module_sampling.py` | 9 对同名函数中，部分逐行重复；用户要求减少腐化 | P1 | 消除重复行，单一实现 | 三重门槛：逐行文本等价 + 影响脚本对现表逐值回归（容差 0）+ 同种子 RNG 调用序列一致；`zmat`/`build_background`/`draw_module` 默认不共享；`_parse_args` 不抽；阶段脚本 00–07 不合并。**诚实说明：净增约 4 个文件，精简体现在重复行而非文件数** |
| 7 | 重写 README 的绝对声明 | `README.md` L72–73 | “No analysis script contains a hard-coded path, accession or threshold. The only constant fixed in source is the random seed.” 已被证伪 | P1 | 声明与实际一致、可 grep 验证 | 先全量 grep 再改写；新措辞限定“正式 FLOW / 产出 results/tables 的路径”，排除 MR 自检与 revision 遗留种子 |
| 8 | 建立机器可读的文稿-结果对照表 | 新建 `revision_20261116/claim_crosscheck.tsv` | 第 4 项要求逐项核对 | P1 | 每个文稿数字可定位到文件与列 | 字段：claim_id、文本位置、声称值、舍入规则、源文件与列、状态 |
| 9 | 统一注释规范 | 正式包全部 `.py/.R` | 第 5 项要求 | P1 | 文件头写用途/输入/输出，行间只注释非显然逻辑 | 正式包当前已是英文注释；中文注释仅在 `revision_20261116/`，按裁决本轮不动 |
| 10 | 分阶段重跑 | 全包 | 第 3 项要求 | P2 | 证明全部代码可运行 | 阶段 A 无数据自检；阶段 B 缺文件即停并报清单；阶段 C 输出到独立目录并 diff，不覆盖现表；禁止模拟补洞 |
| 11 | 单细胞数字补溯源 | `README.md`、`docs/FLOW.md` | 141,316 cells / 29,592 epithelial / z=−3.8 只在文档中，`results/tables/` 无对应表 | P2 | 文稿数字可溯源 | 先找运行日志；找不到则重跑或有依据地降级表述 |
| 12 | 处理旧稿 `skeleton_ebiomedicine.tex` | `06_EC_IR/paper_ms/` | 其中 51.9%（447/862, p=0.15）与 p=0.29 同现表（51.74%, p=0.1616）及 `skeleton_v3.tex`（51.7%, 446/862, p=0.16）冲突 | P2 | 防止误投旧稿 | 本轮只标记，不改文稿；是否归档/改名由用户定 |

## 三、数据真实性审查结论

### 3.1 随机数
- 未发现无种子随机数。全部 13 处随机数均来自 `np.random.default_rng(seed)` 或 R `set.seed(SEED)`。
- 三类用途：
  1. 统计推断（bootstrap/置换/模块抽样）：保留，种子来自 `config.SEED` 或 `--seed`。
  2. 阳性对照/功效模拟（`10_positive_control.py`，仅模拟效应量，位点/等位基因/SE 用真实数据）：保留，输出未混入正式表，执行时标注 simulated。
  3. 合成自检（`01_mr_estimators.py` 的 `__main__`，seed=1）：保留，加横幅隔离。
- `revision_20261116/proxy_uncertainty.py` 使用 `default_rng(20261116)`，有种子，不违规；输出写到该目录 `results/`，不写正式表。

### 3.2 硬编码数字
- 结果数字：1 处（图 1 PTS）、1 处注释（23.28%）。
- 门限/参数：多数在 `config.PARAMS`；但 `05_coloc_abf_engine.py` 的 `P1/P2/P12_DEFAULT`、`06_cross_cancer/03` 的 `COUPLING_FRACTION=0.8`、`02_load_tcga` 的 `MIN_MEAN_EXPRESSION=0.5`、`01_concordance` 的 `DILUTION_NOMINAL=1e-3` 等仍是脚本常量。它们是算法参数而非结果数据；是否统一进 config 由用户决定（统一不会改变数值，但会动多个文件）。
- 硬编码 accession：`06_second_outcome.py` 的 `SECOND_OUTCOME_NAME="GCST006465"`、`06_load_programme.py` 的 `ATLAS_SERIES`、`10_positive_control.py` 的 canonical gene panel。前两者与 config 重复，建议 P2 改为从 config 派生；基因面板属于专家先验，建议保留并注明为常量。

### 3.3 数据来源与中间变量
- 正式分析路径的数据均先读入 DataFrame/变量再传递，未发现把原始矩阵内嵌代码。
- 例外即第 3.2 节列出的结果数字与 accession。

## 四、数据核对（已执行，真实表值）

| 文稿/README 声称 | 结果文件与值 | 状态 |
|---|---|---|
| 503 gene–locus tests | `coloc_all_tests.tsv` 503 行 | match |
| PP.H4>0.8 共 3 | 3 行（Liver 1，Muscle 2） | match |
| PP.H4>0.5 共 5 | 5 行 | match |
| Liver max PP.H4 0.93 | 0.930889 | rounding_match |
| Muscle max PP.H4 0.962 | 0.962160 | match |
| H3>0.8 共 201 | 201 | match |
| 862-gene programme（608 up/254 down） | `metabolic_program_v2.tsv` 862 行 | match |
| rho +0.31/+0.20/+0.10 | 0.3133 / 0.2022 / 0.0976 | rounding_match |
| responsive 61.6% | 61.62%（99 genes） | rounding_match |
| ΔAUC +0.0004，p=0.15，n=145，events=45 | `incremental_value.tsv` 0.000389, 0.15, 145, 45 | rounding_match |
| 95% CI −0.014 to +0.031 | `incremental_boot.tsv` boot_lo −0.01378 / boot_hi 0.03086，但同行 delta=−0.00156 | **源值存在，但与小节配对不当**；需补变号说明 |
| SKAP1 −0.25 (liver/pancreas), −0.24 (spleen), ovary n.s. | −0.2495 / −0.2482 / −0.2394 / +0.0968 (p=0.33) | match |
| 141,316 cells；29,592 epithelial；z=−3.8 | `results/tables/` 无对应表 | source_missing |
| 51.9%（447/862, p=0.15）；p=0.29 | 现表 51.74%（p=0.1616）；`skeleton_v3.tex` 51.7%（446/862, p=0.16） | **mismatch（旧稿）** |

结论：当前主稿 `skeleton_v3.tex` 与结果表一致；README 有一处配对表述需修正；旧稿 `skeleton_ebiomedicine.tex` 不应作为投稿版。

## 五、运行验证（已执行部分）

- Python 全包 `compileall`：通过。
- R 全包 `parse`：5/5 通过。
- 模块导入：36/36 通过（需用 `/opt/anaconda3/bin/python3.12`；系统默认 `python3`=3.10 的 scipy 动态库损坏）。
- 输入检查：`01_check_inputs.py` 报告 **全部原始数据缺失**（GWAS、GTEx、五个 eQTL 组织、三个 bulk 队列）。因此完整重跑当前不可执行；需先取得数据（可由用户手动下载或授权下载）。
- **需披露的副作用**：导入探针触发了 `02_figure1_design.py` 的模块级作图，覆盖了 `results/figures/fig1_design.pdf` 与 `.png`（仅这两个文件；结果表未动）。新图由当前代码确定性生成，但旧版本已被覆盖、无备份。旧 PNG 副本 `results/figures/png/fig1_design.png`（11-08）分辨率不同，不能逐字节还原。请决定：保留新生成图（与当前代码一致），还是由你从备份恢复。

## 六、需要你决定的事项

1. 是否批准上述 P0/P1/P2 计划，按序执行。
2. 是否授权下载缺失的公共数据以完成完整重跑；若否，重跑只能覆盖不依赖原始数据的阶段。
3. 是否允许新建 `code/common/`（净增约 4 个文件）以消除重复，而不是合并阶段脚本。
4. 是否允许修改 `README.md` 的绝对声明与增量 CI 表述（不改任何结果数字）。
5. `fig1_design.pdf/png` 被重新生成，保留还是恢复。
6. `skeleton_ebiomedicine.tex` 旧稿：归档、改名，还是删除。

---

## 七、执行记录（2026-11-17）

批准范围：P0 全部 + 新建 `code/common/` + 核对表。
未执行：README 改写、脚本阈值统一（未获勾选）；完整重跑（等待数据）。

### 已完成的代码改动

| 文件 | 改动 |
|---|---|
| `code/07_figures/02_figure1_design.py` | 三个 rho 与标签改为从 `results/tables/concordance_by_cohort.tsv` 读入；新增显式键映射与黄金断言；加入 `main()` 与 `__main__` 保护；文件头说明同步更新 |
| `code/02_colocalisation/05_coloc_abf_engine.py` | 删除本地 `_logsum`，改用 `common.numerics.logsum` |
| `code/02_colocalisation/08_prior_sensitivity.py` | 同上 |
| `code/02_colocalisation/10_positive_control.py` | 删除本地 `resolve_input`、`remote_size`，改用 `common.io_utils` |
| `code/02_colocalisation/15_locus_screen.py` | 同上；删除重复的本地 `resolve_input` |
| `code/05_mendelian_randomisation/02_select_instruments.py` | 删除注释中的 “23.28%” |
| `code/**`、`config/config.py`、`config/config.R`（共 42 个源文件） | 删除作者姓名与邮箱；R 中 “journal requirement” 改为中性英文 |
| `code/common/__init__.py`、`io_utils.py`、`numerics.py` | 新增；仅容纳通过等价门槛的函数 |

说明：因 `code` 与 Python 标准库 `code` 模块重名，导入方式为把 `code/` 加入
`sys.path` 后 `from common.io_utils import ...`，物理目录仍为 `code/common/`。

### 未抽取的函数（未过门槛）

- `zmat`、`build_background`、`draw_module`：虽 `zmat`/`draw_module` 的 AST（去 docstring）等价，
  但缺少本地数据，无法做“零容差输出回归”，按计划默认不共享。
- `_find_member`、`genes_near_loci`、`_logsum`（原两份）实现不同，保留合并后的单一实现于 common。
- `_parse_args`：保留两份，不抽取。

### 验证结果

- Python `compileall`：通过（全包）。
- R `parse`：5/5 通过。
- 模块导入：39/39 通过；导入图脚本不再产生任何文件写入。
- 图 1 重渲染：PNG 与改前逐字节一致；PDF 仅元数据差异。
- `logsum` 数值回归：对空数组、常数、负大值、正负混合、含 inf 五类输入，
  与改前两份实现逐值完全一致。
- `resolve_input`：相对路径与绝对路径行为与改前一致。
- `results/tables/`：未被本次改动写入；最新文件时间仍为 2026-11-16。

### 新发现的文稿/文档不一致（已写入 claim_crosscheck.tsv）

- C27/C28：README 的两处绝对声明仍不成立（未获准改 README）。
- C31：`data/README.md` 指向不存在的 `code/01_programme/`。
- C32：`01_check_inputs.py` 期望 `GTEx_Analysis_v8_eQTL.tar`，`data/README.md` 列的是 v11 SuSiE tar。
- C16：增量分析 CI 与点估计来自不同 CV 实现，README 未写变号说明。
- C24/C25：旧稿 `skeleton_ebiomedicine.tex` 的 51.9%/447/p=0.15 与 p=0.29 与现表不符。

### 回滚

改前代码与配置已打包在
`revision_20261116/backup_20261117/code_config_before_biohex.tgz`。

### 等待你下载的数据（放入 `ec-metabolic-coloc/data/`）

必需的三个文件：
1. `GCST90454186.h.tsv.gz`（GWAS Catalog harmonised hg38，约 373 MB）
2. `GCST006465_GRCh37.tsv.gz`
3. `GTEx_Analysis_v8_eQTL.tar`（与 `01_check_inputs.py` 一致；若要按 data/README 用 v11 SuSiE，
   需先统一检查脚本与文档，见 C32）

eQTL 五组织可用 `--remote` 走 tabix，不必本地存；bulk 三个队列
（GSE15773、GSE18732、GSE13070）与增量模型所需的 GSE120490、平台注释也需要本地文件。
你下载完成后，我先跑 `01_check_inputs.py` 校验，再按阶段 A→B→C 重跑并逐值 diff。

---

## 八、重跑记录（2026-11-17）

### 数据

- GWAS：`GCST90454186.h.tsv.gz`（md5 与 EBI 一致）、`GCST006465_GRCh37.tsv.gz`（md5 一致）。
- GTEx：`GTEx_Analysis_v8_eQTL.tar`（1.49 GB，流式读取，无需解压）。
- eQTL：五组织通过 `--remote` 读取；远端文件列序与 `OUT_COLS` 不同。
- Bulk：GSE15773/18732/13070/120490/24215/24422/17025/196033 + GPL570/GPL6480 注释。
- MR 暴露：`GCST90002238.tsv.gz`、`GCST006900_bmi.txt.gz`。
- TCGA：从 `/Volumes/tjogzt4T/data/xena/tcgapancan` 复制表达与表型矩阵。
- 单细胞：未取得 `ec_integrated_v1.h5ad`。

### 跑通并复现的层

| 层 | 结果 | 与随包表 |
|---|---|---|
| Colocalisation | 503 检验，PP.H4>0.8 共 3，PP.H3>0.8 共 201；五组织计数 30/57/171/143/102 | 完全一致 |
| Prior sensitivity | p12=1e-6…1e-4 对应 H4>0.8 为 0/3/3/5/7 | 与 README“0 到 7”一致 |
| Incremental（显式 v2） | ΔAUC +0.000389，p=0.15，n=145，events=45；CI −0.01378/+0.03086 | 完全一致 |
| MR | 8 个估计、条件 F、第二结局数值全部一致 | 数值一致，行序不同 |
| Dual outcome | 17q21.32 仍为共享信号：Liver SKAP1 PP.H4=0.922，Muscle 0.874；2 个 >0.8 | 结论一致；检验数 352 vs 随包 294，属远端 eQTL 版本差异 |
| Figures | fig2：3/503、5 个 >0.5、55 个 H3>0.8；fig3：+0.0004/−0.0038 | 一致 |

### 重跑中发现并修复的代码缺陷

1. `03_eqtl_catalogue_reader.py`：远端表头解析；CRLF 尾字符；tabix 非零退出原先被当成空窗口，现改为三次重试后报错。
2. `07_run_colocalisation.py`：`--remote` 且无本地文件时按数字位置取列，列序错位导致 0 条可用记录；改为从远端表头解析列位置；tabix 失败重试。
3. `09_dual_outcome.py`：用 hg19 的 `pos19` 对比 hg38 的 `lead`，chr17 偏移约 1.9 Mb，17q21.32 被整体丢弃；改为按 `lead19` 取窗。另将远端列位置改为从远端表头解析。

第一次 adipose 远端抓取因 tabix 静默失败只得到 122 个检验；修复并重抓后为 143，总数回到 503。

### 未跑通、已定位原因的层

| 层 | 原因 | 证据 |
|---|---|---|
| 862-gene programme 构建 | 正式包 `05_build_programme.py` 只生成 146 基因的胰岛素应答程序；862 基因程序由 `06_EC_IR/w4/w4_program_v2.py` 从 GSE25401 生成，GSE25401 未列入 data/README | C34 |
| Concordance rho 重算 | `DE_GSE13070/18732/15773` 既未随包提供，也无正式包代码生成；生成脚本只在 `06_EC_IR/w4` | C35 |
| Cross-cancer | 提供的 TCGA Xena 矩阵行名是 Entrez ID，`01_map_gene_symbols.R` 要求 HUGO symbol；symbol 版镜像当前不可达 | C41 |
| Single-cell | `ec_integrated_v1.h5ad` 不存在；本地 `atlas_core_scored.h5ad` 为 140,812 细胞，与文稿 141,316 不同 | C42 |
| 组合表 | 无脚本组装 `coloc_all_tests.tsv`；本次由主 agent 从五个 `coloc_<tissue>.tsv` 拼接验证 | C43 |

### 结论

- 已发表的核心数字中，共定位、先验敏感性、增量分析、MR、主图均可由现有代码和公开数据复现。
- 正式包不能自足重建 862 基因程序与三个 concordance DE 表；这两项需要把 `06_EC_IR/w4` 的生成脚本和 GSE25401 纳入正式包，属于新增代码，待你批准。
- 交叉癌与单细胞层需要替换输入格式或补回原始衍生对象，暂不能声称已复现。

---

## 九、5–7 项修复记录（2026-11-17）

### 5. 文档表述
- `README.md`：增量结果补充“交叉验证重复时点估计变号”，明确这是“无证据”而非“等价”。
- `README.md`：删除“唯一固定常量是随机种子”和“无硬编码 path/accession/threshold”两处绝对声明，改为可验证表述：分析路径随机性由 `config.SEED`/`--seed` 控制；结果数字从表读取；自检 fixture 与遗留 revision 种子单独列出；具名算法常量属于参数。
- `data/README.md`：补入 stage 02 必需的 `GTEx_Analysis_v8_eQTL.tar`；保留 v11 SuSiE 作为 fine-mapping 输入；修正不存在的 `code/01_programme/` 指向，并说明 862 表由遗留 `w4_program_v2.py` 生成、正式包只消费锁定表。
- `data/README.md`：补记远端 eQTL 不按版本冻结，第二结局检验数可随 release 变化，但 17q21.32 结论稳定。
- `docs/FLOW.md`：Stage 03 写明默认使用 v2 程序与模块；v1 需显式传参；说明 bootstrap 表的 `delta_auc` 列不是主点估计。

### 6. 组合表
- `07_run_colocalisation.py` 新增 `write_combined()`：每次写完单个组织后自动重建 `coloc_all_tests.tsv`，并打印覆盖组织数。
- 验证：5/5 组织、503 行、3 个 PP.H4>0.8。

### 7. 默认值
- `07_incremental_value.py` 默认改为 `metabolic_program_v2.tsv` + `random_modules_v2.tsv`。
- 验证：默认命令直接得到 `delta_auc_IR=0.000389`、p=0.15，与文稿一致；v1 仍可通过显式参数复现。

核对表中 C16、C27、C28、C31、C32、C33、C43 已标为 `RESOLVED`。

---

## 十、1–4 项处理记录（2026-11-17）

### 1. 862 程序生成器
- 新增 `code/03_burden/02b_state_programme_de.R`：从 GSE15773、GSE18732、GSE13070、GSE25401 重建四个 DE 表。
- 新增 `code/03_burden/05b_build_programme_v2.py`：从 GSE25401 DE 重建 862 程序与 100 套等大随机模块。
- 四个 DE 表与 `06_EC_IR/w4` 原表逐值完全一致（最大绝对差 0）。
- **历史版本差异**：当前 GSE25401 DE 表时间戳为 2026-11-05，晚于锁定程序（2026-11-04）。按同一规则得到 **894** 基因（新增 40、移除 8），不是锁定的 862。`06_EC_IR/w4/dryrun_894.py` 已记录这 48 个基因的变动。
- 生成器已加保护：与锁定表不一致时写 `metabolic_program_v2_rebuild.tsv` / `random_modules_v2_rebuild.tsv`，打印差异并以退出码 2 结束，**不覆盖锁定表**。
- 结论：生成流程已移植；锁定 862 的精确重建需要 2026-11-04 版 DE，该文件已被覆盖，无法恢复。

### 2. Concordance DE 生成器
- 同一 `02b_state_programme_de.R` 生成 `DE_GSE15773_adipose_IR_state.tsv`、`DE_GSE18732_muscle_DM_state.tsv`、`DE_GSE13070_muscle_baseline_IR.tsv`。
- `01_concordance_pipeline.py` 重跑结果与随包表逐值完全一致：rho +0.0976 / +0.2022 / +0.3133，一致率与 CI 全部相同。

### 3. Cross-cancer 映射
- `01_map_gene_symbols.R` 现在从 `org.Hs.eg.db` 生成 `data/ENTREZ2SYMBOL.tsv`。
- `02_load_tcga_expression.py` 在流式读取时把 Entrez 映射为 symbol，并按 symbol 取均值折叠重复。
- 四步全部跑通（load 16,312 基因 × 10,441 样本；30 个癌种、13 个正常组织）。
- 与附录差异：可评分程序 **793/862**（附录 770/862）；肿瘤中位 rho **−0.0285**（附录 −0.029）；>0.4 为 **3/30**（附录 4/30）；正常中位 +0.473（附录 +0.483）。方向与结论一致，属 symbol 矩阵重建的近似。

### 4. 单细胞层
- 精确的 `ec_integrated_v1.h5ad` 已找到并复制到 `data/`：141,316 细胞 × 36,601 基因。
- `01_project_programme.py`：IR 差值 −0.059，MWU p=0.109。
- `02_null_calibration.py`：**z = −3.803**（文稿 −3.8），经验 p=1.000，与文稿一致。
- 上皮细胞数：代码为 **34,625**；README 原写 29,592，已更正为 34,625（与 `skeleton_v3.tex` 一致）。

### 最终状态
- 所有分析层现在均可运行：共定位、先验敏感性、增量、MR、单细胞、交叉癌、图件、DE 与一致性。
- 唯一未闭合的精确复现点：锁定 862 程序无法由当前 DE 重建（历史 DE 丢失）；生成器已安全隔离并记录差异。

---

## 十一、清理与归档（2026-11-17）

- 本地项目由约 30 GB 降至 8.5 GB；系统盘可用由 16 GB 升至 38 GB。
- 约 21 GB 可重下/重复/派生文件已移动到
  `/Volumes/tjogzt4T/archive_gongbing_20261117/`，保留原路径结构，并生成
  `MANIFEST.txt` 与 `MANIFEST_paths.txt`。包括：
  - `06_EC_IR/data/QTD000{266,281,296,346}.all.tsv.gz`（约 11 GB，公开可重下）
  - `06_EC_IR/data/GTEx_Analysis_v8_eQTL.tar`（与正式包重复）
  - `06_EC_IR/data/finngen_R12_*.gz` 及索引（可选，公开）
  - `06_EC_IR/data/GCST90454186.h.tsv.gz`、`GCST006465_GRCh37.tsv.gz`（与正式包重复）
  - `06_EC_IR/a1/data/GCST90002238.tsv.gz`（与正式包重复）
  - `06_EC_IR/w3/results/atlas_core_{qc,scored}.h5ad`（衍生；精确 atlas 已在正式包）
  - `_tmp_ramach_*`、`a4/cache/expr.npz`、GPL `*.soft.gz` 中间文件
- 已删除：全部 `__pycache__`、`.DS_Store`、`/tmp` 临时日志。
- 已归档：`archive_20261117/biohex_runs_20261116.tar.gz`、
  `rerun_logs_20261117.tar.gz`、`download_logs_20261117.tar.gz`。
- 保留：`ec-metabolic-coloc/data/`（约 5 GB，复现所需）、`results/`、代码、文稿、
  `06_EC_IR/w2/data` 与 `w4/data` 的小型输入、`06_EC_IR/data/GTEx_Analysis_v11_eQTL_SuSiE.tar`。
