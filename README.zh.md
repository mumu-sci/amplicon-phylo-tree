# amplicon-phylo-tree（中文说明）

一个给 AI agent 用的技能：**从二代扩增子数据（16S / 18S / ITS / nifH / amoA / COI…）
构建可发表、可辩护的系统发育树**。

它是「工作流文档 + 参数化辅助脚本」，不是黑箱：每个阈值都按本研究确定并落盘，
每道门都有退出码，每一次覆盖都留痕。

> 技能本体在 [`skills/amplicon-phylo-tree/`](skills/amplicon-phylo-tree/)。
> 让 agent 读它里面的 `SKILL.md`。

---

## 它做什么

```
探测环境 → 定类群口径 → 输入审计 → 读码框/链向归一 → 比对 → 比对体检
   → 建树 → 生根 → 验收 → 交付
```

每一步都产出审稿人会要的东西：

| 门 | 检查什么 | 不通过时 |
|---|---|---|
| **G1** | 读码框/链向可判定比例（仅蛋白编码） | **只报告**；不可判定的序列仍入树 |
| **G1b** | 输入审计：方向、同源性、混合扩增区段、嵌合体步骤 | 方向混杂 → 修正后重跑；同源性不足或区段混杂 → **停下报告** |
| **G2** | 比对体检（碎裂、字母表、模糊碱基策略） | 结构错误**硬停**；比例指标默认 FAIL，可带理由覆盖并留痕 |
| **G3** | 模型 vs 数据类型（四层 + fail-closed） | **硬停** |
| **G4** | 树结构与 tip 集合（split 身份、根度数、枝长） | **硬停** |
| **G5** | 可追溯（run manifest、工具版本、SHA256） | 报出并停下 |

## 它**不做**什么

- **不自动删类群**。碎裂的补救是参考骨架或分块比对——两者都保留全部 tip。
- **不做真实性判决**。它不判定哪棵树"对"，只判定你做的分析**是否如你所声明的那样站得住**。
- **不算下游生态指标**。PD / UniFrac / NRI-NTI / βNTI 一概不代算；它给你树和必须写出的声明。
- **不做上游去噪/聚类**。DADA2 之类在它上游。
- **不捆绑二进制或数据库**。参考库请自备。
- **不实现片段插入**（SEPP / EPA-ng 只做路由与文档，这是明确的设计决定）。

## 依赖

- **Python 3.9+**（仅标准库）
- **MAFFT**，以及 **FastTree** 或 **IQ-TREE**
- **HMMER** 与/或 **BLAST+**（编码标志基因的输入审计）
- **R 4.x** + **ape** + **phangorn**（生根，`scripts/root_tree.R`）
- **vsearch**（建议；嵌合体筛查）

`python skills/amplicon-phylo-tree/scripts/probe_tools.py --require mafft,tree`
会报告缺什么并记录版本。

## 安装

把技能文件夹复制进你的 agent 的技能目录：

```bash
cp -r skills/amplicon-phylo-tree <你的 skills 目录>/amplicon-phylo-tree
```

文件夹名**必须**保持 `amplicon-phylo-tree`（要与 `SKILL.md` 的 `name` 一致）。

## 校验你的副本

```bash
cd skills/amplicon-phylo-tree
python scripts/check_consistency.py     # 文档与脚本是否一致
python scripts/regression_tests.py      # 黑箱回归测试（v1.0.1 为 36 条）
python scripts/selftest.py              # 坐标与遗传密码契约
```

## 已知局限（照实写）

- **FastTree 是近似 ML**：它的支持率是 SH-like 局部检验，**不可**与 UFBoot2 / SH-aLRT 横向比较。
  技能会强制你写明这一点。
- **中点法是约定，不是证据**：无外群时技能默认做生根敏感性分析，且**拒绝**让 `NOT RUN` 的报告过验收门。
- **ITS 不是群落系统发育标记**：跨属以上的 ITS 树被拒绝用于指标；属内树需收窄 + 结构约束比对，且只作补充材料。
- **内存启发式只是规划用**：其中 α（B/格）已被实测反驳；文档如实写出，而不是藏起来。
- **实测锚点来自开发者机器**，并且被明确标注为"仅供数量级参考"。请自己跑基准测试。

## 许可

MIT，见 [LICENSE](LICENSE)。

## 版本

1.0.1 —— 见 [CHANGELOG.md](CHANGELOG.md)。
