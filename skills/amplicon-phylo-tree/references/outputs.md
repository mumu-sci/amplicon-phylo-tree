# Deliverables and Document Templates

---

## 1. Delivery checklist

| File | Description |
|---|---|
| `tree.rooted.nwk` | **The single downstream entry point** (rooted). For UniFrac / PD / NRI-NTI / betaNTI, this tree may ONLY be used |
| `tree.rooted.nwk.candN.nwk` | Candidate-root trees, one per candidate root position for rooting sensitivity (written by `root_tree.R`; the first-line comment holds the candidate name and the placement method). **If you want to compare rooting sensitivity with downstream metrics, compute them yourself on this batch of candidate trees** — this skill reports the root split ONLY |
| `tree.unrooted.nwk` | For inspection and figures only. **Support values exist on the unrooted tree alone** |
| `aln.fasta` | The alignment used for tree inference (recommended to archive compressed alongside) |
| `qc/tools.json` | Tool-version record (`probe_tools.py`) |
| `qc/taxa.log` | Taxon scope (rule text verbatim + retained proportion + the three discard items) |
| `qc/audit.json` | Input audit (orientation, homology, mixed regions, `--region-map` summary and SHA256) |
| `qc/aln.json` | Alignment health check (diagnostic matrix, gap-rate definition, **`gate_overrides` override record**) |
| `qc/tree.log` | Full log of the tree-inference tool (G3 evidence) |
| `qc/rooting.txt` | Rooting report (including a **rooting sensitivity** section) |
| `qc/verify.json` | G3 + G4 acceptance report |
| `run_manifest.txt` | See §2 |
| `ops_log.md` | See §5 |
| `methods_zh.md` / `methods_en.md` | See §3 / §4 |

> **Naming is not enforced.** If your project or pipeline has its own contract (for example,
> some workflows require `99_tree.nwk`, or require files under `_docs/`), follow that; this
> skill requires only that the **content is complete** and **traceability** is preserved.

---

## 2. run_manifest (traceability, G5)

No naming convention is enforced, but every run MUST produce a manifest containing:

```
run_id      : <timestamp or an id you assign>
label       : <study label>
date        : <ISO date>
tool_versions:
  mafft     : v7.526
  fasttree  : 2.1.11 Double precision
  iqtree    : 2.4.0
  hmmsearch : 3.4
  R         : 4.5.2 (ape 5.8-1, phangorn 2.12.1)
parameters  : see params.md (D1–D11, each with its rationale)
commands    : <list the exact commands actually executed, in order, with all parameters>
inputs:
  <relative path>  sha256=<...>
outputs:
  <relative path>  sha256=<...>
deviations  : <every place you departed from the default, plus every item overridden by --accept-and-record and its reason>
```

How you generate it is up to you (by hand, `sha256sum`, or a script). **The purpose is: to be able to reconstruct this run completely a year later (or when a reviewer asks).**

---

## 3. Methods template (Chinese -- kept in Chinese on purpose; paste it as-is)

> **系统发育树构建。** 用于建树的 `<标志基因>` 代表序列按以下标准筛选：`<D3 的规则原文>`，
> 共保留 `<N>` 条序列（占全部 ASV 的 `<x>`%，占全部 `<marker>` reads 的 `<y>`%）。
> 序列方向与同源性经审计（反向互补检出 `<n>` 条；未检出/已修正）；嵌合体检测使用
> `<vsearch --uchime3_denovo / 显式跳过并说明>`。
> 〔**仅蛋白编码基因**〕序列的读码框与链向经统一：以 `<HMM 名称与版本>` 模型
> （`hmmsearch --cut_tc`）判定，`<n>` 条获得受信命中并统一到 `<框/链>`。
> 〔**走氨基酸时**〕因类群间核酸同一性较低（域内成对氨基酸同一性中位数 `<p_aa>`），
> 第三位密码子趋于饱和，故在氨基酸层面比对并**回译为密码子比对**
> （`<PAL2NAL / trimAl -backtrans / MACSE>`）后再建树。
>
> 多序列比对使用 **MAFFT `<版本>`**（`<实际算法与参数>`）；比对矩阵为
> `<N>` 条序列 × `<L>` 列，空位比例为 `<g>`%（空位率定义为 gap 字符数 /（序列数 × 列数））。
> 〔**若分块**〕比对分 `<k>` 块独立进行后以空位填充拼接；**跨块枝长不具生物学含义**，
> 可解释关系限于各块内部。
>
> 系统发育树使用 **`<FastTree 版本>` / `<IQ-TREE 版本>`** 构建，替换模型为 `<模型>`。
> 〔**FastTree**〕采用近似最大似然（GTR+CAT 近似、NNI/SPR 局部搜索）；
> 枝支持度为 1,000 次重抽样的 SH-like 局部检验值，**不等价于非参数自举，
> 也不可与 UFBoot2 / SH-aLRT 数值横向比较**。
> 〔**IQ-TREE**〕模型经 ModelFinder 自动选择；枝支持度同时用 UFBoot2（1,000 次）
> 与 SH-aLRT（1,000 次），分别以 ≥ 95 与 ≥ 80 视为强支持；UFBoot2 基于位点似然重抽样，
> 与标准非参数自举不可直接比较。
> 本文所有枝长与支持度均来自同一流程，未与其他软件结果拼接。
>
> 〔**生根**〕`<使用外群 X 生根>` / `<因无可用外群，使用中点法
> （phangorn::midpoint）生根；为评估根的敏感性，我们比较了 <n> 个候选根位置
> （脚本已一并输出各候选根树），其中 <k> 个改变了根两侧的 split，故涉及根的结论
> 按此限制报告>`。
> **注意**：这里只声明 **root split** 的敏感性；UniFrac / βNTI 等下游指标要在同一批
> 候选根树上另算并单独声明（本 skill 不代算）。
>
> **结论边界。** `<按标志基因类别填写，见下方"可解释上限"清单>`
> 〔可选〕为评估类群筛选的敏感性，我们在更松/更严的规则（`<规则>`）下重算，
> 主要结论未改变（补充材料 S`<x>`）。

### 「Interpretation ceiling」selected by class (**exactly one entry is REQUIRED**)

- `rrna`：本数据集扩增区段长度约 `<a>–<b>` bp，对亲缘关系很近的类群（如属内不同种）分辨率有限，
  故系统发育相关性解释主要针对较高分类层级。
- `its`：**本树不适用于跨类群定名或群落系统发育指标**；ITS 的长度异质性与拷贝内变异
  （不完全协同进化）限制了其分辨率与拓扑置信度，本树仅用于展示组内关系。
- `functional-nuclear`：`<marker>` 的不同功能簇（如 `<簇名>`）为古老旁系同源支系；
  在约 `<l>` bp 扩增子与核酸层面的同源性不足以解析跨簇深层关系，此类节点不可作为进化结论。
  本树适用于**本数据集内部**的多样性/群落系统发育指标。
- `organelle`：本区段（`<l>` bp）第三位密码子饱和明显，深层节点不可靠；遗传密码表使用 `<表号>`。

---

## 4. Methods template (English)

> **Phylogenetic tree construction.** Representative sequences of `<marker>` were
> selected using the rule `<D3 rule verbatim>`, retaining `<N>` sequences
> (`<x>`% of ASVs, `<y>`% of `<marker>` reads). Sequence orientation and homology
> were audited (reverse-complemented sequences detected: `<n>`); chimeras were
> screened with `<vsearch --uchime3_denovo / explicitly skipped, with reason>`.
> 〔protein-coding only〕Reading frame and strand were normalised per sequence using
> `<HMM name and version>` with trusted cutoffs (`hmmsearch --cut_tc`), unifying
> `<n>` sequences to `<frame/strand>`.
> 〔amino-acid route〕Because nucleotide identity between clades was low
> (median within-domain amino-acid identity `<p_aa>`), so that third codon positions
> are saturated, sequences were aligned at the amino-acid level and the alignment was
> **back-translated to codons** (`<PAL2NAL / trimAl -backtrans / MACSE>`) before tree
> inference.
>
> Multiple sequence alignment used **MAFFT `<version>`** (`<algorithm and parameters>`);
> the final matrix was `<N>` sequences x `<L>` columns with `<g>`% gaps (gaps defined as
> gap characters / (sequences x columns)). 〔if blocked〕Alignment was performed in
> `<k>` blocks and concatenated with gap padding; **branch lengths between blocks are
> not biologically interpretable**, and interpretable relationships are confined to
> within-block comparisons.
>
> The tree was inferred with **`<FastTree>` / `<IQ-TREE>`** under `<model>`.
> 〔FastTree〕This is **approximate** maximum likelihood (GTR+CAT approximation,
> NNI/SPR local search); support values are SH-like local tests from 1,000 resamples,
> which are **not equivalent to non-parametric bootstrap and cannot be compared with
> UFBoot2/SH-aLRT values**.
> 〔IQ-TREE〕The model was selected by ModelFinder; branch support used both UFBoot2
> (1,000 replicates) and SH-aLRT (1,000 replicates), with UFBoot2 >= 95 and
> SH-aLRT >= 80 taken as strong support. UFBoot2 resamples site likelihoods and is not
> directly comparable to standard non-parametric bootstrap. All branch lengths and
> support values come from a single pipeline; results from different programs were not
> combined.
>
> 〔rooting〕`<The tree was rooted using outgroup X>` / `<As no suitable outgroup was
> available, the tree was rooted by the midpoint method (phangorn::midpoint); to assess
> root sensitivity, <n> candidate root positions were compared (the candidate trees are
> written alongside the main tree), and <k> of them changed the bipartition across the
> root, so root-dependent conclusions are reported with that limitation>`.
> **Note**: this declares the sensitivity of the **root split** only.  Downstream
> metrics (UniFrac / betaNTI) must be computed by you on the same candidate trees and
> declared separately; this skill does not compute them.
>
> **Interpretation limits.** `<one sentence from the class-specific list above>`

---

## 5. ops_log template

```markdown
# <marker> phylogenetic tree -- ops log

Date: <ISO>      run_id: <id>      operator: <>
Study label: <>

## 1. Parameters (see params.md; only non-default items are excerpted here)
| ID | Value | Rationale | Deviation from default? |
|---|---|---|---|

## 2. Execution steps and commands
<List, step by step, the commands actually executed verbatim, their wall time, and key output>

## 3. Gate results
| Gate | Result | Numbers / evidence |
|---|---|---|
| G1 reading frame/strand |  |  |
| G1b input audit |  |  |
| G2 alignment health check |  |  |
| G3 model assertion |  | what the four layers of evidence are |
| G4 tree structure and tips |  | two-way set difference / root degree / zero-length branches |
| G5 traceability |  |  |

## 4. Deviations and overrides (**every single one MUST be written down with its reason**)
<which defaults were changed, which gates were overridden by --accept-and-record, who decided>

## 5. Pitfalls hit and solutions
<symptom / root cause / handling>

## 6. Unresolved / to be confirmed
<e.g. no suitable outgroup available; a reference database version for some taxon still to be locked; some region still suboptimal>
```
