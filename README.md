# COPD髓系细胞：公共单细胞数据再分析

本代码包配套中文 Research Letter《公共肺单细胞数据整合与COPD髓系抗原呈递状态的系统分析》。将本目录的**内容**作为 GitHub 仓库根目录即可，无需上传整个原始工作目录。

## 可以复现什么

- **重建投稿材料**：从现有数值结果生成一份可编辑SVG拼图、中文Word正文、单页原始共定位Excel和配套RIS文献。
- **重跑髓系分析**：从已整合的 `Myeloid/Myeloid.Rds` 开始，依次完成注释与供者聚合、拟时序、疾病比较、GSEA、候选模块及通路轨迹分析。
- **可选前期步骤**：从既有整合全细胞对象与GWAS汇总重跑scPagwas；对Myeloid原有分辨率进行诊断。

本包不从FASTQ重建整合对象，也不重新计算共定位。`coloc.CSV` 是既有输入，候选阈值为 **PP.H4 > 0.6**。主图不展示PP.H4概率，候选–通路热图展示功能背景，不能单凭共表达推断疾病特异性或因果性。

## 安装

只重绘投稿材料需要 Python 3.12、下列Python包及Poppler的 `pdftocairo`：

```bash
python -m pip install -r requirements-render.txt
# Ubuntu/Debian（尚未安装Poppler时）：sudo apt install poppler-utils
```

完整分析另需：

```bash
python -m pip install -r requirements-analysis.txt
Rscript tools/check_r_dependencies.R
```

R环境使用Seurat 5，另需Matrix、ggplot2、ggrepel、limma、statmod、nlme、cluster。缺包时可安装：

```r
install.packages(c("Seurat", "ggplot2", "ggrepel", "statmod", "nlme"))
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install("limma")
```

可选scPagwas阶段还需要scPagwas、bigparallelr、bigreadr、bigstatsr及其参考注释。原运行环境版本见 `resources/software_versions.json`；Python显式版本保存在requirements文件中。这不是覆盖所有传递依赖的环境锁文件。完整分析涉及大型Seurat对象，内存需求明显高于仅重绘。

## 准备数据

把数据保存在仓库之外的 `DATA_ROOT`，维持原项目的相对路径。输入清单及数据来源见 [docs/INPUTS.md](docs/INPUTS.md)。不要把绝对路径写入脚本；统一使用 `--data-root`。

代码在任何当前目录都可调用，例如：

```bash
python /path/to/repository/run.py --data-root /path/to/data --stage render --check
python /path/to/repository/run.py --data-root /path/to/data --stage render --output-dir /path/to/output
```

`render` 只读取已有统计结果，不重新估计模型。输出恰好为：

1. `Myeloid_COPD_Composite.svg`：单份可编辑矢量拼图；A直接嵌入指定PDF的矢量内容，原PDF文字可能为轮廓。
2. `Myeloid_COPD_Letter_CN.docx`：五段正文、可获得性声明、图注及8篇参考文献。
3. `Myeloid_Coloc_Candidate_Pathways.xlsx`：沿用既有文件名，内容仅为一个无配色sheet的原始共定位结果。
4. `Myeloid_COPD_References.ris`：与正文引用对应。

只修改正文/参考文献时编辑 `resources/manuscript_cn.json` 与 `resources/references.json`，运行 `--stage documents`。主图逐面板来源见 [docs/FIGURE_MAP.md](docs/FIGURE_MAP.md)。

## 重新运行分析

先检查输入，再执行相应阶段。以下命令会在 `DATA_ROOT/Myeloid` 下写入结果；如需保留旧结果，请使用单独的数据工作目录。

```bash
python run.py --data-root /path/to/data --stage annotation --check
python run.py --data-root /path/to/data --stage analysis
python run.py --data-root /path/to/data --stage render --output-dir /path/to/output
```

`analysis` 顺序执行 `02_Myeloid.R → 03_Myeloid_trajectory.R → 04_Myeloid_disease_pathways.R`。后者串联候选功能、GSEA、组间关联强度及通路拟时序的辅助脚本。也可分别指定 `annotation`、`trajectory` 或 `disease`。新增的 `diagnostics` 仅核对已有分辨率，不能替代人工marker复核；没有重新优化resolution。

```bash
python run.py --data-root /path/to/data --stage audit
python run.py --data-root /path/to/data --stage diagnostics
SCPAGWAS_CORES=16 python run.py --data-root /path/to/data --stage scpagwas
```

供者归并规则：GSE168191的同一供者rep1/rep2合并；其他来源使用带数据集前缀的orig.ident作为供者单位。本对象为39个样本标签、33个供者单位（COPD 16、对照17）。各GEO完整系列的样本量与本分析实际纳入子集不同。

## 代码核验与研究边界

```bash
python tools/validate_package.py
Rscript tools/check_r_dependencies.R
```

已对所有Python/R脚本进行语法检查，使用现有输入实际运行独立重绘入口并核对SVG、Word、Excel、RIS；没有为本次整理重新计算大型矩阵、拟时序图或全部模型。跨平台、软件版本及随机图算法可能造成数值差别。原始脚本与整理版校验信息见 `resources/source_code_manifest.json`。

State5评分比较和候选关联检验使用既有合并供者模型；来源相关的敏感性分析另行报告。所有统计量保留其既有检验家族，未将探索性筛选改写为确定的疾病机制。原始coloc细胞标签 `immune_Proliferating` 本身不能单独确定髓系归属，沿用既有候选范围时需结合原参考数据注释解释。

## 公开代码与投稿

上传步骤见 [docs/GITHUB_UPLOAD.md](docs/GITHUB_UPLOAD.md)。目前未代为创建或发布远程仓库。上传后将**真实公开仓库地址**传入：

```bash
python run.py --data-root /path/to/data --stage documents --output-dir /path/to/output --repository-url https://github.com/YOUR_ACCOUNT/YOUR_REPOSITORY
```

没有提供地址时，Word会明确保留“待作者上传后填写；尚未公开”标记。投稿前还需给出整合对象、所用GWAS/共定位汇总的可访问来源或存档，而不只提供代码URL。中文版按约1,000英文词的篇幅结构组织；英文正式稿仍须重新核对期刊词数。作者、单位、资助及利益冲突等投稿信息由作者补充。

本包没有替作者指定代码许可或版权归属；创建公开仓库时请选定适合该研究的许可证，并保留第三方工具及原数据的引用要求。
