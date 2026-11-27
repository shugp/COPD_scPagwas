# 输入及数据来源

所有路径均相对于 `--data-root`。代码目录可整体移动，数据无需复制进GitHub仓库。

## 原始研究

| GEO数据集 | 原始论文 | 本对象供者数 | 本对象样本标签数 |
|---|---|---:|---:|
| [GSE168191](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168191) | Basil et al., Nature 2022, [DOI](https://doi.org/10.1038/s41586-022-04552-0) | 9 | 15 |
| [GSE171541](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171541) | Huang et al., Respiratory Research 2022, [DOI](https://doi.org/10.1186/s12931-022-02293-2) | 9 | 9 |
| [GSE173896](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173896) | Watanabe et al., AJRCMB 2022, [DOI](https://doi.org/10.1165/rcmb.2021-0555OC) | 9 | 9 |
| [GSE196638](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196638) | Wang et al., Immunity 2023, [DOI](https://doi.org/10.1016/j.immuni.2023.01.032) | 6 | 6 |

上述数量指本地分析对象，不代表GEO系列全部样本。GSE168191纳入的是外周肺样本，六名供者具有重复文库。GSE196638原研究涉及CD45⁺免疫细胞富集；不同来源的细胞构成不应直接等同于组织中的真实丰度。原始组织和实验策略见原研究。

## 完整髓系分析的起点

- `Myeloid/Myeloid.Rds`：已预处理的Seurat对象，需保留RNA counts/data、UMAP、Harmony，以及 `orig.ident`、`data_source`、`status`、`cell_type_jin`、`RNA_snn_res.0.5` 等元数据。
- `coloc.CSV`：外部共定位汇总，原列包括 `PP.H4.abf`、`gene`、`Gene`、`cell_type`、`population` 等；严格使用 `PP.H4.abf > 0.6`。脚本不从PP.H4推断基因效应方向。
- `pagwas/inte_all_gene_PCC.csv`：先前scPagwas的全细胞基因PCC结果。
- `Myeloid/disease_pathways/MSigDB_Hallmark_2020.gmt` 与 `Reactome_2022.gmt`：原分析使用的Enrichr基因集快照。精确版本校验及来源URL见 `resources/input_manifest.json`。

上游Seurat对象的从头构建、原始GWAS选择和变异级共定位不属于本包的可执行范围。作者需为这些中间输入提供可访问存档及产生方法，才能使第三方完成与本文相同范围的再分析；仅从GEO下载FASTQ不足以直接运行本包。

## 重建主图及文稿

运行 `python run.py --data-root DATA_ROOT --stage render --check` 可列出并检查全部必需文件。主图用到：

- `gPAS/scPagwas.inte_all.bar_positie_nagtive1.pdf`
- `Myeloid/trajectory/cell_pseudotime_state.csv`
- `Myeloid/results_cleaned/umap_coordinates.csv`
- `Myeloid/results_cleaned/annotation_marker_expression.csv`
- `Myeloid/results_cleaned/coloc_overall_gene_tests.csv`
- `Myeloid/disease_pathways/gsea_all_results.csv`
- `Myeloid/disease_pathways/antigen_residual_coexpression.csv`
- `Myeloid/disease_pathways/candidate_module_reactome_enrichment.csv`
- `Myeloid/disease_pathways/Reactome_2022.gmt`
- `Myeloid/disease_pathways/pathway_pseudotime/antigen_pathway_fitted_curves.csv`
- `Myeloid/disease_pathways/pathway_pseudotime/antigen_pathway_disease_pseudotime_tests.csv`
- `Myeloid/disease_pathways/pathway_pseudotime/antigen_source_curve_sensitivity.csv`

Word、Excel及RIS还读取根目录 `coloc.CSV` 和代码包内的两个JSON资源。不需要临时pickle、临时缓存CSV或临时拼图脚本。

可选scPagwas阶段另需 `inte_all.rds`、`copdgwas.txt` 和scPagwas自带参考注释。原输入文件很大，未打包或上传。
