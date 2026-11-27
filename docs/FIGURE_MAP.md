# 主图与分析代码对应

最终拼图由 `tools/render_figure.py` 统一重绘；不会改变统计值。所有面板保留为矢量，除原始PDF可能包含的字形轮廓外，新增标签是SVG文字。

| 面板 | 展示内容 | 产生数值结果的代码/输入 |
|---|---|---|
| A | scPagwas阳性/阴性细胞比例 | 指定原PDF；可选上游 `analysis/01_scPagwas_fast.R` |
| B | 五类Myeloid UMAP | `analysis/02_Myeloid.R`、`myeloid_annotation.R`；沿用输入对象UMAP |
| C | Myeloid marker气泡图 | `analysis/02_Myeloid.R`，输出 `annotation_marker_expression.csv` |
| D | COPD/对照分组扩散拟时序 | `analysis/03_Myeloid_trajectory.R`、`myeloid_trajectory.py`；Harmony前30维，30邻居，图State分辨率0.15 |
| E | State内疾病比较的GSEA | `analysis/04_Myeloid_disease_pathways.R`、`myeloid_gsea.py` |
| F | 控制供者和State的候选–伙伴偏相关 | `analysis/06_Myeloid_result_panels.R` |
| G | 候选共表达模块与抗原呈递通路联系 | `analysis/myeloid_antigen_function.R`；通路成员标记从原Reactome GMT独立读取 |
| H—K | 四条通路的疾病分组拟时序曲线 | `analysis/08_Myeloid_pathway_pseudotime.R`、`myeloid_pathway_source_curve_sensitivity.R` |

G的颜色表示模块过度富集FDR，并非候选基因自身的疾病效应或共定位概率。空心圆是所用Reactome快照的直接成员注释，未收录不能排除已知功能。比如HLA-DQB1属于MHC体系，但在所用快照的相应通路集合中未被列入。

## 正文中的其他统计与敏感性分析

- 总体Myeloid候选表达：`myeloid_coloc.R`（从02调用）。
- State及拟时序表达：`03_Myeloid_trajectory.R`。
- 疾病×拟时序基因交互及State内比较：`04_Myeloid_disease_pathways.R`。
- State5抗原呈递评分、候选评分关联及top100残差模块：`myeloid_antigen_function.R`。
- 疾病组关联强度差异：`07_Myeloid_differential_association.R`、`myeloid_differential_association.py`、`myeloid_differential_association_sensitivity.py`。
- 全通路探索与来源曲线敏感性：`myeloid_all_pathway_mixed_models.R`、`myeloid_pathway_screen_curves.R`、`myeloid_pathway_source_curve_sensitivity.R`。
- State组成及其他附属图：`myeloid_result_panels.py`、`05_Myeloid_plots.py`、`myeloid_pathway_pseudotime.py`。

宽泛Myeloid注释使用resolution 0.5；轨迹图划分State使用Leiden分辨率0.15。这是两个不同图与任务的参数，不应混为同一种细胞类型分辨率。
