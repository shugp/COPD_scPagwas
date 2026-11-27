# COPD髓系细胞：公共单细胞数据再分析

## 安装

 Python 3.12：


R环境使用Seurat 5，另需Matrix、ggplot2、ggrepel、limma、statmod、nlme、cluster

```r
install.packages(c("Seurat", "ggplot2", "ggrepel", "statmod", "nlme"))
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install("limma")
```

可选scPagwas阶段还需要scPagwas、bigparallelr、bigreadr、bigstatsr及其参考注释。

