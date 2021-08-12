
library(tidyverse)

rvg11.binding <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/7_tobias/RVG11_BINDetect_single_output/allTFBS_bound.RVG11.bed',col_names = F)
rvg1.binding <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/7_tobias/RVG1_BINDetect_single_output/allTFBS_bound.RVG1.bed',col_names = F)
rvg4.binding <- read_tsv('/Volumes/BlairPerry.Data1/__Cvv_Venom_Regulation/_newAnalysis_August2021/ATAC/7_tobias/RVG4_BINDetect_single_output/allTFBS_bound.RVG4.bed',col_names = F)

rvg11.min <- min(rvg11.binding$X10)
rvg1.min <- min(rvg1.binding$X10)
rvg4.min <- min(rvg4.binding$X10)

boundThresholds <- tribble(
  ~sample, ~threshold,
  'RVG11', rvg11.min,
  'RVG1', rvg1.min,
  'RVG4', rvg4.min
)

write_tsv(boundThresholds,'analysis/4_atacseq/tobias_footprinting/BoundThresholds_08.11.21.tsv')
