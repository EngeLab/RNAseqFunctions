#requires sp.scRNAseqData package
#run from package root with: source('inst/rawdata/processedCounts.R')

#counts
packages <- c("sp.scRNAseqData", "tidyverse")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

lg <- grepl("^s", colnames(countsSorted2))
rv <- matrixStats::rowMaxs(countsSorted2[, lg])
select <- order(rv, decreasing = TRUE)[1:2000]

pro.counts <- rbind(countsSortedERCC2[, lg], countsSorted2[select, lg]) %>%
  .[, rev(order(colnames(.)))]

#meta
pro.meta <- countsSortedMeta2 %>%
filter(sample %in% colnames(countsSorted2)[lg]) %>%
select(sample, cellTypes)

save(pro.counts, pro.meta, file = "data/processedCounts.rda", compress = "bzip2")
