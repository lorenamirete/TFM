if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

library('biomaRt')

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                  host="grch38.ensembl.org", path="/biomart/martservice", 
                  dataset="hsapiens_gene_ensembl")

attributes = listAttributes(ensembl)
attributes[1:30,]


library(readxl)
Panel <- read_excel("~/Downloads/2_Llista gens Panell IDPs NGS VH.xlsx", 
                    sheet = "Llista completa 323 gens", col_names = FALSE)
View(Panel)


listPanel <-(Panel$...1)

chromosomes = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y')
#'transcript_start', 'transcript_end'

local = getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
  filters = 'hgnc_symbol',
  values = listPanel,
  mart = ensembl
)

local2= getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
                      filters = c('hgnc_symbol', 'chromosome_name'),
                      values = list(listPanel, chromosomes),
                      mart = ensembl
)

print(local2)
View(local2)

setwd('/home/lmirete/Downloads/1000G')
write.csv(local2,"posiciones38.csv")
