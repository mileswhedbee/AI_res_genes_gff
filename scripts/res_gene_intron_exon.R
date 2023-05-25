rm(list = ls())
library(dplyr)
load(file = "output/res_gene_overlap")
# make sep col for gene ids
strsplit(overlap.gtf$Description, split = "; gene_id") -> split_tmp
as.data.frame(split_tmp) -> tmp2
t(tmp2) -> tmp2
as.data.frame(tmp2) -> tmp2
overlap.gtf$Gene_id = tmp2$V2
# number of genes edited
overlap.gtf %>% filter(Locus_type == "gene") -> tmp
unique(tmp$Gene_id) %>% length() -> N_genes
# N transcripts
# number of genes edited
overlap.gtf %>% filter(Locus_type == "transcript") -> tmp
unique(tmp$Description) %>% length() -> N_transcripts
# number of introns edited
overlap.gtf %>% filter(Locus_type == "intron") -> tmp
unique(tmp$Gene_id) %>% length() -> N_introns_genes
unique(tmp$Gene_id) -> introns_genes_ids
# number of exons edited
overlap.gtf %>% filter(Locus_type == "exon") -> tmp
unique(tmp$Gene_id) %>% length() -> N_exons_genes
unique(tmp$Gene_id) -> exons_genes_ids
# N of genes exons are edited and introns
which(is.element(introns_genes_ids, exons_genes_ids)) %>% length() -> N_exon_intron_overlap

library(VennDiagram)
grid.newpage()
draw.pairwise.venn(area1 = N_introns_genes, area2 = N_exons_genes, N_exon_intron_overlap,
                   category = c("Introns","Exons"))
