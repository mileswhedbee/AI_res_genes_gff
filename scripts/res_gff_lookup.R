rm(list = ls())
library(dplyr)
#read.csv(file = "data/both_drains.csv") -> both_drains
#read.delim(file = "data/PRET-male-geneID.annotations.gtf", header = F) -> gtf
#save(both_drains, gtf, file = "data/both_drains_andGTF")
load(file = "data/both_drains_andGTF")


# drop unneeded cols
gtf %>% select(-V2,-V6,-V8) -> gtf
both_drains %>% select(-Start.0base.) -> both_drains
colnames(gtf) = c("Chrom","Locus_type","Start","Stop","Sense","Description")

unique(both_drains$Id) -> uniq_ids

tmp.out = data.frame()
tmp.out2 = data.frame()
for (chrom in unique(both_drains$X.Chrom)) {   # for each chrom 
  both_drains %>% filter(X.Chrom == chrom) -> tmp.chrom.drains # filter into temp all the same chrom 
  gtf %>% filter(Chrom == chrom) -> tmp.chrom.gtf # match the same chrom filter into temp
  message(Sys.time()," processing chrom: ", chrom)
  for (res in unique(tmp.chrom.drains$End.1base.)) {   # iterate through uniq set of res
        tmp.chrom.gtf %>% filter(Start <= res & Stop >= res) -> tmp.overlap
        tmp.chrom.drains %>% filter(End.1base. == res) -> tmp.overlap.res
    if(nrow(tmp.overlap) > 0) {    # if N rows more than zero
      tmp.overlap$Res = res
      tmp.out = rbind(tmp.out, tmp.overlap) # then store tmp.overlap
      tmp.out2 = rbind(tmp.out2, tmp.overlap.res) # and .res
      
    }
    }
  }
gene_overlap.res = tmp.out2
overlap.gtf = tmp.out
#load(file = "output/res_gene_overlap")

# estimate AG rate AG+ and TC-
# estimate error rate - number post-filtered / N pre-filter
gene_overlap.res %>% nrow() -> tmp.error_pre
gene_overlap.res %>% filter(Type == "AG" & Strand == "+") -> tmp
gene_overlap.res %>% filter(Type == "TC" & Strand == "-") -> tmp2
rbind(tmp, tmp2) -> gene_overlap.res
gene_overlap.res %>% nrow() -> tmp.error_post
error_rate = tmp.error_post / tmp.error_pre

tmp.out = data.frame()
for (id in unique(gene_overlap.res$Id)) {
  gene_overlap.res %>% filter(Id == id) -> tmp
  tmp$N_fish = nrow(tmp)    # nrows equal to N fish editing
  tmp$DP_sum = sum(tmp$DP)  # sum depth counts
  tmp$SR_sum = sum(tmp$Supporting_reads) # sum supporting read counts
  tmp$EL_mean = mean(tmp$Edit_level)     # mean of edit levels
  tmp.out = rbind(tmp.out, tmp)          # rbind and store 
}
gene_overlap.res = tmp.out
# calc total length of genes
overlap.gtf %>% filter(Locus_type == "gene") -> tmp
tmp$diff = tmp$Stop - tmp$Start   # calc difference between boundaries of gene (i.e. start/stop)
tmp.out = data.frame()
for (id in unique(tmp$Description)) {
  tmp %>% filter(Description == id) -> tmp2
  if(nrow(tmp2) > 0) {
    tmp2[1,] -> tmp2
    }
  tmp.out = rbind(tmp.out, tmp2)
}
overlap.gtf.dedup = tmp.out
Total_length_nt = sum(tmp$diff)
unique(gene_overlap.res$Id_overlap) %>% length() -> N_uniq_res_genes
# N res per gene
N_uniq_res_genes / nrow(overlap.gtf.dedup) -> res_genes_rate
# N res per 100kb
N_uniq_res_genes / Total_length_nt -> res_per_nt
res_per_nt * 1e05 -> res_per_100kb_genes
# calc for total
unique(both_drains$Id_overlap) %>% length() -> N_uniq_total_res
(N_uniq_total_res - N_uniq_res_genes) / (697.709e06 - Total_length_nt) -> res_per_nt_whole_genome
res_per_nt_whole_genome * 1e05 -> res_per_100kb_whole_genome
res_per_100kb_genes / res_per_100kb_whole_genome -> fold_diff

save(gene_overlap.res, overlap.gtf, gtf, file = "output/res_gene_overlap")

hist(gene_overlap.res$N_fish)
hist(gene_overlap.res$Edit_level)
