
#!/path/to/Rscript
#
# example
# Rscript scanoneqtl_c.R -n Ab -g Ab_gen.txt -p Ab_phe.txt -x 1000
# 
# ------------ RIL input genotype file format -------------
# id	chr	N0060	N0061	N0062	N0063
# 1_1_2146985	1	0	0	2	2
# 1_2146986_2350291	1	2	0	2	2
# 1_2350292_2549669	1	2	0	2	2
# ------------ RIL input phenotype file format -------------
# id	HD	PH
# N0060	0.922409076	0.769724412
# N0061	0.021655533	-0.008497137
# N0062	0.322102747	0.155586838
# N0063	-0.171848196	0.543512724
#-------------------------------------------------
#
library(optparse)
library(qtl)
option_list <- list(
  make_option(c("-n", "--name"), type = "character", default = FALSE, help = "please give a ril name"),
  make_option(c("-g", "--geno"), type = "character", default = FALSE, help = "please give a .txt filename"),
  make_option(c("-p", "--pheno"), type = "character", default = FALSE, help = "please give a .txt filename"),
  make_option(c("-x", "--xperm"), type = "character", default = 1000, help = "please give a number")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# pass parameter
name <- opt$n
g <- opt$g
p <- opt$p
nperm <- opt$x
# name="Ab"
# nperm=10
# g <- "Ab_gen.txt"
# g <- read.table("Ab_gen.txt", header = T)
g <- read.table(g, header = T)
g2 <- g
colnames(g2) <- c("id", "", colnames(g)[3:ncol(g)])
# write.csv(g2, "mapthis_Ab_gen_rot.csv", quote = F, row.names = F)
write.csv(g2, paste("mapthis_",name,"_gen_rot.csv",sep = ""), quote = F, row.names = F)

# p <- "Ab_phe.txt"
# p <- read.table("Ab_phe.txt", header = T)
p <- read.table(p, header = T)

p2 <- t(p)
colnames(p2) <- p2[1,]
p3 <- p2[-1,]
p4 <- data.frame(id=rownames(p3), p3)
write.csv(p4, paste("mapthis_",name,"_phe_rot.csv",sep = ""), quote = F, row.names = F)


aa5 <- read.cross("csvsr", ".",
                  genfile = paste("mapthis_",name,"_gen_rot.csv",sep = ""), 
                  phefile = paste("mapthis_",name,"_phe_rot.csv",sep = ""),
                  crosstype = "riself", 
                  genotypes = c(0,1,2), alleles = c(0,2))

endmap=NULL
for(i in 1:12){
  map1=mm=mmm=NULL
  map1 <- est.map(aa5, chr = i, error.prob = 0.0001, map.function = "morgan")
  mm <- unlist(map1)
  mmm <- data.frame(mm)
  endmap <- rbind(endmap, mmm)
}
write.table(endmap, paste(name,"_geneticmap.txt",sep=""), quote = F, col.names = T, row.names = T, sep = "\t")

g3 <- data.frame(g2[,1:2],endmap[,1],g2[,3:ncol(g2)])
colnames(g3) <- c("id", "", "", colnames(g2)[3:ncol(g2)])
write.csv(g3, paste("mapthis_",name,"_gen_rot2.csv",sep = ""), quote = F, row.names = F)

aa5 <- read.cross("csvsr", ".",
                  genfile = paste("mapthis_",name,"_gen_rot2.csv",sep = ""), 
                  phefile = paste("mapthis_",name,"_phe_rot.csv",sep = ""),
                  crosstype = "riself", 
                  genotypes = c(0,1,2), alleles = c(0,2))
lis <- calc.genoprob(aa5, step=1, error.prob=0.001)

pncol <- ncol(p)
pcolnames <- colnames(p)
for(i in 2:ncol(p)){
  hd=hdperm=threshold=hdp2=hdp3=NULL
  hd <- scanone(cross = lis, pheno.col = i) # id在第一行，则从2开始
  hdperm <- scanone(cross = lis, pheno.col = i, n.perm = nperm)
  threshold <- summary(hdperm, alpha=0.05)[1]
  write.table(hd, paste(name,"_",pcolnames[i],"_scanone.txt", sep = ""), quote = F)
  write.table(threshold, paste(name,"_",pcolnames[i],"_threshold.txt", sep = ""), quote = F)
  hdp2 <- hd[hd$lod>threshold,]
  hdp3 <- hdp2[hdp2$pos%%1>0,]
  write.table(hdp3[,-2], paste(name,"_",pcolnames[i],"_scanoneQTL.txt", sep = ""), quote = F)
}






