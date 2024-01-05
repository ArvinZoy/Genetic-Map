#
# 计算Genetic Map ==== 
# 1. chisqtest >0.05 
# 2. lagranre for all

setwd("")
name3 <- c("Ab") # bias 无问题的
#如果某一条染色体小于5个SNP，则随机增加10个SNP
for(i in c(name3)){
  # i="Az"
  geno=geno1=geno2=x2p=geno3=geno1.3.1=NULL
  geno <- read.table(paste0("binall.",i,".txt"), header = T, stringsAsFactors = F)
  geno1 <- data.frame(SNP=paste0("M",geno[,1],"_",geno[,2],"_",geno[,3]), geno) # SNP名称调整，方便后面比对
  x2p <- apply(geno1[, -c(1:5)], 1, function(x){chisq.test(c(length(which(x==0)), length(which(x==2))))$p.value}) # 卡方测验
  geno1.1 <- cbind(geno1,x2p)
  geno1.2=geno1.3=NULL
  for(m in 1:12){ # 为了后面可以完整插值方法，需要加入每条染色体两头的SNP
    geno1.1.chr =NULL
    geno1.1.chr <- geno1.1[geno1.1[,2]==m, ]
    geno1.2 <- rbind(geno1.2, geno1.1.chr[c(1,nrow(geno1.1.chr)),]) # 染色体两头的SNP
    geno1.3 <- rbind(geno1.3, geno1.1.chr[-c(1,nrow(geno1.1.chr)),]) # 染色体非两头的SNP
  }
  geno2 <- rbind(geno1.2, geno1.3[geno1.3[,ncol(geno1.3)]>0.05, ]) # 用来计算总GD的SNP
  geno2 <- geno2[order(geno2[,2], geno2[,3]), ]
  geno1.3.1 <- geno1.3[geno1.3[,ncol(geno1.3)]<=0.05, ]
  
  for(m in 1:12){
    #m=9
    geno2.chr=geno1.3.1.chr=geno1.3.1.chr.b=NULL
    geno2.chr <- geno2[geno2[,2]==m, ]
    geno1.3.1.chr <- geno1.3.1[geno1.3.1[,2]==m, ]
    geno1.3.1.chr.b <- geno1.3.1[geno1.3.1[,2]!=m, ]
    if(nrow(geno2.chr)<=5){
      cat(m,"\t")
      if(nrow(geno1.3.1.chr)<10){
        geno2 = rbind(geno2, geno1.3.1.chr)
      } else {
        geno2 = rbind(geno2, geno1.3.1.chr[c(sample(1:nrow(geno1.3.1.chr),10,replace = F)), ])
      }
    } else {
      geno2 = geno2
    }
  }
  
  geno2 <- geno2[,-ncol(geno2)]
  geno2 <- geno2[order(geno2[,2], geno2[,3]), ]
  # geno2 <- geno1[x2p>0.05, ] # 筛选符合孟德尔分离定律的bin
  cat(i,"All bin number: ",nrow(geno1),", not bias bin number : ", nrow(geno2), "; rate : ", nrow(geno2)/nrow(geno1), "\n")
  # geno3 <- geno1[x2p<=0.05, ]
  recombinationFreq=geneticDistanceInterval=geneticDistanceAccumulate=NULL
  for(j in 1:12){ # All bin genetic map
    #i=1
    chr=NULL
    chr <- geno1[geno1[,2]==j, ]
    recombinationFreqCHR=geneticDistanceCHRinterval=geneticDistanceCHRaccumulate=0
    for(k in 2:nrow(chr)){
      # k=2
      aa=r=m=NULL
      aa <- (as.numeric(chr[k, 5:ncol(chr)] + chr[k-1, 5:ncol(chr)])) # 以相邻标记相加的方式方便计算geneticmap
      r <- 0.5*length(which(aa==2))/length(which(aa==0 | aa==4)) # RIL # 计算RIL重组率
      if(r<0.38){m <- 25*log((1+2*r)/(1-2*r))} else{m<-50}# Kosambi 方法转换RF 2 cM；其中设定RF最大值0.4
      recombinationFreqCHR <- c(recombinationFreqCHR, r)
      geneticDistanceCHRaccumulate <- c(geneticDistanceCHRaccumulate, m+sum(geneticDistanceCHRinterval))
      geneticDistanceCHRinterval <- c(geneticDistanceCHRinterval, m)
    }
    recombinationFreq <- c(recombinationFreq, recombinationFreqCHR)
    geneticDistanceInterval <- c(geneticDistanceInterval, geneticDistanceCHRinterval)
    geneticDistanceAccumulate <- c(geneticDistanceAccumulate, geneticDistanceCHRaccumulate)
  }
  GD1 <- data.frame(geno1[,1:4], recombinationFreq, geneticDistanceInterval, geneticDistanceAccumulate, Chi_Square=x2p)
  write.table(GD1, paste0("binall.geneticDistance.",i,".chisq0.0"), quote = F, col.names = T, row.names = F, sep = "\t")
  cat(i,"done for all bin", "\n")
  
  
  recombinationFreq=geneticDistanceInterval=geneticDistanceAccumulate=NULL
  for(j in 1:12){ # not bias bin genetic map
    #j=12
    chr=NULL
    chr <- geno2[geno2[,2]==j, ]
    recombinationFreqCHR=geneticDistanceCHRinterval=geneticDistanceCHRaccumulate=0
    for(k in 2:nrow(chr)){
      # k=2
      # cat(k,"\t")
      aa=r=m=NULL
      aa <- (as.numeric(chr[k, 5:ncol(chr)] + chr[k-1, 5:ncol(chr)])) # 以相邻标记相加的方式方便计算geneticmap
      r <- 0.5*length(which(aa==2))/length(which(aa==0 | aa==4)) # RIL # 计算RIL重组率
      if(r<0.38){m <- 25*log((1+2*r)/(1-2*r))} else{m<-50}# Kosambi 方法转换RF 2 cM；其中设定RF最大值0.4
      recombinationFreqCHR <- c(recombinationFreqCHR, r)
      geneticDistanceCHRaccumulate <- c(geneticDistanceCHRaccumulate, m+sum(geneticDistanceCHRinterval))
      geneticDistanceCHRinterval <- c(geneticDistanceCHRinterval, m)
    }
    recombinationFreq <- c(recombinationFreq, recombinationFreqCHR)
    geneticDistanceInterval <- c(geneticDistanceInterval, geneticDistanceCHRinterval)
    geneticDistanceAccumulate <- c(geneticDistanceAccumulate, geneticDistanceCHRaccumulate)
  }
  tt=NULL
  tt <- merge(geno2[,1:4], geno1.1[,c(1,ncol(geno1.1))], all.x = T,by = "SNP")
  tt <- tt[order(tt[,2],tt[,3]), ]
  GD2 <- data.frame(geno2[,1:4], recombinationFreq, geneticDistanceInterval, geneticDistanceAccumulate, Chi_Square=tt$x2p)
  write.table(GD2, paste0("binall.geneticDistance.",i,".chisq0.05"), quote = F, col.names = T, row.names = F, sep = "\t")
  cat(i,"done for not bias bin 0.05", "\n")
}

# 插值法插入其他剩余bin
name4 <- c("Ab") # bias 无问题的

lagrange <- function(x, xk, yk){ # 拉格朗日插入法
  n=lagr=NULL
  n <- length(xk)
  lagr <- 0
  for(i in 1:n){
    Li <- 1
    for(j in 1:n){
      if(i != j) Li <- Li*(x-xk[j])/(xk[i]-xk[j])
    }
    lagr <- Li*yk[i]+lagr
  }
  return(lagr)
}

# name5 <- c("Bl", "Bu")
for(k in name4){
  # k="Az"
  gd=gd1=gd2=gd3=NULL
  gd <- read.table(paste0("binall.geneticDistance.",k,".chisq0.05"), header = T, stringsAsFactors = F)
  gd1 <- read.table(paste0("binall.geneticDistance.",k,".chisq0.0"), header = T, stringsAsFactors = F)
  gd2 <- merge(gd1[,c(1:4)], gd[,c(1,5:8)],  by.x = "SNP", by.y = "SNP", all.x = T)
  gd2 <- gd2[order(gd2[,2], gd2[,3]), ]
  gd3 <- gd2[which(is.na(gd2[,5])), ]
  gd3.all=i=j=NULL
  for(i in 1:12){
    #i=8
    gdchr=gd3chr=NULL
    gdchr <- gd[gd[,2]==i, ]
    gd3chr <- gd3[gd3[,2]==i, ]
    if(nrow(gd3chr)>0){ # 判断gd3chr某染色体没有snp的
      for(j in 1:nrow(gd3chr)){
        #j=1
        leftGD=rightGD=leftPos=rightPos=NULL
        leftPos <- max(gdchr[which(gdchr[,4] < gd3chr[j,4]), 4])
        rightPos <- min(gdchr[which(gdchr[,4] > gd3chr[j,4]), 4])
        leftGD <- max(gdchr[which(gdchr[,4] < gd3chr[j,4]), 7])
        rightGD <- min(gdchr[which(gdchr[,4] > gd3chr[j,4]), 7])
        gd3chr[j, 7] <- lagrange(x = gd3chr[j,4], xk = c(leftPos, rightPos), yk = c(leftGD,rightGD))
      }
      cat(i,"\t")
      gd3.all <- rbind(gd3.all, gd3chr)
    }
  }
  gd.all=NULL
  gd.all <- rbind(gd, gd3.all)
  gd.all <- gd.all[order(gd.all[,2], gd.all[,3]), ]
  
  gd.all.2=i=j=NULL
  for(i in 1:12){
    gd.all.chr=NULL
    gd.all.chr <- gd.all[gd.all[,2]==i, ]
    for(j in 2:nrow(gd.all.chr)){ # 调整interval
      gd.all.chr[j, 6] <- gd.all.chr[j,7]-gd.all.chr[j-1,7]
    }
    gd.all.2 <- rbind(gd.all.2, gd.all.chr)
  }
  write.table(gd.all.2, paste0("binall.geneticDistance.",k,".ALL"), quote = F, sep = "\t", row.names = F, col.names = T)
  cat(k," done", "\n")
}
#
# 检测结果准确性 =====
for(i in name4){
  #i="Ab"
  aa=NULL
  aa <- read.table(paste0("binall.geneticDistance.",i,".ALL"), header = T)
  ifelse(length(which(is.na(aa[,6])))>0, print(paste0(i, " interval with na")), print(paste0(i, " done")))
  ifelse(length(which(is.na(aa[,7])))>0, print(paste0(i, " Accumulate with na")), print(paste0(i, " done")))
}
# ====







