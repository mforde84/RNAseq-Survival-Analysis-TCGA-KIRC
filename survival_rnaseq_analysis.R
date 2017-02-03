library(survival)
library(limma)

# input counts, and filter genes whose expression is zero in more than half the samples
rna <- as.matrix(read.table("KIRC.normalized.data.txt", header=T, row.names=1, sep="\t"))
x <- t(apply(rna,1,as.numeric))
r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
remove <- which(r > dim(rna)[2]*0.5)
rna <- rna[-remove,]

# get the index of the normal/control samples
n_index <- which(substr(colnames(rna),14,14) == "1")
t_index <- which(substr(colnames(rna),14,14) == "0")

# voom normalization
cond <- factor(ifelse(seq(1, dim(rna)[2],1) %in% t_index, 1,  0))
d <- model.matrix(~1 + cond)
x <- t(apply(rna, 1, as.numeric))
rna_vm <- voom(x, d, plot=FALSE)$E
colnames(rna_vm) <- gsub("\\.","-",substr(colnames(rna),1,12))

# zscore scaling against normals
mean_n <- rowMeans(rna_vm[, n_index])
sd_n <- apply(rna_vm[, n_index], 1, sd)
tumor <- rna_vm[,t_index]
z_rna <- matrix(nrow=nrow(tumor), ncol=ncol(tumor))
colnames(z_rna) <- colnames(rna_vm[,t_index])
rownames(z_rna) <- rownames(rna_vm[,t_index])
rownames(z_rna) <- sapply(rownames(z_rna), function(x) unlist(strsplit(x,"\\|"))[[1]])
for(i in 1:nrow(tumor)){
 for(j in 1:ncol(tumor)){
  z_rna[i,j] <- (tumor[i, j] - mean_n[i]) / sd_n[i]
 }
}

# input clinical information
clinical <- read.table("KIRC.clinical.csv", header=TRUE, sep="\t")
all_clin <- data.frame(cbind(clinical[,2],clinical[,3],clinical[,3]))
colnames(all_clin) <- c("new_tumor_days", "death_days", "followUp_days")
rownames(all_clin) <- clinical$IDs

# time to tumor 
all_clin$new_time <- c()
for (i in 1:length(as.numeric(as.character(all_clin$new_tumor_days)))){
 all_clin$new_time[i] <- ifelse(is.na(as.numeric(as.character(all_clin$new_tumor_days))[i]), 
  as.numeric(as.character(all_clin$followUp_days))[i],
  as.numeric(as.character(all_clin$new_tumor_days))[i])
}

# time to death 
all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
 all_clin$new_death[i] <- ifelse(is.na(as.numeric(as.character(all_clin$death_days))[i]),
  as.numeric(as.character(all_clin$followUp_days))[i],
  as.numeric(as.character(all_clin$death_days))[i])
}

# death censor event
all_clin$death_event <- ifelse(clinical$vital == "alive", 0, 1)

# create event vector for RNASeq data
# where signficantly differentially enriched gene > 1.96 zscore
event_rna <- t(apply(z_rna, 1, function(x) ifelse(abs(x) > 1.96,1,0)))

# indices for matched samples
ind_tum <- which(unique(colnames(z_rna)) %in% rownames(all_clin))
ind_clin <- which(rownames(all_clin) %in% colnames(z_rna))

# gene of interest
ind_gene <- which(rownames(z_rna) == "TP53")

# survival analysis
s <- survfit(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin], 
 all_clin$death_event[ind_clin]) ~ event_rna[ind_gene, ind_tum])
s1 <- tryCatch(
 survdiff(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin], 
  all_clin$death_event[ind_clin]) ~ event_rna[ind_gene, ind_tum]), 
  error = function(e) 
  return(NA)
)

# extract p.value
pv <- ifelse(is.na(s1), next, (round(1 - pchisq(s1$chisq, length(s1$n) - 1),3)))[[1]]

# plot curves
plot(s, col=c(1:3), frame=F, lwd=2,main=paste("KIRK",rownames(z_rna)[ind_gene],sep="\n"))

# add lines for the median survival
x1 <- ifelse(is.na(as.numeric(summary(s)$table[,'median'][1])), "NA", 
 as.numeric(summary(s)$table[,'median'][1]))
x2 <- as.numeric(summary(s)$table[,'median'][2])
if(x1 != "NA" || x2 != "NA"){
 lines(c(0, x1), c(0.5, 0.5), col="blue")
 lines(c(x1, x1), c(0, 0.5), col="black")
 lines(c(x2, x2), c(0, 0.5), col="red")
}

# add legend
legend(1800, 0.995, 
 legend=paste("p.value = ", pv[[1]], sep=""), 
 bty="n", 
 cex=1.4)
legend(max(as.numeric(as.character(all_clin$death_days)[ind_clin]), na.rm = T) * 0.7, 0.94, 
 legend=c(paste("NotAltered=",x1), paste("Altered=",x2)), 
 bty="n",
 cex=1.3, 
 lwd=3, 
 col=c("black","red"))
