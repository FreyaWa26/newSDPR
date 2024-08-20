
geno = read.table("tmp.vcf")
geno1 = sapply(geno[,10:ncol(geno)], function(x) substr(x, 1, 1))
geno2 = sapply(geno[,10:ncol(geno)], function(x) substr(x, 3, 3))

a = read.table("tmp.tsv", header=F)
b = as.matrix(a[,7:ncol(a)])
a[nrow(a), 6] = a[nrow(a), 6] + 1

res = list()
for (i in 1:nrow(a)) {
        res[[i]] = do.call(rbind, replicate(a[i,6], b[i,], simplify=F))
}
anc = do.call(rbind, res)
print(dim(geno1))
print(dim(anc))

idx_hap1 = seq(1, ncol(anc), 2)
idx_hap2 = seq(2, ncol(anc), 2)

anc1 = anc[, idx_hap1]
anc2 = anc[, idx_hap2]

geno_pop1 = matrix(0, nrow(geno), ncol(geno1))
geno_pop2 = matrix(0, nrow(geno), ncol(geno2))

for (i in 1:nrow(geno1)) {
    idx = which(anc1[i,] == "0")
    geno_pop1[i, idx] = as.numeric(geno1[i, idx])
    geno_pop2[i, -idx] = as.numeric(geno1[i, -idx])

    idx = which(anc2[i,] == "0")
    geno_pop1[i, idx] = geno_pop1[i, idx] + as.numeric(geno2[i, idx])
    geno_pop2[i, -idx] = geno_pop2[i, -idx] + as.numeric(geno2[i, -idx])
}

#write.table(geno_pop1, file="X_afr1.txt", row.names=F, col.names=F, quote=F, sep=" ", append=F)
#write.table(geno_pop2, file="X_eur1.txt", row.names=F, col.names=F, quote=F, sep=" ", append=F)

X_afr = read.table("X_afr.txt", header=F)
X_eur = read.table("X_eur.txt", header=F)
X_afr = apply(X_afr, 2, as.numeric)
X_eur = apply(X_eur, 2, as.numeric)
pheno = read.table("../phenotype/SDPR_admix/height.txt", header=F)
idx = which(!is.na(pheno$V3))
all.equal(X_afr, geno_pop1[,idx], check.attributes=F)
all.equal(X_eur, geno_pop2[,idx], check.attributes=F)








