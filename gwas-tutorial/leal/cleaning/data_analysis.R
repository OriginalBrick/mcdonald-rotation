# Sex Check
sexcheck = read.table("GWAS_sex_checking.sexcheck", header=T)
names(sexcheck)
sex_problem = sexcheck[which(sexcheck$STATUS=="PROBLEM"),]
sex_problem

# Duplicates
dups = read.table("duplicates.genome", header = T)
problem_pairs = dups[which(dups$PI_HAT > 0.4),]
problem_pairs