# read.EIMMo_mmu.R
# To read and preprocess mouse data from EIMMo.
# 
# Status: 9/2017 - Has NOT been updated to use current annotations available through Bioconductor
#
######
#
# 1. Read data from downloaded file
#
######

org = 'mmu'
eimmo_mmu.file = "./EIMMo/v5_Jan2011/mm_targets_FullList_flat.tab"
eimmo_mmu.data = read.delim(eimmo_mmu.file)

library(data.table)

# MySQL parameters
source("build.multiMiR_DB.parameters.R")

# functions
source("multiMiR.R")

######
#
# 2. Preprocess data
#
######

dim(eimmo_mmu.data)
#[1] 2087811       9

# 2.1 split column 'Annot' into RefSeq ID & microRNA ID
eimmo_mmu.data2 = strsplit(as.character(eimmo_mmu.data$Annot), ":")
eimmo_mmu.data2 <- do.call("rbind", eimmo_mmu.data2)
eimmo_mmu.data2 = cbind(eimmo_mmu.data2, eimmo_mmu.data$p)	# 2087811 rows
colnames(eimmo_mmu.data2) = c("RefSeq","miRNA","p")
eimmo_mmu.data2 = unique(eimmo_mmu.data2)			# 1913360 rows

# 2.2 update mature miRNA IDs
eimmo_mmu.matureID = unique(eimmo_mmu.data2[,2])
eimmo_mmu.newID = update.matureID(eimmo_mmu.matureID)
# in case some new IDs == NA
eimmo_mmu.combineID = eimmo_mmu.newID
eimmo_mmu.combineID[is.na(eimmo_mmu.combineID)] = eimmo_mmu.matureID[is.na(eimmo_mmu.combineID)]
m = match(eimmo_mmu.data2[,2], eimmo_mmu.matureID)
eimmo_mmu.data2[,2] = eimmo_mmu.combineID[m]

# 2.3 get mature miRNA accession numbers
eimmo_mmu.matureACC = matureID2matureACC(eimmo_mmu.combineID)
eimmo_mmu.data2 = cbind(mature_mirna_acc=eimmo_mmu.matureACC[m], eimmo_mmu.data2)

# 2.4 get gene symbols, Entrez & Ensembl gene IDs
refseq = unique(eimmo_mmu.data2[,2])
length(refseq)
#[1] 24909
eimmo_mmu.IDs = convert.gene.IDs(org=org, ID=refseq, ID.type='refseq')
eimmo_mmu.IDs = unique(eimmo_mmu.IDs)
dim(eimmo_mmu.IDs)
#[1] 25879     4

# 2.5 combine eimmo_mmu.IDs with eimmo_mmu.data2
eimmo_mmu.data2.DT = data.table(eimmo_mmu.data2)
setkey(eimmo_mmu.data2.DT, RefSeq)
eimmo_mmu.IDs.DT = data.table(eimmo_mmu.IDs)
setkey(eimmo_mmu.IDs.DT, refseq_acc)
eimmo_mmu.data3.DT = eimmo_mmu.IDs.DT[eimmo_mmu.data2.DT, allow.cartesian=TRUE]
eimmo_mmu.data3.DT = eimmo_mmu.data3.DT[,list(mature_mirna_acc,miRNA,gene_symbol,entrez_gene_id,ensembl_gene_id,refseq_acc,p)]
eimmo_mmu.data3.DT = unique(eimmo_mmu.data3.DT)
dim(eimmo_mmu.data3.DT)
#[1] 2341389       7

eimmo_mmu = data.frame(eimmo_mmu.data3.DT)

# 2.6 replace NA by '' in miRNA and gene ID columns
for (i in 1:6) {
    tmp = as.character(eimmo_mmu[,i])
    tmp[is.na(tmp)] = ''
    eimmo_mmu[,i] = tmp
}

# 2.7 remove column refseq_acc
# NOTE: a big reduction from 2341389 to 1457665
eimmo_mmu = unique(eimmo_mmu[,-6])
dim(eimmo_mmu)
#[1] 1457665       6

# 2.8 remove records with '' in gene_symbol, entrez_gene_id & ensembl_gene_id
remove = which(eimmo_mmu$gene_symbol == '' & eimmo_mmu$entrez_gene_id == '' & eimmo_mmu$ensembl_gene_id == '')
if (length(remove) > 0) {eimmo_mmu = eimmo_mmu[-remove,]}
dim(eimmo_mmu)
#[1] 1449133       6

######
#
# 3. Populate MySQL tables mirna, target & eimmo
#
######

# add org to the table
eimmo_mmu = cbind(org=org, eimmo_mmu)

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', eimmo_mmu[,1:3])

# update target table
target.uid = update.target.table(con, 'target', eimmo_mmu[,c(1,4:6)])

# then insert eimmo_mmu into table eimmo
eimmo_mmu = data.frame(mirna.uid, target.uid, eimmo_mmu[,7])
fields = dbListFields(con, 'eimmo')
colnames(eimmo_mmu) = fields
dim(eimmo_mmu)
#[1] 1449133       3
dim(unique(eimmo_mmu))
#[1] 1449133       3
dbWriteTable(con, 'eimmo', eimmo_mmu, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 4. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.EIMMo_mmu.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(m, refseq, i, tmp, remove, drv, con)
save.image("read.EIMMo_mmu.RData")

