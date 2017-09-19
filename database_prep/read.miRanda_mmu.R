# read.miRanda_mmu.R
# To read and preprocess mouse data from miRanda.
# 
# Status: 9/2017 - Has NOT been updated to use current annotations available through Bioconductor
#
######
#
# 1. Read data from downloaded file
#
######

org = 'mmu'
miranda_mmu.file = "./miRanda/2010_Aug/mouse_predictions_S_C_aug2010.txt"
miranda_mmu.data = read.delim(miranda_mmu.file)

library(data.table)

# MySQL parameters
source("build.multiMiR_DB.parameters.R")

# functions
source("multiMiR.R")

MIRANDA_CONSERVE = 0.5

######
#
# 2. Preprocess data
#
######

# 2.0 remove some columns
dim(miranda_mmu.data)
#[1] 819078     19
dim(unique(miranda_mmu.data))
#[1] 819078     19
miranda_mmu.data2 = miranda_mmu.data[,c(1:4,6,15:19)]
#miranda_mmu.data2 = miranda_mmu.data[,c(1:4,6,10:13,15:19)]
miranda_mmu.data2 = unique(miranda_mmu.data2)
dim(miranda_mmu.data2)
#[1] 817512     10

# 2.1 update mature miRNA IDs
miranda_mmu.matureID = unique(as.character(miranda_mmu.data2$mirna_name))
miranda_mmu.newID = update.matureID(miranda_mmu.matureID) # in case some new IDs == NA
miranda_mmu.combineID = miranda_mmu.newID
miranda_mmu.combineID[is.na(miranda_mmu.combineID)] = miranda_mmu.matureID[is.na(miranda_mmu.combineID)]
m = match(miranda_mmu.data2$mirna_name, miranda_mmu.matureID)
miranda_mmu.data2$mirna_name = miranda_mmu.combineID[m]

# 2.2 get Ensembl gene IDs
# NOTE: Gene symbols and Ensembl gene IDs are retrieved using Entrez gene IDs.
#	It is not checked whether the retrieved gene symbols are the same as
#	the original gene symbols in miranda.data.
#	The retrieved gene symbols are used eventually.
miranda_mmu.IDs = NULL
entrez = unique(miranda_mmu.data2$gene_id)
miranda_mmu.IDs = convert.gene.IDs(org=org, ID=entrez, ID.type='entrez')
miranda_mmu.IDs = unique(miranda_mmu.IDs)

# 2.3 combine miranda_mmu.IDs with miranda_mmu.data2
miranda_mmu.data2.DT = data.table(miranda_mmu.data2)
setkey(miranda_mmu.data2.DT, gene_id)
miranda_mmu.IDs.DT = data.table(miranda_mmu.IDs)
setkey(miranda_mmu.IDs.DT, entrez_gene_id)
miranda_mmu.data3.DT = miranda_mmu.IDs.DT[miranda_mmu.data2.DT, allow.cartesian=TRUE]
# NOTE: column 'gene_symbol.1' in miranda_mmu.data3.DT is the same as column
# 'gene_symbol' in miranda_mmu.data2
miranda_mmu.data3.DT = miranda_mmu.data3.DT[,list(X.mirbase_acc,mirna_name,gene_symbol,gene_symbol.1,entrez_gene_id,ensembl_gene_id,ext_transcript_id,conservation,align_score,seed_cat,energy,mirsvr_score)]
miranda_mmu.data3.DT = unique(miranda_mmu.data3.DT)

miranda_mmu.data3 = data.frame(miranda_mmu.data3.DT)
dim(miranda_mmu.data3)
#[1] 846817     12

# Use column 'gene_symbol', instead of 'gene_symbol.1', but replace NA's in
# 'gene_symbol' with corresponding records in 'gene_symbol.1'
gene_symbol = as.character(miranda_mmu.data3$gene_symbol)
gene_symbol[is.na(gene_symbol)] = as.character(miranda_mmu.data3$gene_symbol.1[is.na(gene_symbol)])
miranda_mmu.data3$gene_symbol = gene_symbol
miranda_mmu = unique(miranda_mmu.data3[,-4])
dim(miranda_mmu)
#[1] 846817     11

# 2.4 replace NA by '' in miRNA and gene ID columns
for (i in 1:6) {
    tmp = as.character(miranda_mmu[,i])
    tmp[is.na(tmp)] = ''
    miranda_mmu[,i] = tmp
}

# 2.5 only keep records with conservation >= MIRANDA_CONSERVE
if (!is.null(MIRANDA_CONSERVE)) {
    sum(miranda_mmu$conservation >= MIRANDA_CONSERVE) / nrow(miranda_mmu)	# 0.9180827
    miranda_mmu = miranda_mmu[miranda_mmu$conservation >= MIRANDA_CONSERVE,]
}
dim(miranda_mmu)
#[1] 777448     11

# 2.6 remove a few columns (ext_transcript_id, align_score, seed_cat, energy)
miranda_mmu = unique(miranda_mmu[,c(1:5,7,11)])
dim(miranda_mmu)
#[1] 776936      7

######
#
# 3. Populate MySQL tables mirna, target & miranda
#
######

# add org to the table
miranda_mmu = cbind(org=org, miranda_mmu)

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', miranda_mmu[,1:3])

# update target table
target.uid = update.target.table(con, 'target', miranda_mmu[,c(1,4:6)])

# then insert miranda_mmu into table miranda
miranda_mmu = cbind(mirna.uid, target.uid, miranda_mmu[,7:8])
fields = dbListFields(con, 'miranda')
colnames(miranda_mmu) = fields
dim(miranda_mmu)
#[1] 776936      4
dim(unique(miranda_mmu))
#[1] 776936      4
dbWriteTable(con, 'miranda', miranda_mmu, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 4. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.miRanda_mmu.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(m, entrez, gene_symbol, i, tmp, drv, con)
save.image("read.miRanda_mmu.RData")

