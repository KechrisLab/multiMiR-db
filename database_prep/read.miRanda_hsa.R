# read.miRanda_hsa.R
# To read and preprocess human data from miRanda.
# 
# Status: 9/2017 - Has NOT been updated to use current annotations available through Bioconductor
#
######
#
# 1. Read data from downloaded file
#
######

org = 'hsa'
miranda_hsa.file = "./miRanda/2010_Aug/human_predictions_S_C_aug2010.txt"
miranda_hsa.data = read.delim(miranda_hsa.file)

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
dim(miranda_hsa.data)
#[1] 1097064      19
dim(unique(miranda_hsa.data))
#[1] 1097064      19
miranda_hsa.data2 = miranda_hsa.data[,c(1:4,6,15:19)]
#miranda_hsa.data2 = miranda_hsa.data[,c(1:4,6,10:13,15:19)]
miranda_hsa.data2 = unique(miranda_hsa.data2)
dim(miranda_hsa.data2)
#[1] 1084841      10

# 2.1 update mature miRNA IDs
miranda_hsa.matureID = unique(as.character(miranda_hsa.data2$mirna_name))
miranda_hsa.newID = update.matureID(miranda_hsa.matureID) # in case some new IDs == NA
miranda_hsa.combineID = miranda_hsa.newID
miranda_hsa.combineID[is.na(miranda_hsa.combineID)] = miranda_hsa.matureID[is.na(miranda_hsa.combineID)]
m = match(miranda_hsa.data2$mirna_name, miranda_hsa.matureID)
miranda_hsa.data2$mirna_name = miranda_hsa.combineID[m]

# 2.2 get Ensembl gene IDs
# NOTE: Gene symbols and Ensembl gene IDs are retrieved using Entrez gene IDs.
#	It is not checked whether the retrieved gene symbols are the same as
#	the original gene symbols in miranda.data.
#	The retrieved gene symbols are used eventually.
miranda_hsa.IDs = NULL
entrez = unique(miranda_hsa.data2$gene_id)
miranda_hsa.IDs = convert.gene.IDs(org=org, ID=entrez, ID.type='entrez')
miranda_hsa.IDs = unique(miranda_hsa.IDs)

# 2.3 combine miranda_hsa.IDs with miranda_hsa.data2
miranda_hsa.data2.DT = data.table(miranda_hsa.data2)
setkey(miranda_hsa.data2.DT, gene_id)
miranda_hsa.IDs.DT = data.table(miranda_hsa.IDs)
setkey(miranda_hsa.IDs.DT, entrez_gene_id)
miranda_hsa.data3.DT = miranda_hsa.IDs.DT[miranda_hsa.data2.DT, allow.cartesian=TRUE]
# NOTE: column 'gene_symbol.1' in miranda_hsa.data3.DT is the same as column
# 'gene_symbol' in miranda_hsa.data2
miranda_hsa.data3.DT = miranda_hsa.data3.DT[,list(X.mirbase_acc,mirna_name,gene_symbol,gene_symbol.1,entrez_gene_id,ensembl_gene_id,ext_transcript_id,conservation,align_score,seed_cat,energy,mirsvr_score)]
miranda_hsa.data3.DT = unique(miranda_hsa.data3.DT)

miranda_hsa.data3 = data.frame(miranda_hsa.data3.DT)
dim(miranda_hsa.data3)
#[1] 1328687      12

# Use column 'gene_symbol', instead of 'gene_symbol.1', but replace NA's in
# 'gene_symbol' with corresponding records in 'gene_symbol.1'
gene_symbol = as.character(miranda_hsa.data3$gene_symbol)
gene_symbol[is.na(gene_symbol)] = as.character(miranda_hsa.data3$gene_symbol.1[is.na(gene_symbol)])
miranda_hsa.data3$gene_symbol = gene_symbol
miranda_hsa = unique(miranda_hsa.data3[,-4])
dim(miranda_hsa)
#[1] 1328687      11

# 2.4 replace NA by '' in miRNA and gene ID columns
for (i in 1:6) {
    tmp = as.character(miranda_hsa[,i])
    tmp[is.na(tmp)] = ''
    miranda_hsa[,i] = tmp
}

# 2.5 only keep records with conservation >= MIRANDA_CONSERVE
if (!is.null(MIRANDA_CONSERVE)) {
    sum(miranda_hsa$conservation >= MIRANDA_CONSERVE) / nrow(miranda_hsa)	#[1] 0.9383038
    miranda_hsa = miranda_hsa[miranda_hsa$conservation >= MIRANDA_CONSERVE,]
}
dim(miranda_hsa)
#[1] 1246712      11

# 2.6 remove a few columns (ext_transcript_id, align_score, seed_cat, energy)
miranda_hsa = unique(miranda_hsa[,c(1:5,7,11)])
dim(miranda_hsa)
#[1] 1236550       7

######
#
# 3. Populate MySQL tables mirna, target & miranda
#
######

# add org to the table
miranda_hsa = cbind(org=org, miranda_hsa)

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', miranda_hsa[,1:3])

# update target table
target.uid = update.target.table(con, 'target', miranda_hsa[,c(1,4:6)])

# then insert miranda_hsa into table miranda
miranda_hsa = cbind(mirna.uid, target.uid, miranda_hsa[,7:8])
fields = dbListFields(con, 'miranda')
colnames(miranda_hsa) = fields
dim(miranda_hsa)
#[1] 1236550       4
dim(unique(miranda_hsa))
#[1] 1236550       4
dbWriteTable(con, 'miranda', miranda_hsa, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 4. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.miRanda_hsa.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(m, entrez, gene_symbol, i, tmp, drv, con)
save.image("read.miRanda_hsa.RData")

