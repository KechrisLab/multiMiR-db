# read.EIMMo_hsa.R
# To read and preprocess human data from EIMMo.
# 
# Status: 9/2017 - Has NOT been updated to use current annotations available through Bioconductor
#
######
#
# 1. Read data from downloaded file
#
######

org = 'hsa'
eimmo_hsa.file = "./EIMMo/v5_Jan2011/hg_targets_FullList_flat.tab"
eimmo_hsa.data = read.delim(eimmo_hsa.file)

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

dim(eimmo_hsa.data)
#[1] 7042277       9

# 2.1 split column 'Annot' into RefSeq ID & microRNA ID
eimmo_hsa.data2 = strsplit(as.character(eimmo_hsa.data$Annot), ":")
eimmo_hsa.data2 <- do.call("rbind", eimmo_hsa.data2)
eimmo_hsa.data2 = cbind(eimmo_hsa.data2, eimmo_hsa.data$p)	# 7042277 rows
colnames(eimmo_hsa.data2) = c("RefSeq","miRNA","p")
eimmo_hsa.data2 = unique(eimmo_hsa.data2)			# 6412704 rows

# 2.2 update mature miRNA IDs
eimmo_hsa.matureID = unique(eimmo_hsa.data2[,2])
eimmo_hsa.newID = update.matureID(eimmo_hsa.matureID)
# in case some new IDs == NA
eimmo_hsa.combineID = eimmo_hsa.newID
eimmo_hsa.combineID[is.na(eimmo_hsa.combineID)] = eimmo_hsa.matureID[is.na(eimmo_hsa.combineID)]
m = match(eimmo_hsa.data2[,2], eimmo_hsa.matureID)
eimmo_hsa.data2[,2] = eimmo_hsa.combineID[m]

# 2.3 get mature miRNA accession numbers
eimmo_hsa.matureACC = matureID2matureACC(eimmo_hsa.combineID)
eimmo_hsa.data2 = cbind(mature_mirna_acc=eimmo_hsa.matureACC[m], eimmo_hsa.data2)

# 2.4 get gene symbols, Entrez & Ensembl gene IDs
refseq = unique(eimmo_hsa.data2[,2])
length(refseq)
#[1] 31303
eimmo_hsa.IDs = convert.gene.IDs(org=org, ID=refseq, ID.type='refseq')
eimmo_hsa.IDs = unique(eimmo_hsa.IDs)
dim(eimmo_hsa.IDs)
#[1] 31408     4

# 2.5 combine eimmo_hsa.IDs with eimmo_hsa.data2
eimmo_hsa.data2.DT = data.table(eimmo_hsa.data2)
setkey(eimmo_hsa.data2.DT, RefSeq)
eimmo_hsa.IDs.DT = data.table(eimmo_hsa.IDs)
setkey(eimmo_hsa.IDs.DT, refseq_acc)
eimmo_hsa.data3.DT = eimmo_hsa.IDs.DT[eimmo_hsa.data2.DT, allow.cartesian=TRUE]
eimmo_hsa.data3.DT = eimmo_hsa.data3.DT[,list(mature_mirna_acc,miRNA,gene_symbol,entrez_gene_id,ensembl_gene_id,refseq_acc,p)]
eimmo_hsa.data3.DT = unique(eimmo_hsa.data3.DT)
dim(eimmo_hsa.data3.DT)
#[1] 6624480       7

eimmo_hsa = data.frame(eimmo_hsa.data3.DT)

# 2.6 replace NA by '' in miRNA and gene ID columns
for (i in 1:6) {
    tmp = as.character(eimmo_hsa[,i])
    tmp[is.na(tmp)] = ''
    eimmo_hsa[,i] = tmp
}

# 2.7 remove column refseq_acc
# NOTE: a big reduction from 6624480 to 3973770
eimmo_hsa = unique(eimmo_hsa[,-6])
dim(eimmo_hsa)
#[1] 3973770       6

# 2.8 remove records with '' in gene_symbol, entrez_gene_id & ensembl_gene_id
remove = which(eimmo_hsa$gene_symbol == '' & eimmo_hsa$entrez_gene_id == '' & eimmo_hsa$ensembl_gene_id == '')
if (length(remove) > 0) {eimmo_hsa = eimmo_hsa[-remove,]}
dim(eimmo_hsa)
#[1] 3959112       6

######
#
# 3. Populate MySQL tables mirna, target & eimmo
#
######

# add org to the table
eimmo_hsa = cbind(org=org, eimmo_hsa)

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', eimmo_hsa[,1:3])

# update target table
target.uid = update.target.table(con, 'target', eimmo_hsa[,c(1,4:6)])

# then insert eimmo_hsa into table eimmo
eimmo_hsa = data.frame(mirna.uid, target.uid, eimmo_hsa[,7])
fields = dbListFields(con, 'eimmo')
colnames(eimmo_hsa) = fields
dim(eimmo_hsa)
#[1] 3959112       3
dim(unique(eimmo_hsa))
#[1] 3959112       3
dbWriteTable(con, 'eimmo', eimmo_hsa, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 4. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.EIMMo_hsa.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(m, refseq, i, tmp, remove, drv, con)
save.image("read.EIMMo_hsa.RData")

