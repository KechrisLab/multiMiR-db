# read.PicTar_hsa.R
# To read and preprocess human data from PicTar.
# 
# Status: 9/2017 - Has NOT been updated to use current annotations available through Bioconductor
#
######
#
# 1. Read data from downloaded file
#
######

org = 'hsa'
pictar_hsa.file = "./PicTar/v2/pictar_hg19_mammals.bulk_download.csv"
pictar_hsa.data = read.csv(pictar_hsa.file)

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

# 2.1 update mature miRNA IDs
pictar_hsa.matureID = unique(as.character(pictar_hsa.data$miRNA))
pictar_hsa.newID = update.matureID(pictar_hsa.matureID)
# in case some new IDs == NA
pictar_hsa.combineID = pictar_hsa.newID
pictar_hsa.combineID[is.na(pictar_hsa.combineID)] = pictar_hsa.matureID[is.na(pictar_hsa.combineID)]
m = match(pictar_hsa.data$miRNA, pictar_hsa.matureID)
pictar_hsa.data$miRNA = pictar_hsa.combineID[m]

# 2.2 get mature miRNA accession numbers
pictar_hsa.matureACC = matureID2matureACC(pictar_hsa.combineID)
pictar_hsa.data = cbind(pictar_hsa.data, mature_mirna_acc=pictar_hsa.matureACC[m])

# 2.3 get gene symbols, Entrez & Ensembl gene IDs using RefSeq IDs
refseq = unique(as.character(pictar_hsa.data$X.id))
pictar_hsa.IDs = convert.gene.IDs(org=org, ID=refseq, ID.type='refseq')
pictar_hsa.IDs = unique(pictar_hsa.IDs)

# 2.4 combine pictar_hsa.IDs with pictar_hsa.data
pictar_hsa.data.DT = data.table(pictar_hsa.data)
setkey(pictar_hsa.data.DT, X.id)
pictar_hsa.IDs.DT = data.table(pictar_hsa.IDs)
setkey(pictar_hsa.IDs.DT, refseq_acc)
pictar_hsa.data2.DT = pictar_hsa.IDs.DT[pictar_hsa.data.DT, allow.cartesian=TRUE]
pictar_hsa.data2.DT = pictar_hsa.data2.DT[,list(mature_mirna_acc,miRNA,gene_symbol,entrez_gene_id,ensembl_gene_id,refseq_acc,score)]
pictar_hsa.data2.DT = unique(pictar_hsa.data2.DT)

pictar_hsa = data.frame(pictar_hsa.data2.DT)

# 2.5 replace NA by '' in miRNA and gene ID columns
for (i in 1:6) {
    tmp = as.character(pictar_hsa[,i])
    tmp[is.na(tmp)] = ''
    pictar_hsa[,i] = tmp
}

######
#
# 3. Remove column refseq_acc in pictar_hsa --
#    there is no refseq_acc in the MySQL target table. However, refseq_acc is
#    the only gene ID in the original PicTar data file. Removing refseq_acc will
#    remove the interactions involving genes with only refseq_acc available.
#    A better solution in the future: add refseq_acc to the MySQL target table.
#
######

# there are 12301 (1.7%) interactions involving genes with refseq_acc only.
length(which(pictar_hsa$gene_symbol=='' & pictar_hsa$entrez_gene_id=='' & pictar_hsa$ensembl_gene_id==''))
#[1] 12301
length(which(pictar_hsa$gene_symbol=='' & pictar_hsa$entrez_gene_id=='' & pictar_hsa$ensembl_gene_id=='')) / nrow(pictar_hsa)
#[1] 0.01703978

# a big drop from 721899 to 413206
dim(pictar_hsa)
#[1] 721899      7
pictar_hsa = pictar_hsa[,-6]
pictar_hsa = unique(pictar_hsa)
dim(pictar_hsa)
#[1] 413206      6

remove = which(pictar_hsa$gene_symbol=='' & pictar_hsa$entrez_gene_id=='' & pictar_hsa$ensembl_gene_id=='')
if (length(remove) > 0) {pictar_hsa = pictar_hsa[-remove,]}
dim(pictar_hsa)
#[1] 404066      6

######
#
# 4. Populate MySQL tables mirna, target & pictar
#
######

# add org to the table
pictar_hsa = cbind(org=org, pictar_hsa)

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', pictar_hsa[,1:3])

# update target table
target.uid = update.target.table(con, 'target', pictar_hsa[,c(1,4:6)])

# then insert pictar_hsa into table pictar
pictar_hsa = data.frame(mirna.uid, target.uid, pictar_hsa[,7])
fields = dbListFields(con, 'pictar')
colnames(pictar_hsa) = fields
dbWriteTable(con, 'pictar', pictar_hsa, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 5. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.PicTar_hsa.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(m, refseq, i, tmp, remove, drv, con)
save.image("read.PicTar_hsa.RData")


