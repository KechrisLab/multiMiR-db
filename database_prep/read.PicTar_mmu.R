# read.PicTar_mmu.R
# To read and preprocess mouse data from PicTar.
# 
# Status: 9/2017 - Has NOT been updated to use current annotations available through Bioconductor
#
######
#
# 1. Read data from downloaded file
#
######

org = 'mmu'
pictar_mmu.file = "./PicTar/v2/pictar_mm9_mammals.bulk_download.csv"
pictar_mmu.data = read.csv(pictar_mmu.file)

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
pictar_mmu.matureID = unique(as.character(pictar_mmu.data$miRNA))
pictar_mmu.newID = update.matureID(pictar_mmu.matureID)
# in case some new IDs == NA
pictar_mmu.combineID = pictar_mmu.newID
pictar_mmu.combineID[is.na(pictar_mmu.combineID)] = pictar_mmu.matureID[is.na(pictar_mmu.combineID)]
m = match(pictar_mmu.data$miRNA, pictar_mmu.matureID)
pictar_mmu.data$miRNA = pictar_mmu.combineID[m]

# 2.2 get mature miRNA accession numbers
pictar_mmu.matureACC = matureID2matureACC(pictar_mmu.combineID)
pictar_mmu.data = cbind(pictar_mmu.data, mature_mirna_acc=pictar_mmu.matureACC[m])

# 2.3 get gene symbols, Entrez & Ensembl gene IDs using RefSeq IDs
refseq = unique(as.character(pictar_mmu.data$X.id))
pictar_mmu.IDs = convert.gene.IDs(org=org, ID=refseq, ID.type='refseq')
pictar_mmu.IDs = unique(pictar_mmu.IDs)

# 2.4 combine pictar_mmu.IDs with pictar_mmu.data
pictar_mmu.data.DT = data.table(pictar_mmu.data)
setkey(pictar_mmu.data.DT, X.id)
pictar_mmu.IDs.DT = data.table(pictar_mmu.IDs)
setkey(pictar_mmu.IDs.DT, refseq_acc)
pictar_mmu.data2.DT = pictar_mmu.IDs.DT[pictar_mmu.data.DT, allow.cartesian=TRUE]
pictar_mmu.data2.DT = pictar_mmu.data2.DT[,list(mature_mirna_acc,miRNA,gene_symbol,entrez_gene_id,ensembl_gene_id,refseq_acc,score)]
pictar_mmu.data2.DT = unique(pictar_mmu.data2.DT)

pictar_mmu = data.frame(pictar_mmu.data2.DT)

# 2.5 replace NA by '' in miRNA and gene ID columns
for (i in 1:6) {
    tmp = as.character(pictar_mmu[,i])
    tmp[is.na(tmp)] = ''
    pictar_mmu[,i] = tmp
}

######
#
# 3. Remove column refseq_acc in pictar_mmu --
#    there is no refseq_acc in the MySQL target table. However, refseq_acc is
#    the only gene ID in the original PicTar data file. Removing refseq_acc will
#    remove the interactions involving genes with only refseq_acc available.
#    A better solution in the future: add refseq_acc to the MySQL target table.
#
######

# there are 3982 (0.8%) interactions involving genes with refseq_acc only.
length(which(pictar_mmu$gene_symbol=='' & pictar_mmu$entrez_gene_id=='' & pictar_mmu$ensembl_gene_id==''))
#[1] 3982
length(which(pictar_mmu$gene_symbol=='' & pictar_mmu$entrez_gene_id=='' & pictar_mmu$ensembl_gene_id=='')) / nrow(pictar_mmu)
#[1] 0.007951945

# a big drop from 500758 to 304834
dim(pictar_mmu)
#[1] 500758      7
pictar_mmu = pictar_mmu[,-6]
pictar_mmu = unique(pictar_mmu)
dim(pictar_mmu)
#[1] 304834      6

remove = which(pictar_mmu$gene_symbol=='' & pictar_mmu$entrez_gene_id=='' & pictar_mmu$ensembl_gene_id=='')
if (length(remove) > 0) {pictar_mmu = pictar_mmu[-remove,]}
dim(pictar_mmu)
#[1] 302236      6

######
#
# 4. Populate MySQL tables mirna, target & mirecords
#
######

# add org to the table
pictar_mmu = cbind(org=org, pictar_mmu)

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', pictar_mmu[,1:3])

# update target table
target.uid = update.target.table(con, 'target', pictar_mmu[,c(1,4:6)])

# then insert pictar_mmu into table pictar_mmu
pictar_mmu = data.frame(mirna.uid, target.uid, pictar_mmu[,7])
fields = dbListFields(con, 'pictar')
colnames(pictar_mmu) = fields
dbWriteTable(con, 'pictar', pictar_mmu, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 5. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.PicTar_mmu.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(m, refseq, i, tmp, remove, drv, con)
save.image("read.PicTar_mmu.RData")


