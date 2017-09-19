# read.miRDB.R
# To read and preprocess data from miRDB.
# 
# Status: 9/2017 - Has been updated to use current annotations available through Bioconductor
#

######
#
# 1. Load data, parameters & libraries
#
######

mirdb.file = "./miRDB_v5.0_prediction_result.txt"
mirdb.data = read.delim(mirdb.file, header=F)
colnames(mirdb.data) = c('mirna','refseq_acc','score')

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
mirdb.matureID = unique(as.character(mirdb.data$mirna))
mirdb.newID = update.matureID(mirdb.matureID)
# in case some new IDs == NA
mirdb.combineID = mirdb.newID
mirdb.combineID[is.na(mirdb.combineID)] = mirdb.matureID[is.na(mirdb.combineID)]
m = match(mirdb.data$mirna, mirdb.matureID)
mirdb.data$mirna = mirdb.combineID[m]

# 2.2 get mature miRNA accession numbers
mirdb.matureACC = matureID2matureACC(mirdb.combineID)
mirdb.data = cbind(mature_mirna_acc=mirdb.matureACC[m], mirdb.data)

# 2.3 get gene symbols, Entrez & Ensembl gene IDs using RefSeq IDs
mirdb.org = unique(substr(mirdb.data$mirna, 1, 3))
mirdb.IDs = NULL
for (org in mirdb.org) {
    m = grep(org, mirdb.data$mirna)
    refseq = unique(as.character(mirdb.data$refseq_acc[m]))
    ID = convert.gene.IDs(org=org, ID=refseq, ID.type='refseq')
    mirdb.IDs = rbind(mirdb.IDs, ID)
}
mirdb.IDs = unique(mirdb.IDs)
dim(mirdb.IDs)
#[1] 92811     4

# 2.6 combine mirdb.IDs with mirdb.data
mirdb.data.DT = data.table(mirdb.data)
setkey(mirdb.data.DT, refseq_acc)
mirdb.IDs.DT = data.table(mirdb.IDs)
setkey(mirdb.IDs.DT, refseq_acc)
mirdb.data2.DT = mirdb.IDs.DT[mirdb.data.DT, allow.cartesian=TRUE]
mirdb.data2.DT = mirdb.data2.DT[,list(mature_mirna_acc,mirna,gene_symbol,entrez_gene_id,ensembl_gene_id,refseq_acc,score)]
mirdb.data2.DT = unique(mirdb.data2.DT)

mirdb.all = data.frame(mirdb.data2.DT)

# 2.7 species-specific predictions -- only keep human & mouse for now
#mirdb_cfa = mirdb[grep('cfa', mirdb$mirna),]
#mirdb_cfa = cbind(org='Canis lupus familiaris', mirdb_cfa)
#mirdb_gga = mirdb[grep('gga', mirdb$mirna),]
#mirdb_gga = cbind(org='Gallus gallus', mirdb_gga)
mirdb_rno = mirdb.all[grep('rno', mirdb.all$mirna),]
mirdb_rno = cbind(org='rno', mirdb_rno)
mirdb_mmu = mirdb.all[grep('mmu', mirdb.all$mirna),]
mirdb_mmu = cbind(org='mmu', mirdb_mmu)
mirdb_hsa = mirdb.all[grep('hsa', mirdb.all$mirna),]
mirdb_hsa = cbind(org='hsa', mirdb_hsa)
#mirdb = rbind(mirdb_hsa, mirdb_mmu, mirdb_rno, mirdb_gga, mirdb_cfa)
mirdb = rbind(mirdb_hsa, mirdb_mmu,mirdb_rno)

# 2.8 replace NA by '' in miRNA and gene ID columns
for (i in 2:7) {
    tmp = as.character(mirdb[,i])
    tmp[is.na(tmp)] = ''
    mirdb[,i] = tmp
}

######
#
# 3. Remove column refseq_acc in mirdb --
#    there is no refseq_acc in the MySQL target table. However, refseq_acc is
#    the only gene ID in the original miRDB data file. Removing refseq_acc will
#    remove the interactions involving genes with only refseq_acc available.
#    A better solution in the future: add refseq_acc to the MySQL target table.
#
######

# there are 42032 (2.4%) interactions involving genes with refseq_acc only.
length(which(mirdb$gene_symbol=='' & mirdb$entrez_gene_id=='' & mirdb$ensembl_gene_id==''))
#[1] 42032
length(which(mirdb$gene_symbol=='' & mirdb$entrez_gene_id=='' & mirdb$ensembl_gene_id=='')) / nrow(mirdb)
#[1] 0.023847

# a big drop from 1762570 to 1100327
dim(mirdb)
#[1] 1762570       8
mirdb = mirdb[,-7]
mirdb = unique(mirdb)
dim(mirdb)
#[1] 1100327       7

remove = which(mirdb$gene_symbol=='' & mirdb$entrez_gene_id=='' & mirdb$ensembl_gene_id=='')
if (length(remove) > 0) {mirdb = mirdb[-remove,]}
dim(mirdb)
#[1] 1071047       7

######
#
# 4. Populate MySQL tables mirna, target & mirdb
#
######

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', mirdb[,1:3])

# update target table
target.uid = update.target.table(con, 'target', mirdb[,c(1,4:6)])

# then insert mirdb into table mirdb
mirdb = data.frame(mirna.uid, target.uid, mirdb[,7])
fields = dbListFields(con, 'mirdb')
colnames(mirdb) = fields
dbWriteTable(con, 'mirdb', mirdb, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 5. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.miRDB.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(m, org, refseq, ID, i, tmp, remove, drv, con)
save.image("./new.read.miRDB.RData")

