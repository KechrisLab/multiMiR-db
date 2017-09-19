# read.PITA_mmu.R
# To read and preprocess mouse data from PITA.
# 
# Status: 9/2017 - Has NOT been updated to use current annotations available through Bioconductor
#
######
#
# 1. Read data from downloaded file
#
######

org = 'mmu'
pita_mmu.file = "./PITA/v6/PITA_sites_hg18_0_0_ALL.tab"
pita_mmu.data = read.delim(pita_mmu.file)

PITA_CONSERVE = 0.5
#PITA_CONSERVE = NULL

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

# 2.1 remove some columns
pita_mmu.data = pita_mmu.data[,c(1:3,5,8,9)]

# 2.2 update mature miRNA IDs
pita_mmu.matureID = unique(as.character(pita_mmu.data$microRNA))
pita_mmu.newID = update.matureID(pita_mmu.matureID)
# in case some new IDs == NA
pita_mmu.combineID = pita_mmu.newID
pita_mmu.combineID[is.na(pita_mmu.combineID)] = pita_mmu.matureID[is.na(pita_mmu.combineID)]
m = match(pita_mmu.data$microRNA, pita_mmu.matureID)
pita_mmu.data$microRNA = pita_mmu.combineID[m]

# 2.3 get mature miRNA accession numbers
pita_mmu.matureACC = matureID2matureACC(pita_mmu.combineID)
pita_mmu.data = cbind(mature_mirna_acc=pita_mmu.matureACC[m], pita_mmu.data)

# 2.4 NOTE: in column RefSeq, some have multiple accession numbers separated by
#     ';', "NM_000018;NM_001033859" for example.
#     ONLY keep the first accession number.
length(grep(";", pita_mmu.data$RefSeq))
#[1] 557923
pita_mmu.data$RefSeq = gsub(";\\w*", "", as.character(pita_mmu.data$RefSeq), perl=TRUE)

# 2.5 get gene symbols, Entrez & Ensembl gene IDs using RefSeq IDs
refseq = unique(as.character(pita_mmu.data$RefSeq))
pita_mmu.IDs = convert.gene.IDs(org=org, ID=refseq, ID.type='refseq')
pita_mmu.IDs = unique(pita_mmu.IDs)

# 2.6 combine pita_mmu.IDs with pita_mmu.data
pita_mmu.data.DT = data.table(pita_mmu.data)
setkey(pita_mmu.data.DT, RefSeq)
pita_mmu.IDs.DT = data.table(pita_mmu.IDs)
setkey(pita_mmu.IDs.DT, refseq_acc)
pita_mmu.data2.DT = pita_mmu.IDs.DT[pita_mmu.data.DT, allow.cartesian=TRUE]
pita_mmu.data2.DT = pita_mmu.data2.DT[,list(mature_mirna_acc,microRNA,Name,gene_symbol,entrez_gene_id,ensembl_gene_id,refseq_acc,Seed,ddG,Conservation)]
pita_mmu.data2.DT = unique(pita_mmu.data2.DT)
pita_mmu.data2 = data.frame(pita_mmu.data2.DT)
dim(pita_mmu.data2)
#[1] 5163467      10

# 2.7 Use column 'gene_symbol', instead of 'Name', but replace NA's in
#     'gene_symbol' with corresponding records in 'Name'
sum(is.na(pita_mmu.data2$gene_symbol))
#[1] 67002
gene_symbol = as.character(pita_mmu.data2$gene_symbol)
gene_symbol[is.na(gene_symbol)] = as.character(pita_mmu.data2$Name[is.na(gene_symbol)])
pita_mmu.data2$gene_symbol = gene_symbol
pita_mmu = unique(pita_mmu.data2[,-3])
dim(pita_mmu)
#[1] 5163467       9

# 2.8 replace NA by '' in miRNA and gene ID columns
for (i in 1:6) {
    tmp = as.character(pita_mmu[,i])
    tmp[is.na(tmp)] = ''
    pita_mmu[,i] = tmp
}

# 2.9 remove columns 'refseq_acc' & 'Seed'
pita_mmu = unique(pita_mmu[,-c(6,7)])
dim(pita_mmu)
#[1] 5163153       7

# 2.10 remove records with '' in gene_symbol, entrez_gene_id & ensembl_gene_id
remove = which(pita_mmu$gene_symbol == '' & pita_mmu$entrez_gene_id == '' & pita_mmu$ensembl_gene_id == '')
if (length(remove) > 0) {pita_mmu = pita_mmu[-remove,]}
dim(pita_mmu)
#[1] 5163153       7

# 2.11 remove records with conservation < PITA_CONSERVE
if (!is.null(PITA_CONSERVE)) {
    pita_mmu = pita_mmu[pita_mmu$Conservation >= PITA_CONSERVE,]
}
dim(pita_mmu)
#[1] 1321249       7

######
#
# 3. Populate MySQL tables mirna, target & pita
#
######

# add org to the table
pita_mmu = cbind(org=org, pita_mmu)

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', pita_mmu[,1:3])

# update target table
target.uid = update.target.table(con, 'target', pita_mmu[,c(1,4:6)])

# then insert pita_mmu into table pita
pita_mmu = data.frame(mirna.uid, target.uid, pita_mmu[,7:8])
fields = dbListFields(con, 'pita')
colnames(pita_mmu) = fields
dim(pita_mmu)
#[1] 1321249       4
dim(unique(pita_mmu))
#[1] 1321249       4
dbWriteTable(con, 'pita', pita_mmu, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 4. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.PITA_mmu.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(m, refseq, gene_symbol, i, tmp, remove, drv, con)
save.image("read.PITA_mmu.RData")

