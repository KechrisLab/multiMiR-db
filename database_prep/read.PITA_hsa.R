# read.PITA_hsa.R
# To read and preprocess human data from PITA.
# 
# Status: 9/2017 - Has NOT been updated to use current annotations available through Bioconductor
#
######
#
# 1. Read data from downloaded file
#
######

org = 'hsa'
pita_hsa.file = "./PITA/v6/PITA_sites_hg18_0_0_ALL.tab"
pita_hsa.data = read.delim(pita_hsa.file)

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
pita_hsa.data = pita_hsa.data[,c(1:3,5,8,9)]

# 2.2 update mature miRNA IDs
pita_hsa.matureID = unique(as.character(pita_hsa.data$microRNA))
pita_hsa.newID = update.matureID(pita_hsa.matureID)
# in case some new IDs == NA
pita_hsa.combineID = pita_hsa.newID
pita_hsa.combineID[is.na(pita_hsa.combineID)] = pita_hsa.matureID[is.na(pita_hsa.combineID)]
m = match(pita_hsa.data$microRNA, pita_hsa.matureID)
pita_hsa.data$microRNA = pita_hsa.combineID[m]

# 2.3 get mature miRNA accession numbers
pita_hsa.matureACC = matureID2matureACC(pita_hsa.combineID)
pita_hsa.data = cbind(mature_mirna_acc=pita_hsa.matureACC[m], pita_hsa.data)

# 2.4 NOTE: in column RefSeq, some have multiple accession numbers separated by
#     ';', "NM_000018;NM_001033859" for example.
#     ONLY keep the first accession number.
pita_hsa.data$RefSeq = gsub(";\\w*", "", as.character(pita_hsa.data$RefSeq), perl=TRUE)

# 2.5 get gene symbols, Entrez & Ensembl gene IDs using RefSeq IDs
refseq = unique(as.character(pita_hsa.data$RefSeq))
pita_hsa.IDs = convert.gene.IDs(org=org, ID=refseq, ID.type='refseq')
pita_hsa.IDs = unique(pita_hsa.IDs)

# 2.6 combine pita_hsa.IDs with pita_hsa.data
pita_hsa.data.DT = data.table(pita_hsa.data)
setkey(pita_hsa.data.DT, RefSeq)
pita_hsa.IDs.DT = data.table(pita_hsa.IDs)
setkey(pita_hsa.IDs.DT, refseq_acc)
pita_hsa.data2.DT = pita_hsa.IDs.DT[pita_hsa.data.DT, allow.cartesian=TRUE]
pita_hsa.data2.DT = pita_hsa.data2.DT[,list(mature_mirna_acc,microRNA,Name,gene_symbol,entrez_gene_id,ensembl_gene_id,refseq_acc,Seed,ddG,Conservation)]
pita_hsa.data2.DT = unique(pita_hsa.data2.DT)
pita_hsa.data2 = data.frame(pita_hsa.data2.DT)

# 2.7 Use column 'gene_symbol', instead of 'Name', but replace NA's in
#     'gene_symbol' with corresponding records in 'Name'
sum(is.na(pita_hsa.data2$gene_symbol))
#[1] 180340
gene_symbol = as.character(pita_hsa.data2$gene_symbol)
gene_symbol[is.na(gene_symbol)] = as.character(pita_hsa.data2$Name[is.na(gene_symbol)])
pita_hsa.data2$gene_symbol = gene_symbol
pita_hsa = unique(pita_hsa.data2[,-3])
dim(pita_hsa)
#[1] 7711997       9

# 2.8 replace NA by '' in miRNA and gene ID columns
for (i in 1:6) {
    tmp = as.character(pita_hsa[,i])
    tmp[is.na(tmp)] = ''
    pita_hsa[,i] = tmp
}

# 2.9 remove columns 'refseq_acc' & 'Seed'
pita_hsa = unique(pita_hsa[,-c(6,7)])
dim(pita_hsa)
#[1] 7710936       7

# 2.10 remove records with '' in gene_symbol, entrez_gene_id & ensembl_gene_id
remove = which(pita_hsa$gene_symbol == '' & pita_hsa$entrez_gene_id == '' & pita_hsa$ensembl_gene_id == '')
if (length(remove) > 0) {pita_hsa = pita_hsa[-remove,]}
dim(pita_hsa)
#[1] 7710936       7

# 2.11 remove records with conservation < PITA_CONSERVE
if (!is.null(PITA_CONSERVE)) {
    pita_hsa = pita_hsa[pita_hsa$Conservation >= PITA_CONSERVE,]
}
dim(pita_hsa)
#[1] 1874890       7

######
#
# 3. Populate MySQL tables mirna, target & pita
#
######

# add org to the table
pita_hsa = cbind(org=org, pita_hsa)

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', pita_hsa[,1:3])

# update target table
target.uid = update.target.table(con, 'target', pita_hsa[,c(1,4:6)])

# then insert pita_hsa into table pita
pita_hsa = data.frame(mirna.uid, target.uid, pita_hsa[,7:8])
fields = dbListFields(con, 'pita')
colnames(pita_hsa) = fields
dim(pita_hsa)
#[1] 1874890       4
dim(unique(pita_hsa))
#[1] 1874890       4
dbWriteTable(con, 'pita', pita_hsa, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 4. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.PITA_hsa.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(m, refseq, gene_symbol, i, tmp, remove, drv, con)
save.image("read.PITA_hsa.RData")

