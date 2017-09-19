# read.TargetScan_hsa_conserved.R
# To read and preprocess human data from TargetScan.
# 
# Status: 9/2017 - Has been updated to use current annotations available through Bioconductor
#
######
#
# 1. Read data from downloaded file
#
######

org = 'hsa'
targetscan_hsa.file = "./TargetScanv7.1/hsa/Conserved_Site_Context_Scores.txt"
targetscan_hsa.data = read.delim(targetscan_hsa.file)	# 1450123 rows

library(data.table)

# MySQL parameters
source("build.multiMiR_DB.parameters.R")


# functions
source("multiMiR.R")

SAME_ORG = TRUE

######
#
# 2. Preprocess data
#
######

dim(targetscan_hsa.data)
#[1] 1450123      12

# 2.1 only human miRNAs?
summary(as.factor(substr(targetscan_hsa.data$miRNA, 1, 4)))
#bta-   cfa-   gga-   hsa-   mdo-   mml-   mmu-   ptr-   rno-   xtr- 
#  172529 155011  46638 257532  77463 200443 193357 187484 138299  21367 
if (SAME_ORG) {
    targetscan_hsa.data = targetscan_hsa.data[grep(org, targetscan_hsa.data$miRNA),]
    dim(targetscan_hsa.data)
#[1] 257532     12
}

# 2.2 update mature miRNA IDs
targetscan_hsa.matureID = unique(as.character(targetscan_hsa.data$miRNA))
targetscan_hsa.newID = update.matureID(targetscan_hsa.matureID)
# in case some new IDs == NA
targetscan_hsa.combineID = targetscan_hsa.newID
targetscan_hsa.combineID[is.na(targetscan_hsa.combineID)] = targetscan_hsa.matureID[is.na(targetscan_hsa.combineID)]
m = match(targetscan_hsa.data$miRNA, targetscan_hsa.matureID)
targetscan_hsa.data$miRNA = targetscan_hsa.combineID[m]

# 2.3 get mature miRNA accession numbers
targetscan_hsa.matureACC = matureID2matureACC(targetscan_hsa.combineID)
targetscan_hsa.data = cbind(mature_mirna_acc=targetscan_hsa.matureACC[m], targetscan_hsa.data)

# 2.4 get Ensembl gene IDs
# NOTE: Gene symbols and Ensembl gene IDs are retrieved using Entrez gene IDs.
#	It is not checked whether the retrieved gene symbols are the same as
#	the original gene symbols in targetscan_hsa.data.
#	The retrieved gene symbols are used eventually.

#Remove .1 suffix of ENSGxxxxx.1
targetscan_hsa.data$Gene.ID=gsub("\\..*","",targetscan_hsa.data$Gene.ID)
targetscan_hsa.data$Transcript.ID=gsub("\\..*","",targetscan_hsa.data$Transcript.ID)

targetscan_hsa.IDs = NULL
ensembl = unique(targetscan_hsa.data$Gene.ID)
targetscan_hsa.IDs = convert.gene.IDs(org=org, ID=ensembl, ID.type='ensembl')
targetscan_hsa.IDs = unique(targetscan_hsa.IDs)

# 2.5 combine targetscan_hsa.IDs with targetscan_hsa.data
targetscan_hsa.data.DT = data.table(targetscan_hsa.data)
setkey(targetscan_hsa.data.DT, Gene.ID)
targetscan_hsa.IDs.DT = data.table(targetscan_hsa.IDs)
setkey(targetscan_hsa.IDs.DT, ensembl_gene_id)
targetscan_hsa.data2.DT = targetscan_hsa.IDs.DT[targetscan_hsa.data.DT, allow.cartesian=TRUE]

targetscan_hsa.data2.DT = targetscan_hsa.data2.DT[,list(mature_mirna_acc,miRNA,gene_symbol,Gene.Symbol,entrez_gene_id,ensembl_gene_id,Transcript.ID,Site.Type,context...score,context...score.percentile)]
targetscan_hsa.data2.DT = unique(targetscan_hsa.data2.DT)

targetscan_hsa.data2 = data.frame(targetscan_hsa.data2.DT)
dim(targetscan_hsa.data2)
#[1] 257527     10

# Use column 'gene_symbol' (the new one), instead of 'Gene.Symbol', but replace
# NA's in 'gene_symbol' with corresponding records in 'Gene.Symbol'
gene_symbol = as.character(targetscan_hsa.data2$gene_symbol)
gene_symbol[is.na(gene_symbol)] = as.character(targetscan_hsa.data2$Gene.Symbol[is.na(gene_symbol)])
targetscan_hsa.data2$gene_symbol = gene_symbol
targetscan_hsa = unique(targetscan_hsa.data2[,-4])
dim(targetscan_hsa)
#[1] 257527      9

# 2.6 replace NA by '' in miRNA and gene ID columns
for (i in 1:6) {
    tmp = as.character(targetscan_hsa[,i])
    tmp[is.na(tmp)] = ''
    targetscan_hsa[,i] = tmp
}

# 2.7 Remove some columns (Transcript.ID & context..score.percentile)
targetscan_hsa = unique(targetscan_hsa[,c(1:5,7,8)])
dim(targetscan_hsa)
#[1] 257527      7

# 2.8 Remove pairs with context+ score 'NULL'
m = targetscan_hsa$context...score != 'NULL'
targetscan_hsa = targetscan_hsa[m,]
dim(targetscan_hsa)
#[1] 257527      7

######
#
# 3. Populate MySQL tables mirna, target & targetscan
#
######

# add org to the table
targetscan_hsa = cbind(org=org, targetscan_hsa, conserved_site='Y')

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', targetscan_hsa[,1:3])

# update target table
target.uid = update.target.table(con, 'target', targetscan_hsa[,c(1,4:6)])

# then insert targetscan_hsa into table targetscan
targetscan_hsa = cbind(mirna.uid, target.uid, targetscan_hsa[,7:9])
fields = dbListFields(con, 'targetscan')
colnames(targetscan_hsa) = fields
dim(targetscan_hsa)
#[1] 621218      5
dim(unique(targetscan_hsa))
#[1] 621218      5
dbWriteTable(con, 'targetscan', targetscan_hsa, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 4. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.TargetScan_hsa_conserved.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(m, entrez, gene_symbol, i, tmp, drv, con)
save.image("read.TargetScan_hsa_conserved.RData")


