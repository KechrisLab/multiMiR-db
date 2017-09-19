# read.TargetScan_mmu_nonconserved_site.R
# To read and preprocess mouse data from TargetScan.
# 
# Status: 9/2017 - Has been updated to use current annotations available through Bioconductor
#
######
#
# 1. Read data from downloaded file
#
######

org = 'mmu'
targetscan_mmu.file = "./TargetScanv7.1/mmu/Nonconserved_Site_Context_Scores.txt"
targetscan_mmu.data = read.delim(targetscan_mmu.file)	# 17075529 rows

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

dim(targetscan_mmu.data)
#[1] 35381113       12

# 2.1 only mouse miRNAs?
summary(as.factor(substr(targetscan_mmu.data$miRNA, 1, 4)))
#     bta-     cfa-     gga-     hsa-     mdo-     mml-     mmu-     ptr-     rno-     xtr- 
#2354375  1315272   729242 10147597  1349598  3106867 10775035  1967728  3573934    61465 

if (SAME_ORG) {
    targetscan_mmu.data = targetscan_mmu.data[grep(org, targetscan_mmu.data$miRNA),]
    dim(targetscan_mmu.data)
#[1] 10775035       12
}

# 2.2 update mature miRNA IDs
targetscan_mmu.matureID = unique(as.character(targetscan_mmu.data$miRNA))
targetscan_mmu.newID = update.matureID(targetscan_mmu.matureID)
# in case some new IDs == NA
targetscan_mmu.combineID = targetscan_mmu.newID
targetscan_mmu.combineID[is.na(targetscan_mmu.combineID)] = targetscan_mmu.matureID[is.na(targetscan_mmu.combineID)]
m = match(targetscan_mmu.data$miRNA, targetscan_mmu.matureID)
targetscan_mmu.data$miRNA = targetscan_mmu.combineID[m]

# 2.3 get mature miRNA accession numbers
targetscan_mmu.matureACC = matureID2matureACC(targetscan_mmu.combineID)
targetscan_mmu.data = cbind(mature_mirna_acc=targetscan_mmu.matureACC[m], targetscan_mmu.data)

# 2.4 get Ensembl gene IDs
# NOTE: Gene symbols and Ensembl gene IDs are retrieved using Entrez gene IDs.
#	It is not checked whether the retrieved gene symbols are the same as
#	the original gene symbols in targetscan_mmu.data.
#	The retrieved gene symbols are used eventually.
targetscan_mmu.data$Gene.ID=gsub("\\..*","",targetscan_mmu.data$Gene.ID)
targetscan_mmu.data$Transcript.ID=gsub("\\..*","",targetscan_mmu.data$Transcript.ID)
targetscan_mmu.IDs = NULL
ensembl = unique(targetscan_mmu.data$Gene.ID)
targetscan_mmu.IDs = convert.gene.IDs(org=org, ID=ensembl, ID.type='ensembl')
targetscan_mmu.IDs = unique(targetscan_mmu.IDs)

# 2.5 combine targetscan_mmu.IDs with targetscan_mmu.data
targetscan_mmu.data.DT = data.table(targetscan_mmu.data)
setkey(targetscan_mmu.data.DT, Gene.ID)
targetscan_mmu.IDs.DT = data.table(targetscan_mmu.IDs)
setkey(targetscan_mmu.IDs.DT, ensembl_gene_id)
targetscan_mmu.data2.DT = targetscan_mmu.IDs.DT[targetscan_mmu.data.DT, allow.cartesian=TRUE]

targetscan_mmu.data2.DT = targetscan_mmu.data2.DT[,list(mature_mirna_acc,miRNA,gene_symbol,Gene.Symbol,entrez_gene_id,ensembl_gene_id,Transcript.ID,Site.Type,context...score,context...score.percentile)]
targetscan_mmu.data2.DT = unique(targetscan_mmu.data2.DT)

targetscan_mmu.data2 = data.frame(targetscan_mmu.data2.DT)
dim(targetscan_mmu.data2)
#[1] 10348189       10

# Use column 'gene_symbol' (the new one), instead of 'Gene.Symbol', but replace
# NA's in 'gene_symbol' with corresponding records in 'Gene.Symbol'
gene_symbol = as.character(targetscan_mmu.data2$gene_symbol)
gene_symbol[is.na(gene_symbol)] = as.character(targetscan_mmu.data2$Gene.Symbol[is.na(gene_symbol)])
targetscan_mmu.data2$gene_symbol = gene_symbol
targetscan_mmu = unique(targetscan_mmu.data2[,-4])
dim(targetscan_mmu)
#[1] 10348189        9

# 2.6 replace NA by '' in miRNA and gene ID columns
for (i in 1:6) {
    tmp = as.character(targetscan_mmu[,i])
    tmp[is.na(tmp)] = ''
    if(i>=5 ){
      tmp[!startsWith(tmp, "ENS")] = ''
    }
    targetscan_mmu[,i] = tmp
}

# 2.7 Remove some columns (Transcript.ID & context..score.percentile)
targetscan_mmu = unique(targetscan_mmu[,c(1:5,7,8)])
dim(targetscan_mmu)
#[1] 10348189        7

# 2.8 Remove pairs with context+ score 'NULL'
m = targetscan_mmu$context...score != 'NULL'
targetscan_mmu = targetscan_mmu[m,]
dim(targetscan_mmu)
#[1] 10254744        7

######
#
# 3. Populate MySQL tables mirna, target & targetscan
#
######

# add org to the table
targetscan_mmu = cbind(org=org, targetscan_mmu, conserved_site='N')

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', targetscan_mmu[,1:3])

# update target table
target.uid = update.target.table(con, 'target', targetscan_mmu[,c(1,4:6)])

# then insert targetscan_mmu into table targetscan
targetscan_mmu = cbind(mirna.uid, target.uid, targetscan_mmu[,7:9])
fields = dbListFields(con, 'targetscan')
colnames(targetscan_mmu) = fields
dim(targetscan_mmu)
#[1] 10254744        5
dim(unique(targetscan_mmu))
#[1] 10254744        5
dbWriteTable(con, 'targetscan', targetscan_mmu, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 4. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.TargetScan_mmu_nonconserved_site.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(m, entrez, gene_symbol, i, tmp, drv, con)
save.image("read.TargetScan_mmu_nonconserved_site.RData")


