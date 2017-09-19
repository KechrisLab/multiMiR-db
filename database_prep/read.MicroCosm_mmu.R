# read.MicroCosm_mmu.R
# To read and preprocess mouse data from MicroCosm.
# 
# Status: 9/2017 - Has NOT been updated to use current annotations available through Bioconductor
#
######
#
# 1. Read data from downloaded file
#
######

org = 'mmu'
microcosm_mmu.file = "./MicroCosm/v5/v5.txt.mus_musculus"
microcosm_mmu.data = read.delim(microcosm_mmu.file, skip=4)

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
# (of which X..GROUP, METHOD, FEATURE & PHASE are constant)
microcosm_mmu.data = microcosm_mmu.data[,c(2,10:13)]

# some duplicates after removing these columns
dim(microcosm_mmu.data)
#[1] 750375      5
dim(unique(microcosm_mmu.data))
#[1] 745596      5

# 2.2 only keep microRNAs from mouse
summary(as.factor(substr(microcosm_mmu.data$SEQ, 1, 3)))
#   bta    gga    hsa    mml    mmu    rno 
#  1814   2741 204663    941 534885   5331 
microcosm_mmu.data = microcosm_mmu.data[grep(org, microcosm_mmu.data$SEQ),]

# 2.3 update mature miRNA IDs
microcosm_mmu.matureID = unique(as.character(microcosm_mmu.data$SEQ))
microcosm_mmu.newID = update.matureID(microcosm_mmu.matureID)
# in case some new IDs == NA
microcosm_mmu.combineID = microcosm_mmu.newID
microcosm_mmu.combineID[is.na(microcosm_mmu.combineID)] = microcosm_mmu.matureID[is.na(microcosm_mmu.combineID)]
m = match(microcosm_mmu.data$SEQ, microcosm_mmu.matureID)
microcosm_mmu.data$SEQ = microcosm_mmu.combineID[m]

# 2.4 get mature miRNA accession numbers
microcosm_mmu.matureACC = matureID2matureACC(microcosm_mmu.combineID)
microcosm_mmu.data = cbind(mature_mirna_acc=microcosm_mmu.matureACC[m], microcosm_mmu.data)

# 2.5 get gene symbols, Entrez & Ensembl gene IDs using Ensembl transcript IDs
# (column TRANSCRIPT_ID has no empty IDs but column EXTERNAL_NAME (gene symbol)
# has some empty IDs)
sum(is.na(microcosm_mmu.data$TRANSCRIPT_ID))
#[1] 0
sum(is.na(microcosm_mmu.data$EXTERNAL_NAME))
#[1] 0
sum(microcosm_mmu.data$TRANSCRIPT_ID == '')
#[1] 0
sum(microcosm_mmu.data$EXTERNAL_NAME == '')
#[1] 26559
ensembl = unique(as.character(microcosm_mmu.data$TRANSCRIPT_ID))
microcosm_mmu.IDs = convert.gene.IDs(org=org, ID=ensembl, ID.type='ensembl_transcript')
microcosm_mmu.IDs = unique(microcosm_mmu.IDs)

# 2.6 combine microcosm_mmu.IDs with microcosm_mmu.data
microcosm_mmu.data.DT = data.table(microcosm_mmu.data)
setkey(microcosm_mmu.data.DT, TRANSCRIPT_ID)
microcosm_mmu.IDs.DT = data.table(microcosm_mmu.IDs)
setkey(microcosm_mmu.IDs.DT, ensembl_transcript_id)
microcosm_mmu.data2.DT = microcosm_mmu.IDs.DT[microcosm_mmu.data.DT, allow.cartesian=TRUE]

# 2.7 For Ensembl transcript IDs whose gene symbols could not be found by
#     biomaRt, use the original gene symbols (column EXTERNAL_NAME) to find
#     Entrez & Ensembl gene IDs
microcosm_mmu.data2.DT1 = microcosm_mmu.data2.DT[!is.na(microcosm_mmu.data2.DT$gene_symbol),]
microcosm_mmu.data2.DT2 = microcosm_mmu.data2.DT[is.na(microcosm_mmu.data2.DT$gene_symbol),]
symbol = unique(as.character(microcosm_mmu.data2.DT2$EXTERNAL_NAME))
symbol = setdiff(symbol, '')
microcosm_mmu.IDs2 = convert.gene.IDs(org=org, ID=symbol, ID.type='symbol')
microcosm_mmu.IDs2 = unique(microcosm_mmu.IDs2)

# 2.8 combine microcosm_mmu.IDs2 with microcosm_mmu.data2.DT2
setkey(microcosm_mmu.data2.DT2, EXTERNAL_NAME)
microcosm_mmu.IDs2.DT = data.table(microcosm_mmu.IDs2)
setkey(microcosm_mmu.IDs2.DT, gene_symbol)
microcosm_mmu.data2.DT3 = microcosm_mmu.IDs2.DT[microcosm_mmu.data2.DT2, allow.cartesian=TRUE]
microcosm_mmu.data2.DT3 = microcosm_mmu.data2.DT3[,list(mature_mirna_acc,SEQ,gene_symbol,entrez_gene_id,ensembl_gene_id,ensembl_transcript_id,SCORE,PVALUE_OG)]

# 2.9 combine microcosm_mmu.data2.DT3 with microcosm_mmu.data2.DT1
microcosm_mmu.data2.DT1 = microcosm_mmu.data2.DT1[,list(mature_mirna_acc,SEQ,gene_symbol,entrez_gene_id,ensembl_gene_id,ensembl_transcript_id,SCORE,PVALUE_OG)]
microcosm_mmu.data3.DT = rbind(microcosm_mmu.data2.DT1, microcosm_mmu.data2.DT3)
microcosm_mmu.data3.DT = unique(microcosm_mmu.data3.DT)

microcosm_mmu = data.frame(microcosm_mmu.data3.DT)

# 2.10 replace NA by '' in miRNA and gene ID columns
for (i in 1:6) {
    tmp = as.character(microcosm_mmu[,i])
    tmp[is.na(tmp)] = ''
    microcosm_mmu[,i] = tmp
}

######
#
# 3. Remove column ensembl_transcript_id in microcosm_mmu --
#    there is no ensembl_transcript_id in the MySQL target table.
#    However, removing ensembl_transcript_id will remove the interactions
#    involving genes with only ensembl_transcript_id available.
#
######

# there are 17416 (3.1%) interactions involving genes with ensembl_transcript_id only.
length(which(microcosm_mmu$gene_symbol=='' & microcosm_mmu$entrez_gene_id=='' & microcosm_mmu$ensembl_gene_id==''))
#[1] 17416
length(which(microcosm_mmu$gene_symbol=='' & microcosm_mmu$entrez_gene_id=='' & microcosm_mmu$ensembl_gene_id=='')) / nrow(microcosm_mmu)
#[1] 0.03116684

dim(microcosm_mmu)
#[1] 558799      8
microcosm_mmu = microcosm_mmu[,-6]
microcosm_mmu = unique(microcosm_mmu)
dim(microcosm_mmu)
#[1] 550280      7

remove = which(microcosm_mmu$gene_symbol=='' & microcosm_mmu$entrez_gene_id=='' & microcosm_mmu$ensembl_gene_id=='')
if (length(remove) > 0) {microcosm_mmu = microcosm_mmu[-remove,]}
dim(microcosm_mmu)
#[1] 534735      7

######
#
# 4. Populate MySQL tables mirna, target & microcosm
#
######

# add org to the table
microcosm_mmu = cbind(org=org, microcosm_mmu)

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', microcosm_mmu[,1:3])

# update target table
target.uid = update.target.table(con, 'target', microcosm_mmu[,c(1,4:6)])

# then insert microcosm_mmu into table microcosm
# NOTE: some duplicates after removing column #8 (PVALUE_OG) -- same score &
#       different p-values for multiple sites in the same gene?
# Alternative: keep column #8 (PVALUE_OG)?
microcosm_mmu = data.frame(mirna.uid, target.uid, microcosm_mmu[,7])
fields = dbListFields(con, 'microcosm')
colnames(microcosm_mmu) = fields
dim(microcosm_mmu)
#[1] 534735      3
dim(unique(microcosm_mmu))
#[1] 532285      3
dbWriteTable(con, 'microcosm', microcosm_mmu, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 5. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.MicroCosm_mmu.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(m, ensembl, symbol, i, tmp, remove, drv, con)
save.image("read.MicroCosm_mmu.RData")

