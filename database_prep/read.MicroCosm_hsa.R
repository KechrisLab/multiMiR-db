# read.MicroCosm_hsa.R
# To read and preprocess human data from MicroCosm.
# 
# Status: 9/2017 - Has NOT been updated to use current annotations available through Bioconductor
#
######
#
# 1. Read data from downloaded file
#
######

org = 'hsa'
microcosm_hsa.file = "./MicroCosm/v5/v5.txt.homo_sapiens"
microcosm_hsa.data = read.delim(microcosm_hsa.file, skip=4)

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
microcosm_hsa.data = microcosm_hsa.data[,c(2,10:13)]

# some duplicates after removing these columns
dim(microcosm_hsa.data)
#[1] 879049      5
dim(unique(microcosm_hsa.data))
#[1] 872980      5

# 2.2 only keep microRNAs from human
summary(as.factor(substr(microcosm_hsa.data$SEQ, 1, 3)))
#   gga    hsa    mml    mmu    rno 
#  2948 728288   1088 140867   5858 
microcosm_hsa.data = microcosm_hsa.data[grep(org, microcosm_hsa.data$SEQ),]

# 2.3 update mature miRNA IDs
microcosm_hsa.matureID = unique(as.character(microcosm_hsa.data$SEQ))
microcosm_hsa.newID = update.matureID(microcosm_hsa.matureID)
# in case some new IDs == NA
microcosm_hsa.combineID = microcosm_hsa.newID
microcosm_hsa.combineID[is.na(microcosm_hsa.combineID)] = microcosm_hsa.matureID[is.na(microcosm_hsa.combineID)]
m = match(microcosm_hsa.data$SEQ, microcosm_hsa.matureID)
microcosm_hsa.data$SEQ = microcosm_hsa.combineID[m]

# 2.4 get mature miRNA accession numbers
microcosm_hsa.matureACC = matureID2matureACC(microcosm_hsa.combineID)
microcosm_hsa.data = cbind(mature_mirna_acc=microcosm_hsa.matureACC[m], microcosm_hsa.data)

# 2.5 get gene symbols, Entrez & Ensembl gene IDs using Ensembl transcript IDs
# (column TRANSCRIPT_ID has no empty IDs but column EXTERNAL_NAME (gene symbol)
# has some empty IDs)
sum(is.na(microcosm_hsa.data$TRANSCRIPT_ID))
#[1] 0
sum(is.na(microcosm_hsa.data$EXTERNAL_NAME))
#[1] 0
sum(microcosm_hsa.data$TRANSCRIPT_ID == '')
#[1] 0
sum(microcosm_hsa.data$EXTERNAL_NAME == '')
#[1] 14994
ensembl = unique(as.character(microcosm_hsa.data$TRANSCRIPT_ID))
microcosm_hsa.IDs = convert.gene.IDs(org=org, ID=ensembl, ID.type='ensembl_transcript')
microcosm_hsa.IDs = unique(microcosm_hsa.IDs)

# 2.6 combine microcosm_hsa.IDs with microcosm_hsa.data
microcosm_hsa.data.DT = data.table(microcosm_hsa.data)
setkey(microcosm_hsa.data.DT, TRANSCRIPT_ID)
microcosm_hsa.IDs.DT = data.table(microcosm_hsa.IDs)
setkey(microcosm_hsa.IDs.DT, ensembl_transcript_id)
microcosm_hsa.data2.DT = microcosm_hsa.IDs.DT[microcosm_hsa.data.DT, allow.cartesian=TRUE]

# 2.7 For Ensembl transcript IDs whose gene symbols could not be found by
#     biomaRt, use the original gene symbols (column EXTERNAL_NAME) to find
#     Entrez & Ensembl gene IDs
microcosm_hsa.data2.DT1 = microcosm_hsa.data2.DT[!is.na(microcosm_hsa.data2.DT$gene_symbol),]
microcosm_hsa.data2.DT2 = microcosm_hsa.data2.DT[is.na(microcosm_hsa.data2.DT$gene_symbol),]
symbol = unique(as.character(microcosm_hsa.data2.DT2$EXTERNAL_NAME))
symbol = setdiff(symbol, '')
microcosm_hsa.IDs2 = convert.gene.IDs(org=org, ID=symbol, ID.type='symbol')
microcosm_hsa.IDs2 = unique(microcosm_hsa.IDs2)

# 2.8 combine microcosm_hsa.IDs2 with microcosm_hsa.data2.DT2
setkey(microcosm_hsa.data2.DT2, EXTERNAL_NAME)
microcosm_hsa.IDs2.DT = data.table(microcosm_hsa.IDs2)
setkey(microcosm_hsa.IDs2.DT, gene_symbol)
microcosm_hsa.data2.DT3 = microcosm_hsa.IDs2.DT[microcosm_hsa.data2.DT2, allow.cartesian=TRUE]
microcosm_hsa.data2.DT3 = microcosm_hsa.data2.DT3[,list(mature_mirna_acc,SEQ,gene_symbol,entrez_gene_id,ensembl_gene_id,ensembl_transcript_id,SCORE,PVALUE_OG)]

# 2.9 combine microcosm_hsa.data2.DT3 with microcosm_hsa.data2.DT1
microcosm_hsa.data2.DT1 = microcosm_hsa.data2.DT1[,list(mature_mirna_acc,SEQ,gene_symbol,entrez_gene_id,ensembl_gene_id,ensembl_transcript_id,SCORE,PVALUE_OG)]
microcosm_hsa.data3.DT = rbind(microcosm_hsa.data2.DT1, microcosm_hsa.data2.DT3)
microcosm_hsa.data3.DT = unique(microcosm_hsa.data3.DT)

microcosm_hsa = data.frame(microcosm_hsa.data3.DT)

# 2.10 replace NA by '' in miRNA and gene ID columns
for (i in 1:6) {
    tmp = as.character(microcosm_hsa[,i])
    tmp[is.na(tmp)] = ''
    microcosm_hsa[,i] = tmp
}

######
#
# 3. Remove column ensembl_transcript_id in microcosm_hsa --
#    there is no ensembl_transcript_id in the MySQL target table.
#    However, removing ensembl_transcript_id will remove the interactions
#    involving genes with only ensembl_transcript_id available.
#
######

# there are 10579 (1.3%) interactions involving genes with ensembl_transcript_id only.
length(which(microcosm_hsa$gene_symbol=='' & microcosm_hsa$entrez_gene_id=='' & microcosm_hsa$ensembl_gene_id==''))
#[1] 10579
length(which(microcosm_hsa$gene_symbol=='' & microcosm_hsa$entrez_gene_id=='' & microcosm_hsa$ensembl_gene_id=='')) / nrow(microcosm_hsa)
#[1] 0.01335361

dim(microcosm_hsa)
#[1] 792220      8
microcosm_hsa = microcosm_hsa[,-6]
microcosm_hsa = unique(microcosm_hsa)
dim(microcosm_hsa)
#[1] 772924      7

remove = which(microcosm_hsa$gene_symbol=='' & microcosm_hsa$entrez_gene_id=='' & microcosm_hsa$ensembl_gene_id=='')
if (length(remove) > 0) {microcosm_hsa = microcosm_hsa[-remove,]}
dim(microcosm_hsa)
#[1] 762987      7

######
#
# 4. Populate MySQL tables mirna, target & microcosm
#
######

# add org to the table
microcosm_hsa = cbind(org=org, microcosm_hsa)

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', microcosm_hsa[,1:3])

# update target table
target.uid = update.target.table(con, 'target', microcosm_hsa[,c(1,4:6)])

# then insert microcosm_hsa into table microcosm
# NOTE: some duplicates after removing column #8 (PVALUE_OG) -- same score &
#       different p-values for multiple sites in the same gene?
# Alternative: keep column #8 (PVALUE_OG)?
microcosm_hsa = data.frame(mirna.uid, target.uid, microcosm_hsa[,7])
fields = dbListFields(con, 'microcosm')
colnames(microcosm_hsa) = fields
dim(microcosm_hsa)
#[1] 762987      3
dim(unique(microcosm_hsa))
#[1] 755987      3
dbWriteTable(con, 'microcosm', microcosm_hsa, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 5. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.MicroCosm_hsa.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(m, ensembl, symbol, i, tmp, remove, drv, con)
save.image("read.MicroCosm_hsa.RData")

