# read.DIANA-microT.R
# To read and preprocess data from DIANA-microT.
# The data include 4 species (cel, dme, hsa & mmu).
# 
# Status: 9/2017 - Has NOT been updated to use current annotations available through Bioconductor
#
######
#
# 1. Read data from downloaded file
#
######

diana_microt.file = "./DIANA-microT/microT_CDS_data.out"
diana_microt.data = read.csv(diana_microt.file, sep=",")

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

dim(diana_microt.data)
#[1] 11598429        7

# 2.1 Keep miRNAs from human and mouse only -- no mouse
summary(as.factor(substr(diana_microt.data$miRNA, 1, 4)))
#   cel-    dme-    hsa-    mmu- 
# 254289  445878 7337705 3560557 
keep = substr(diana_microt.data$miRNA, 1, 3) == 'hsa' | substr(diana_microt.data$miRNA, 1, 3) == 'mmu'
diana_microt.data = diana_microt.data[keep,]
dim(diana_microt.data)
#[1] 10898262        7

# 2.2 update mature miRNA IDs
diana_microt.matureID = unique(as.character(diana_microt.data$miRNA))
diana_microt.newID = update.matureID(diana_microt.matureID)
# in case some new IDs == NA
diana_microt.combineID = diana_microt.newID
diana_microt.combineID[is.na(diana_microt.combineID)] = diana_microt.matureID[is.na(diana_microt.combineID)]
m = match(diana_microt.data$miRNA, diana_microt.matureID)
diana_microt.data$miRNA = diana_microt.combineID[m]

# 2.3 get mature miRNA accession numbers
diana_microt.matureACC = matureID2matureACC(diana_microt.combineID)
diana_microt.data = cbind(mature_mirna_acc=diana_microt.matureACC[m], diana_microt.data)

# 2.4 get gene symbol & Entrez gene IDs
diana_microt.org = unique(substr(diana_microt.data$miRNA, 1, 3))
diana_microt.IDs = NULL
for (org in diana_microt.org) {
    m = grep(org, diana_microt.data$miRNA)
    ensembl = unique(diana_microt.data$GeneID[m])
    ID = convert.gene.IDs(org=org, ID=ensembl, ID.type='ensembl')
    diana_microt.IDs = rbind(diana_microt.IDs, ID)
}
diana_microt.IDs = unique(diana_microt.IDs)

# 2.5 combine diana_microt.IDs with diana_microt.data
diana_microt.data.DT = data.table(diana_microt.data)
setkey(diana_microt.data.DT, GeneID)
diana_microt.IDs.DT = data.table(diana_microt.IDs)
setkey(diana_microt.IDs.DT, ensembl_gene_id)
diana_microt.data2.DT = diana_microt.IDs.DT[diana_microt.data.DT, allow.cartesian=TRUE]
# remove column 'TranscriptID'
diana_microt.data2.DT = diana_microt.data2.DT[,list(mature_mirna_acc,miRNA,gene_symbol,GeneName,entrez_gene_id,ensembl_gene_id,miTG_score,UTR3_hit,CDS_hit)]
diana_microt.data2.DT = unique(diana_microt.data2.DT)
diana_microt.data2 = data.frame(diana_microt.data2.DT)
dim(diana_microt.data2)
#[1] 11411773        9

# NOTE: in diana_microt.data2.DT, column 'gene_symbol' is from diana_microt.IDs
# and column 'GeneName' is from diana_microt.data.
# Use column 'gene_symbol', instead of 'GeneName', but replace NA's in
# 'gene_symbol' with corresponding records in 'GeneName'
gene_symbol = as.character(diana_microt.data2$gene_symbol)
gene_symbol[is.na(gene_symbol)] = as.character(diana_microt.data2$GeneName[is.na(gene_symbol)])
diana_microt.data2$gene_symbol = gene_symbol
diana_microt = unique(diana_microt.data2[,-4])
dim(diana_microt)
#[1] 11411773        8

# 2.6 replace NA by '' in miRNA and gene ID columns
for (i in 1:5) {
    tmp = as.character(diana_microt[,i])
    tmp[is.na(tmp)] = ''
    diana_microt[,i] = tmp
}

sum(diana_microt$gene_symbol == '' & diana_microt$entrez_gene_id == '' & diana_microt$ensembl_gene_id == '')
#[1] 0
sum(diana_microt$mature_mirna_acc == '' & diana_microt$miRNA == '')
#[1] 0

######
#
# 3. Populate MySQL tables mirna, target & diana_microt
#
######

# add org to the table
org = substr(diana_microt$miRNA, 1, 3)
summary(as.factor(org))
#    hsa     mmu 
#7664602 3747171 
diana_microt = cbind(org=org, diana_microt)

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', diana_microt[,1:3])

# update target table
target.uid = update.target.table(con, 'target', diana_microt[,c(1,4:6)])

# then insert diana_microt into table diana_microt
diana_microt = cbind(mirna.uid, target.uid, diana_microt[,7:9])
fields = dbListFields(con, 'diana_microt')
colnames(diana_microt) = fields
dim(diana_microt)
#[1] 11411773        5
dim(unique(diana_microt))
#[1] 11411773        5
dbWriteTable(con, 'diana_microt', diana_microt, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 4. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.DIANA-microT.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(keep, m, org, ensembl, ID, gene_symbol, i, tmp, drv, con)
save.image("read.DIANA-microT.RData")

