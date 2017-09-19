# read.miRTarBase.R
# To read and preprocess data from miRTarBase.
# 
# Status: 9/2017 - Has been updated to use current annotations available through Bioconductor
#

#
# 1. Read data from downloaded file
#

mirtarbase.file = "miRTarBase_MTI.csv"
mirtarbase.data = read.delim(mirtarbase.file, sep=",", quote="")

# MySQL parameters
source("build.multiMiR_DB.parameters.R")

# functions
source("multiMiR.R")

#
# 2. Preprocess data
#

# 2.1 Keep miRNA-target pairs in human only, or in mouse only
# NOTE: For other species (e.g. Arabidopsis), need to update file
#       multiMiR.species.csv to specify the correct bioMart database to use.
keep = (mirtarbase.data[,3] == 'Homo sapiens' & mirtarbase.data[,6] == 'Homo sapiens') | (mirtarbase.data[,3] == 'Mus musculus' & mirtarbase.data[,6] == 'Mus musculus') | (mirtarbase.data[,3] == 'Rattus norvegicus' & mirtarbase.data[,6] == 'Rattus norvegicus')
dim(mirtarbase.data)
#[1] 34025     9
mirtarbase.data = mirtarbase.data[keep,]
dim(mirtarbase.data)
#[1] 30154     9

# 2.2 remove records without miRNA IDs
remove = mirtarbase.data$miRNA == '' | is.na(mirtarbase.data$miRNA)
mirtarbase.data = mirtarbase.data[!remove,]
dim(mirtarbase.data)
#[1] 30109     9



# 2.3 update mature miRNA IDs
mirtarbase.matureID = unique(as.character(mirtarbase.data$miRNA))
mirtarbase.newID = update.matureID(mirtarbase.matureID)
# some new IDs == NA
mirtarbase.combineID = mirtarbase.newID
mirtarbase.combineID[is.na(mirtarbase.combineID)] = mirtarbase.matureID[is.na(mirtarbase.combineID)]

# 2.4 add mature miRNA accession numbers
mirtarbase.matureACC = matureID2matureACC(mirtarbase.combineID)
# some ACCs == NA
m = match(mirtarbase.data$miRNA, mirtarbase.matureID)
mirtarbase.data2 = cbind(mirtarbase.data[,1], mirtarbase.matureACC[m], mirtarbase.combineID[m], mirtarbase.data[3:ncol(mirtarbase.data)])
mirtarbase.data2 = unique(mirtarbase.data2)
dim(mirtarbase.data2)
#[1] 30109    10

# 2.5 remove records without target gene information
#remove =  (mirtarbase.data2[,5] == '' | is.na(mirtarbase.data2[,5])) & (mirtarbase.data2[,6] == '' | is.na(mirtarbase.data2[,6]))
remove = mirtarbase.data2[,6] == '' | is.na(mirtarbase.data2[,6])
mirtarbase.data2 = mirtarbase.data2[!remove,]
dim(mirtarbase.data2)
#[1] 30100    10

# 2.6 add Ensembl gene IDs
# NOTE: Gene symbols and Ensembl gene IDs are retrieved using Entrez gene IDs.
#	It is not checked whether the retrieved gene symbols are the same as
#	the original gene symbols in mirtarbase.data.
#	The retrieved gene symbols are used eventually.
mirtarbase.org = unique(as.character(mirtarbase.data$Species..Target.Gene.))
mirtarbase.IDs = NULL
for (org in mirtarbase.org) {
    m = grep(org, mirtarbase.data$Species..Target.Gene.)
    entrez = unique(mirtarbase.data$Target.Gene..Entrez.Gene.ID.[m])
    ID = convert.gene.IDs(org=org, ID=entrez, ID.type='entrez')
    mirtarbase.IDs = rbind(mirtarbase.IDs, ID)
}
mirtarbase.IDs = unique(mirtarbase.IDs)
dim(mirtarbase.IDs)
#[1] 15518     3

# 2.7 combine mirtarbase.IDs with mirtarbase.data2 using Entrez gene IDs as keys
colnames(mirtarbase.data2)[5:6] = colnames(mirtarbase.IDs)[1:2]
mirtarbase.IDs.DT = data.table(mirtarbase.IDs)
setkeyv(mirtarbase.IDs.DT, colnames(mirtarbase.IDs)[2])
mirtarbase.data2.DT = data.table(mirtarbase.data2)
setkeyv(mirtarbase.data2.DT, colnames(mirtarbase.IDs)[2])
#mirtarbase.data3 = mirtarbase.IDs.DT[mirtarbase.data2.DT]	# error
mirtarbase.data3 = mirtarbase.IDs.DT[mirtarbase.data2.DT, allow.cartesian=TRUE]
dim(mirtarbase.data3)
#[1] 33205    12

# 2.8 What if new gene symbol != miRTarBase gene symbol?
# (1) new gene symbol is NA: replaced with miRTarBase gene symbol
sum(is.na(mirtarbase.data3$gene_symbol))	# 160
sum(is.na(mirtarbase.data3$i.gene_symbol))	# 0
mirtarbase.data3$gene_symbol[is.na(mirtarbase.data3$gene_symbol)] = as.character(mirtarbase.data3$i.gene_symbol[is.na(mirtarbase.data3$gene_symbol)])
# (2) new gene symbol == '': replaced with miRTarBase gene symbol
sum(mirtarbase.data3$gene_symbol == '', na.rm=T)	# 238
sum(mirtarbase.data3$i.gene_symbol == '', na.rm=T)	# 0
mirtarbase.data3$gene_symbol[mirtarbase.data3$gene_symbol == ''] = as.character(mirtarbase.data3$i.gene_symbol[mirtarbase.data3$gene_symbol == ''])
# (3) new gene symbol != miRTarBase gene symbol: use the new gene symbol
mirtarbase = unique(data.frame(mirtarbase.data3))
sum(mirtarbase$gene_symbol != mirtarbase$i.gene_symbol)	# 847
mirtarbase = mirtarbase[,-8]
dim(mirtarbase)
#[1] 33205    11

# 2.9 species columns
mirtarbase[,7] = as.character(mirtarbase[,7])
mirtarbase[,8] = as.character(mirtarbase[,8])
mirtarbase[mirtarbase[,7] == 'Homo sapiens', 7] = 'hsa'
mirtarbase[mirtarbase[,8] == 'Homo sapiens', 8] = 'hsa'
mirtarbase[mirtarbase[,7] == 'Mus musculus', 7] = 'mmu'
mirtarbase[mirtarbase[,8] == 'Mus musculus', 8] = 'mmu'
mirtarbase[mirtarbase[,7] == 'Rattus norvegicus', 7] = 'rno'
mirtarbase[mirtarbase[,8] == 'Rattus norvegicus', 8] = 'rno'

# 2.10 replace NA by '' in miRNA and gene ID columns
for (i in c(1:3,5:8)) {
    tmp = as.character(mirtarbase[,i])
    tmp[is.na(tmp)] = ''
    mirtarbase[,i] = tmp
}

######
#
# 3. Populate MySQL tables mirna, target & mirtarbase
#
######

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', mirtarbase[,c(7,5,6)])
# update target table
# NOTE: no RefSeq numbers
target.uid = update.target.table(con, 'target', mirtarbase[,c(8,1,2,3)])
# then insert mirtarbase into table mirtarbase
mirtarbase = cbind(mirna.uid, target.uid, mirtarbase[,9:11])
fields = dbListFields(con, 'mirtarbase')
colnames(mirtarbase) = fields
dbWriteTable(con, 'mirtarbase', mirtarbase, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 4. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.miRTarBase.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(keep, remove, m, org, ID, entrez, i, tmp, con, drv)
save.image("new.read.miRTarBase.RData")

