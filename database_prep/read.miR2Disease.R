# read.miR2Disease.R
# To read and preprocess data from miR2Disease.
# 
# Status: 9/2017 - Has NOT been updated to use current annotations available through Bioconductor
#
######
#
# 1. Read data from downloaded file
#
######

mir2disease.file = "./miR2Disease/2011/AllEntries.txt"
mir2disease.data = read.delim(mir2disease.file, header=F)
colnames(mir2disease.data) = c('miRNA','disease','miRNA_regulation','experiment','year','title')

# MySQL parameters
source("build.multiMiR_DB.parameters.R")

# functions
source("multiMiR.R")

######
#
# 2. Preprocess data
#
######

dim(mir2disease.data)
#[1] 2877    6

# 2.1 clean up column 'miRNA_regulation'
mir2disease.data$miRNA_regulation[mir2disease.data$miRNA_regulation == 'hepatocellular carcinoma (HCC)'] = ''

# 2.2 clean up column 'miRNA'
mir2disease.data$miRNA = gsub(' ', '', mir2disease.data$miRNA)
mir2disease.data$miRNA = gsub('--', '-', mir2disease.data$miRNA)
mir2disease.data$miRNA = gsub('miR-BART', 'ebv-miR-BART', mir2disease.data$miRNA)

# 2.3 update mature miRNA IDs
mir2disease.matureID = unique(as.character(mir2disease.data$miRNA))
mir2disease.newID = update.matureID(mir2disease.matureID)	# some new IDs == NA
mir2disease.combineID = mir2disease.newID
mir2disease.combineID[is.na(mir2disease.combineID)] = mir2disease.matureID[is.na(mir2disease.combineID)]
m = match(mir2disease.data$miRNA, mir2disease.matureID)
mir2disease.data$miRNA = mir2disease.combineID[m]

# 2.4 add mature miRNA accession numbers
mir2disease.matureACC = matureID2matureACC(mir2disease.combineID)
mir2disease.data = cbind(mature_mirna_acc=mir2disease.matureACC[m], mir2disease.data)
mir2disease.data = unique(mir2disease.data)
dim(mir2disease.data)
#[1] 2877    7

# 2.5 add miRNA organism
mir2disease.org = substr(mir2disease.data$miRNA, 1, 4)
mir2disease.org = gsub('-','',mir2disease.org)
summary(as.factor(mir2disease.org))
# ebv  hsa 
#   2 2875 
mir2disease.data = cbind(mir2disease.data[,1:2], mirna_org=mir2disease.org, mir2disease.data[,3:7])

# 2.6 keep miRNAs in human or mouse only
mir2disease = mir2disease.data
mir2disease = mir2disease[mir2disease$mirna_org == 'hsa',]
dim(mir2disease)
#[1] 2875    8

# 2.7 replace NA by '' in miRNA ID columns
for (i in 1:3) {
    tmp = as.character(mir2disease[,i])
    tmp[is.na(tmp)] = ''
    mir2disease[,i] = tmp
}

# NOTE: there are different spelling of the same disease in mir2disease$disease
# use sort(unique(as.character(mir2disease$disease))) to see

######
#
# 3. Populate MySQL tables mirna, target & mir2disease
#
######

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', mir2disease[,c(3,1,2)])

# update target table - no target info

# then insert mir2disease into table mir2disease
mir2disease = cbind(mirna.uid, mir2disease[,4:8])
fields = dbListFields(con, 'mir2disease')
colnames(mir2disease) = fields
dbWriteTable(con, 'mir2disease', mir2disease, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 4. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.miR2Disease.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(m, i, tmp, drv, con)
save.image("read.miR2Disease.RData")

