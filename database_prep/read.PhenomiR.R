# read.PhenomiR.R
# To read and preprocess data from Pharmaco-miR VerSe (Verified Sets).
# 
# Status: 9/2017 - Has NOT been updated to use current annotations available through Bioconductor
#
######
#
# 1. Read data from downloaded file
#
######

phenomir.file = "./PhenomiR/2011_v2.0/phenomir-2.0.tbl"
phenomir.data = read.delim(phenomir.file)

# MySQL parameters
source("build.multiMiR_DB.parameters.R")

# functions
source("multiMiR.R")

######
#
# 2. Preprocess data
#
######

dim(phenomir.data)
#[1] 11144    13
phenomir.data = unique(phenomir.data)
dim(phenomir.data)
#[1] 11130    13

# 2.1 in columns 'miRNA' & 'accession', most are precursor miRNAs (with
#     accession 'MI00*'), but some are mature miRNAs (with accession '')
# SOLUTION: make these 2 columns into 4 columns:
#     2 for precursor miRNAs (pre_mirna, pre_accession) and
#     2 for mature miRNAs (mature_mirna, mature_accession).

# 2.1.1 first, split the data into 2 (ones with accession 'MI00*' and others
# with accession ''
phenomir.data.1 = phenomir.data[phenomir.data$accession != '',]
phenomir.data.2 = phenomir.data[phenomir.data$accession == '',]
colnames(phenomir.data.1)[6:7] = c('pre_mirna','pre_accession')
colnames(phenomir.data.2)[6:7] = c('mature_mirna','mature_accession')

# 2.1.2 update precursor miRNA IDs in phenomir.data.1
phenomir.pre_mirna = unique(as.character(phenomir.data.1$pre_mirna))
phenomir.new.pre_mirna = update.precursorID(phenomir.pre_mirna)
phenomir.combine.pre_mirna = phenomir.new.pre_mirna
phenomir.combine.pre_mirna[is.na(phenomir.combine.pre_mirna)] = phenomir.pre_mirna[is.na(phenomir.combine.pre_mirna)]
m = match(phenomir.data.1$pre_mirna, phenomir.pre_mirna)
phenomir.data.1$pre_mirna = phenomir.combine.pre_mirna[m]

# 2.1.3 find mature miRNA IDs & accessions using precursor accessions
# (in phenomir.data.1)
phenomir.pre_acc = unique(as.character(phenomir.data.1$pre_accession))
phenomir.IDs = precursorACC2matureACCID(phenomir.pre_acc)
# combine phenomir.IDs with phenomir.data.1
phenomir.data.1.DT = data.table(phenomir.data.1)
setkey(phenomir.data.1.DT, pre_accession)
phenomir.IDs.DT = data.table(phenomir.IDs)
setkey(phenomir.IDs.DT, mirna_acc)
phenomir.data.1.DT.2 = phenomir.IDs.DT[phenomir.data.1.DT, allow.cartesian=TRUE]
phenomir.data.1.DT.2 = phenomir.data.1.DT.2[,list(mirna_acc,pre_mirna,mature_acc,mature_name,disease,class,expression,name,method,pmid)]
setnames(phenomir.data.1.DT.2, 'pre_mirna', 'mirna_id')

# 2.1.4 update mature miRNA IDs in phenomir.data.2
phenomir.matureID = unique(as.character(phenomir.data.2$mature_mirna))
phenomir.newID = update.matureID(phenomir.matureID)	# some new IDs == NA
phenomir.combineID = phenomir.newID
phenomir.combineID[is.na(phenomir.combineID)] = phenomir.matureID[is.na(phenomir.combineID)]
m = match(phenomir.data.2$mature_mirna, phenomir.matureID)
phenomir.data.2$mature_mirna = phenomir.combineID[m]

# 2.1.5 find mature accessions, precursor miRNA IDs & accessions using
# mature miRNA IDs (in phenomir.data.2)
phenomir.mature_mirna = unique(as.character(phenomir.data.2$mature_mirna))
phenomir.IDs = matureID2precursorACCID(phenomir.mature_mirna)
# combine phenomir.IDs with phenomir.data.2
phenomir.data.2.DT = data.table(phenomir.data.2)
setkey(phenomir.data.2.DT, mature_mirna)
phenomir.IDs.DT = data.table(phenomir.IDs)
setkey(phenomir.IDs.DT, mature_name)
phenomir.data.2.DT.2 = phenomir.IDs.DT[phenomir.data.2.DT, allow.cartesian=TRUE]
phenomir.data.2.DT.2 = phenomir.data.2.DT.2[,list(mirna_acc,mirna_id,mature_acc,mature_name,disease,class,expression,name,method,pmid)]

# 2.1.6 combine phenomir.data.1.DT.2 & phenomir.data.2.DT.2
phenomir.data.3 = rbind(phenomir.data.1.DT.2, phenomir.data.2.DT.2)
dim(phenomir.data.3)
#[1] 20268    10
phenomir.data.3 = data.frame(unique(phenomir.data.3))
dim(phenomir.data.3)
#[1] 15629    10

# 2.2 add miRNA organism
phenomir.org = substr(phenomir.data.3$mirna_id, 1, 4)
phenomir.org = gsub('-','',phenomir.org)
na = is.na(phenomir.org)
phenomir.org[na] = substr(phenomir.data.3$mature_name[na], 1, 4)
phenomir.org[na] = gsub('-','',phenomir.org[na])
phenomir.data.3 = cbind(mirna_org=phenomir.org, phenomir.data.3)

phenomir = unique(phenomir.data.3)

# 2.3 replace NA by '' in miRNA ID columns
for (i in 1:5) {
    tmp = as.character(phenomir[,i])
    tmp[is.na(tmp)] = ''
    phenomir[,i] = tmp
}

######
#
# 3. Populate MySQL tables mirna, target & phenomir
#
######

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', phenomir[,c(1,4,5)])

# update target table
# NO targets in phenomir

# then insert phenomir into table phenomir
phenomir = cbind(mirna.uid, phenomir[,c(2,3,6:11)])
fields = dbListFields(con, 'phenomir')
colnames(phenomir) = fields
dbWriteTable(con, 'phenomir', phenomir, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 4. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.PhenomiR.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(m, i, tmp, con, drv)
save.image("read.PhenomiR.RData")

