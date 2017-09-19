# read.Pharmaco-miR_VerSe.R
# To read and preprocess data from Pharmaco-miR VerSe (Verified Sets).
# 
# Status: 9/2017 - Has NOT been updated to use current annotations available through Bioconductor
#
######
#
# 1. Read data from downloaded file
#
######

pharmaco_mir.file = "./Pharmaco-miR_VerSe/pharmacomir_VERSE_DB.csv"
pharmaco_mir.data = read.csv(pharmaco_mir.file)

# MySQL parameters
source("build.multiMiR_DB.parameters.R")

# functions
source("multiMiR.R")

######
#
# 2. Preprocess data
#
######

# 2.1 try to guess organisms:
#     If title contains 'mouse' or 'mice', organism = mouse;
#     otherwise, organism = human
pharmaco_mir.org = rep('hsa', nrow(pharmaco_mir.data))
pharmaco_mir.org[c(grep('mice', pharmaco_mir.data$Paper), grep('mouse', pharmaco_mir.data$Paper))] = 'mmu'

# 2.2 add species abbreviations to mature miRNA IDs
pharmaco_mir.mir = paste(pharmaco_mir.org, as.character(pharmaco_mir.data$miRNA), sep="-")
pharmaco_mir.data$miRNA = pharmaco_mir.mir

# 2.3 update mature miRNA IDs
pharmaco_mir.matureID = unique(pharmaco_mir.mir)
pharmaco_mir.newID = update.matureID(pharmaco_mir.matureID)	# some new IDs == NA
pharmaco_mir.combineID = pharmaco_mir.newID
pharmaco_mir.combineID[is.na(pharmaco_mir.combineID)] = pharmaco_mir.matureID[is.na(pharmaco_mir.combineID)]
m = match(pharmaco_mir.data$miRNA, pharmaco_mir.matureID)
pharmaco_mir.data$miRNA = pharmaco_mir.combineID[m]

# 2.4 add mature miRNA accession numbers
pharmaco_mir.matureACC = matureID2matureACC(pharmaco_mir.combineID)
pharmaco_mir.data = cbind(mature_mirna_acc=pharmaco_mir.matureACC[m], pharmaco_mir.data)

# 2.5 add miRNA organism
pharmaco_mir.data = cbind(org=pharmaco_mir.org, pharmaco_mir.data)

# 2.6 retrieve Entrez & Ensembl gene IDs using gene symbols
pharmaco_mir.IDs = NULL
for (org in unique(pharmaco_mir.org)) {
    m = grep(org, pharmaco_mir.org)
    symbol = unique(pharmaco_mir.data$Gene[m])
    ID = convert.gene.IDs(org=org, ID=symbol, ID.type='symbol')
    if (!is.null(ID)) {
        if (nrow(ID) > 0) {
            ID = cbind(ID, org=org)
            pharmaco_mir.IDs = rbind(pharmaco_mir.IDs, ID)
        }
    }
}
pharmaco_mir.IDs = unique(pharmaco_mir.IDs)
pharmaco_mir.IDs$gene_symbol = toupper(pharmaco_mir.IDs$gene_symbol)

# 2.7 combine pharmaco_mir.IDs with pharmaco_mir.data
pharmaco_mir.data.DT = data.table(pharmaco_mir.data)
setkey(pharmaco_mir.data.DT, org, Gene)
pharmaco_mir.IDs.DT = data.table(pharmaco_mir.IDs)
setkey(pharmaco_mir.IDs.DT, org, gene_symbol)
pharmaco_mir.data2.DT = pharmaco_mir.IDs.DT[pharmaco_mir.data.DT, allow.cartesian=TRUE]
pharmaco_mir.data2.DT = pharmaco_mir.data2.DT[,list(org,mature_mirna_acc,miRNA,gene_symbol,entrez_gene_id,ensembl_gene_id,Drug,PubMed.Id)]

pharmaco_mir = unique(data.frame(pharmaco_mir.data2.DT))

# 2.8 replace NA by '' in miRNA and gene ID columns
for (i in 1:6) {
    tmp = as.character(pharmaco_mir[,i])
    tmp[is.na(tmp)] = ''
    pharmaco_mir[,i] = tmp
}

######
#
# 3. Populate MySQL tables mirna, target & pharmaco_mir
#
######

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', pharmaco_mir[,1:3])

# update target table
target.uid = update.target.table(con, 'target', pharmaco_mir[,c(1,4:6)])

# then insert pharmaco_mir into table pharmaco_mir
pharmaco_mir = cbind(mirna.uid, target.uid, pharmaco_mir[,7:8])
fields = dbListFields(con, 'pharmaco_mir')
colnames(pharmaco_mir) = fields
dbWriteTable(con, 'pharmaco_mir', pharmaco_mir, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 4. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.Pharmaco-miR_VerSe.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(m, org, symbol, ID, i, tmp, drv, con)
save.image("read.Pharmaco-miR_VerSe.RData")

