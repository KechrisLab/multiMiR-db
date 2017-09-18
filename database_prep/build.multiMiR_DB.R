# build.multiMiR_DB.R
# To build the multiMiR MySQL database.

######
#
# 1. Load libraries, set default parameters
#
######

# use RMySQL to create our database 
library(RMySQL)

# use mirbase.db for microRNA information
#library(mirbase.db)

# use biomaRt for gene information
library(biomaRt)

if (!"ALT_BIOMART_HOST" %in% ls(pat = "ALT_BIOMART_HOST")) ALT_BIOMART_HOST = FALSE
if (ALT_BIOMART_HOST) {
    source("multiMiR.alt_biomaRt_host.R")
}else {
    source("multiMiR.R")
}

# MySQL parameters
source("build.multiMiR_DB.parameters.R")

######
#
# 2. Create user, database and tables
#
######

# NOTE: Need to create user 'multimir' first
if (CREATE.USER) {
    cat("Creating MySQL user", DB.USER, "...\n")
    command = paste("CREATE USER '", DB.USER, "'@'localhost' IDENTIFIED BY '", DB.PASS, "';", sep="")
    write(command, file=CREATE.USER.FILE, append=FALSE)
    command = paste("CREATE USER '", DB.USER, "'@'%' IDENTIFIED BY '", DB.PASS, "';", sep="")
    write(command, file=CREATE.USER.FILE, append=TRUE)
    command = paste("GRANT ALL PRIVILEGES ON multimir.* TO '", DB.USER, "'@'localhost';", sep="")
    write(command, file=CREATE.USER.FILE, append=TRUE)
    command = paste("GRANT ALL PRIVILEGES ON multimir.* TO '", DB.USER, "'@'%';", sep="")
    write(command, file=CREATE.USER.FILE, append=TRUE)
    command = paste("mysql -h localhost -u root -pUbu@Think1 < ", CREATE.USER.FILE, sep="")
    x = system(command, intern=TRUE)
}

# Create MySQL database DB
if (CREATE.DB) {
    cat("Creating MySQL database", DB, "...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("echo 'CREATE DATABASE ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
}

# Create tables in database DB
if (CREATE.TABLES) {
    cat("Creating tables in MySQL database", DB, "...\n")
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " ", DB, " < ", DB.SCHEMA.FILE, sep="")
    x = system(command, intern=TRUE)
}

######
#
# 3. Read and process miRNA-target data and populate tables
#
######

# read and process miRecords - table mirecords
if (PROCESS.mirecords) {
    cat("Processing miRecords ...\n")
    x = system("R -f read.miRecords.R", intern=TRUE)
}else {
    cat("miRecords: populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.miRecords.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process miRTarBase - table mirtarbase
if (PROCESS.mirtarbase) {
    cat("Processing miRTarBase ...\n")
    x = system("R -f read.miRTarBase.R", intern=TRUE)
}else {
    cat("miRTarBase: populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.miRTarBase.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process TarBase - table tarbase
if (PROCESS.tarbase) {
    cat("Processing TarBase ...\n")
    x = system("R -f read.TarBase.R", intern=TRUE)
}else {
    cat("TarBase: populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.TarBase.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process miR2Disease - table mir2disease
if (PROCESS.mir2disease) {
    cat("Processing miR2Disease ...\n")
    x = system("R -f read.miR2Disease.R", intern=TRUE)
}else {
    cat("miR2Disease: populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.miR2Disease.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process PhenomiR - table phenomir
if (PROCESS.phenomir) {
    cat("Processing PhenomiR ...\n")
    x = system("R -f read.PhenomiR.R", intern=TRUE)
}else {
    cat("PhenomiR: populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.PhenomiR.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process Pharmaco-miR Verified Set - table pharmaco_mir
if (PROCESS.pharmaco_mir) {
    cat("Processing Pharmaco-miR ...\n")
    x = system("R -f read.Pharmaco-miR_VerSe.R", intern=TRUE)
}else {
    cat("Pharmaco-miR: populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.Pharmaco-miR_VerSe.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process miRDB - table mirdb
if (PROCESS.mirdb) {
    cat("Processing miRDB ...\n")
    x = system("R -f read.miRDB.R", intern=TRUE)
}else {
    cat("miRDB: populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.miRDB.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process PicTar (human) - table pictar
if (PROCESS.pictar_hsa) {
    cat("Processing PicTar (human) ...\n")
    x = system("R -f read.PicTar_hsa.R", intern=TRUE)
}else {
    cat("PicTar (human): populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.PicTar_hsa.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process PicTar (mouse) - table pictar
if (PROCESS.pictar_mmu) {
    cat("Processing PicTar (mouse) ...\n")
    x = system("R -f read.PicTar_mmu.R", intern=TRUE)
}else {
    cat("PicTar (mouse): populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.PicTar_mmu.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process MicroCosm (human) - table microcosm
if (PROCESS.microcosm_hsa) {
    cat("Processing MicroCosm (human) ...\n")
    x = system("R -f read.MicroCosm_hsa.R", intern=TRUE)
}else {
    cat("MicroCosm (human): populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.MicroCosm_hsa.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process MicroCosm (mouse) - table microcosm
if (PROCESS.microcosm_mmu) {
    cat("Processing MicroCosm (mouse) ...\n")
    x = system("R -f read.MicroCosm_mmu.R", intern=TRUE)
}else {
    cat("MicroCosm (mouse): populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.MicroCosm_mmu.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process miRanda (human) - table miranda
if (PROCESS.miranda_hsa) {
    cat("Processing miRanda (human) ...\n")
    x = system("R -f read.miRanda_hsa.R", intern=TRUE)
}else {
    cat("miRanda (human): populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.miRanda_hsa.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process miRanda (mouse) - table miranda
if (PROCESS.miranda_mmu) {
    cat("Processing miRanda (mouse) ...\n")
    x = system("R -f read.miRanda_mmu.R", intern=TRUE)
}else {
    cat("miRanda (mouse): populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.miRanda_mmu.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process TargetScan (human) - table targetscan
if (PROCESS.targetscan_hsa) {
    cat("Processing TargetScan (human) ...\n")
    x = system("R -f read.TargetScan_hsa_conserved.R", intern=TRUE)
}else {
    cat("TargetScan (human): populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.TargetScan_hsa_conserved.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process TargetScan (mouse) - table targetscan
if (PROCESS.targetscan_mmu) {
    cat("Processing TargetScan (mouse) ...\n")
    x = system("R -f read.TargetScan_mmu_conserved.R", intern=TRUE)
}else {
    cat("TargetScan (mouse): populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.TargetScan_mmu_conserved.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process DIANA-microT - table diana_microt
if (PROCESS.diana_microt) {
    cat("Processing DIANA-microT ...\n")
    x = system("R -f read.DIANA-microT.R", intern=TRUE)
}else {
    cat("DIANA-microT: populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.DIANA-microT.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process EIMMo (human) - table eimmo
if (PROCESS.eimmo_hsa) {
    cat("Processing EIMMo (human) ...\n")
    x = system("R -f read.EIMMo_hsa.R", intern=TRUE)
}else {
    cat("EIMMo (human): populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.EIMMo_hsa.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process EIMMo (mouse) - table eimmo
if (PROCESS.eimmo_mmu) {
    cat("Processing EIMMo (mouse) ...\n")
    x = system("R -f read.EIMMo_mmu.R", intern=TRUE)
}else {
    cat("EIMMo (mouse): populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.EIMMo_mmu.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process PITA (human) - table pita
if (PROCESS.pita_hsa) {
    cat("Processing PITA (human) ...\n")
    x = system("R -f read.PITA_hsa.R", intern=TRUE)
}else {
    cat("PITA (human): populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.PITA_hsa.R.sql", sep="")
    x = system(command, intern=TRUE)
}

# read and process PITA (mouse) - table pita
if (PROCESS.pita_mmu) {
    cat("Processing PITA (mouse) ...\n")
    x = system("R -f read.PITA_mmu.R", intern=TRUE)
}else {
    cat("PITA (mouse): populating tables ...\n")
    command = paste("echo 'DROP DATABASE IF EXISTS ", DB, "' | mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, sep="")
    x = system(command, intern=TRUE)
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " < ./SQL/multimir_mysqldump_after_read.PITA_mmu.R.sql", sep="")
    x = system(command, intern=TRUE)
}

######
#
# 4. Read and process metadata
#
######

if (RELOAD.META.SCHEMA) {
    cat("Re-creating meta tables in MySQL database", DB, "...\n")
    command = paste("mysql -h ", DB.HOST, " -u ", DB.USER, " -p", DB.PASS, " ", DB, " < ", META.SCHEMA.FILE, sep="")
    x = system(command, intern=TRUE)
}

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# table metadata
cat("Processing metadata ...\n")
metadata.file = "multiMiR.metadata.csv"
metadata = read.delim(metadata.file)
dbWriteTable(con, 'metadata', metadata[1:4,], row.names=FALSE, append=TRUE)
if (dbGetQuery(con, "SELECT count(*) FROM metadata") != nrow(metadata[1:4,])) {
    stop("Table metadata doesn't have the correct number of records!")
}

# table map_metadata
cat("Processing map_metadata ...\n")
map_metadata = data.frame(matrix(metadata[-(1:4),2], ncol=4, byrow=TRUE))
map_metadata = data.frame(TABLES, map_metadata[,c(1,3,4,2)])
colnames(map_metadata) = dbListFields(con, 'map_metadata')
dbWriteTable(con, 'map_metadata', map_metadata, row.names=FALSE, append=TRUE)
if (dbGetQuery(con, "SELECT count(*) FROM map_metadata") != nrow(map_metadata)) {
    stop("Table map_metadata doesn't have the correct number of records!")
}

# table map_counts
cat("Processing map_counts ...\n")
map_counts = NULL
for (table in TABLES) {
  if (table %in% c('mir2disease','phenomir')) {
    query = paste("SELECT count(*) FROM ", table, " AS i INNER JOIN mirna AS m ON (m.mature_mirna_uid=i.mature_mirna_uid) WHERE m.org='hsa'", sep='')
  }else {
    query = paste("SELECT count(*) FROM ", table, " AS i INNER JOIN mirna AS m INNER JOIN target AS t ON (m.mature_mirna_uid=i.mature_mirna_uid AND i.target_uid=t.target_uid) WHERE m.org='hsa' AND t.org='hsa'", sep='')
  }
  hsa.n = dbGetQuery(con, query)
  if (table %in% c('mir2disease','phenomir')) {
    query = paste("SELECT count(*) FROM ", table, " AS i INNER JOIN mirna AS m ON (m.mature_mirna_uid=i.mature_mirna_uid) WHERE m.org='mmu'", sep='')
  }else {
    query = paste("SELECT count(*) FROM ", table, " AS i INNER JOIN mirna AS m INNER JOIN target AS t ON (m.mature_mirna_uid=i.mature_mirna_uid AND i.target_uid=t.target_uid) WHERE m.org='mmu' AND t.org='mmu'", sep='')
  }
  mmu.n = dbGetQuery(con, query)
  if (table %in% c('mir2disease','phenomir')) {
    query = paste("SELECT count(*) FROM ", table, " AS i INNER JOIN mirna AS m ON (m.mature_mirna_uid=i.mature_mirna_uid) WHERE m.org='rno'", sep='')
  }else {
    query = paste("SELECT count(*) FROM ", table, " AS i INNER JOIN mirna AS m INNER JOIN target AS t ON (m.mature_mirna_uid=i.mature_mirna_uid AND i.target_uid=t.target_uid) WHERE m.org='rno' AND t.org='rno'", sep='')
  }
  rno.n = dbGetQuery(con, query)
  query = paste("SELECT count(*) FROM ", table, sep='')
  total.n = dbGetQuery(con, query)
  map_counts = rbind(map_counts, c(table, hsa.n[[1]], mmu.n[[1]], rno.n[[1]], total.n[[1]]))
}
colnames(map_counts) = dbListFields(con, 'map_counts')
map_counts = data.frame(map_counts)
dbWriteTable(con, 'map_counts', map_counts, row.names=FALSE, append=TRUE)
if (dbGetQuery(con, "SELECT count(*) FROM map_counts") != nrow(map_counts)) {
  stop("Table map_counts doesn't have the correct number of records!")
}

######
#
# 5. calculate preset cutoffs
#
######

cat("Updating multimir_cutoffs.rda ...\n")
multimir_cutoffs <- get.multimir.cutoff(con=con)
save(multimir_cutoffs, file="./multimir_cutoffs.rda")

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 6. Dump the final version of the database
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_FINAL.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

save.image("build.multiMiR_DB.RData")

