# read.TarBase.R
# To read and preprocess data from TarBase.
# 
# Status: 9/2017 - Has NOT been updated to use current annotations available through Bioconductor
#
#####
# 1. Read data from downloaded file
#

tarbase.file = "./TarBase/v6.0/tarbase_data.csv"
tarbase.data = read.csv(tarbase.file)
dim(tarbase.data)
#[1] 66011    13
tarbase.data = unique(tarbase.data)	# a lot of duplicated rows
dim(tarbase.data)
#[1] 36032    13

# MySQL parameters
source("build.multiMiR_DB.parameters.R")

# functions
source("multiMiR.R")

library(data.table)

#
# 2. Preprocess data
#

# 2.1 there is 1 row without target gene information
#     (ensgid == "\\N" & gene.name == "") -- remove it for now
remove = tarbase.data$ensgid == "\\N" & tarbase.data$gene.name == ""
tarbase.data = tarbase.data[!remove,]

# 2.2 change "\\N" to NA in Ensembl gene IDs (column ensgid)
tarbase.data$ensgid[tarbase.data$ensgid == "\\N"] = NA

# 2.3 update mature miRNA IDs
# NOTE: In column miRNA.name, there are 3 entries that are not named in a
#	traditional way. They are lna_let-7b, bantam, & edited-hsa-mir-376a-5p.
#	There are other entries with '/', such as hsa-mir-221/222cluster,
#	hsa-mir-17/20, & hsa-mir-132/mir-212.
#	Leave them as they are now.
tarbase.matureID = unique(as.character(tarbase.data$miRNA.name))
tarbase.newID = update.matureID(tarbase.matureID)       # some new IDs == NA
tarbase.combineID = tarbase.newID
tarbase.combineID[is.na(tarbase.combineID)] = tarbase.matureID[is.na(tarbase.combineID)]
m = match(tarbase.data$miRNA.name, tarbase.matureID)
tarbase.data2 = cbind(miRNA.mimat=tarbase.data[,1], miRNA.name=tarbase.combineID[m], tarbase.data[3:ncol(tarbase.data)])
tarbase.data2 = unique(tarbase.data2)
dim(tarbase.data2)
#[1] 36031    13

# 2.4 get species information from miRNA names
# NOTE: assume species is same for both miRNA & target since there is no other
#	species information
tarbase.org = substr(tarbase.data2$miRNA.name,1,4)
tarbase.org = gsub('-','',tarbase.org)
# process miRNA names with 'edit'
if (sum(tarbase.org == 'edit') > 0) {
    edit = strsplit(as.character(tarbase.data2$miRNA.name[tarbase.org == 'edit']), split='-')
    edit.org = NULL
    for (i in 1:length(edit)) {edit.org[i] = edit[[i]][2]}
}
tarbase.org[tarbase.org == 'edit'] = edit.org
tarbase.org.unique = unique(tarbase.org)
# NOTE: 'lna_' & 'bant' are not species names in tarbase.org

# 2.5 correct some errors in gene symbols (column gene.name)
#     delete 'd' at the end of some Arabidopsis gene names (as in AT1G63230d)
tarbase.gene.name = as.character(tarbase.data2$gene.name)
m = grep('AT.G.....d', tarbase.gene.name, ignore.case=TRUE)
tarbase.gene.name[m] = gsub('d','',tarbase.gene.name[m])
# change Os0611310 to Os06g11310
tarbase.gene.name[tarbase.gene.name == 'Os0611310'] = 'Os06g11310'
# add 'LOC_' to the beginning of rice gene names
m = grep('OS..G.....', tarbase.gene.name, ignore.case=TRUE)
tarbase.gene.name[m] = paste('LOC_', tarbase.gene.name[m], sep='')
tarbase.data2 = cbind(tarbase.data2[,1:3], gene.name=tarbase.gene.name, tarbase.data2[,5:ncol(tarbase.data2)])

# 2.6 add Entrez Gene IDs
# NOTE: Entrez and Ensembl gene IDs are retrieved using gene symbols.
#	It is not checked whether the retrieved Ensembl gene IDs are the same as
#	the original Ensembl gene IDs in tarbase.data.
#	The retrieved Ensembl gene IDs are used eventually.
# On 10/29/2013: errors ("Incorrect BioMart name") for some plant species
# - to fix.
# Right now only keep human and mouse entries.
keep = tarbase.org == 'hsa' | tarbase.org == 'mmu'
tarbase.org = tarbase.org[keep]
tarbase.data2 = tarbase.data2[keep,]
dim(tarbase.data2)
#[1] 30471    13
tarbase.IDs = NULL
for (org in c('hsa','mmu')) {
    m = grep(org, tarbase.data2$miRNA.name)
    symbol = unique(tarbase.data2$gene.name[m])
    ID = convert.gene.IDs(org=org, ID=symbol, ID.type='mixed')
    if (!is.null(ID)) {
	if (nrow(ID) > 0) {
	    ID = cbind(ID, org=org)
	    tarbase.IDs = rbind(tarbase.IDs, ID)
	}
    }
}
tarbase.IDs = unique(tarbase.IDs)
dim(tarbase.IDs)
#[1] 14251     4

# 2.7 combine tarbase.IDs with tarbase.data2
tarbase.data2.DT = cbind(gene_symbol=tarbase.data2[,4], org=tarbase.org, ensembl_gene_id=tarbase.data2[,3], tarbase.data2[,c(1,2,5:ncol(tarbase.data2))])
tarbase.data2.DT = data.table(tarbase.data2.DT)
setkey(tarbase.data2.DT, gene_symbol, org)
tarbase.IDs.DT = data.table(tarbase.IDs[,c(1,4,3,2)])
setkey(tarbase.IDs.DT, gene_symbol, org)
tarbase.data3.DT = tarbase.IDs.DT[tarbase.data2.DT, allow.cartesian=TRUE]
# NOTE: column 'ensembl_gene_id.1' in tarbase.data3.DT is the same as column
# 'ensgid' in tarbase.data2

tarbase.data3.DT = tarbase.data3.DT[,list(miRNA.mimat,miRNA.name,gene_symbol,entrez_gene_id,ensembl_gene_id,reporter_gene,nothern_blot,western_blot,qPCR,proteomics,microarray,sequencing,degradome_seq,other,org,ensembl_gene_id.1)]
tarbase.data3.DT = unique(tarbase.data3.DT)
tarbase.data3 = data.frame(tarbase.data3.DT)
# if new ensembl_gene_id is NA, use the old ensembl_gene_id
tarbase.data3$ensembl_gene_id[is.na(tarbase.data3$ensembl_gene_id)] = as.character(tarbase.data3$ensembl_gene_id.1[is.na(tarbase.data3$ensembl_gene_id)])

# One of the examples that tarbase.data3 may miss the original ensembl_gene_id
#		gene_symbol	ensembl_gene_id	entrez_gene_id
# tarbase.data2	AC010724.6-6	ENSG00000205270	NA
# tarbase.data3	AC010724.6-6	NA	NA
# tarbase.IDs	no 'AC010724.6-6' since neither ensembl_gene_id nor entrez_gene_id was found
# If using the method in ../mirtar.validated.db/read.TarBase.R, tarbase.data2
# & tarbase.data3 will be same for this example. But the method here is much
# faster (<0.1 second vs. ~2 hours).

# One of the examples that tarbase.data3 may miss the new entrez_gene_id
#		gene_symbol	ensembl_gene_id	entrez_gene_id
# tarbase.data2	Glyma08g01450	NA	NA
# tarbase.data3	Glyma08g01450	NA	NA
# tarbase.IDs	''	GLYMA08G01450	100805328

# 2.8 combine experiments (it runs for 45 min - quicker way?)
tarbase.pair = unique(tarbase.data3[,1:5])
experiments = c('Reporter assay','Northern blot','Western blot','qRT-PCR','Proteomics','Microarray','Sequencing','Degradome sequencing','Other')
tarbase.data4 = NULL
for (i in 1:nrow(tarbase.pair)) {
    m = rep(TRUE, nrow(tarbase.data3))
    for (j in 1:ncol(tarbase.pair)) {
	if (!is.na(tarbase.pair[i,j])) {
	    m = m & tarbase.data3[,j] == tarbase.pair[i,j]
	}
    }
    m = which(m)
    positive = NULL
    negative = NULL
    unknown = NULL
    for (j in m) {
	for (k in 6:14) {
	    if (tarbase.data3[j,k] == 'POSITIVE') {
		positive = paste(positive, experiments[k-5], sep="//")
	    }else if (tarbase.data3[j,k] == 'NEGATIVE') {
		negative = paste(negative, experiments[k-5], sep="//")
	    }else if (tarbase.data3[j,k] == 'UNKNOWN') {
		unknown = paste(unknown, experiments[k-5], sep="//")
	    }
	}
    }
    if (!is.null(positive)) {
	positive = sub("//","",positive)
	combine = cbind(tarbase.pair[i,], experiment=positive, support_type='positive', org=tarbase.data3[j,15])
	tarbase.data4 = rbind(tarbase.data4, combine)
    }
    if (!is.null(negative)) {
	negative = sub("//","",negative)
	combine = cbind(tarbase.pair[i,], experiment=negative, support_type='negative', org=tarbase.data3[j,15])
	tarbase.data4 = rbind(tarbase.data4, combine)
    }
    if (!is.null(unknown)) {
	unknown = sub("//","",unknown)
	combine = cbind(tarbase.pair[i,], experiment=unknown, support_type='unknown', org=tarbase.data3[j,15])
	tarbase.data4 = rbind(tarbase.data4, combine)
    }
    if (is.null(positive) & is.null(negative) & is.null(unknown)) {
	warning("No experiments for ", tarbase.pair[i,1], ' ', tarbase.pair[i,2], ' ', tarbase.pair[i,3], ' ', tarbase.pair[i,4], ' ', tarbase.pair[i,5])
    }
}
tarbase.data4 = unique(tarbase.data4)

# 2.9 add pubmed_id column
tarbase = cbind(tarbase.data4, pubmed_id=NA)

# 2.10 replace NA by '' in miRNA and gene ID columns
for (i in c(1:5,8)) {
    tmp = as.character(tarbase[,i])
    tmp[is.na(tmp)] = ''
    tarbase[,i] = tmp
}

# 2.11 change "'CUSTOM'" to "CUSTOM" in column miRNA.mimat
tarbase$miRNA.mimat[tarbase$miRNA.mimat == "'CUSTOM'"] = "CUSTOM"

######
#
# 3. Populate MySQL tables mirna, target & tarbase
#
######

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', tarbase[,c(8,1,2)])

# update target table
# NOTE: no RefSeq numbers
target.uid = update.target.table(con, 'target', tarbase[,c(8,3,4,5)])

# then insert tarbase into table tarbase
tarbase = cbind(mirna.uid, target.uid, tarbase[,c(6,7,9)])
fields = dbListFields(con, 'tarbase')
colnames(tarbase) = fields
dbWriteTable(con, 'tarbase', tarbase, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 4. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.TarBase.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(remove, m, edit, edit.org, i, keep, org, symbol, ID, j, k, positive, negative, unknown, combine, tmp, drv, con)
save.image("read.TarBase.RData")


