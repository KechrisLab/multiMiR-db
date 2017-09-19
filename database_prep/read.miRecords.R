# read.miRecords.R
# To read and preprocess data from miRecords.
# 
# Status: 9/2017 - Has NOT been updated to use current annotations available through Bioconductor
#

######
#
# 1. Load data, parameters & libraries
#
######

# the newest version is version 4
mirecords.file = "./miRecords/2013_Apr/miRecords_version4.csv"
mirecords.data = read.delim(mirecords.file)

# MySQL parameters
source("build.multiMiR_DB.parameters.R")

# functions
source("multiMiR.R")

######
#
# 2. Preprocess data
#
######

# 2.1 change 'unk...' to '' in mirecords.data$Target.site_position
#mirecords.data$Target.site_position[mirecords.data$Target.site_position == ''] = NA
#mirecords.data$Target.site_position[grep('unk', mirecords.data$Target.site_position)] = NA
mirecords.data$Target.site_position[grep('unk', mirecords.data$Target.site_position)] = ''

# 2.2 change NA to '' in mirecords.data$Target.site_number
mirecords.data$Target.site_number[is.na(mirecords.data$Target.site_number)] = ''

# 2.3 combine information of supporting experiments from columns
#     Test_method_inter, Test_method_inter_site, & Post.mutation_method_site.
mirecords.exp = paste(mirecords.data$Test_method_inter, mirecords.data$Test_method_inter_site, mirecords.data$Post.mutation_method_site, sep="//")
mirecords.exp = gsub("N/A", "", mirecords.exp)
mirecords.exp = gsub("unknown", "", mirecords.exp)
mirecords.exp = gsub("\\}\\{", "//", mirecords.exp)
mirecords.exp = gsub("\\{", "", mirecords.exp)
mirecords.exp = gsub("\\}", "", mirecords.exp)
mirecords.exp = gsub("////", "//", mirecords.exp)
mirecords.exp = gsub("^//", "", mirecords.exp)
mirecords.exp = gsub("//$", "", mirecords.exp)
mirecords.exp = gsub("^ ", "", mirecords.exp)
mirecords.exp = gsub(" $", "", mirecords.exp)
mirecords.exp = gsub("NA//", "", mirecords.exp)
mirecords.exp = gsub("//NA", "", mirecords.exp)
mirecords.exp = gsub("^NA$", "", mirecords.exp)
# remove duplicates
for (i in 1:length(mirecords.exp)) {
    e = unique(strsplit(mirecords.exp[i], "//")[[1]])
    if (length(e) == 0) {
#	mirecords.exp[i] = NA
	mirecords.exp[i] = ''
    }else if (length(e) >= 1) {
	mirecords.exp[i] = paste(e, collapse = "//")
    }
}
mirecords.data2 = cbind(mirecords.data[,c(7,6,3,4,2,5,17,1)], mirecords.exp)
mirecords.data2 = unique(mirecords.data2)

# 2.4 process species columns
mirecords.data2$miRNA_species = gsub("^ ", "", mirecords.data2$miRNA_species)
mirecords.data2$miRNA_species = gsub(" $", "", mirecords.data2$miRNA_species)
mirecords.data2$miRNA_species[mirecords.data2$miRNA_species == 'Mus musculus miR-24'] = 'Mus musculus'
mirecords.data2$miRNA_species[mirecords.data2$miRNA_species == 'us musculus'] = 'Mus musculus'
mirecords.data2$miRNA_species[mirecords.data2$miRNA_species == 'Xenopus. tropicalis'] = 'Xenopus tropicalis'

mirecords.data2$Target.gene_species_scientific = gsub("^ ", "", mirecords.data2$Target.gene_species_scientific)
mirecords.data2$Target.gene_species_scientific = gsub(" $", "", mirecords.data2$Target.gene_species_scientific)
mirecords.data2$Target.gene_species_scientific = gsub("  ", " ", mirecords.data2$Target.gene_species_scientific)
mirecords.data2$Target.gene_species_scientific = gsub("\\. ", " ", mirecords.data2$Target.gene_species_scientific)
mirecords.data2$Target.gene_species_scientific = gsub('NM_213232.1', 'Danio rerio', mirecords.data2$Target.gene_species_scientific)
mirecords.data2$Target.gene_species_scientific[mirecords.data2$Target.gene_species_scientific == 'human'] = 'Homo sapiens'
mirecords.data2$Target.gene_species_scientific[mirecords.data2$Target.gene_species_scientific == 'Xenopus'] = 'Xenopus tropicalis'

# 2.5 remove records without target gene information
# there are 3 lines to remove
#mirecords.data2[mirecords.data2$Target.gene_name == '',]
#    miRNA_mature_ID miRNA_species Target.gene_name Target.gene_Refseq_acc
#946      mmu-miR-93  Mus musculus                               NC_001560
#947    mmu-miR-146a  Mus musculus                               NC_001560
#948    mmu-miR-378*  Mus musculus                               NC_001560
#        Target.gene_species_scientific Target.site_number Target.site_position
#946 Vesicular stomatitis Indiana virus                  1                 1873
#947 Vesicular stomatitis Indiana virus                  1                 2073
#948 Vesicular stomatitis Indiana virus                  1                 1462
#    Pubmed_id mirecords.exp
#946  17613256              
#947  17613256              
#948  17613256              
mirecords.data2 = mirecords.data2[mirecords.data2$Target.gene_name != '',]

# 2.6 process column miRNA_mature_ID
mirecords.data2$miRNA_mature_ID = gsub("has", "hsa", mirecords.data2$miRNA_mature_ID)
mirecords.data2$miRNA_mature_ID = gsub("\\[", "", mirecords.data2$miRNA_mature_ID)
mirecords.data2$miRNA_mature_ID = gsub("\\]", "", mirecords.data2$miRNA_mature_ID)
mirecords.data2$miRNA_mature_ID[grep('P-27', mirecords.data2$miRNA_mature_ID)] = 'P-27-5p'
# add species abbreviations to the mature miRNA IDs if they don't have any
mirecords.data2$miRNA_mature_ID = add.org.to.mature.miRNA.ID(ID=mirecords.data2$miRNA_mature_ID, org=mirecords.data2$miRNA_species)

# 2.7 update mature miRNA IDs
mirecords.data2$miRNA_mature_ID = gsub(" ", "", mirecords.data2$miRNA_mature_ID)
mirecords.matureID = unique(as.character(mirecords.data2$miRNA_mature_ID))
mirecords.newID = update.matureID(mirecords.matureID)	# some new IDs == NA
mirecords.combineID = mirecords.newID
mirecords.combineID[is.na(mirecords.combineID)] = mirecords.matureID[is.na(mirecords.combineID)]
m = match(mirecords.data2$miRNA_mature_ID, mirecords.matureID)
mirecords.data2$miRNA_mature_ID = mirecords.combineID[m]

# 2.8 add mature miRNA accession numbers
mirecords.matureACC = matureID2matureACC(mirecords.combineID)
mirecords.data2 = cbind(mature_mirna_acc=mirecords.matureACC[m], mirecords.data2)
mirecords.data2 = unique(mirecords.data2)

# 2.9 combine experiments for the same record
mirecords.pair = unique(mirecords.data2[,1:9])
if (nrow(mirecords.pair) < nrow(mirecords.data2)) {
    mirecords.data3 = NULL
    for (i in 1:nrow(mirecords.pair)) {
	m = rep(TRUE, nrow(mirecords.data2))
	for (j in 2:ncol(mirecords.pair)) {
	    if (!is.na(mirecords.pair[i,j])) {
		m = m & mirecords.data2[,j] == mirecords.pair[i,j]
	    }
	}
	m = which(m)
	m = m[mirecords.data2$mirecords.exp[m] != '']
	exp = paste(mirecords.data2$mirecords.exp[m], collapse="//")
	mirecords.data3 = rbind(mirecords.data3, cbind(mirecords.pair[i,], mirecords.exp=exp))
    }
}else {
    mirecords.data3 = mirecords.data2
}

# 2.10 add Entrez & Ensembl gene IDs
# NOTE: Entrez & Ensembl gene IDs are retrieved using gene symbols.
mirecords.org.unique = unique(as.character(mirecords.data3$Target.gene_species_scientific))
mirecords.IDs = NULL
for (org in mirecords.org.unique) {
    m = mirecords.data3$Target.gene_species_scientific == org
    symbol = unique(as.character(mirecords.data3$Target.gene_name[m]))
    if (org == "[Virus]") {
	next
    }
    ID = convert.gene.IDs(org=org, ID=symbol, ID.type='mixed')
    if (!is.null(ID)) {
	if (nrow(ID) > 0) {
	    ID = cbind(ID, org=org)
            mirecords.IDs = rbind(mirecords.IDs, ID)
	}
    }
}
mirecords.IDs = unique(mirecords.IDs)

# 2.11 combine mirecords.IDs & mirecords.data3
mirecords.data4 = NULL
n = ncol(mirecords.data3)
for (i in 1:nrow(mirecords.data3)) {
    m = which(mirecords.IDs$gene_symbol == as.character(mirecords.data3[i,4]) & mirecords.IDs$org == as.character(mirecords.data3[i,6]))
    if (length(m) == 0) {
	combine = cbind(mirecords.data3[i,1:4], entrez_gene_id=NA, ensembl_gene_id=NA, mirecords.data3[i,5:n])
    }else {
	combine = cbind(mirecords.data3[i,1:4], entrez_gene_id=mirecords.IDs[m,2], ensembl_gene_id=mirecords.IDs[m,3], mirecords.data3[i,5:n])
    }
    mirecords.data4 = rbind(mirecords.data4, combine)
}

mirecords = unique(mirecords.data4)
mirecords = cbind(mirecords[,1:10], experiment=mirecords[,12], support_type=NA, pubmed_id=mirecords[,11])
dim(mirecords)
#[1] 3469   13

# Keep miRNA-target pairs in human only, or in mouse only
keep = (mirecords[,3] == 'Homo sapiens' & mirecords[,8] == 'Homo sapiens') | (mirecords[,3] == 'Mus musculus' & mirecords[,8] == 'Mus musculus')
mirecords = mirecords[keep,]
dim(mirecords)
#[1] 2874   13
mirecords[mirecords[,3] == 'Homo sapiens', 3] = 'hsa'
mirecords[mirecords[,8] == 'Homo sapiens', 8] = 'hsa'
mirecords[mirecords[,3] == 'Mus musculus', 3] = 'mmu'
mirecords[mirecords[,8] == 'Mus musculus', 8] = 'mmu'

# replace NA by '' in miRNA and gene ID columns
for (i in 1:8) {
    tmp = as.character(mirecords[,i])
    tmp[is.na(tmp)] = ''
    mirecords[,i] = tmp
}

######
#
# 3. Populate MySQL tables mirna, target & mirecords
#
######

drv = dbDriver("MySQL")
con = dbConnect(drv, host=DB.HOST, user=DB.USER, password=DB.PASS, dbname=DB)

# update mirna table
mirna.uid = update.mirna.table(con, 'mirna', mirecords[,c(3,1,2)])

# update target table
target.uid = update.target.table(con, 'target', mirecords[,c(8,4:6)])

# then insert mirecords into table mirecords
mirecords = cbind(mirna.uid, target.uid, mirecords[,9:13])
fields = dbListFields(con, 'mirecords')
colnames(mirecords) = fields
dbWriteTable(con, 'mirecords', mirecords, row.names=FALSE, append=TRUE)

dbDisconnect(con)
dbUnloadDriver(drv)

######
#
# 4. Back up the MySQL database and save RData
#
######

if (DUMP.DB) {
    dump.file = "./SQL/multimir_mysqldump_after_read.miRecords.R.sql"
    command = paste("mysqldump --host=", DB.HOST, " --user=", DB.USER, " --password=", DB.PASS, " --databases ", DB, " --add-drop-database > ", dump.file, sep="")
    system(command, intern=TRUE)
}

rm(i, e, m, org, symbol, ID, n, combine, keep, tmp, con, drv)
save.image("read.miRecords.RData")

