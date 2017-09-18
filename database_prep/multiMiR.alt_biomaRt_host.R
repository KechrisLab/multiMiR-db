# multiMiR.alt_biomaRt_host.R
# This file uses alternative biomaRt server
# http://www.ensembl.org/biomart/martservice/ instead of http://www.biomart.org


library(mirbase.db)
library(biomaRt)


# Convert a list of mature miRNA IDs to a list of mature miRNA accession numbers
matureID2matureACC <- function (matureID) {
    matureACC <- NULL
    for (i in 1:length(matureID)) {
	query <- paste("SELECT mature_acc FROM mirna_mature WHERE mature_name LIKE '", matureID[i], "' LIMIT 1", sep='')
	result <- dbGetQuery(mirbase_dbconn(), query)
	matureACC[i] <- ifelse(length(result[[1]]) == 1, result[[1]], NA)
#	if (length(result[[1]]) == 1) {
#	    matureACC[i] <- result[[1]]
#	}else {
#	    matureACC[i] <- NA
#	}
    }
    return (matureACC)
}


# Update mature miRNA IDs based on miRBase
update.matureID <- function (matureID) {
    newID <- NULL
    for (i in 1:length(matureID)) {
	query1 <- paste("SELECT count(*) FROM mirna_mature WHERE mature_name LIKE '", matureID[i], "'", sep='')
	result1 <- dbGetQuery(mirbase_dbconn(), query1)
	if (result1[[1]] > 0) {
	# match to mature miRNA ID in miRBase
	    newID[i] = matureID[i]
	}else {
	# no match, search in previous mature miRNA IDs
	    query2 <- paste("SELECT mature_name FROM mirna_mature WHERE previous_mature_id LIKE '", matureID[i], "' OR previous_mature_id LIKE '", matureID[i], ";%' OR previous_mature_id LIKE '%;", matureID[i], ";%' OR previous_mature_id LIKE '%;", matureID[i], "' LIMIT 1", sep='')
	    result2 <- dbGetQuery(mirbase_dbconn(), query2)
	    if (length(result2[[1]]) == 1) {
		newID[i] <- result2[[1]]
	    }else {
		newID[i] <- NA
	    }
	}
    }
    return (newID)
}


# Update precursor miRNA IDs based on miRBase
update.precursorID <- function (precursorID) {
    newID <- NULL
    for (i in 1:length(precursorID)) {
	query1 <- paste("SELECT count(*) FROM mirna WHERE mirna_id LIKE '", precursorID[i], "'", sep='')
	result1 <- dbGetQuery(mirbase_dbconn(), query1)
	if (result1[[1]] > 0) {
	# match to precursor miRNA ID in miRBase
	    newID[i] = precursorID[i]
	}else {
	# no match, search in previous precursor miRNA IDs
	    query2 <- paste("SELECT mirna_id FROM mirna WHERE previous_mirna_id LIKE '", precursorID[i], "' OR previous_mirna_id LIKE '", precursorID[i], ";%' OR previous_mirna_id LIKE '%;", precursorID[i], ";%' OR previous_mirna_id LIKE '%;", precursorID[i], "' LIMIT 1", sep='')
	    result2 <- dbGetQuery(mirbase_dbconn(), query2)
	    if (length(result2[[1]]) == 1) {
		newID[i] <- result2[[1]]
	    }else {
		newID[i] <- NA
	    }
	}
    }
    return (newID)
}


# Convert a list of precursor miRNA accession numbers to mature miRNA IDs &
# accession numbers
precursorACC2matureACCID <- function (precursorACC) {
    precursorACC <- paste(precursorACC, collapse="','")
    precursorACC <- paste("('", precursorACC, "')", sep='')

    query <- paste("SELECT mirna.mirna_acc, mirna.mirna_id, mirna_mature.mature_acc, mirna_mature.mature_name FROM mirna, mirna_pre_mature, mirna_mature WHERE mirna._id = mirna_pre_mature._id AND mirna_pre_mature.auto_mature = mirna_mature.auto_mature AND mirna.mirna_acc IN ", precursorACC, sep='')
    result <- dbGetQuery(mirbase_dbconn(), query)
    return (result)
}


# Convert a list of mature miRNA IDs to precursor miRNA IDs & accession numbers
matureID2precursorACCID <- function (matureID) {
    matureID <- paste(matureID, collapse="','")
    matureID <- paste("('", matureID, "')", sep='')

    query <- paste("SELECT mirna.mirna_acc, mirna.mirna_id, mirna_mature.mature_acc, mirna_mature.mature_name FROM mirna, mirna_pre_mature, mirna_mature WHERE mirna._id = mirna_pre_mature._id AND mirna_pre_mature.auto_mature = mirna_mature.auto_mature AND mirna_mature.mature_name IN ", matureID, sep='')
    result <- dbGetQuery(mirbase_dbconn(), query)
    return (result)
}


# convert between gene symbols, Ensembl & Entrez gene IDs using biomaRt
# ID.type could be:
#	'symbol', to convert from gene symbols to Ensembl & Entrez gene IDs;
#	'ensembl', to convert from Ensembl gene IDs to gene symbols & Entrez gene IDs;
#	'entrez', to convert from Entrez gene IDs to gene symbols & Ensembl gene IDs;
#	'mixed', the input is a mix of at least 2 of the 3 types of IDs;
#	'refseq', to convert from RefSeq accessions to gene symbols, Ensembl & Entrez gene IDs.
convert.gene.IDs <- function (org='hsa', ID=NULL, ID.type=NULL) {
    # check Internet connection (needed for biomaRt)
    # -- TO DO

    # find organism in the organism table
    org.table <- read.delim("multiMiR.species.csv")
    org.info <- find.org(org, org.table)

    if (is.na(org.info[1])) {
	warning("Species ", org, " could not be found!")
	return (NULL)
    }

    # convert gene IDs
    converted <- NULL
    if (ID.type == 'symbol') {
	converted <- Symbol2EnsemblEntrez(org=org, org.info=org.info, ID)
    }else if (ID.type == 'ensembl') {
	converted <- Ensembl2EntrezSymbol(org=org, org.info=org.info, ID)
    }else if (ID.type == 'ensembl_transcript') {
	converted <- EnsemblTranscript2Others(org=org, org.info=org.info, ID)
    }else if (ID.type == 'entrez') {
	converted <- Entrez2EnsemblSymbol(org=org, org.info=org.info, ID)
    }else if (ID.type == 'mixed') {	# a mix of gene symbols, Entrez & Ensembl gene IDs
	converted.1 <- Symbol2EnsemblEntrez(org=org, org.info=org.info, ID)
	converted.2 <- Ensembl2EntrezSymbol(org=org, org.info=org.info, ID)
	converted.3 <- Entrez2EnsemblSymbol(org=org, org.info=org.info, ID)
	converted <- rbind(converted.1, converted.2, converted.3)
	converted <- unique(converted)
    }else if (ID.type == 'refseq') {
	converted <- RefSeq2Others(org=org, org.info=org.info, ID)
    }

    # remove the rows with only 1 ID
    NA.count <- apply(converted, 1, function(x) {sum(is.na(x))})
    keep <- NA.count < (ncol(converted) - 1)
    converted <- converted[keep,]

    return (converted)
}


# Convert gene symbols to Ensembl & Entrez gene IDs using biomaRt
Symbol2EnsemblEntrez <- function (org='hsa', org.info, symbol) {
    IDs <- NULL
    if (!is.na(org.info[,4])) {
	org.data.type <- 'ENSEMBL_MART_ENSEMBL'
	org.data <- useMart(org.data.type, dataset=as.character(org.info[,4]), host="www.ensembl.org")
	symbol.column <- as.character(org.info[,6])
	if (is.na(symbol.column)) {symbol.column <- 'wikigene_name'}
	IDs <- getBM(attributes=c(symbol.column,'entrezgene','ensembl_gene_id'), filters=symbol.column, values=symbol, mart=org.data)
	# for rnorvegicus_gene_ensembl: mgi_symbol or rgd_symbol?
#    }else if (!is.na(org.info[,5])) {
#	org.data.type <- 'plants_mart_18'
#	org.data <- useMart(org.data.type, dataset=as.character(org.info[,5]))
#	symbol.column <- as.character(org.info[,6])
#	if (is.na(symbol.column)) {symbol.column <- 'wikigene_name'}
#	IDs <- getBM(attributes=c(symbol.column,'entrezgene','ensembl_gene_id'), filters=symbol.column, values=symbol, mart=org.data)
    }else {
	IDs <- cbind(symbol, rep(NA, length(symbol)), rep(NA, length(symbol)))
    }
    colnames(IDs) = c('gene_symbol', 'entrez_gene_id', 'ensembl_gene_id')
    return (IDs)
}


# Convert Entrez gene IDs to gene symbols & Ensembl gene IDs using biomaRt
Entrez2EnsemblSymbol <- function (org='hsa', org.info, entrez) {
    IDs <- NULL
    if (!is.na(org.info[,4])) {
	org.data.type <- 'ENSEMBL_MART_ENSEMBL'
	org.data <- useMart(org.data.type, dataset=as.character(org.info[,4]), host="www.ensembl.org")
	symbol.column <- as.character(org.info[,6])
	if (is.na(symbol.column)) {symbol.column <- 'wikigene_name'}
	IDs <- getBM(attributes=c(symbol.column,'entrezgene','ensembl_gene_id'), filters='entrezgene', values=entrez, mart=org.data)
#    }else if (!is.na(org.info[,5])) {
#	org.data.type <- 'plants_mart_18'
#	org.data <- useMart(org.data.type, dataset=as.character(org.info[,5]))
#	symbol.column <- as.character(org.info[,6])
#	if (is.na(symbol.column)) {symbol.column <- 'wikigene_name'}
#	IDs <- getBM(attributes=c(symbol.column,'entrezgene','ensembl_gene_id'), filters='entrezgene', values=entrez, mart=org.data)
    }else {
	IDs <- cbind(rep(NA, length(entrez)), entrez, rep(NA, length(entrez)))
    }
    colnames(IDs) = c('gene_symbol', 'entrez_gene_id', 'ensembl_gene_id')
    return (IDs)
}


# Convert Ensembl gene IDs to gene symbols & Entrez gene IDs using biomaRt
Ensembl2EntrezSymbol <- function (org='hsa', org.info, ensembl) {
    IDs <- NULL
    if (!is.na(org.info[,4])) {
	org.data.type <- 'ENSEMBL_MART_ENSEMBL'
	org.data <- useMart(org.data.type, dataset=as.character(org.info[,4]), host="www.ensembl.org")
	symbol.column <- as.character(org.info[,6])
	if (is.na(symbol.column)) {symbol.column <- 'wikigene_name'}
	IDs <- getBM(attributes=c(symbol.column,'entrezgene','ensembl_gene_id'), filters='ensembl_gene_id', values=ensembl, mart=org.data)
#    }else if (!is.na(org.info[,5])) {
#	org.data.type <- 'plants_mart_18'
#	org.data <- useMart(org.data.type, dataset=as.character(org.info[,5]))
#	symbol.column <- as.character(org.info[,6])
#	if (is.na(symbol.column)) {symbol.column <- 'wikigene_name'}
#	IDs <- getBM(attributes=c(symbol.column,'entrezgene','ensembl_gene_id'), filters='ensembl_gene_id', values=ensembl, mart=org.data)
    }else {
	IDs <- cbind(rep(NA, length(ensembl)), rep(NA, length(ensembl)), ensembl)
    }
    colnames(IDs) = c('gene_symbol', 'entrez_gene_id', 'ensembl_gene_id')
    return (IDs)
}


# Convert Ensembl transcript IDs to gene symbols, Entrez & Ensembl gene IDs using biomaRt
EnsemblTranscript2Others <- function (org='hsa', org.info, ensembl) {
    IDs <- NULL
    if (!is.na(org.info[,4])) {
	org.data.type <- 'ENSEMBL_MART_ENSEMBL'
	org.data <- useMart(org.data.type, dataset=as.character(org.info[,4]), host="www.ensembl.org")
	symbol.column <- as.character(org.info[,6])
	if (is.na(symbol.column)) {symbol.column <- 'wikigene_name'}
	IDs <- getBM(attributes=c(symbol.column,'entrezgene','ensembl_gene_id','ensembl_transcript_id'), filters='ensembl_transcript_id', values=ensembl, mart=org.data)
#    }else if (!is.na(org.info[,5])) {
#	org.data.type <- 'plants_mart_18'
#	org.data <- useMart(org.data.type, dataset=as.character(org.info[,5]))
#	symbol.column <- as.character(org.info[,6])
#	if (is.na(symbol.column)) {symbol.column <- 'wikigene_name'}
#	IDs <- getBM(attributes=c(symbol.column,'entrezgene','ensembl_gene_id','ensembl_transcript_id'), filters='ensembl_transcript_id', values=ensembl, mart=org.data)
    }else {
	IDs <- cbind(rep(NA, length(ensembl)), rep(NA, length(ensembl)), rep(NA, length(ensembl)), ensembl)
    }
    colnames(IDs) = c('gene_symbol', 'entrez_gene_id', 'ensembl_gene_id', 'ensembl_transcript_id')
    return (IDs)
}


# Convert RefSeq IDs to gene symbols, Entrez & Ensembl gene IDs
RefSeq2Others <- function (org='hsa', org.info, refseq) {
    IDs <- NULL
    if (!is.na(org.info[,4]) & !is.na(org.info[,7])) {
	org.data.type <- 'ENSEMBL_MART_ENSEMBL'
	org.data <- useMart(org.data.type, dataset=as.character(org.info[,4]), host="www.ensembl.org")
	symbol.column <- as.character(org.info[,6])
	refseq.column <- as.character(org.info[,7])
	if (is.na(symbol.column)) {symbol.column <- 'wikigene_name'}
	IDs <- getBM(attributes=c(symbol.column,'entrezgene','ensembl_gene_id',refseq.column), filters=refseq.column, values=refseq, mart=org.data)
	if (length(grep('predicted',refseq.column)) == 0) {
	    refseq.column.2 <- paste(refseq.column, "_predicted", sep='')
	    IDs.2 <- getBM(attributes=c(symbol.column,'entrezgene','ensembl_gene_id',refseq.column.2), filters=refseq.column.2, values=refseq, mart=org.data)
	    if (nrow(IDs.2) > 0) {
		colnames(IDs.2)[4] <- refseq.column
		IDs <- rbind(IDs, IDs.2)
		IDs <- unique(IDs)
	    }
	}
#    }else if (!is.na(org.info[,5]) & !is.na(org.info[,7])) {
#	org.data.type <- 'plants_mart_18'
#	org.data <- useMart(org.data.type, dataset=as.character(org.info[,5]))
#	symbol.column <- as.character(org.info[,6])
#	refseq.column <- as.character(org.info[,7])
#	if (is.na(symbol.column)) {symbol.column <- 'wikigene_name'}
#	IDs <- getBM(attributes=c(symbol.column,'entrezgene','ensembl_gene_id',refseq.column), filters=refseq.column, values=refseq, mart=org.data)
    }else {
	IDs <- cbind(rep(NA, length(refseq)), rep(NA, length(refseq)), rep(NA, length(refseq)), refseq)
    }
    colnames(IDs) = c('gene_symbol', 'entrez_gene_id', 'ensembl_gene_id', 'refseq_acc')
    return (IDs)
}


# Find the organism in an organism table
find.org <- function (org='hsa', org.table) {
    m1 = match(org, org.table[,1])	# organism abbreviation
    if (is.na(m1)) {
	m2 = grep(org, org.table[,2])	# organism name
#	org.info <- ifelse(length(m2) == 0, NA, org.table[m2,])
	if (length(m2) == 0) {
	    org.info <- NA
	}else {
	    org.info <- org.table[m2,]
	}
    }else {
	org.info <- org.table[m1,]
    }
    return (org.info)
}


# Some mature miRNA IDs don't have species abbreviation at the beginning -
# add species abbreviations to these IDs
add.org.to.mature.miRNA.ID <- function (ID, org) {
    first4 <- substr(ID, 1, 4)
    add.list <- first4 %in% c('mir-', 'miR-', 'let-')

    # find organism in the organism table
    org.table <- read.delim("multiMiR.species.csv")
    add.org = NULL
    for (i in 1:sum(add.list)) {
	org.info <- find.org(org=org[add.list][i], org.table)
	add.org = c(add.org, as.character(org.info[,1]))
    }

    ID[add.list] = paste(add.org, ID[add.list], sep='-')

    return (ID)
}


# Convert between organism abbreviation and organism name
org.conversion <- function (org, type='abbr2name') {
    org.table <- read.delim("multiMiR.species.csv")
    if (type == 'abbr2name') {
	m <- match(org, org.table[,1])
	converted <- as.character(org.table[m,2])
    }else if (type == 'name2abbr') {
	m <- match(org, org.table[,2])
	converted <- as.character(org.table[m,1])
    }
    return (converted)
}


# Calculate preset cutoffs
# top.p, the top percentile (i.e. 0.1 is top 10th percentile)
# top.n, the top nth miRNA-target interactions (i.e. 1000 is the top 1000 interactions)
calculate.cutoffs <- function (table=c('diana_microt','eimmo_hsa','microcosm_hsa','miranda_hsa','pita_hsa','targetscan_hsa'), top.p=seq(0.1,0.6,0.1), top.n=1:6*100000) {
    cutoffs = list()
    for (t in table) {

    }
    return (cutoffs)
}


# get score cutoff in a table.
# If no parameters are provided, the default is to calculate score cutoffs for
# top 1% to top 100% and for top 10,000 up to all predictions in each of the 8
# predicted tables in human and mouse.
get.multimir.cutoff <- function(db=NULL, table=c('diana_microt','eimmo','microcosm','miranda','mirdb','pictar','pita','targetscan'), org=c('hsa','mmu'), cutoff=NULL) {
    if (is.null(db)) return (NULL)
    CUTOFFS = list()
    if (is.null(cutoff)) {
	for (t in table) {
	    for (o in org) {
		# prepare query to retrieve all scores from the table
		if (t == 'diana_microt') {
		    q <- paste("SELECT miTG_score FROM ", t, " WHERE org = '", o, "' ORDER BY miTG_score DESC", sep="")
		}else if (t == 'eimmo') {
		    q <- paste("SELECT p FROM ", t, " WHERE org = '", o, "' ORDER BY p DESC", sep="")
		}else if (t == 'microcosm' | t == 'mirdb' | t == 'pictar') {
		    q <- paste("SELECT score FROM ", t, " WHERE org = '", o, "' ORDER BY score DESC", sep="")
		}else if (t == 'miranda') {
		    q <- paste("SELECT mirsvr_score FROM ", t, " WHERE org = '", o, "' ORDER BY mirsvr_score", sep="")
		}else if (t == 'pita') {
		    q <- paste("SELECT ddG FROM ", t, " WHERE org = '", o, "' ORDER BY ddG", sep="")
		}else if (t == 'targetscan') {
		    q <- paste("SELECT context_plus_score FROM ", t, " WHERE org = '", o, "' ORDER BY context_plus_score", sep="")
		}
		# retrieve scores
		scores <- dbGetQuery(db, q)
		scores <- scores[,1]

		# calculate cutoffs for top 1% to top 99%
		name <- paste(t, o, sep=".")
		if (t %in% c('miranda','pita','targetscan')) {
		    CUTOFFS[[name]] = quantile(scores, probs=1:99/100, names=FALSE)
		}else {
		    CUTOFFS[[name]] = quantile(scores, probs=99:1/100, names=FALSE)
		}
		names(CUTOFFS[[name]]) = paste(1:99, "%", sep="")

		# calculate cutoffs for top 10,000 to all predictions
		n <- as.integer(length(scores)/10000)
		for (i in 1:n) {
		    CUTOFFS[[name]][[paste(i, '0000', sep='')]] = scores[i*10000]
		}
		CUTOFFS[[name]][['count']] = length(scores)
	    }
	}
    }else {

    }
    return (CUTOFFS)
}

