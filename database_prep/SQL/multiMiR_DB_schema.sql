
--
-- Table structure for table `mirna`
--

DROP TABLE IF EXISTS `mirna`;
CREATE TABLE `mirna` (
  mature_mirna_uid INTEGER UNSIGNED AUTO_INCREMENT,  -- mature miRNA unique ID
  org VARCHAR(4) NOT NULL,			-- organism abbreviation
  mature_mirna_acc VARCHAR(20) default NULL,	-- mature miRNA accession
  mature_mirna_id VARCHAR(20) default NULL,	-- mature miRNA ID/name
  PRIMARY KEY (mature_mirna_uid),
  KEY org (org),
  KEY mature_mirna_acc (mature_mirna_acc),
  KEY mature_mirna_id (mature_mirna_id)
);

--
-- Table structure for table `target`
--

DROP TABLE IF EXISTS `target`;
CREATE TABLE `target` (
  target_uid INTEGER UNSIGNED AUTO_INCREMENT,	-- target gene unique ID
  org VARCHAR(4) NOT NULL,			-- organism abbreviation
  target_symbol VARCHAR(80) default NULL,	-- target gene symbol
  target_entrez VARCHAR(10) default NULL,	-- target gene Entrez gene ID
  target_ensembl VARCHAR(20) default NULL,	-- target gene Ensembl gene ID
  PRIMARY KEY (target_uid),
  KEY org (org),
  KEY target_symbol (target_symbol),
  KEY target_entrez (target_entrez),
  KEY target_ensembl (target_ensembl)
);

--
-- Table structure for table `mirecords`
--

DROP TABLE IF EXISTS `mirecords`;
CREATE TABLE `mirecords` (
  mature_mirna_uid INTEGER UNSIGNED NOT NULL,	-- mature miRNA unique ID
  target_uid INTEGER UNSIGNED NOT NULL,		-- target gene unique ID
  target_site_number INT(10) default NULL,	-- target site number
  target_site_position INT(10) default NULL,	-- target site position
  experiment VARCHAR(160) default NULL,		-- supporting experiment
  support_type VARCHAR(40) default NULL,	-- type of supporting experiment
  pubmed_id VARCHAR(10) default NULL,		-- PubMed ID
  FOREIGN KEY (mature_mirna_uid)
    REFERENCES mirna(mature_mirna_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  FOREIGN KEY (target_uid)
    REFERENCES target(target_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT
);

--
-- Table structure for table `mirtarbase`
--

DROP TABLE IF EXISTS `mirtarbase`;
CREATE TABLE `mirtarbase` (
  mature_mirna_uid INTEGER UNSIGNED NOT NULL,	-- mature miRNA unique ID
  target_uid INTEGER UNSIGNED NOT NULL,		-- target gene unique ID
  experiment VARCHAR(160) NOT NULL,		-- supporting experiment
  support_type VARCHAR(40) NOT NULL,		-- type of supporting experiment
  pubmed_id VARCHAR(10) default NULL,		-- PubMed ID
  FOREIGN KEY (mature_mirna_uid)
    REFERENCES mirna(mature_mirna_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  FOREIGN KEY (target_uid)
    REFERENCES target(target_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT
);

--
-- Table structure for table `tarbase`
--

DROP TABLE IF EXISTS `tarbase`;
CREATE TABLE `tarbase` (
  mature_mirna_uid INTEGER UNSIGNED NOT NULL,	-- mature miRNA unique ID
  target_uid INTEGER UNSIGNED NOT NULL,		-- target gene unique ID
  experiment VARCHAR(160) NOT NULL,		-- supporting experiment
  support_type VARCHAR(40) NOT NULL,		-- type of supporting experiment
  pubmed_id VARCHAR(10) default NULL,		-- PubMed ID
  FOREIGN KEY (mature_mirna_uid)
    REFERENCES mirna(mature_mirna_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  FOREIGN KEY (target_uid)
    REFERENCES target(target_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT
);

--
-- Table structure for table `miranda`
--

DROP TABLE IF EXISTS `miranda`;
CREATE TABLE `miranda` (
  mature_mirna_uid INTEGER UNSIGNED NOT NULL,	-- mature miRNA unique ID
  target_uid INTEGER UNSIGNED NOT NULL,		-- target gene unique ID
  conservation REAL NOT NULL, 			-- conservation score
  mirsvr_score REAL NOT NULL,			-- mirSVR downregulation score
  FOREIGN KEY (mature_mirna_uid)
    REFERENCES mirna(mature_mirna_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  FOREIGN KEY (target_uid)
    REFERENCES target(target_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  KEY mirsvr_score (mirsvr_score)
);

--
-- Table structure for table `targetscan`
--

DROP TABLE IF EXISTS `targetscan`;
CREATE TABLE `targetscan` (
  mature_mirna_uid INTEGER UNSIGNED NOT NULL,	-- mature miRNA unique ID
  target_uid INTEGER UNSIGNED NOT NULL,		-- target gene unique ID
  site_type INTEGER UNSIGNED NOT NULL,		-- site type (see http://www.targetscan.org/faqs.html)
  context_plus_score REAL NOT NULL,		-- context+ score
  FOREIGN KEY (mature_mirna_uid)
    REFERENCES mirna(mature_mirna_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  FOREIGN KEY (target_uid)
    REFERENCES target(target_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  KEY context_plus_score (context_plus_score)
);

--
-- Table structure for table `diana_microt`
--

DROP TABLE IF EXISTS `diana_microt`;
CREATE TABLE `diana_microt` (
  mature_mirna_uid INTEGER UNSIGNED NOT NULL,	-- mature miRNA unique ID
  target_uid INTEGER UNSIGNED NOT NULL,		-- target gene unique ID
  miTG_score REAL NOT NULL,			-- miRNA target gene score
  UTR3_hit INTEGER UNSIGNED NOT NULL,		-- number of 3'-UTR binding sites
  CDS_hit INTEGER UNSIGNED NOT NULL,		-- number of CDS binding sites
  FOREIGN KEY (mature_mirna_uid)
    REFERENCES mirna(mature_mirna_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  FOREIGN KEY (target_uid)
    REFERENCES target(target_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  KEY miTG_score (miTG_score)
);

--
-- Table structure for table `eimmo`
--

DROP TABLE IF EXISTS `elmmo`;
CREATE TABLE `elmmo` (
  mature_mirna_uid INTEGER UNSIGNED NOT NULL,	-- mature miRNA unique ID
  target_uid INTEGER UNSIGNED NOT NULL,		-- target gene unique ID
  p REAL NOT NULL,				-- "the posterior probability that the site is under evolutionnary selective pressure"
  FOREIGN KEY (mature_mirna_uid)
    REFERENCES mirna(mature_mirna_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  FOREIGN KEY (target_uid)
    REFERENCES target(target_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  KEY p (p)
);

--
-- Table structure for table `pita`
--

DROP TABLE IF EXISTS `pita`;
CREATE TABLE `pita` (
  mature_mirna_uid INTEGER UNSIGNED NOT NULL,	-- mature miRNA unique ID
  target_uid INTEGER UNSIGNED NOT NULL,		-- target gene unique ID
  ddG REAL NOT NULL,				-- ddG = dGduplex (microRNA-target hybridization energy) - dGopen (energy required to make the target site accessible)
  conservation REAL NOT NULL,			-- site conservation (range 0~1)
  FOREIGN KEY (mature_mirna_uid)
    REFERENCES mirna(mature_mirna_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  FOREIGN KEY (target_uid)
    REFERENCES target(target_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  KEY ddG (ddG)
);

--
-- Table structure for table `microcosm`
--

DROP TABLE IF EXISTS `microcosm`;
CREATE TABLE `microcosm` (
  mature_mirna_uid INTEGER UNSIGNED NOT NULL,	-- mature miRNA unique ID
  target_uid INTEGER UNSIGNED NOT NULL,		-- target gene unique ID
  score REAL NOT NULL,				-- miRanda score
  FOREIGN KEY (mature_mirna_uid)
    REFERENCES mirna(mature_mirna_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  FOREIGN KEY (target_uid)
    REFERENCES target(target_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  KEY score (score)
);

--
-- Table structure for table `pictar`
--

DROP TABLE IF EXISTS `pictar`;
CREATE TABLE `pictar` (
  mature_mirna_uid INTEGER UNSIGNED NOT NULL,	-- mature miRNA unique ID
  target_uid INTEGER UNSIGNED NOT NULL,		-- target gene unique ID
  score REAL NOT NULL,				-- score
  FOREIGN KEY (mature_mirna_uid)
    REFERENCES mirna(mature_mirna_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  FOREIGN KEY (target_uid)
    REFERENCES target(target_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  KEY score (score)
);

--
-- Table structure for table `mirdb`
--

DROP TABLE IF EXISTS `mirdb`;
CREATE TABLE `mirdb` (
  mature_mirna_uid INTEGER UNSIGNED NOT NULL,	-- mature miRNA unique ID
  target_uid INTEGER UNSIGNED NOT NULL,		-- target gene unique ID
  score REAL NOT NULL,				-- score
  FOREIGN KEY (mature_mirna_uid)
    REFERENCES mirna(mature_mirna_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  FOREIGN KEY (target_uid)
    REFERENCES target(target_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  KEY score (score)
);

--
-- Table structure for table `mir2disease`
--

DROP TABLE IF EXISTS `mir2disease`;
CREATE TABLE `mir2disease` (
  mature_mirna_uid INTEGER UNSIGNED NOT NULL,	-- mature miRNA unique ID
  disease VARCHAR(100) NOT NULL,		-- disease
  mirna_regulation VARCHAR(20) NOT NULL,	-- how miRNA is regulated in the disease (up-regulated, down-regulated or normal)
  experiment VARCHAR(40) NOT NULL,		-- supporting experiment
  year INTEGER UNSIGNED NOT NULL,		-- year of the paper
  title TEXT NOT NULL,				-- title of the paper
  FOREIGN KEY (mature_mirna_uid)
    REFERENCES mirna(mature_mirna_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  KEY disease (disease)
);

--
-- Table structure for table `pharmaco_mir`
--

DROP TABLE IF EXISTS `pharmaco_mir`;
CREATE TABLE `pharmaco_mir` (
  mature_mirna_uid INTEGER UNSIGNED NOT NULL,	-- mature miRNA unique ID
  target_uid INTEGER UNSIGNED NOT NULL,		-- target gene unique ID
  drug VARCHAR(40) NOT NULL,			-- disease
  pubmed_id VARCHAR(10) NOT NULL,		-- PubMed ID
  FOREIGN KEY (mature_mirna_uid)
    REFERENCES mirna(mature_mirna_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  FOREIGN KEY (target_uid)
    REFERENCES target(target_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  KEY drug (drug)
);

--
-- Table structure for table `phenomir`
--

DROP TABLE IF EXISTS `phenomir`;
CREATE TABLE `phenomir` (
  mature_mirna_uid INTEGER UNSIGNED NOT NULL,	-- mature miRNA unique ID
  pre_mirna_acc VARCHAR(20) default NULL,	-- precursor miRNA accession
  pre_mirna_id VARCHAR(20) default NULL,	-- precursor miRNA ID
  disease VARCHAR(60) NOT NULL,			-- disease
  disease_class VARCHAR(20) NOT NULL,		-- disease class
  mirna_expression VARCHAR(30) NOT NULL,	-- how miRNA is expressed in the disease
  study VARCHAR(40) NOT NULL,			-- type of study (in cells, patients, etc)
  experiment VARCHAR(40) NOT NULL,		-- supporting experiment
  pubmed_id VARCHAR(10) NOT NULL,		-- PubMed ID
  FOREIGN KEY (mature_mirna_uid)
    REFERENCES mirna(mature_mirna_uid)
    ON UPDATE CASCADE ON DELETE RESTRICT,
  KEY disease (disease),
  KEY disease_class (disease_class)
);

--
-- Table structure for table `metadata`
--

DROP TABLE IF EXISTS `metadata`;
CREATE TABLE metadata (
  name VARCHAR(80) PRIMARY KEY,
  value VARCHAR(255)
);

--
-- Table structure for table `map_metadata`
--

DROP TABLE IF EXISTS `map_metadata`;
CREATE TABLE map_metadata (
  map_name VARCHAR(80) PRIMARY KEY,
  source_name VARCHAR(80) NOT NULL,
  source_version VARCHAR(10),
  source_date VARCHAR(20),
  source_url VARCHAR(255) NOT NULL
);

--
-- Table structure for table `map_counts`
--

DROP TABLE IF EXISTS `map_counts`;
CREATE TABLE map_counts (
  map_name VARCHAR(80) PRIMARY KEY,
  count INTEGER UNSIGNED NOT NULL
);

