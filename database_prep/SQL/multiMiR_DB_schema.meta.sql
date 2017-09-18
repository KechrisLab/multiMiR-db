-- multiMiR_DB_schema.meta.sql

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

