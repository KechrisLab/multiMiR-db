# build.multiMiR_DB.parameters.R

CREATE.USER = FALSE
CREATE.USER.FILE = "./SQL/create_multiMiR_user.sql"
CREATE.DB = TRUE
CREATE.TABLES = TRUE
DB.HOST = 'localhost'
DB.USER = 'user'
DB.PASS = 'password'
DB = 'multimir'
DB.SCHEMA.FILE = "./SQL/multiMiR_DB_schema.sql"
DUMP.DB = TRUE

META.SCHEMA.FILE = "./SQL/multiMiR_DB_schema.meta.sql"
RELOAD.META.SCHEMA = TRUE

PROCESS.mirecords = FALSE
PROCESS.mirtarbase = FALSE
PROCESS.tarbase = FALSE
PROCESS.mir2disease = FALSE
PROCESS.phenomir = FALSE
PROCESS.pharmaco_mir = FALSE
PROCESS.mirdb = FALSE
PROCESS.pictar_hsa = FALSE
PROCESS.pictar_mmu = FALSE
PROCESS.microcosm_hsa = FALSE
PROCESS.microcosm_mmu = FALSE
PROCESS.miranda_hsa = FALSE
PROCESS.miranda_mmu = FALSE
PROCESS.targetscan_hsa = FALSE
PROCESS.targetscan_mmu = FALSE
PROCESS.diana_microt = FALSE
PROCESS.eimmo_hsa = FALSE
PROCESS.eimmo_mmu = FALSE
PROCESS.pita_hsa = FALSE
PROCESS.pita_mmu = FALSE

TABLES = c(
'diana_microt',
'elmmo',
'microcosm',
'mir2disease',
'miranda',
'mirdb',
'mirecords',
'mirtarbase',
'pharmaco_mir',
'phenomir',
'pictar',
'pita',
'tarbase',
'targetscan')

