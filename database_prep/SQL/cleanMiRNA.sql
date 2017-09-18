select * from mirna where 
mature_mirna_uid not in (select mature_mirna_uid from mirecords)
and mature_mirna_uid not in (select mature_mirna_uid from diana_microt)
and mature_mirna_uid not in (select mature_mirna_uid from targetscan)
and mature_mirna_uid not in (select mature_mirna_uid from pita)
and mature_mirna_uid not in (select mature_mirna_uid from pictar)
and mature_mirna_uid not in (select mature_mirna_uid from miranda)
and mature_mirna_uid not in (select mature_mirna_uid from microcosm)
and mature_mirna_uid not in (select mature_mirna_uid from elmmo)
and mature_mirna_uid not in (select mature_mirna_uid from tarbase)
and mature_mirna_uid not in (select mature_mirna_uid from mir2disease)
and mature_mirna_uid not in (select mature_mirna_uid from phenomir)
and mature_mirna_uid not in (select mature_mirna_uid from pharmaco_mir);