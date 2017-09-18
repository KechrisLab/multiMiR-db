select * from target where 
target_uid not in (select target_uid from mirecords) 
and target_uid not in (select target_uid from diana_microt)
and target_uid not in (select target_uid from targetscan)
and target_uid not in (select target_uid from pita)
and target_uid not in (select target_uid from pictar)
and target_uid not in (select target_uid from miranda)
and target_uid not in (select target_uid from microcosm)
and target_uid not in (select target_uid from elmmo)
and target_uid not in (select target_uid from tarbase)
and target_uid not in (select target_uid from pharmaco_mir);