#!/Library/Frameworks/Python.framework/Versions/3.5/bin/python3

##########
# Spencer Mahaffey
# 12/12/2016
#
# Used to migrate existing tables to new versions of the database.  
# Comment out updated tables to skip them.
# Will join the mirna and target tables from the source DB into the source table
# then will check the destination database for the corresponding ID.
# If the ID exists uses the existing ID in the row otherwise inserts the miRNA or target
# uses the new ID in the row as it is inserted.



import mysql.connector

def insertTarget(conn,target):
  add_row = ("INSERT INTO target ( org, target_symbol, target_entrez,target_ensembl) VALUES ( %s, %s, %s, %s)")
  data_row=(target['organism'],target['sym'],target['entrez'],target['ensembl'])
  insCursor=conn.cursor()
  insCursor.execute(add_row, data_row)
  targetID = insCursor.lastrowid
  conn.commit()
  insCursor.close()
  return targetID

def insertMiRNA(conn,mirna):
  add_row = ("INSERT INTO mirna ( org, mature_mirna_acc, mature_mirna_id) VALUES ( %s, %s, %s)")
  data_row=(mirna['organism'],mirna['acc'],mirna['mid'])
  insCursor=conn.cursor()
  insCursor.execute(add_row, data_row)
  mID = insCursor.lastrowid
  conn.commit()
  insCursor.close()
  return mID

def getNewMiRNA(cache,mrna,destConn):
  newMUid=0
  key=mrna['organism']+"_"+mrna['acc']+"_"+mrna['mid']
  if(key in cache):
    newMUid=cache[key]
  else:
    checkQ=("Select mature_mirna_uid from mirna where org = %s and mature_mirna_acc= %s and mature_mirna_id=%s ")
    mcursor = destConn.cursor()
    mcursor.execute(checkQ,(mrna['organism'],mrna['acc'],mrna['mid']))
    mcount=0
    for (destMUid) in mcursor:
  	  newMUid=destMUid[0]
  	  cache[key]=destMUid[0]
  	  mcount+=1
    if(mcount==0):
      newMUid=insertMiRNA(destConn,mrna)
      cache[key]=newMUid
    elif(mcount>1):
  	  print("error wrong number of miRNAs:",mcount)
  	  newMUid=0
    mcursor.close()
  return newMUid

def getNewTarget(cache,tar,destConn):
  newTUid=0
  key=tar['organism']+"_"+tar['sym']+"_"+tar['entrez']+"_"+tar['ensembl']
  if(key in cache):
    newTUid=cache[key]
  else:
    checkTQ=("Select target_uid from target where org = %s and target_symbol= %s and target_entrez=%s and target_ensembl=%s")
    tcursor = destConn.cursor()
    tcursor.execute(checkTQ,(tar['organism'],tar['sym'],tar['entrez'],tar['ensembl']))
    tcount=0
    for (destTUid) in tcursor:
  	  newTUid=destTUid[0]
  	  cache[key]=destTUid[0]
  	  tcount+=1
    if(tcount==0):
  	  newTUid=insertTarget(destConn,srcTarget[tuid])
  	  cache[key]=newTUid
    elif(tcount>1):
      print("error wrong number of targets:",tcount)
      newTUid=0
    tcursor.close()
  return newTUid

def updateDianaMicroT(srcConn,destConn,srcMiRNA,srcTarget):
  query = ("SELECT * FROM diana_microt ")
  cursor = srcConn.cursor()
  cursor.execute(query)
  table=[]
  for (muid,tuid,score,utr,cds) in cursor:
  	tmp={'muid':muid,'tuid':tuid,'score':score,'utr':utr,'cds':cds}
  	table.append(tmp)
  cursor.close()

  print ("DianaMicroT:",len(table))
  count=0
  cachedMUID={}
  cachedTUID={}
  for x in table:
    mrna=srcMiRNA[x['muid']]
    tar=srcTarget[x['tuid']]
    newMUid=getNewMiRNA(cachedMUID,mrna,destConn)
    newTUid=getNewTarget(cachedTUID,tar,destConn)
    #print("insert ",newMUid,"\t",newTUid,x['score'],x['utr'],x['cds'],tar['sym'],mrna['acc'])
    add_row = ("INSERT INTO diana_microt (mature_mirna_uid, target_uid, miTG_score, UTR3_hit, CDS_hit) VALUES (%s, %s, %s, %s, %s)")
    data_row=(newMUid,newTUid,x['score'],x['utr'],x['cds'])
    insCursor=destConn.cursor()
    insCursor.execute(add_row, data_row)
    destConn.commit()
    insCursor.close()
    count+=1
    if(count%50000==0):
      print("processed ",count)

def updateEimmo(srcConn,destConn,srcMiRNA,srcTarget):
  query = ("SELECT * FROM elmmo ")
  cursor = srcConn.cursor()
  cursor.execute(query)
  table={}
  for (muid,tuid,p) in cursor:
  	table={'muid':muid,'tuid':tuid,'p':p}
  cursor.close()

  print ("Eimmo:",len(table))
  count=0
  cachedMUID={}
  cachedTUID={}
  for x in table:
    mrna=srcMiRNA[x['muid']]
    tar=srcTarget[x['tuid']]
    newMUid=getNewMiRNA(cachedMUID,mrna,destConn)
    newTUid=getNewTarget(cachedTUID,tar,destConn)
    add_row = ("INSERT INTO elmmo (mature_mirna_uid, target_uid, p) VALUES (%s, %s, %s)")
    data_row=(newMUid,newTUid,x['p'])
    insCursor=destConn.cursor()
    insCursor.execute(add_row, data_row)
    destConn.commit()
    insCursor.close()
    count+=1
    if(count%10000==0):
      print("processed ",count)

def updateMicrocosom(srcConn,destConn,srcMiRNA,srcTarget):
  query = ("SELECT * FROM microcosom ")
  cursor = srcConn.cursor()
  cursor.execute(query)
  table={}
  for (muid,tuid,score) in cursor:
  	table={'muid':muid,'tuid':tuid,'score':score}
  cursor.close()

  print ("Microcosom:",len(table))
  count=0
  cachedMUID={}
  cachedTUID={}
  for x in table:
    mrna=srcMiRNA[x['muid']]
    tar=srcTarget[x['tuid']]
    newMUid=getNewMiRNA(cachedMUID,mrna,destConn)
    newTUid=getNewTarget(cachedTUID,tar,destConn)
    add_row = ("INSERT INTO microcosom (mature_mirna_uid, target_uid, score) VALUES (%s, %s, %s)")
    data_row=(newMUid,newTUid,x['score'])
    insCursor=destConn.cursor()
    insCursor.execute(add_row, data_row)
    destConn.commit()
    insCursor.close()
    count+=1
    if(count%10000==0):
      print("processed ",count)

def updateMiranda(srcConn,destConn,srcMiRNA,srcTarget):
  query = ("SELECT * FROM miranda ")
  cursor = srcConn.cursor()
  cursor.execute(query)
  table={}
  for (muid,tuid,conservation,mirsvr_score) in cursor:
  	table={'muid':muid,'tuid':tuid,'conservation':conservation,'score':mirsvr_score}
  cursor.close()

  print ("Miranda:",len(table))
  count=0
  cachedMUID={}
  cachedTUID={}
  for x in table:
    mrna=srcMiRNA[x['muid']]
    tar=srcTarget[x['tuid']]
    newMUid=getNewMiRNA(cachedMUID,mrna,destConn)
    newTUid=getNewTarget(cachedTUID,tar,destConn)
    add_row = ("INSERT INTO miranda (mature_mirna_uid, target_uid, conservation, mirsvr_score) VALUES (%s, %s, %s ,%s)")
    data_row=(newMUid,newTUid,x['conservation'],x['score'])
    insCursor=destConn.cursor()
    insCursor.execute(add_row, data_row)
    destConn.commit()
    insCursor.close()
    count+=1
    if(count%10000==0):
      print("processed ",count)

def updateMir2Disease(srcConn,destConn,srcMiRNA):
  query = ("SELECT * FROM mir2disease ")
  cursor = srcConn.cursor()
  cursor.execute(query)
  table={}
  for (muid,disease,mirnaReg,exp,year,title) in cursor:
  	table={'muid':muid,'disease':disease,'reg':mirnaReg,'exp':exp,'year':year,'title':title}
  cursor.close()

  print ("Mir2disease:",len(table))
  count=0
  cachedMUID={}
  for x in table:
    mrna=srcMiRNA[x['muid']]
    newMUid=getNewMiRNA(cachedMUID,mrna,destConn)
    add_row = ("INSERT INTO mir2disease (mature_mirna_uid, disease, mirna_regulation, experiment, year, title) VALUES (%s, %s, %s,%s, %s, %s, %s)")
    data_row=(newMUid,x['disease'],x['reg'],x['exp'],x['year'],x['title'])
    insCursor=destConn.cursor()
    insCursor.execute(add_row, data_row)
    destConn.commit()
    insCursor.close()
    count+=1
    if(count%10000==0):
      print("processed ",count)

def updateMirecords(srcConn,destConn,srcMiRNA,srcTarget):
  query = ("SELECT * FROM mirecords ")
  cursor = srcConn.cursor()
  cursor.execute(query)
  table={}
  for (muid,tuid,tSitNum,tSitePos,exp,sup,pub) in cursor:
  	table={'muid':muid,'tuid':tuid,'tsite':tSitNum,'tsitepos':tSitePos,'exp':exp,'support':sup,'pubmed':pub}
  cursor.close()

  print ("Mirecords:",len(table))
  count=0
  cachedMUID={}
  cachedTUID={}
  for x in table:
    mrna=srcMiRNA[x['muid']]
    tar=srcTarget[x['tuid']]
    newMUid=getNewMiRNA(cachedMUID,mrna,destConn)
    newTUid=getNewTarget(cachedTUID,tar,destConn)
    add_row = ("INSERT INTO mirecords (mature_mirna_uid, target_uid, target_site_number, target_site_positon, experiment, support_type, pubmed_id) VALUES (%s, %s, %s,%s, %s, %s, %s)")
    data_row=(newMUid,newTUid,x['tsite'],x['tsitepos'],x['exp'],x['support'],x['pubmed'])
    insCursor=destConn.cursor()
    insCursor.execute(add_row, data_row)
    destConn.commit()
    insCursor.close()
    count+=1
    if(count%10000==0):
      print("processed ",count)

def updatePharmacoMir(srcConn,destConn,srcMiRNA,srcTarget):
  query = ("SELECT * FROM pharmaco_mir ")
  cursor = srcConn.cursor()
  cursor.execute(query)
  table={}
  for (muid,tuid,drug,pub) in cursor:
  	table={'muid':muid,'tuid':tuid,'drug':drug, 'pubmed_id':pub}
  cursor.close()

  print ("PharmacoMir:",len(table))
  count=0
  cachedMUID={}
  cachedTUID={}
  for x in table:
    mrna=srcMiRNA[x['muid']]
    tar=srcTarget[x['tuid']]
    newMUid=getNewMiRNA(cachedMUID,mrna,destConn)
    newTUid=getNewTarget(cachedTUID,tar,destConn)
    add_row = ("INSERT INTO pharmaco_mir (mature_mirna_uid, target_uid, drug, pubmed_id) VALUES (%s, %s, %s, %s)")
    data_row=(newMUid,newTUid,x['drug'],x['pubmed_id'])
    insCursor=destConn.cursor()
    insCursor.execute(add_row, data_row)
    destConn.commit()
    insCursor.close()
    count+=1
    if(count%10000==0):
      print("processed ",count)

def updatePhenoMir(srcConn,destConn,srcMiRNA):
  query = ("SELECT * FROM phenomir ")
  cursor = srcConn.cursor()
  cursor.execute(query)
  table={}
  for (muid,pmirnaAcc,pmirnaID,disease,disease_class,expr,study,exp,pub) in cursor:
  	table={'muid':muid,'preMiRnaACC':pmirnaAcc,'preMiRnaID':pmirnaID,'disease':disease,'diseaseClass':diseaseClass,'expr':expr,'study':study,'exp':exp,'pub':pub}
  cursor.close()

  print ("Phenomir:",len(table))
  count=0
  cachedMUID={}
  cachedTUID={}
  for x in table:
    mrna=srcMiRNA[x['muid']]
    tar=srcTarget[x['tuid']]
    newMUid=getNewMiRNA(cachedMUID,mrna,destConn)
    newTUid=getNewTarget(cachedTUID,tar,destConn)
    add_row = ("INSERT INTO phenomir (mature_mirna_uid, target_uid, pre_mirna_acc, pre_mirna_id, disease, disease_class,mirna_expression,study,experiment,pubmed_id) VALUES (%s, %s, %s,%s, %s, %s,%s, %s, %s)")
    data_row=(newMUid,newTUid,x['preMiRnaACC'],x['preMiRnaID'],x['disease'],x['diseaseClass'],x['expr'],x['study'],x['exp'],x['pub'])
    insCursor=destConn.cursor()
    insCursor.execute(add_row, data_row)
    destConn.commit()
    insCursor.close()
    count+=1
    if(count%10000==0):
      print("processed ",count)


def updatePictar(srcConn,destConn,srcMiRNA,srcTarget):
  query = ("SELECT * FROM pictar ")
  cursor = srcConn.cursor()
  cursor.execute(query)
  table={}
  for (muid,tuid,score) in cursor:
  	table={'muid':muid,'tuid':tuid,'score':score}
  cursor.close()

  print ("Pictar:",len(table))
  count=0
  cachedMUID={}
  cachedTUID={}
  for x in table:
    mrna=srcMiRNA[x['muid']]
    tar=srcTarget[x['tuid']]
    newMUid=getNewMiRNA(cachedMUID,mrna,destConn)
    newTUid=getNewTarget(cachedTUID,tar,destConn)
    add_row = ("INSERT INTO pictar (mature_mirna_uid, target_uid, score) VALUES (%s, %s, %s)")
    data_row=(newMUid,newTUid,x['score'])
    insCursor=destConn.cursor()
    insCursor.execute(add_row, data_row)
    destConn.commit()
    insCursor.close()
    count+=1
    if(count%10000==0):
      print("processed ",count)

def updatePITA(srcConn,destConn,srcMiRNA,srcTarget):
  query = ("SELECT * FROM pita ")
  cursor = srcConn.cursor()
  cursor.execute(query)
  table={}
  for (muid,tuid,ddg,conservation) in cursor:
  	table={'muid':muid,'tuid':tuid,'ddg':ddg,'conservation':conservation}
  cursor.close()

  print ("PITA:",len(table))
  count=0
  cachedMUID={}
  cachedTUID={}
  for x in table:
    mrna=srcMiRNA[x['muid']]
    tar=srcTarget[x['tuid']]
    newMUid=getNewMiRNA(cachedMUID,mrna,destConn)
    newTUid=getNewTarget(cachedTUID,tar,destConn)
    add_row = ("INSERT INTO pita (mature_mirna_uid, target_uid, ddG, conservation) VALUES (%s, %s, %s, %s)")
    data_row=(newMUid,newTUid,x['ddg'],x['conservation'])
    insCursor=destConn.cursor()
    insCursor.execute(add_row, data_row)
    destConn.commit()
    insCursor.close()
    count+=1
    if(count%10000==0):
      print("processed ",count)

def updateTarBase(srcConn,destConn,srcMiRNA,srcTarget):
  query = ("SELECT * FROM tarbase ")
  cursor = srcConn.cursor()
  cursor.execute(query)
  table={}
  for (muid,tuid,exp,sup,pub) in cursor:
  	table={'muid':muid,'tuid':tuid,'exp':exp,'support':sup,'pubmed':pub}
  cursor.close()

  print ("TarBase:",len(table))
  count=0
  cachedMUID={}
  cachedTUID={}
  for x in table:
    mrna=srcMiRNA[x['muid']]
    tar=srcTarget[x['tuid']]
    newMUid=getNewMiRNA(cachedMUID,mrna,destConn)
    newTUid=getNewTarget(cachedTUID,tar,destConn)
    add_row = ("INSERT INTO tarbase (mature_mirna_uid, target_uid, experiment, support_type, pubmed_id) VALUES (%s, %s, %s ,%s, %s)")
    data_row=(newMUid,newTUid,x['exp'],x['support'],x['pubmed'])
    insCursor=destConn.cursor()
    insCursor.execute(add_row, data_row)
    destConn.commit()
    insCursor.close()
    count+=1
    if(count%10000==0):
      print("processed ",count)


def updateTargetScan(srcConn,destConn,srcMiRNA,srcTarget):
  query = ("SELECT * FROM targetscan ")
  cursor = srcConn.cursor()
  cursor.execute(query)
  table={}
  for (muid,tuid,site,score) in cursor:
  	table={'muid':muid,'tuid':tuid,'site':site,'score':score}
  cursor.close()

  print ("TargetScan:",len(table))
  count=0
  cachedMUID={}
  cachedTUID={}
  for x in table:
    mrna=srcMiRNA[x['muid']]
    tar=srcTarget[x['tuid']]
    newMUid=getNewMiRNA(cachedMUID,mrna,destConn)
    newTUid=getNewTarget(cachedTUID,tar,destConn)
    add_row = ("INSERT INTO targetscan (mature_mirna_uid, target_uid, site_type, context_plus_score) VALUES (%s, %s, %s, %s)")
    data_row=(newMUid,newTUid,x['site'],x['score'])
    insCursor=destConn.cursor()
    insCursor.execute(add_row, data_row)
    destConn.commit()
    insCursor.close()
    count+=1
    if(count%10000==0):
      print("processed ",count)
srcConn = mysql.connector.connect(user='user', password='password',host='localhost',database='multimir')
destConn = mysql.connector.connect(user='user', password='password',host='localhost',database='multimir_new')

srcMiRNA={}
srcTarget={}

#Get Source MiRNAs
query = ("SELECT * FROM mirna ")
cursor = srcConn.cursor()
cursor.execute(query)
for (uid, org, mirnaACC, mirnaID) in cursor:
  srcMiRNA[uid]={'organism': org, 'acc': mirnaACC, 'mid':mirnaID}
cursor.close()

#Get Source Targets
query = ("SELECT * FROM target ")
cursor = srcConn.cursor()
cursor.execute(query)
for (tuid, org, symbol, entrez, ensembl) in cursor:
  srcTarget[tuid]={'organism': org, 'sym': symbol, 'entrez':entrez, 'ensembl':ensembl}
cursor.close()

#Migrate diana_microt
updateDianaMicroT(srcConn,destConn,srcMiRNA,srcTarget)

#Migrate elmmo
updateEimmo(srcConn,destConn,srcMiRNA,srcTarget)

#Migrate microcosom
updateMicrocosom(srcConn,destConn,srcMiRNA,srcTarget)

#Migrate miranda
updateMiranda(srcConn,destConn,srcMiRNA,srcTarget)

#Migrate mir2disease
updateMir2Disease(srcConn,destConn,srcMiRNA)

#Migrate mirdb
updateMirDB(srcConn,destConn,srcMiRNA,srcTarget)


#Migrate mirecords
updateMirecords(srcConn,destConn,srcMiRNA,srcTarget)

#Migrate mirtarbase
updateMirTarBase()

#Migrate pharmaco_mir
updatePharmacoMir(srcConn,destConn,srcMiRNA,srcTarget)
#Migrate phenomir
updatePhenoMir(srcConn,destConn,srcMiRNA)
#Migrate pictar
updatePictar(srcConn,destConn,srcMiRNA,srcTarget)
#Migrate pita
updatePITA(srcConn,destConn,srcMiRNA,srcTarget)
#Migrate tarbase
updateTarBase(srcConn,destConn,srcMiRNA,srcTarget)
#Migrate targetscan
updateTargetScan(srcConn,destConn,srcMiRNA,srcTarget)

destConn.close()
srcConn.close()