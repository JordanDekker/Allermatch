#!/usr/bin/env python

import os,re,string,sys, httplib, urllib

conn = httplib.HTTPConnection('apbioinf004.wurnet.nl')
BaseParams = {
    'database'      :	'All',
    'cutOff'	    :	'35.0',
    'wordlength'    :	'6',
    'against'       :   '',
    'allAlignments' :	'0',
    'rawOutput'     :   '1',
    'method'	    :	'window',
    'seq'           :	'aiscgqvasaiapcisyargqgsgpsagccsgvrslnnaarttadrraacnclknaaagvsglnagnaasipskcgvsipytiststdcsrvnspam' }

def RUNRUNRUN(params):
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain"}
    Params = urllib.urlencode(BaseParams)
    conn.request("POST",
                 "http://apbioinf004.wurnet.nl/~mf/Allermatch/allermatch.py/search",
                 Params, headers)
    response = conn.getresponse()
    if response.reason != "OK":
        raise "ERROR FROM SERVER!!!!!!!!!!" + str(params)
    
    data = response.read()
    conn.close()
    lines = map(string.strip,data.split("\n"))
    results = []
    for line in lines:
        record = {}
        if not line: continue
        for item in line.split("\t"):
            key,val =  item.split(":",1)
            record[key] = val
        results.append(record)
    return results

method = sys.argv[1]
BaseParams['method'] = method
db = sys.argv[2]
database_files = {
    'SwissProt'              : 'SwissProt.fasta',
    'WHO-IUIS'               : 'WHO_All.fasta',
    'SwissProt and WHO-IUIS' : 'All.fasta' }
for x in database_files.keys():
    if database_files[x] == db:
        dbName = x
        print "database",x
        break
    #print database_files
else: raise 'invalid database', db
BaseParams['database'] = dbName

if method == 'wordmatch' and len(sys.argv) > 3:
    wordsize = int(sys.argv[3])
else: wordsize = 6
BaseParams['wordlength'] = str(wordsize)
#read db
seqs = []; firstRun = 1
for line in map(string.strip, open(db).readlines()):
    if not line: continue
    if (line[0] == '>')  and (firstRun):
        firstRun = 0
        #print '1'
        head = line
        seq = ''
    elif line[0] == '>':
        #print '2'
        seqs.append((head,seq))
        head = line
        seq = ''
    else:
        #print '3'
        seq += line
seqs.append((head,seq))
print "loaded", len(seqs), "sequences"

for rec in seqs:
    headLine,seq = rec
    head = headLine[1:].split("\t")
    AMId = head[0]
    sys.stderr.write("%10s %-15s %-15s : " % (method, dbName, AMId) )
    BaseParams['seq'] = seq
    results = RUNRUNRUN(BaseParams)
    resIds = map(lambda x: x['Allergen id'], results)
    if not AMId in resIds:
        print "\t".join([AMId, '-1', str(len(seq))])
        sys.stderr.write("%8s %8s\n" % (-1, len(seq)))
        continue
    #kick out self
    resIds.remove(AMId)
    print "\t".join([AMId, str(len(resIds)), str(len(seq))])
    sys.stderr.write("%8s %8s\n" % (len(resIds), len(seq)))

    

