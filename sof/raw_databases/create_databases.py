#!/usr/bin/env python

import string,sys,os,re

seqs = {}
dbCount = {}
for l in map(string.strip, open(sys.argv[1]).readlines()):
    if not l: continue
    if l[0] == '#': continue
    spl = l.split("\t")
    if len(spl) != 11:
        print 'Error in line '
        print '   ', l
        sys.exit()
    ID, DB = spl[0], spl[1]
    dbCount[DB] = dbCount.get(DB,0) + 1
    if seqs.has_key(ID): raise 'DUPLICATE KEY' + ID
    seqs[ID] = l

print 'Imported', len(seqs), 'sequences'
print 'Per database'
for k in dbCount.keys():
    print " * %s : %d " % (k, dbCount[k])

def export(fn, ss, dblist, dbl_pre, dbl_db):
    doubles = []
    c = 0
    F = open(fn, 'w')
    for seq in ss.values():

        spl = seq.split("\t")
        DB = spl[1]
        #print DB
        #print dblist
        if not DB in dblist: continue

        ID = spl[0]
        export_id = ID
        REALID = ID[3:]
        if REALID in doubles: continue
        

        #find doubles:
        for dblSeq in ss.values():
            dblSpl = dblSeq.split("\t")
            if not dblSpl[1] in dblist: continue
            dblId = dblSpl[0]
            if (dblId[3:] == REALID) and (dblId != ID):
                #this sequence has a double:
                doubles.append(REALID)
                export_id = dbl_pre + REALID
                spl[1] = dbl_db
#                print " #xprt %20s #curr %20s #dbl %20s " % (export_id, ID, dblId)
        c+=1
        F.write(">%s\t%s\n%s\n" %(export_id, "\t".join(spl[1:10]), spl[10]))
    print 'Exported %d sequences to "%s" ' % (c,fn)
    F.close()
    
#output databases
extensie = '.fasta'
#export('SwissProt'+extensie, seqs, ['SwissProt'], 'sp_', 'SwissProt')
export('SwissProt'+extensie, seqs, ['UniProt'], 'sp_', 'SwissProt')
export('WHO_All'  +extensie, seqs, ['WHO-IUIS'], 'wa_', 'WHO-IUIS Allergen')
export('All'      +extensie, seqs, ['WHO-IUIS',
                                    'UniProt'], 'al_', 'SwissProt and WHO-IUIS')
