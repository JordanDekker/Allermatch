#!/usr/bin/env python2.3

import sys,string,re,os

if len(sys.argv) > 1: dr = sys.argv[1]
else: dr = '.'
fileList = os.listdir(dr)
fileList.sort()
for File in fileList:
    if not File: continue
    if File[-8:] != '.selfLen': continue

    _lf = File.split('.')
    Db = _lf[1]
    Met = _lf[0]
    WL  = int(_lf[3])
    HL = {}
    for line in map(string.strip, open(File).readlines()[2:]):
        Id, Hits, Len = line.split()
        Hits = int(Hits)
        Len = int(Len)
        HL[Hits] = HL.get(Hits,0) + 1
    All = sum(HL.values())
    noH = HL.get(0,0)
    noS = HL.get(-1,0)


    
    print "%20s %10s %2d All %6d 1: %6d 0: %6d %6.2f " % (
        Db, Met, WL, All, noH, noS, float(noS+noH)/All*100)

    
    
