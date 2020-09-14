import  random, time, os, sys, string, re, copy
sys.path.append(os.path.dirname(__file__))
from allermatch_settings import *

def extractSeq(seq):
    s = ''
    for l in map(string.strip, seq.split("\n")):
        if not l: continue
        if l[0] == '>': continue
        s += l
    return s

def formatSeq(sq):
    r = ''; s = sq
    while len(s) > 0:
      r += s[:60] + "\n"
      s = s[60:]
    return r

def formatDetails(rec):
    #push rec-record in local namespace
    for x in rec.keys(): locals()[x.replace(' ', '_')] = rec[x]
    prc = '%'; sqf =  formatSeq(rec["Sequence"])
    return """<table cellpadding='5' cellspacing='0' width='100%(prc)s'>
        <tr><td valign='top'><b>Allergen Id</b></td>
            <td valign='top'>%(Allergen_id)s</td></tr>
        <tr bgcolor='#dfdeda'><td valign='top'><b>Allermatch<sup>tm</sup> Database</b>
            </td><td> %(Database_Name)s</td></tr>
        <tr><td valign='top'><b>Allergen Name</b>
            </td><td> %(Allergen_Name)s</td></tr>
        <tr bgcolor='#dfdeda'><td  valign='top'><b>Source database</b></td>
            <td valign='top'> %(Source_db)s</td></tr>
        <tr><td  valign='top'><b>Accession Id</b></td>
            <td valign='top'> %(Accession_id)s</td></tr>
        <tr bgcolor='#dfdeda'><td  valign='top'><b>External link</b></td>
             <td valign='top'><a href='%(Hyperlink)s'>
             %(Hyperlink)s</a></td></tr>
        <tr><td  valign='top'><b>Scientific Name</b></td>
            <td valign='top'>%(Species_name)s</td></tr>
        <tr bgcolor='#dfdeda'><td  valign='top'><b>Common Name</b></td>
            <td valign='top'>%(English_name)s</td></tr>
        <tr><td  valign='top'><b>Description</b></td>
            <td valign='top'>%(Remark)s</td></tr>
        <tr bgcolor='#dfdeda'><td  valign='top'><b>Size mature protein</b></td>
            <td valign='top'>%(Size_mature_protein)s</td></tr>
        <tr><td  valign='top'><b>Sequence</b></td>
            <td valign='top'><pre>%(sqf)s</pre></td></tr>
      </table> """ % vars()

def formatDetailsRow(rec, rowNo, showSeq=0, link = 0, bgcolor='#dfdeda', database='AllergenDB'):
    if rowNo % 2 == 0: bgcolor=bgcolor
    else: bgcolor='#FFFFFF'
    if link: l = '<a name="%s"></a>' % (rec["Allergen id"])
    else: l = ''
    x  = "<tr bgcolor='%s'>" % (bgcolor)
    x += "<td valign='top'>%s</td>"   %(l + rec["Allergen id"])
    x += "<td valign='top'>%s</td>"   %(l + rec["Database Name"])
    x += "<td valign='top'>%s</td>"   %(rec["Allergen Name"])
    x += "<td valign='top'>%s</td>"   %(rec["Source db"])
    x += "<td valign='top'><a href='%s'>%s</a></td>"      %(rec["Hyperlink"],rec["Accession id"])
    x += "<td valign='top'>%s</td>"   %(rec["Species name"])
    x += "<td valign='top'>%s</td>"   %(rec["English name"])
    x += "<td valign='top'>%s</td>"   %(rec["Remark"])
    x += "<td valign='top'>%s</td>"   %(rec["Size mature protein"])
    if showSeq:
      x += "</td></tr><tr bgcolor='%s'><td></td><td colspan='7' valign='top'>" % (bgcolor)
      x += "<pre>%s</pre></td></tr>\n"   %(formatSeq(rec["Sequence"]))
    x += "</tr>"
    return x


def getRandomFileName():
    return os.path.join('/tmp/', str(time.time())+'.'+str(random.randint(0,99999))+'.aa')

def saveSeqToRandomFile(s):
    rf = getRandomFileName()
    O = open(rf, 'w')
    O.write(">query\n")
    O.write(string.strip(s))
    O.close()
    return rf

def formatAlignment(ali,cod, startlines = 2):
    ls = ali.split("\n")
    r = ''
    for l in ls[:startlines]: r += l.strip() + "\n"
    for l in ls[startlines:]:
      if l[:6] == cod[:6]:
        l = '%-10s%s' % (cod[:10], l[6:])
      else:
        l = '%-10s%s' % (l[:6],l[6:])
    r += l + "\n"
    return r

def maskIt(start,end,mask):
    for x in range(start,end):
      mask[x] = mask[x]+1
    return mask

#def send_header():
 #   req.content_type = "text/html"
  #  req.send_http_header()


###################################################################################################
# List all proteins in the database ###DOES NOT WORK <FIX MARTIJN> possibly link to DB missing
###################################################################################################
def list(req, database):
    db = database_loaded[database]
    keys = db.keys()
    keys.sort()
    rv = """
        %s
        <table cellspacing='0' cellpadding = '3'>
           <tr bgcolor='#dfdeda'>
              <td valign='top'> <b> Allergen <br> Id </b> </td>
              <td valign='top'> <b> Allermatch <br> Database  </b> </td>
              <td valign='top'> <b> Allergen <br> Name   </b> </td>
              <td valign='top'> <b> Source <br> Database  </b> </td>
              <td valign='top'> <b> Source <br> Accession Id    </b> </td>
              <td valign='top'> <b> Scientific <br> Name </b> </td>
              <td valign='top'> <b> Common <br> Name </b> </td>
              <td valign='top'> <b> Remark </b> </td>
              <td valign='top'> <b> Mature Protein <br> Length </b> </td>
           </tr>
           <tr align='right' bgcolor='#dfdeda'><td><small>*1</small></td>
              <td> <small> *2</small></td>
              <td> <small> *3</small></td>
              <td> <small> *4</small></td>
              <td> <small> *5</small></td>
              <td> <small> *6</small></td>
              <td> <small> *7</small></td>
              <td> <small> *8</small></td>
              <td> <small> *9</small></td>
           </tr>
        %s
       </table>
       <small><ol>
          <li> Allermatch<sup>tm</sup> id
          <li> Allermatch<sup>tm</sup> database
          <li> Allergen name
          <li> Source database name
          <li> External accession identifier with a link to the external database
          <li> The Latin species name
          <li> Species name in English
          <li> Remark, for example which signalling proteins have been removed
          <li> Lenght of the mature protein (as present in the database)
        %s
        """ % (
            header("Allermatch allergen finder: List all database sequences"),
            "\n".join(map(lambda x: formatDetailsRow(db[keys[x]],x), range(len(keys)))),
            footer())
#    for k in kl:
#        rv += "<li>%(k)s" % vars()
    return rv


###################################################################################################
#start of the full fasta alignemnt
###################################################################################################
def CalcFull(seq, database):
    db = database_loaded[database]
    rv = []
    #perform the full alignment... easiest.
    tmpFile = saveSeqToRandomFile(seq)
    commandline = "%s %s %s %s" % (fastaEx, fastaPmFll, tmpFile, database_files[database])
    fastaPipe = os.popen(commandline)
    fastaRes = fastaPipe.read()
    try:
        fastaPipe.close()
    except:
        #this does go wrong at times.. try to ignore now
        pass

    rv.append("<h4>output of the fasta analysis:\n</h4><pre>")
    sumModus = 0; allmodus = 0; firstalignment = 0
    curAlign = ''; curId = ''
    for line in fastaRes.split("\n"):
        el = line.strip()
        if el[:1] == '#': continue
        if el[:20] == 'The best scores are:':
            rv.append(line); sumModus = 1
        elif sumModus and not el:
            rv.append(line); sumModus = 0
        elif sumModus:
            eid = el.split()[0]
            pxa = '<a href="#%s">' % (eid)
            try:
                pxb = el[:el.index('http')]
            except:
                pxb = el[:50]
            pxc = '</a>'
            rv.append("%s%-50s%s%s" % (pxa,pxb,pxc, el[50:]))
        elif el[0:2] == '>>' and not allmodus:
            curId = el[2:].split()[0]
            allmodus = 1; firstalignment = 1
            rv.append("<table cellspacing='0' width='100%' cellpadding='5'>")
        elif allmodus and el.find('residues') > -1 and el.find('query') > -1:
            #end of alignment
            allmodus = 0
            rv.append(formatAlignment(curAlign, curId))
            rv.append('</pre></td></tr></table><pre>')
        elif allmodus and el[0:2] != '>>':
            curAlign += line + "\n"
        elif el[0:2] != '>>':
            rv.append( line)

        if allmodus and el[0:2]=='>>':
            if not firstalignment:
                #close previous cell
                rv.append( formatAlignment(curAlign, curId))
                rv.append( "</pre></td></tr>")
            firstalignment = 0
            curAlign = ''
            curId = el[2:].split()[0]
            rv.append( formatDetailsRow(db[curId],0, showSeq=0, link=1, bgcolor='dfdeda', database=database) )
            rv.append( "<tr><td colspan='7'><pre>" )
    rv.append("</pre></br>")
    return rv

###################################################################################################
#start of sliding 80aa window method
###################################################################################################
def CalcWindow(seq, database, cutOff, rawOutput):
    #print(seq)
    db = database_loaded[database]
    _database_list = {}
    rv = []
    # do an initial Fasta... To display global alignment score
    globalIdent = {}
    tmpFile = saveSeqToRandomFile(seq)
    commandline = "%s %s %s %s" % (fastaEx, fastaPmAll, tmpFile, database_files[database])
    #print(commandline)
    fastaPipe = os.popen(commandline)
    fastaRes = fastaPipe.read()
    O = open("/home/name/Documents/verwijder.txt", 'w')
    O.write(fastaRes)
    #print(fastaRes)
    #print("______")
    try:
        fastaPipe.close()
    except:
        #gaat vaak mis.. ignore dus maar, wie weet werkt dat
        pass
    for align in fastaRes.split('>>')[2:]:
        ident, idlen = re.search(
            '([0-9\.]+)\% *identity \(.+?\) in (\d*) aa overlap', align).groups()
        #print(align)
        #print(ident, idlen)
        name = align[:align.find("\n")].split()[0]
        globalIdent[name] = (float(ident), int(idlen))

    noWindows = 0
    hitList  = {}
    for pos in range(0, max(1,len(seq)-79)):
        noWindows += 1
        lseq = seq[pos:pos+80]
        #store seq
        tmpFile = saveSeqToRandomFile(lseq)

        commandline = "%s %s %s %s" % (fastaEx, fastaPmWin, tmpFile, database_files[database])
        fastaPipe = os.popen(commandline)
        fastaRes = fastaPipe.read()
        try:
            fastaPipe.close()
        except:
            #this does go wrong at times.. try to ignore now
            pass

        for align in fastaRes.split('>>')[2:]: #first two are nonsense
            #get identity
            name = align[:align.find("\n")].split()[0]
            identity, idlen = re.search(' ([0-9\.]+)\% *identity \(.+?\) in (\d*) aa overlap',
                                        align).groups()
            identity = float(identity)
            idlen = float(idlen)

            ###Correct identity
            if idlen < 80:
                identity = identity * (idlen/80.0)
            if identity < cutOff: continue

            #add to hitlist
            idList = hitList.get(name, [])
            idList.append(identity)
            hitList[name] = idList

    #ok, now do the output
    #Sort on quality of output:
    sortList = []
    for key in hitList.keys():
        scores = hitList[key]
        sortList.append((max(scores), len(scores), key))
    sortList.sort()
    sortList.reverse()
    prc = '%'
    if not rawOutput: rv.append("""
        <table cellpadding='5' cellspacing='0' width='100 %(prc)s'>
           <tr bgcolor='#ababa9'>
              <td valign='top'> <b> Hit <br> No </b> </td>
              <td valign='top'> <b>Db </b> </td>
              <td valign='top'> <b> Description </b> </td>
              <td valign='top'> <b> Best hit <br>  Identity  </b> </td>
              <td valign='top'> <b> No of hits <br> ident &gt; %(cutOff)2.2f </b> </td>
              <td valign='top'> <b> %(prc)s of hits <br> ident &gt; %(cutOff)2.2f </b> </td>
              <td valign='top'> <b> Full <br> Identity </b> </td>
              <td valign='top'> <b> External <br>link </b> </td>
              <td valign='top'> <b> Scientific Name </b> </td>
              <td valign='top'> <b> Detailed <br> Information </b></td>
           </tr>
           <tr align='right' bgcolor='#ababa9'><td><small>*1</small></td>
              <td> <small> *2</small></td>
              <td> <small> *3</small></td>
              <td> <small> *4</small></td>
              <td> <small> *5</small></td>
              <td> <small> *6</small></td>
              <td> <small> *7</small></td>
              <td> <small> *8</small></td>
              <td> <small> *9</small></td>
              <td> <small> *10</small></td>
           </tr>""" % vars() )
    itemNo = 0
    for item in sortList:
        itemNo+=1
        key = string.strip(item[2])

        if itemNo % 2 == 0:
            bgcolor = "#dfdeda"
        else:
            bgcolor = "#FFFFFF"

        simiBestHit = item[0]
        noHits      = item[1]
        percHits    = float(item[1])/noWindows * 100
        globId      = globalIdent[key][0]
        globIdLen   = globalIdent[key][1]
        try:
            swpLink     = db[key]["Hyperlink"]
        except:
            raise key + "\n\n" + str(db) + "\n\n".join(db.keys())
        AccId       = db[key]["Accession id"]
        specName    = db[key]["Species name"]
	description	= db[key]["Remark"]

        _si = '?'
        _seqDb = db[key]["Database Name"]
        if _seqDb == 'WHO-IUIS Allergen': _si = 'WA'
	elif _seqDb == 'AllergenDB': _si = 'AL'
	elif _seqDb == 'UniProt': _si = 'UniProt'
	elif _seqDb == 'GenBank': _si = 'GenBank'
        else: _si = '?'

        _seqSr = db[key]["Source db"]
        _ss = '?'
        if _seqSr == 'UniProt': _ss = 'U'
	elif _seqSr == 'GenBank': _ss = 'G'

        _script     = CGI_SCRIPT

        if not rawOutput: rv.append("""
            <tr bgcolor='%(bgcolor)s'>
              <td valign='top'> %(itemNo)d </td>
              <td valign='top'> %(_si)s </td>
              <td valign='top'> %(description)s </td>
              <td valign='top'> %(simiBestHit)2.2f </td>
              <td valign='top'> %(noHits)d </td>
              <td valign='top'> %(percHits)2.2f </td>
              <td valign='top'> %(globId)2.2f / %(globIdLen)d </td>
              <td valign='top'> <a href='%(swpLink)s'> %(AccId)s<sup>%(_ss)s</sup></small> </a></td>
              <td valign='top'> %(specName)s </td>
              <td valign='top' align='middle'>
                 <form action='/allermatchsearch/search' method='POST'
                       enctype='application/x-www-form-urlencoded'>
                   <input name='Go' value='Go' type='submit'>
                   <input type='hidden' name='database' value='%(database)s'>
                   <input type='hidden' name='against' value='%(key)s'>
                   <input type='hidden' name='cutOff' value='%(cutOff)s'>
                   <input type='hidden' name='wordlength' value='6'>
                   <input type='hidden' name='allAlignments' value='0'>
                   <input type='hidden' name='method' value='windowSingle'>
                   <input type='hidden' name='seq' value='%(seq)s'>
                </form>
              </td>
            </tr> """ % locals())
        else:
            xpr = {}
            xpr['itemNo'] = itemNo
            xpr['_si'] = _si
            xpr['key'] = key
            xpr['simiBestHit'] = simiBestHit
            xpr['noHits'] = noHits
            xpr['percHits'] = percHits
            xpr['globId'] = globId
            xpr['swpLink'] = swpLink
            xpr['AccId'] = AccId
            xpr['_ss'] = _ss
            xpr['specName'] = specName
            rv.append("\t".join(map(lambda X: "%s:%s" % X, xpr.items())))

    if not rawOutput: rv.append("""</table>
       <br> Analyzed %(noWindows)d windows
       <small><ol>
       <li> Number of the hit, best hit comes first
       <li> External Database:
          <ul>
	    <li>UniProt : UniProt Protein Knowledgebase
            <li>GenBank : GenBank NCBI (RefSeqProtein)
          </ul>
       <li> Description of the sequence
       <li> Identity of the best hit (percent identical amino acids in the
            aligned 80-amino-acid sliding window)
       <li> The number of hits the input sequence had with this allergen
       <li> The percentage of windows analysed for this input sequence hitting
            this allergen
       <li> Results of a fasta alignment of the complete input sequence against
            this database sequence. The first number is the percentage of
            identity. The second number is the length of sequence over which fasta aligned
       <li> External database accession id linking to this database, the superscript ids
            indicate which database this is:
              <ul>
              <li> U : UniProt
              <li> G : GenBank NCBI (RefSeqProtein)
              </ul>
       <li> Scientific name of the species.
       <li> Links to a page with specific details on this database sequence, the
            complete fasta alignment and the part of the input sequence aligning
            to the database sequence.
       """ % locals())
    #print(rv)
    return rv


###################################################################################################
#Start of wordmatch method
###################################################################################################
def CalcWordmatch(seq, database, wordlength, rawOutput):
    db = database_loaded[database]
    rv = []
    #for each 6-mer in the incoming sequence
    hits = {}
    nowindows = 0
    for pos in range(0,len(seq)-wordlength+1):
        nowindows += 1
        word = string.lower(seq[pos:pos+wordlength])
        for k in db.keys():
            sq = db[k]['Sequence']
            c = sq.count(word)
            hits[k] = hits.get(k,0)+ min(c,1)
    if not rawOutput: rv.append("""<table width='100%' cellpadding='3' cellspacing='0'>
          <tr bgcolor='#ababa9'>
            <td valign='top'><b>No</b></td>
            <td valign='top'><b>Db</b> </td>
          <td valign='top'><b>Description</b></td>
          <td valign='top'><b>No of exact<br>wordmatches</b></td>
          <td valign='top'><b>% of exact<br>wordmatches</b></td>
          <td valign='top'><b>External<br>db</b></td>
          <td valign='top'><b>Scientific Name</b></td>
          <td valign='top'><b>Detailed<br>Information</b></td>
          </tr>
          <tr align='right' bgcolor='#ababa9'><td><small>*1</small></td>
            <td><small>*2</small></td>
            <td><small>*3</small></td>
            <td><small>*4</small></td>
            <td><small>*5</small></td>
            <td><small>*6</small></td>
            <td><small>*7</small></td>
            <td><small>*8</small></td>
          </tr> """)

    hitList = []
    for k in  hits.keys():
        hitList.append((hits[k],k))
        hitList.sort()
        hitList.reverse()

    itemNo = 0
    for hitem in hitList:
        itemNo += 1
        if itemNo % 2 == 0: bgcolor = "#dfdeda"
        else: bgcolor = "#FFFFFF"
        hits,key = hitem
        if hits>0:
            _script = CGI_SCRIPT #workaround to ge this variable local
            link="""
              <form action='/allermatchsearch/search' method='POST'>
                 <input name='Go' value='Go' type='submit'>
                 <input type='hidden' name='against' value='%(key)s'>
                 <input type='hidden' name='method' value='wordmatchSingle'>
                 <input type='hidden' name='seq' value='%(seq)s'>
                 <input type='hidden' name='database' value='%(database)s'>
                 <input type='hidden' name='cutOff' value='35'>
                 <input type='hidden' name='wordlength' value='%(wordlength)d'>
              </form>""" % vars()
            Remark =   db[key]['Remark'] ##Remark
            PercHit =   float(hits) / float(len(seq) - wordlength + 1) * 100
            hyperlink = db[key]["Hyperlink"]
            swissacc =  db[key]["Accession id"]
            specName =  db[key]['Species name']

            _seqDb = db[key]["Database Name"]
            if _seqDb == 'WHO-IUIS Allergen': _si = 'WA'
            elif _seqDb == 'AllergenDB': _si = 'AL'
	    elif _seqDb == 'UniProt': _si = 'UniProt'
	    elif _seqDb == 'GenBank': _si = 'GenBank'
            else: _si = '?'
            _seqSr = db[key]["Source db"]
            if _seqSr == 'UniProt': _ss = 'U'
            elif _seqSr == 'GenBank': _ss = 'G'
            else: _ss = '?'


            if not rawOutput: rv.append("""
               <tr bgcolor='%(bgcolor)s'>
                 <td valign='top'> %(itemNo)d </td>
                 <td valign='top'> %(_si)s </td>
                 <td valign='top'> %(Remark)s </td>
                 <td valign='top'> %(hits)s </td>
                 <td valign='top'> %(PercHit)2.2f </td>
                 <td valign='top'><a href='%(hyperlink)s'>
                    %(swissacc)s<sup>%(_ss)s</sup> </a></td>
                 <td valign='top'> %(specName)s </td>
                 <td  valign='top' align='middle'> %(link)s </td></tr>""" % vars())
            else:
                rv.append("\t".join(map(lambda X: "%s:%s" % X, db[key].items())))

    if not rawOutput: rv.append("""</table>
        <br> Analyzed %(nowindows)s windows
        <small><ol>
          <li> Number of the hit, best hit comes first
          <li> External Database:
          <ul>
	    <li>UniProt : UniProt Protein Knowledgebase
            <li>GenBank : GenBank NCBI (RefSeqProtein)
          </ul>
          <li> Description of the sequence
          <li> The number of exact %(wordlength)d aa hits the input
               sequence had with this allergen
          <li> The Percentage of exact hits the input sequence is found
               to hit this allergen sequence
          <li> External database accession id linking to this database, the superscript ids
            indicate which database this is:
              <ul>
              <li> U : UniProt
              <li> G : GenBank NCBI (RefSeqProtein)
              </ul>
          <li> Species name of the allergen
          <li> Links to a page with specific details on this database
               sequence and the part of the input sequence aligning to the
               database sequence. """ % vars())
    return rv


###################################################################################
#Start of wordmatch method against a single sequence
###################################################################################
def CalcWordmatchSingle(seq, database, wordlength, against):
    db = database_loaded[database]
    rv = []
    sq = db[against]["Sequence"]
    rv.append(formatDetails(db[against]))

    sm = ' ' * len(sq)
    sm2 = [0 for x in range(len(sq))]
    sqm = ' ' * len(seq)
    sqm2 = [0 for x in range(len(seq))]

    #recalculate results
    for pos in range(0,len(seq)-wordlength+1):
        word = string.lower(seq[pos:pos+wordlength])

        if sq.count(word) > 0:
            sqm = sqm[:pos] + "#"*(len(word)) + sqm[pos+wordlength:]
            sqm2 = maskIt(pos, pos+len(word), sqm2)
            #mask also the against sequence
            p = -1
            while 1:
                p = sq.find(word,p+1)
                if p > -1:
                    sm = sm[:p] + "#"*(len(word))  + sm[p+len(word):]
                    sm2 = maskIt(p, p+len(word), sm2)
                else: break


    rv.append("""
         <h4>Input Sequence</h4><small>
         The '#' denotes which parts of the
         input sequence have an exact %(wordlength)d aa match against this
         specific database sequence</small><pre>""" %vars() )

    s = seq; m = sqm; m2=''.join(map(str,sqm2)).replace('0',' '); pos = 0
    while len(s) > 0:
        rv.append("%5d : %s\n%5d : %s\n" % (
            pos, s[:60], pos, m[:60]))
        s = s[60:]; m = m[60:]; pos += 60; m2 = m2[60:]

    rv.append("""
       </pre>
       <h4>Database Sequence</h4><small>The '#' denotes which parts of the
       database sequence have an exact %(wordlength)d aa match against the
       input sequence</small><pre>""" % vars() )

    s = sq; m = sm; m2=''.join(map(str,sm2)).replace('0',' ');
    pos = 0
    while len(s) > 0:
        rv.append("%5d : %s\n%5d : %s\n" % (
            pos, s[:60], pos, m[:60]))
        s = s[60:]; m = m[60:]; pos += 60; m2=m2[60:]
    print "</pre>"

    return rv

###################################################################################################
#Start of window method against a single sequence
###################################################################################################
def CalcWindowSingle(seq, database, cutOff, against, allAlignments):
    db = database_loaded[database]
    rv = []
    #perform a global alignment and display the results
    if allAlignments:
        nextAlign = 0
        nextAlignText = 'Hide all alignments'
    else:
        nextAlign = 1
        nextAlignText = 'Show all alignments'

    rv.append("""<form action='/allermatchsearch/search' method='POST'
            enctype='application/x-www-form-urlencoded'>
         <input name='Go' value='%(nextAlignText)s' type='submit'>
         <input type='hidden' name='against' value='%(against)s'>
         <input type='hidden' name='cutOff' value='%(cutOff)s'>
         <input type='hidden' name='wordlength' value='6'>
         <input type='hidden' name='database' value='%(database)s'>
         <input type='hidden' name='allAlignments' value='%(nextAlign)s'>
         <input type='hidden' name='method' value='windowSingle'>
         <input type='hidden' name='seq' value='%(seq)s'>
      </form> """ % vars() )

    tmpFile = saveSeqToRandomFile(seq)
    commandline = "%s %s %s %s" % (fastaEx, fastaPmAll, tmpFile, database_files[database])
    fastaPipe = os.popen(commandline)
    fastaRes = fastaPipe.read()
    try:
        fastaPipe.close()
    except:
        #this does go wrong at times.. try to ignore now
        pass

    for hit in fastaRes.split('>>')[2:]:
        if hit.find(against) > -1:
            globalAlign = hit

    rv.append(formatDetails(db[against]))

    aga = db[against]['Sequence']
    seqmask = ' ' * len(seq)
    agamask = ' ' * len(aga)
    allAlign = []

    noWindows = 0
    for pos in range(0, max(1,len(seq)-79)):
        noWindows += 1
        lseq = seq[pos:pos+80]
        #print "<li>"+lseq
        #store seq
        tmpFile = saveSeqToRandomFile(lseq)

        commandline = "%s %s %s %s" % (fastaEx, fastaPmWin, tmpFile, database_files[database])
        fastaPipe = os.popen(commandline)
        fastaRes = fastaPipe.read()
        try:
            fastaPipe.close()
        except:
            #this does go wrong at times.. try to ignore now
            pass


        for align in fastaRes.split('>>')[2:]: #first two are nonsense
            name = align[:align.find("\n")].split()[0]
            if  name != against:
                continue
            identity, idlen = re.search(
                ' ([0-9\.]+)\% *identity \(.+?\) in (\d*) aa overlap',
                align).groups()
            identity = float(identity)
            idlen = int(idlen)

            ###Correct identity
            if idlen < 80:
                identity = identity * (idlen/80.0)
            if identity < cutOff:
                continue

            if allAlignments:
                a = align.split("\n")[1:]
                a[0] = string.strip(a[0])
                allAlign.append("corrected identity %f \n%s" %( identity, "\n".join(a)))

           # print "</pre><br>accepted<br>"
            qal = ''
            dal = ''
            Aal = align.split("\n")[4:]
            doff = 0; doffSet = 0
            for lno in range(len(Aal)):
                l = Aal[lno]
                if l[:5] == 'query': qal += l[7:]
                if l[:5] == against[:5]:
                    dal += l[7:]
                    if not doffSet:
                        doffSet = 1
                        spaces,firstcoord = re.search("^(\s+)(\d+)\s",Aal[lno+1]).groups()
                        doff = -1 * ( len(spaces) - 7 - int(firstcoord) + len(firstcoord) )
            dal = dal.rstrip()
            qal = qal.rstrip()
            sdal = dal.replace('-','')
            agaStart = max(0,len(qal) - len(qal.lstrip()) + doff) #(off by one???)
            agaEnd   = doff + len(sdal) - len(dal) + len(qal)

            #mask incoming sequence:
            seqmask = seqmask[:pos] + '#'*80 + seqmask[pos+80:]
            agamask = agamask[:agaStart] + '#'*(agaEnd-agaStart) + agamask[agaEnd:]
            break

    #start creating output
    rv.append("""
       <h4>Input Sequence</h4><small>The '#' denotes which parts of the
       input sequence having an 80 aa hit >= %(cutOff)2.2f%% identity against this
       specific database sequence</small>
       <pre>""" % vars() )
    s = seq; p = 0; m = seqmask[:len(s)]
    while len(s) > 0:
        rv.append("%5d : %s\n%5d : %s\n" % (p, s[:60],p, m[:60]))
        s = s[60:]; m = m[60:]
        p += 60

    rv.append(""" </pre><h4>Database Sequence</h4><small>
        The '#' denotes which parts of the
        database sequence having an 80 aa hit >= %(cutOff)2.2f%%
        identity  against the input sequence</small>
        <pre>""" % vars())

    s = aga; p = 0; m = agamask[:len(s)]
    while len(s) > 0:
        rv.append("%5d : %s\n%5d : %s\n" % (p, s[:60],p, m[:60]))
        m = m[60:]; s = s[60:]; p += 60
    rv.append("""
       </pre>
       <h4>Full Alignment</h4>
       <pre>%s</pre>""" % (formatAlignment(globalAlign, against, 3)))
    if allAlignments:
        formAli = "</pre>\n<hr>\n<pre>".join(
            map(lambda x: formatAlignment(x, against, 3), allAlign))

        rv.append("""<h4>All 80 aa window alignments with an identity
             above a cut-off value of %(cutOff)2.2f%% </h4>
             <pre> %(formAli)s </pre> """ % vars())
    return rv



