from Bio import SeqIO
from allermatch_settings import *
import  random, time, os, sys, string, re, copy, itertools



def extract_sequences_from_fasta(seq):
    """
    When a fasta-file was used as query this function extracts the sequences and the id's
    from the file.

    :param seq: A string containing sequences and id's in fasta-format

    :return sequences: A list contatining all the dna sequences from the fasta-file
    :return identifiers: A list containing all the id's from the fasta-file
    """
    fasta_sequences = SeqIO.parse(seq,'fasta')
    sequences = []
    identifiers = []
    for fasta in fasta_sequences:
        #name, sequence = fasta.id, str(fasta.seq)
        sequences.append(str(fasta.seq))
        identifiers.append(str(fasta.id))
    return sequences[1:], identifiers[1:]


def saveSeqToRandomFile2(sequence, id):
    rf = getRandomFileName()
    O = open(rf, 'w')
    O.write(">"+id+"\n")
    O.write(string.strip(sequence))
    O.close()
    return rf

def saveSeqToRandomFile3(sequences ,ids):
    """
    Saves the sequences and the id's to a tmp file. This file will be used for
    FASTA36

    :param sequences: A list contatining all the dna sequences from the fasta-file
    :param identifiers: A list containing all the id's from the fasta-file

    :returns rf: File path to the newly created file
    """
    rf = getRandomFileName()
    O = open(rf, 'w')
    for (s, i)  in zip(sequences, ids):
        O.write(">"+i+"\n")
        O.write(string.strip(s))
        O.write("\n")
    O.close()
    return rf


def getRandomFileName():
    """
    Creates a random filepath in the tmp folder:

    returns: A file path
    """
    return os.path.join('/tmp/', str(time.time())+'.'+str(random.randint(0,99999))+'.aa')


def formatSeq(sq):
    r = ''; s = sq
    while len(s) > 0:
      r += s[:60] + "\n"
      s = s[60:]
    return r

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

def calc_window_multiple_sequences(sequences, ids, database, cutOff, rawOutput):
    """
    Calculates possible allergen hits with Fasta36 of multiple query sequences 
    using the sliding window approach.

    :param sequences: A list containing the query sequences
    :param ids: A list containing the query ID's
    :param database: A string for which database to use for the library sequences
    :param cutOff: A integer for the minimal % identity when considered to be a hit
    :param rawOutput: A integer

    :return rv: A list containing HTML strings for generating a HTML page
    """
    db = database_loaded[database]
    _database_list = {}
    rv = []
    # do an initial Fasta... To display global alignment score
    globalIdent = {}
    for (seq, i) in zip(sequences, ids):
        tmpFile = saveSeqToRandomFile2(seq, i)
        commandline = "%s %s %s %s" % (fastaEx, fastaPmAll, tmpFile, database_files[database])
        fastaPipe = os.popen(commandline)
        fastaRes = fastaPipe.read()
        try:
            fastaPipe.close()
        except:
            #gaat vaak mis.. ignore dus maar, wie weet werkt dat
            pass
        for align in fastaRes.split('>>')[2:]:
            ident, idlen = re.search(
                '([0-9\.]+)\% *identity \(.+?\) in (\d*) aa overlap', align).groups()
            name = align[:align.find("\n")].split()[0]
            globalIdent[name] = (float(ident), int(idlen))

        noWindows = 0
        hitList  = {}
        for pos in range(0, max(1,len(seq)-79)):
            noWindows += 1
            lseq = seq[pos:pos+80]
            #store seq
            tmpFile = saveSeqToRandomFile2(lseq, i)

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
            <button type='button' class='collapsible'>>%(i)s</button>
            <div class='content'>
            <p>
            <table cellpadding='5' cellspacing='0' width='100%(prc)s'>
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
        rv.append("</table></p></div><br><br>")
    if not rawOutput: rv.append("""
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
    return rv


def calc_wordmatch_multiple_sequences(sequences, ids, database, wordlength, rawOutput):
    """
    Calculates possible allergen hits with Fasta36 of multiple query sequences 
    using exact wordmatch method.

    :param sequences: A list containing the query sequences
    :param ids: A list containing the query ID's
    :param database: A string for which database to use for the library sequences
    :param wordlength: A integer for the minimal exact length a sequence need to be a hit
    :param rawOutput: A integer

    :return rv: A list containing HTML strings for generating a HTML page
    """
    db = database_loaded[database]
    rv = []
    for (seq, ID) in zip(sequences, ids):
        rv.append("<button type='button' class='collapsible'>>"+ID+"</button>")
        rv.append(" <div class='content'>\n<p>")
        rv.append("<table cellspacing='0' width='100%' cellpadding='0'>")
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
        rv.append("</table></p></div>")
        rv.append("<br><br>")

    if not rawOutput: rv.append("""
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


def calc_full_multiple_sequences(seq, ids, database):
    """
    Calculates possible allergen hits of multiple query sequences 
    with a FASTA full alignment.

    :param sequences: A list containing the query sequences
    :param ids: A list containing the query ID's
    :param database: A string for which database to use for the library sequences

    :return rv: A list containing HTML strings for generating a HTML page
    """
    db = database_loaded[database]
    rv = []
    tmpFile = saveSeqToRandomFile3(seq, ids)
    commandline = "%s %s %s %s" % (fastaEx, fastaPmFll, tmpFile, database_files[database])
    fastaPipe = os.popen(commandline)
    fastaRes = fastaPipe.read()
    try:
        fastaPipe.close()
    except:
        #this does go wrong at times.. try to ignore now
        pass
    rv.append("<h4>output of the fasta analysis:\n</h4><pre>")
    test = fastaRes.split(">>>")
    for query in test[1:]:
        query_id = query.split("\n")[0].split("-")[0]
        query_length = query.split("\n")[0].split("-")[1]
        rv.append("<button type='button' class='collapsible'>>"+query_id+"  "+query_length+"</button>")
        rv.append(" <div class='content'>\n<p>")
        rv.append("<table cellspacing='0' width='100%' cellpadding='5'>")
        stats = query.split(">>")[0].split("\n\n")[1]
        rv.append(stats)
        best_scores = query.split(">>")[0].split("\n\n")[2]
        rv.append(best_scores.split("\n")[0])
        for line in best_scores.split("\n")[1:]:
            eid = line.split()[0]
            pxa = '<a href="#%s">' % (eid)
            try:
                pxb = line[:line.index('http')]
            except:
                pxb = line[:50]
            pxc = '</a>'
            rv.append("%s%-50s%s%s" % (pxa,pxb,pxc, line[50:]))
        
        rv.append("<br>")

        for match in query.split('>>')[1:]:
            lines = match.split("\n")
            testt = lines[1:]
            curId = lines[0].split()[0]
            rv.append( formatDetailsRow(db[curId],0, showSeq=0, link=1, bgcolor='dfdeda', database=database) )
            rv.append(" <tr><td colspan='7'><pre>")
            rv.append(testt[0])
            rv.append(testt[1])
            rv.append("</pre></td></tr>")

        rv.append("</table>\n</p>\n</div>")

    return rv