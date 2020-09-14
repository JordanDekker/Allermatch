
import string,re,sys,os

#define variables

URL = '137.224.17.7'
#URL = 'http://www.allermatch.org/'
LOC = '/home/name/Downloads/allermatch'

CGI_SCRIPT = URL + '/allermatchsearch'

header_file = LOC + '/htdocs/templates/header.html'
header_file_form = LOC + '/htdocs/templates/header2.html'
footer_file = LOC + '/htdocs/templates/footer.html'
footer_file_form = LOC + '/htdocs/templates/footer2.html'

header_text = open(header_file).read()
header_text_form = open(header_file_form).read()
footer_text = open(footer_file).read()
footer_text_form = open(footer_file_form).read()

TITLE = 'Allermatch'
SCRIPT = 'search'
verbose = 0
fastaEx = '/opt/fasta36-tool/bin/fasta36'
fastaPmAll = ' -Q -E 1000 -d 1000 '
fastaPmWin = ' -Q -H -w 200 -E 10 '
fastaPmFll = ' -Q -H -d 200 -E 100 '

database_files = {
    'AllergenDB_original_sequences'             : LOC + '/sof/AllergenDB_original_sequences.fasta',
    'AllergenDB_propeptides_removed'             : LOC + '/sof/AllergenDB_propeptides_removed.fasta'}

#standard functions
def recify(l,s):
    return {
        "Allergen id"          : l[0], "Database Name"       : l[1],
        "Allergen Name"	       : l[2], "Source db"           : l[3],
        "Accession id"         : l[4], "Hyperlink"           : l[5],
        "Species name"         : l[6], "English name"        : l[7],
        "Remark"               : l[8], "Size mature protein" : l[9],
        "Sequence"             : s }

def importDatabase(db):
    _db = {}
    F = open(db,'r')
    lines = map(string.strip, F.readlines())
    F.close()

    head = ''; seq = ''; firstSeq = 1
    for line in lines:
        if not line: continue
        if (line[0] == '>') and (firstSeq):
            head = line[1:].split("\t"); seq = ''
            firstSeq = 0
        elif line[0] != '>':
            seq += line
        elif (line[0] == '>') and (not firstSeq):
            _db[head[0]] = recify(head,seq)
            head = line[1:].split("\t"); seq = ''

    #last record
    _db[head[0]] = recify(head,seq)

    return _db

#open and load database
database_loaded = {}
for k in database_files.keys():
    database_loaded[k] = importDatabase(database_files[k])
database_default = 'AllergenDB_propeptides_removed'


def header(form = False):
    if form:
        return header_text_form
    else:
        return header_text

def footer(form = False):
    if form:
        return footer_text_form
    else:
        return footer_text

def rawHeader():
    return ""

def rawFooter():
    return ""
