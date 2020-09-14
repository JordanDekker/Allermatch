from flask import Flask, render_template, request
import sys
from allermatch_settings import *
from allermatch_functions import *
from allermatch_functions2 import *
#sys.path.insert(0, "/var/www/html/allermatch/htdocs")

app = Flask(__name__)

@app.route("/")
def home():
    return render_template("index.html")

@app.route("/example.html")
def example():
    return render_template("example.html")

@app.route("/index.html")
def index():
    return render_template("index.html")

@app.route("/database.html")
def database():
    return render_template("database.html")

@app.route("/help.html")
def help():
    return render_template("help.html")

@app.route("/publication.html")
def publication():
    return render_template("publication.html")

@app.route("/about.html")
def about():
    return render_template("about.html")

@app.route("/feedback.html")
def feedback():
    return render_template("feedback.html")

@app.route("/disclaimer.html")
def disclaimer():
    return render_template("disclaimer.html")

@app.route("/copyright.html")
def copyright():
    return render_template("copyright.html")

@app.route("/thanks.html")
def thanks():
    return render_template("thanks.html")

@app.route("/references.html")
def references():
    return render_template("references.html")

@app.route('/allermatchsearch/form')
def form():
    #send_header(req)
    rv = header()
    rv += """
      <h2>Allermatch allergen finder: Input Form</h2>
       This webpage features the following three ways of analysis to identify a relationship
       between your input sequence and an allergen from the database:
	<ul>
         <li><b>80-amino-acid sliding window:</b> The input sequence is
             chopped up in 80-amino-acid windows. For each 80-
             amino acid window, the program counts which allergen it hits
             (with a specific identity).
         <li><b>Full Alignment:</b> Use FASTA to perform a full alignment.
         <li><b>Wordmatch:</b> Look for an exact hit of 6 or more contiguous amino acids in a sequence in the database.
	</ul>

	As explained under <a href="../database.html">Databases</a>, comparisons can be run against either of the following two databases:
	<ol>
	<li>The <b>AllergenDB_propeptides_removed</b> databases contains allergen sequences from which the signal- and propeptide sequences are removed when post-translational modifications (PTMs) were predicted by UniProtKB.
	<li>The <b>AllergenDB_original_sequences</b> database contains the non-processed allergen sequences with PTMs.
	</ol>

<br></br>

         <table border='1' style='border-style: solid; border-width: 1;'
                            cellspacing='0' cellpadding='0'><tr><td>
             <form name='allerform'  method='POST' action = '/allermatchsearch/search'
                   enctype='application/x-www-form-urlencoded'>
             <input type='hidden' name='against' value=''>
             <table celspacing='0' cellpadding='5' width='100%'>
                  <tr><td colspan='2'><b>
                        Copy-paste your amino acid sequence here:</b>
                    </td> </tr>
                  <tr><td colspan='2'>
                      <textarea name='seq' rows='4' cols='80'></textarea>
                    </td> </tr>
                    <tr><td colspan='2'><b>Algorithm:</b></td></tr>
                  <tr><td colspan='1'>
                     <input type='radio' name='method' value='window' checked>
                      Do an 80-amino-acid sliding window alignment
                    </td> <td colspan='1'>
                      <input type='text' name='cutOff' size='10' maxlength='8' value='35'>
                      Cut-off Percentage (only applicable to the 80-amino-acid sliding window)
                    </td> </tr>
                  <tr><td colspan='1'>
                      <input type='radio' name='method' value='wordmatch'>
                      Look for a small exact wordmatch
                    </td> <td colspan='1'>
                      <input type='text' name='wordlength' size='10' maxlength='8' value='6'>
                      Wordlength (only applicable to the exact wordmatch search)
                    </td> </tr>
                  <tr><td colspan='1'>
                      <input type='radio' name='method' value='full'>
                      Do a full FASTA alignment
                    </td> </tr>
                  <tr><td>
                     <b>Select a database:</b>
                     <td> <select name='database'> """
    for db in ["AllergenDB_propeptides_removed", "AllergenDB_original_sequences"]:
        if db == database_default: selected = 'selected'
        else: selected = ''
        rv +=  "%s<option value='%s' %s>%s</option>\n" % (" " * 20, db,selected,db)

    rv += """             </select>
                  </tr></td>
                  <tr><td colspan='2'>
                      <input style='padding-left: 25px; padding-right: 25px;
                                    padding-top:10px; padding-bottom:10px;'
                             name=' Go ' value='Go' type='submit'>
                    </td> </tr>
             </table>
             </form> <p>
             </td></tr>
          </table> <p>

          """
    rv += footer()
    return rv

@app.route('/allermatchsearch/search',  methods=['POST'])
def search():
    database = request.values.get('database')
    seq = request.values.get('seq')
    method = request.values.get('method')
    cutOff = request.values.get('cutOff')
    wordlength = request.values.get('wordlength')
    against = request.values.get('against')
    wordlength = int(wordlength)
    cutOff = float(cutOff)
    
    if ">" in str(seq):
        fasta_file = seq.split("\n")
        sequences ,ids = extract_sequences_from_fasta(fasta_file)
    else:
        seq = extractSeq(seq)
    Form = request.form

    if Form.has_key('rawOutput'): rawOutput = 1
    else: rawOutput = 0

    if not rawOutput: rv = [header(True),'']
    else: rv = [rawHeader(),'']

    if not rawOutput: rv.append('<h4>Database : %s</h4> ' %(database))

    if method == 'full':
        rv[1]=("<h2>Full Alignment</h2>")
        if "sequences" in locals():
            rv.extend(calc_full_multiple_sequences(sequences, ids, database))
        else:
            rv.extend(CalcFull(seq, database))
    elif method == 'window':
        if not rawOutput:  rv[1]=("<h2>80-amino-acid sliding window</h2>")
        if "sequences" in locals():
            rv.extend(calc_window_multiple_sequences(sequences, ids, database, cutOff, rawOutput))
        else:
            rv.extend(CalcWindow(seq, database, cutOff, rawOutput))
    elif method == 'windowSingle':
        rv[1]=("<h2>80-amino-acid sliding window against %(against)s</h2>" % vars())
        if Form.has_key('allAlignments'): allAlignments = int(Form['allAlignments'])
        else: allAlignments = 0
#        allAlign = str(Form['allAlignments'])
        rv.extend(CalcWindowSingle(seq, database, cutOff, against, allAlignments))
    elif method == 'wordmatch':
        if not rawOutput: rv[1]=("<h2>%(wordlength)d Amino Acid Wordmatch</h2>" % vars())
        if "sequences" in locals():
            rv.extend(calc_wordmatch_multiple_sequences(sequences, ids, database, wordlength, rawOutput))
        else:
            rv.extend(CalcWordmatch(seq, database, wordlength, rawOutput))
    elif method == 'wordmatchSingle':
        rv[1]=("<h2>%(wordlength)d Amino Acid Wordmatch against %(against)s</h2>" % vars())
        rv.extend(CalcWordmatchSingle(seq, database, wordlength, against))

    if not rawOutput: rv.append(footer(True))
    return "\n".join(rv)


if __name__ == '__main__':
        app.run(debug=True)

