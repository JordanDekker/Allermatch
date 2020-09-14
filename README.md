# ALLERMATCH

Allermatch<sup>tm</sup> is a unique webtool where you can compare the amino acid sequence of a protein of interest with sequences of allergenic proteins.
This webtool carries out automatically the procedures for predicting the potential allergenicity of proteins by bioinformatics approaches as recommended by the Codex alimentarius and the FAO/WHO Expert consultation on allergenicity of foods derived through modern biotechnology.
The unique features of the Allermatch<sup>tm</sup> webtool allow the user in a user-friendly and time-saving manner to enter the input sequence and retrieve, with a few mouse-clicks, the outcomes of interest in an accurate, concise, and comprehensible format.

## Installation

This is a complete walkthrough from scratch, feel free to skip parts if they are already fulfilled on your server.

### Install apache webserver and needed modules

```shell
apt install aptitude
aptitude update
aptitude install libapache2-mod-python
a2enmod include
a2enmod python
service apache2 restart
```

### Checkout the code, and place it on your server

```shell
cd /var/www/html
git clone https://git.wur.nl/haars001/allermatch.git
```

### Fix index.html (add header and footer)
```shell
cd /var/www/html/allermatch/htdocs
cat header.html index.html footer.html > index.htm
rm index.html
mv index.htm index.html
```

### Adapt configuration for allermatch.org
#### Python settings
```shell
cd /var/www/html/allermatch/htdocs
sed -i 's|^LOC =.*|LOC = '"'"'/var/www/html/allermatch'"'"'|' allermatch_settings.py
sed -i 's|^URL=.*|URL= '"'"'http://www.allermatch.org'"'"'|' allermatch_settings.py
```
(the quotes are strange to be able to add a single quote to the string)

#### Webserver settings
```shell
cd /etc/apache2/sites-available
cat 000-default.conf
<VirtualHost *:80>
        ServerAdmin webmaster@allermatch.org
        DocumentRoot /var/www/html/allermatch/htdocs
        AddHandler mod_python .py
        PythonHandler mod_python.publisher
        PythonDebug On
        ErrorLog ${APACHE_LOG_DIR}/error.log
        CustomLog ${APACHE_LOG_DIR}/access.log combined
</VirtualHost>
```
As the application consists of two parts, the webserver should not point to the root of the application, but to the _htdocs_ directory.

### Load the updated configuration

After setting it all up, running ```apache2ctl graceful``` should give you a running allermatch website.
## Testing

If you enter _aiscgqvasaiapcisyargqgsgpsagccsgvrslnnaarttadrraacnclknaaagvsglnagnaasipskcgvsipytiststdcsrvn_ as query in the form at http://www.allermatch.org/allermatch.py/form , a number of hits against the database should appear, showing that the backend works.
