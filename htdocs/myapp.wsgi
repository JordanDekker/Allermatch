#def application(environ, start_response):
#	status = '200 OK'
#	output = b'Hello World!'
#	response_headers = [('Content-type', 'text/plain'),
#			('Content-Length', str(len(output)))]
#	start_response(status, response_headers)
#	return [output]

import logging
import sys
logging.basicConfig(stream=sys.stderr)
sys.path.insert(0, '/var/www/html/allermatch/htdocs/')
from allermatch import app as application
application.secret_key = 'anything you wish'

