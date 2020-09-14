#!/usr/bin/python
import sys
from django.core.handlers import wsgi

def application(environ, start_response):
    status = '200 OK'

    output = u''
    output += u'sys.version = %s\n' % repr(sys.version)
    output += u'sys.prefix = %s\n' % repr(sys.prefix)

    response_headers = [('Content-type', 'text/plain'),
                        ('Content-Length', str(len(output)))]
    start_response(status, response_headers)

    return [output.encode('UTF-8')]
