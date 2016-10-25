import BaseHTTPServer
import re, os, sys
from byo.io.genome_accessor import GenomeCache
import logging
from optparse import OptionParser
usage = """
usage: %prog [options] <genome_path>
"""

USAGE="""
<html><body>
<h1>Malformed request</h1>
Requests must be of the following form:
<p>
    '/(?P&lt;species&gt;\w+)/(?P&lt;chrom&gt;\S+)\:(?P&lt;start&gt;\d+)\-(?P&lt;end&gt;\d+)(?P&lt;strand&gt;[\+,\-])'
</p>
Example:
<p>
    'http://heartofgold:8000/hg19/chrX:1002345-1002349+'
</p>
Note: Coordinates are zero-based, half open (as used in C, Python, etc.).

</body></html>
"""

genome_cache = None

class GenomeRequestHandler(BaseHTTPServer.BaseHTTPRequestHandler):
    
    def _set_headers(self, code=200):
        self.send_response(code)
        self.send_header('Content-type', 'text/html')
        self.end_headers()

    def do_GET(self):
        M = re.match(r'/(?P<species>\w+)/(?P<chrom>\S+)\:(?P<start>\d+)\-(?P<end>\d+)(?P<strand>[\+,\-])',self.path)
        if not M:
            self._set_headers(400)
            self.wfile.write(USAGE)
        else:
            params = M.groupdict()
            
            species = params['species']
            chrom = params['chrom']
            start = int(params['start'])
            end = int(params['end'])
            strand = params['strand']

            self._set_headers(200)
            seq = genome_cache[species].get_oriented(chrom, start, end, strand)
            self.wfile.write(seq)

    def do_HEAD(self):
        self._set_headers()
        
    def do_POST(self):
        # Doesn't do anything with posted data
        self._set_headers()
        self.wfile.write("<html><body><h1>POST!</h1></body></html>")
   
def run(server_class=BaseHTTPServer.HTTPServer,
        handler_class=GenomeRequestHandler, host='', port=8000):

    server_address = (host, port)
    httpd = server_class(server_address, handler_class)
    httpd.serve_forever()
    

if __name__ == "__main__":
    parser = OptionParser(usage=usage)
    parser.add_option("","--host",dest="host",default='127.0.0.1',help="hostname (default=127.0.0.1)")
    parser.add_option("","--port",dest="port",default=8000,type=int,help="port number to use (default=8000)")
    parser.add_option("","--debug",dest="debug",default=False, action="store_true",help="SWITCH: activate debug output")
    parser.add_option("","--logfile",dest="logfile",default="/dev/stderr", help="log file (default: use stderr)")
    options,args = parser.parse_args()

    if options.debug:
        lvl = logging.DEBUG
    else:
        lvl = logging.INFO
    FORMAT = '%(asctime)-20s\t%(levelname)s\t%(name)s\t%(message)s'
    logging.basicConfig(level=lvl, format=FORMAT, filename=options.logfile, filemode='w')    

    genome_cache = GenomeCache(args[0])
    run(host=options.host, port=options.port)
