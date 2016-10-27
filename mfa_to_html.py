#!/usr/bin/env python

import os,sys,re
from byo.io import fasta_chunks

head = """<html><head>
<style>
table, th, td {
    border: 0px;
    padding: 0px;
    cellpading : 0px;
    cellspacing : 0px;
    font-family : monospace;
}

th, td {
    padding: 0px;
    white-space : nowrap;
}

span { background-color : White; margin : 0px; padding : 0px; }
span.a { background-color : SkyBlue; }
span.c { background-color : Tomato; }
span.g { background-color : Plum; }
span.t { background-color : Khaki; }
span.n { background-color : Silver; }
span.gap { background-color : WhiteSmoke; }
span.w { background-color : White; }

</style></head>
<body>
<h1>Multiple species alignment</h1>
<p>built with 'mfa_to_html.py %s'</p>
\n<table>\n""" % os.path.basename(" ".join(sys.argv))

foot = "</table></body></html>"

def colorize(seq):
    classmap = {
        't' : 't',
        'T' : 't',
        'u' : 't',
        'U' : 't',
        'a' : 'a',
        'A' : 'a',
        'c' : 'c',
        'C' : 'c',
        'g' : 'g',
        'G' : 'g',
        'n' : 'n',
        'N' : 'n',
        '-' : 'gap',
    }
    spans = []
    runlen = 0
    last = ""
    for s in seq:
        if s != last:
            if runlen:
                spans.append('<span class="%s">%s</span>' % (classmap.get(last,'w'),last*runlen))
                runlen = 0
        last = s
        runlen += 1

    if runlen:
        spans.append('<span class="%s">%s</span>' % (classmap.get(last,'w'),last*runlen))

    return "".join(spans)


def highlight_patterns(seq, matches, mapping):
    colormap = {}
    colors = ['#b2e2e2','Silver','#66c2a4','Silver','#2ca25f','Silver','#006d2c']
    
    for m in matches:
        n_groups = len(m.groups())
        for i in range(n_groups):
            start, end = m.span(i+1)
            for j in range(start,end):
                colormap[mapping[j]] = colors[i % len(colors)]
        
    spans = []
    for i, s in enumerate(seq):
        if i in colormap:
            spans.append('<span style="background-color : {0};">{1}</span>'.format(colormap[i],s))
        else:
            spans.append(s)

    return "".join(spans)

        
   
from optparse import *
usage = """
usage: %prog [options] <input_file.msf>

Produces a pretty HTML rendering of a multiple-species FASTA (MSF) sequence alignment file.
"""

parser = OptionParser(usage=usage)
parser.add_option("-w","--wrap",dest="wrap",type=int,default=100,help="line wrap after this many columns. Set to 0 to disable (default=100)")
parser.add_option("-x","--extra",dest="extra",action="store_true",default=False,help="Switch: activate rendering of additional information in the MSF FASTA headers")
parser.add_option("-p","--pattern",dest="pattern",default="",help="highlight occurrences of this regexp pattern (default=off). Does not work with wrapping, yet")
parser.add_option("-d","--destination",dest="outfile",default="",help="destination folder for HTML files)")
parser.add_option("-O","--stdout",dest="stdout",default=False, action="store_true",help="Switch: instead of writing files, write to stdout")
options,args = parser.parse_args()

if options.outfile:
    out = file(options.outfile,'w')
else:
    base, ext = os.path.splitext(args[0])
    out = file(base+'.html','w')

chunks = list(fasta_chunks(file(args[0])))
fa_id_lens = [len(fa_id.split()[0]) for fa_id,seq in chunks]

longest_fa_id = max(*fa_id_lens)+3
parts = [fa_id.split() for fa_id,seq in chunks]
short_headers = [p[0] for p in parts]
extra = [" ".join(p[1:]) if len(p) > 1 else "" for p in parts]
seqs = [seq for fa_id,seq in chunks]

if options.pattern:
    from PyMAF import Alignment
    aln = Alignment(file(args[0]).read())
    matches_list = [[M for M in re.finditer(options.pattern,aln.nogaps[spc])] for spc in aln.species]


out.write(head)
if options.extra:
    for short,ex in zip(short_headers,extra):
        out.write("<tr><td>{short}</td><td>{ex}</td></tr>\n".format(**locals()))
    out.write("</table><table>\n")

def wrap(headers,seqs,max_len=100):
    if not max_len:
        yield headers,0,seqs
    else:
        start = 0
        L = max(*[len(s) for s in seqs])
        #print seqs[:2]
        while L >= max_len:
            yield headers,start,[s[start:start+max_len] for s in seqs]
            start += max_len
            L -= max_len

        if L > 0:
            yield headers,start,[s[start:start+max_len] for s in seqs]

    
for headers, start, wrapped_seqs in wrap(short_headers, seqs, options.wrap):
    out.write('<tr><td style="padding : 1em; width : 3em;">+{start}</td></tr>\n'.format(**locals()))
    
    if options.pattern:
        for short, seq, matches in zip(headers,wrapped_seqs, matches_list):
            colored = highlight_patterns(seq, matches, aln.spc_to_mfa[short])
            out.write('<tr><td style="width : {longest_fa_id}; padding-right : 1em;">{short}</td>\n<td>{colored}</td></tr>\n'.format(**locals()))
    else:
        for short,seq in zip(headers,wrapped_seqs):
            if short.startswith('>'):
                colored = seq
            else:
                colored = colorize(seq)
            out.write('<tr><td style="width : {longest_fa_id}; padding-right : 1em;">{short}</td>\n<td>{colored}</td></tr>\n'.format(**locals()))

out.write(foot)
