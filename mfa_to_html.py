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
<p>built with mfa_to_html.py from '%s'</p>
\n<table>""" % os.path.basename(sys.argv[1])

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


   
chunks = list(fasta_chunks(file(sys.argv[1])))
fa_id_lens = [len(fa_id.split()[0]) for fa_id,seq in chunks]

longest_fa_id = max(*fa_id_lens)+3
parts = [fa_id.split() for fa_id,seq in chunks]
short_headers = [p[0] for p in parts]
extra = [" ".join(p[1:]) if len(p) > 1 else "" for p in parts]
seqs = [seq for fa_id,seq in chunks]

print head

for short,ex in zip(short_headers,extra):
    print "<tr><td>{short}</td><td>{ex}</td></tr>".format(**locals())
print "</table><table>"

def wrap(headers,seqs,max_len=100):
    start = 0
    L = max(*[len(s) for s in seqs])
    #print seqs[:2]
    while L >= max_len:
        yield headers,start,[s[start:start+max_len] for s in seqs]
        start += max_len
        L -= max_len

    if L > 0:
        yield headers,start,[s[start:start+max_len] for s in seqs]

    
for headers,start,wrapped_seqs in wrap(short_headers,seqs):
    print '<tr><td style="padding : 1em; width : 3em;">+{start}</td></tr>'.format(**locals())
    for short,seq in zip(headers,wrapped_seqs):
        if short.startswith('>'):
            colored = seq
        else:
            colored = colorize(seq)
        print '<tr><td style="width : {longest_fa_id}; padding-right : 1em;">{short}</td>\n<td>{colored}</td></tr>'.format(**locals())
    

print foot
