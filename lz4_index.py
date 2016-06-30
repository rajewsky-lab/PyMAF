#!/usr/bin/env python
import time, os, sys
from byo.io.lz import LZFile
import logging

def test_stitching(fname):
    uf = file(fname).read()
    lf = LZFile(fname)
    
    cs = lf.chunk_size
    L = lf.L
    
    tests = [
        (0,10),
        (cs - 10, cs),
        (cs, cs + 10),
        (L-10,L),
        (3*cs - 10, 3*cs + 10),
    ]
    for s,e in tests:
        print "UNCOMPRESSED '{0}'".format(uf[s:e])
        print "COMPRESSED   '{0}'".format(lf[s:e])

def test_speed(fname, n=1000000):
    from time import time
    from numpy.random import randint
    
    uf = file(fname).read()
    lf = LZFile(fname)
    
    rs = 10 # read-size
    starts = [randint(lf.L-rs) for x in xrange(n)]
    
    def benchmark(f):
        t0 = time()
        for s in starts:
            seq = f[s:s+rs]
        
        t1 = time()
        return t1-t0

    print "uncompressed",benchmark(uf)
    print "compressed",benchmark(lf)

def test_random(fname, n=100000):
    from numpy.random import randint
    
    uf = file(fname).read()
    lf = LZFile(fname)
    
    rs = 10 # read-size
    starts = [randint(lf.L-rs) for x in xrange(n)]
    
    for s in starts:
        assert uf[s:s+rs] == lf[s:s+rs]
    
    print len(lf.chunk_cache)
    
    #test_stitching("sacCer3.fa")
    #test_speed("sacCer3.fa")
    #test_random("sacCer3.fa")
    #LZFile.compress_file("test.maf")
    #test_random("test.maf")
    #test_speed("test.maf")


if __name__ == '__main__':
    from optparse import *
    usage = """
    usage: %prog [options] 
    """

    parser = OptionParser(usage=usage)
    parser.add_option("","--debug",dest="debug",default=False,action="store_true",help="activate extensive debug output (default=Off)")
    parser.add_option("-d","--decompress", dest="decompress", default=False, action="store_true", help='decompress an LZ4 compressed file (.lzot/.lzoc suffixes will be removed from filename)')
    options,args = parser.parse_args()

    # prepare logging system
    FORMAT = '%(asctime)-20s\t%(levelname)s\t%(name)-25s\t%(message)s'
    if options.debug:
        logging.basicConfig(level=logging.DEBUG,format=FORMAT)
    else:
        logging.basicConfig(level=logging.INFO,format=FORMAT)
    logger = logging.getLogger('lz4_index.main()')

    if options.decompress:
        for line in LZFile(args[0].replace('.lzot','').replace('.lzoc','')):
            print line.rstrip()
    else:
        LZFile.compress_file(args[0])
        
