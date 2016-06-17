import lz4framed as Z
import os,sys
from time import time

def chunks(fname, chunksize):
    f = file(fname,'rb')
    chunk = f.read(chunksize)
    while chunk:
        yield chunk
        chunk = f.read(chunksize)

MB = 1024.**2

def compress_file(fname, chunksize=1*1024*1024):
    tab_file = file(fname + '.lzot','w')
    comp_file = file(fname + '.lzoc','wb')
    comp_base = 0
    cum_size = 0
    t0 = time()
    
    tab_file.write('{0}\n'.format(chunksize))

    for chunk in chunks(fname, chunksize):
        uncomp_size = len(chunk)

        t1 = time()
        comp = Z.compress(chunk, level=2)
        comp_size = len(comp)
        comp_file.write(comp)
        ratio = 100. * float(comp_size) / uncomp_size 
        t2 = time()
        throughput = cum_size / (t2-t0)

        tab_file.write('{0}\n'.format(comp_base))
        comp_base += comp_size
        cum_size += uncomp_size

        print "compressed {0}MB ({1:.1f}%) in {2:.1f} sec, {3:.2f} MB/s  sec".format(chunksize/MB, ratio, t2-t1, throughput/MB)

    tab_file.write('{0}\n'.format(comp_base))
    tab_file.write('{0}\n'.format(cum_size))

class LZFile(object):
    
    def __init__(self, fname, max_cached = 1000):
        self.basename = fname
        self.load_index(fname + '.lzot')
        self.lz_file = file(fname + '.lzoc','rb')
        self.chunk_cache = {}
        self.max_cached = max_cached
        
    def load_index(self, idxname):
        lines = file(idxname).readlines()
        self.chunk_size = int(lines[0])
        self.L = int(lines[-1])
        self.chunk_starts = [int(b) for b in lines[1:-1]]
        
    def get_chunk(self, i):
        self.lz_file.seek(self.chunk_starts[i])
        comp = self.lz_file.read(self.chunk_starts[i+1] - self.chunk_starts[i])
        return Z.decompress(comp)
    
    def get_chunk_cached(self, i):
        if not i in self.chunk_cache:
            self.chunk_cache[i] = self.get_chunk(i)
            #self.cached_items.append(i)

        if len(self.chunk_cache) > self.max_cached:
            pass
            # not implemented yet: efficient way to discard least used chunks
            
        return self.chunk_cache[i]
            
    def __getslice__(self, start, end):
        cs = self.chunk_size
        out = []
        for chunk_i in range(start / cs, (end / cs) + 1):

            chunk_start = chunk_i * cs
            
            c_start = max(start - chunk_start, 0)
            c_end = min(end - chunk_start, cs)
            #print "CHUNK_I", chunk_i, c_start, c_end, cs
            
            out.append(self.get_chunk_cached(chunk_i)[c_start:c_end])
            
        return "".join(out)
    

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
    
if __name__ == "__main__":
    #test_stitching("sacCer3.fa")
    #test_speed("sacCer3.fa")
    #test_random("sacCer3.fa")
    #compress_file("test.maf")
    #test_random("test.maf")
    #test_speed("test.maf")
    
    compress_file(sys.argv[1])
    
    #compress_file("mm10.fa")
    
    
    
        
    