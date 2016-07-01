#!/usr/bin/env python
import numpy as np
import os,re,sys
import optparse
import logging
from byo.io.lz import LZFile


def pileup_loci(loci_list):
    """
    returns all distinct combinations of overlapping (or single) loci with start 
    and end-coordinates of the block.
    Works by sorting loci by start, then end positions. Then, iterates over all loci, keeping a 
    sorted list of end-coordinates to detect which elements to kick out and which to keep.
    Should be near optimal, probably ~O(N log(N)), dominated by the initial sorting.
    """
    from collections import namedtuple
    L = namedtuple("locus","start,end,annotation_id")

    from bisect import bisect_left
    
    # will keep 'short' elements on top
    loci = []
    ends = []
    last_start = 0

    for l in sorted(loci_list):
        l = L(*l)
        if l.start > last_start:
            # find everything that ends before this new block and yield the combinations
            while ends and (ends[0] <= l.start):
                yield loci,last_start,ends[0]
                # advance the "cursor" until the end of the smallest element
                last_start = ends[0]
                # and kick it out. (Could be more than one element of same size)
                while ends and ends[0] == last_start:
                    ends.pop(0)
                    loci.pop(0)

            # maybe there is still some space, yield until the new locus begins
            if last_start < l.start:
                if loci: yield loci,last_start,l.start

        # now we are at the start position of the new block
        last_start = l.start

        # find appropriate position by end-position
        pos = bisect_left(ends,l.end)
        loci.insert(pos,l)
        ends.insert(pos,l.end)

    # we're through, yield the remaining blocks away
    while ends:
        yield loci,last_start,ends[0]
        loci.pop(0)
        ends.pop(0)

def maf_to_loci(path,ref, compress=True):
    maf_list = []
    prefix = "s %s." % ref
    i = 0
    prev_i = 0
    total = 0

    # preparing the MAF file
    if path.endswith('.gz'):
        # directly read from GZ compressed input
        logger.info("directly reading from GZip compressed file")
        from gzip import GzipFile
        lz = LZFile(
            path.replace('.gz',''), 
            compress_on_open=compress, 
            alt_src = GzipFile(path, 'rb'),
            chunksize = options.chunksize * 1024**2
        )
    else:
        # read from uncompressed input
        logger.info("reading from uncompressed file")
        lz = LZFile(path, compress_on_open=compress, chunksize = options.chunksize * 1024**2)

    for maf_line in lz:
        if maf_line.startswith(prefix):
            parts = re.split(r'\s+',maf_line)
            start,size,strand,total = parts[2:6]
            start = int(start)
            end = start + int(size)
            # record offset of *previous line*, including the score!
            maf_list.append( (start,end,prev_i) )
        prev_i = i
        i += len(maf_line)
    if not total:
        print "reference '%s' never occured in this MAF file. Are you sure you are doing this right?" % ref
    return maf_list,int(total)


if __name__ == "__main__":
    from optparse import OptionParser
    usage = """
    usage: %prog [options] <MAF_file.maf|MAF_file.maf.gz> <reference_species>
    """

    parser = OptionParser(usage=usage)
    parser.add_option("","--debug",dest="debug",default=False,action="store_true",help="activate extensive debug output (default=Off)")
    parser.add_option("-u","--uncompressed", dest="uncompressed", default=False, action="store_true", help='do not automatically compress the MAF with LZ4 (default is compression)')
    parser.add_option("", "--chunksize", dest="chunksize", type=int, default=10, help="chunk-size in MB (default=10)")
    parser.add_option("", "--lz-level", dest="lz_level", type=int, default=2, help="LZ4 compression level. Balance speed vs space (default=2)")
    options,args = parser.parse_args()

    # prepare logging system
    FORMAT = '%(asctime)-20s\t%(levelname)s\t%(name)-25s\t%(message)s'
    if options.debug:
        logging.basicConfig(level=logging.DEBUG,format=FORMAT)
    else:
        logging.basicConfig(level=logging.INFO,format=FORMAT)
    logger = logging.getLogger('maf_sparse_index.main()')

    fname, ref = args

    path = os.path.dirname(fname)
    name = os.path.basename(fname)

    comb_fname = os.path.join(path,name.replace(".maf",".comb_idx").replace('.gz',''))
    bin_fname = os.path.join(path,name.replace(".maf",".comb_bin").replace('.gz',''))

    logger.info("collecting maf-blocks from '{0}' for reference '{1}'".format(fname,ref))
    maf_list,chrom_size = maf_to_loci(fname,ref, compress=not options.uncompressed)
    
    logger.info("collected {0} blocks.".format(len(maf_list)) )

    logger.info("writing combinations to '{0}' and lookup to binary sparse file '{1}'".format(comb_fname,bin_fname))
    comb = file(comb_fname,"w")
    lkup = np.memmap(bin_fname,dtype=np.uint32,mode="w+",shape=(chrom_size,))

    for j,(combination,start,end) in enumerate(pileup_loci(maf_list)):
        lkup[start:end] = j+1
        
        comb_str = ",".join(["%x" % maf_line for start,end,maf_line in combination]) + '\n'
        comb.write(comb_str)
        
        if j and not j % 100000:
            perc = 100. * start / float(chrom_size)
            logger.debug("{0:.2f} % of {1}".format(perc,name))
            
    logger.info("done.")