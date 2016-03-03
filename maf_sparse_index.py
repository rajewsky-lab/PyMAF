import numpy as np
import os,re,sys
import optparse
import logging

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger('maf_sparse_index.py')

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

def maf_to_loci(path,ref):
    maf_list = []
    prefix = "s %s." % ref
    i = 0
    prev_i = 0
    for maf_line in file(path):
        if maf_line.startswith(prefix):
            parts = re.split(r'\s+',maf_line)
            start,size,strand,total = parts[2:6]
            start = int(start)
            end = start + int(size)
            # record offset of *previous line*, including the score!
            maf_list.append( (start,end,prev_i) )
        prev_i = i
        i += len(maf_line)

    return maf_list,int(total)

fname = sys.argv[1]
ref = sys.argv[2]

path = os.path.dirname(fname)
name = os.path.basename(fname)

comb_fname = os.path.join(path,name.replace(".maf",".comb_idx"))
bin_fname = os.path.join(path,name.replace(".maf",".comb_bin"))

log.info("collecting maf-blocks from '%s' for reference '%s'" % (fname,ref))
maf_list,chrom_size = maf_to_loci(fname,"dm6")
log.info("collected %d blocks." % len(maf_list) )

log.info("writing combinations to '%s' and lookup to binary sparse file '%s'" % (comb_fname,bin_fname))
comb = file(comb_fname,"w")
lkup = np.memmap(bin_fname,dtype=np.uint32,mode="w+",shape=(chrom_size,))

for j,(combination,start,end) in enumerate(pileup_loci(maf_list)):
    lkup[start:end] = j+1
    
    comb_str = ",".join(["%x" % maf_line for start,end,maf_line in combination]) + '\n'
    comb.write(comb_str)
    
    if j and not j % 1000:
        perc = 100. * start / float(chrom_size)
        log.debug("%.2f %% of %s" % (perc,name))