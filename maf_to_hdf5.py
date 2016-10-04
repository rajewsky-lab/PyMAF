#!/usr/bin/env python
import numpy as np
import os,re,sys
import optparse
import logging
from PyMAF.maf_hdf5 import MAFBlockDB

if __name__ == "__main__":
    from optparse import OptionParser
    usage = """
    usage: %prog [options] <MAF_file.maf.gz>
    """

    parser = OptionParser(usage=usage)
    parser.add_option("","--debug",dest="debug",default=False,action="store_true",help="activate extensive debug output (default=Off)")
    parser.add_option("","--flush",dest="flush",default=10000,type=int,help="number of blocks to process before flushing to disk (default=10000)")
    parser.add_option("","--threads",dest="threads",default=8,type=int,help="max. number of threads to use for BLOSC compressor (default=8)")
    parser.add_option("","--complevel",dest="complevel",default=2,type=int,help="BLOSC compression level (0=off, 10=max, default=2)")
    options,args = parser.parse_args()

    # prepare logging system
    FORMAT = '%(asctime)-20s\t%(levelname)s\t%(name)-25s\t%(message)s'
    if options.debug:
        logging.basicConfig(level=logging.DEBUG,format=FORMAT)
    else:
        logging.basicConfig(level=logging.INFO,format=FORMAT)

    logger = logging.getLogger('maf_to_hdf5.main()')
    hdf5 = MAFBlockDB()
    if not args:
        parser.error("path to MAF required")
    else:
        hdf5.create_from_gz(args[0], regular_flushes=options.flush, blosc_max_threads=options.threads, complevel=options.complevel)
    
