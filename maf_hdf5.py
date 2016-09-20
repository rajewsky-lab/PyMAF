"""
Part of the PyMAF package. This file holds the HDF5 wrapping of the data 
stored in MAF files as well as code to query the information
"""
__author__ = "Marvin Jens"
__copyright__ = "Copyright 2016, Massachusetts Institute of Technology"
__credits__ = ["Marvin Jens"]
__license__ = "MIT"
__version__ = "0.7"
__email__ = "mjens@mit.edu"
__status__ = "beta"

import tables as tbl
import numpy as np
import os, sys, re
import logging
from time import time

logging.basicConfig(level=logging.DEBUG)

class MAFBlockScores(tbl.IsDescription):
    """
    Keeps the alignment score assigned to a MAF block. Currently unused but may be 
    useful to break ties between overlapping MAF blocks
    """
    block_id = tbl.UInt32Col()
    score    = tbl.Float32Col()


class MAFBlockSpecies(tbl.IsDescription):
    """
    As there are no foreign keys in HDF5 we use the block_id to hold together
    different records for the same MAF block. This holds the coordinates and is 
    thus the main table for PyMAF.
    """
    block_id = tbl.UInt32Col()
    species  = tbl.StringCol(32)

    chrom = tbl.StringCol(32)
    start = tbl.UInt32Col()
    end   = tbl.UInt32Col()
    minus = tbl.BoolCol()

    # sequences are variable length and therefore need to be stored in a separate 
    # VLStringAtom array
    seq_id = tbl.UInt64Col()

    
class MAFBlockPads(tbl.IsDescription):
    """
    Wraps the padding information on how a block connects to neighboring blocks in
    each species (large gaps, etc.). currently this data is not used.
    """
    block_id = tbl.UInt32Col()
    species  = tbl.StringCol(32)

    pre_code  = tbl.StringCol(1)
    pre_num   = tbl.UInt32Col()
    post_code = tbl.StringCol(1)
    post_num  = tbl.UInt32Col()
    

class MAFRecord(object):
    """
    Thin wrapper around each species record in a MAF block.
    """
    def __init__(self, block, seq, block_attrs = ['chrom','species','start','end']):
        self.chrom = block['chrom']
        self.species = block['species']
        self.start = block['start']
        self.end = block['end']
        
        # would do the same but more cryptic.
        #for a in block_attrs:
            #setattr(self,a,block.__getitem__(a))
            
        if block['minus']:
            self.sense = '-'
        else:
            self.sense = '+'

        self.seq = seq
        
    def __str__(self):
        return "MAFRecord({self.species}  {self.chrom}:{self.start}-{self.end}{self.sense}) = '{self.seq}')".format(self=self)

class MAFBlockDB(object):
    """
    The interface between maf.gz (raw MAF), maf.h5 (HDF5) files on the one
    hand, and higher level on the other.
    """
    
    def __init__(self, fname=""):
        if fname:
            self.logger = logging.getLogger("MAFBlockDB({0})".format(fname))
            self.load(fname)
        else:
            self.logger = logging.getLogger("MAFBlockDB()")
    
    def create_from_gz(self, path, regular_flushes=10000, blosc_max_threads=8):
        """
        Converts a raw, gzip compressed MAF file to HDF5 format.
        
        :param path: path to maf.gz file
        :param regular_flushes: number of MAF blocks to process before pushing
            data to disk by PyTables.table.flush()
        :param blosc_max_threads: number of threads for parallel 
            (de-)compression by BLOSC (currently has no effect?)
        """

        from gzip import GzipFile
        tbl.set_blosc_max_threads(blosc_max_threads) # TODO: figure out why BLOSC is still single threaded?
        
        dirname, basename = os.path.split(path)
        name_parts = basename.split('.')[:-1]
        h5path = os.path.join(dirname, ".".join( name_parts + ['h5'] ))
        
        self.logger.info("creating hdf5 table '{0}' from '{1}'".format(h5path, path) )
        filters = tbl.Filters(complevel=5, complib='blosc')
        self.h5 = tbl.open_file(h5path, mode = "w", title = "MAF {0}".format(path), filters=filters)
        
        group = self.h5.create_group("/", 'maf', 'MAF records', filters=filters)
        
        block_tbl = self.h5.create_table(group, 'blocks', MAFBlockSpecies, "MAF block data for each species (a lines)")
        score_tbl = self.h5.create_table(group, 'scores', MAFBlockScores, "scores for each MAF block (s lines)")
        pad_tbl = self.h5.create_table(group, 'pads', MAFBlockPads, "pad codes for each MAF block (i lines)")
        
        blocks = block_tbl.row
        scores = score_tbl.row
        pads   = pad_tbl.row

        seqs = self.h5.create_vlarray(self.h5.root, 'seqs', atom = tbl.VLStringAtom(), title = "all sequences in the MAF blocks", filters=filters)
        n_seqs = 0
        
        maf_block_id = 0
        T0 = time()
        T_last = T0
        for maf_line in GzipFile(path):
            if not maf_line.strip():
                continue
            
            if maf_line.startswith('a'):
                maf_block_id += 1 # a new block!
                scores['block_id'] = maf_block_id
                scores['score'] = float(maf_line.split('=')[1])
                scores.append()

                if maf_block_id and (maf_block_id % regular_flushes == 0):
                    score_tbl.flush()
                    block_tbl.flush()
                    pad_tbl.flush()
                    seqs.flush()
                    
                    T = time()
                    dT = T - T_last
                    T_last = T
                    kbps = (regular_flushes / 1000.) / dT
                    
                    self.logger.debug("processed {2} {0:.0f}k MAF blocks ({1:.1f}k blocks per second)".format(maf_block_id/1000., kbps, maf_block_id) )
                
            elif maf_line.startswith('s'):
                parts = re.split(r'\s+',maf_line)
                loc,start,size,strand,total,seq = parts[1:7]
            
                species,chrom = loc.split('.',1)
                start = int(start)
                size = int(size)
                total = int(total)
                if strand == '+':
                    end = start + size
                else:
                    start, end = total - (start + size), total - start
                    # if you think this is sick, go tell the evil master MAF and his
                    # sidekick Dr. minus, the inventor of the minus strand, to their 
                    # faces. ;)
            
                blocks['block_id'] = maf_block_id
                blocks['species'] = species
                blocks['chrom'] = chrom
                blocks['start'] = start
                blocks['end'] = end
                blocks['minus'] = (strand == '-')
                blocks['seq_id'] = n_seqs
                
                seqs.append(seq.encode('ascii'))
                n_seqs += 1

                blocks.append()

            elif maf_line.startswith('i'):
                parts = re.split(r'\s+',maf_line)
                loc,pre_code,pre_num,post_code,post_num = parts[1:6]
                species,chrom = loc.split('.',1)
                
                pads['block_id'] = maf_block_id
                pads['species'] = species
                pads['pre_code'] = pre_code
                pads['pre_num'] = pre_num
                pads['post_code'] = post_code
                pads['post_num'] = post_num
                pads.append()
            elif maf_line.startswith('e'):
                continue # currently ignored
            else:
                print "ignoring unknown MAF line '{0}'".format(maf_line.strip())
                

        self.logger.info("done processing MAF blocks. Building indices")
        self.h5.root.maf.blocks.cols.chrom.create_index()
        self.h5.root.maf.blocks.cols.species.create_index()
        self.h5.root.maf.blocks.cols.block_id.create_index()    
        self.h5.root.maf.blocks.cols.start.create_index()
        self.h5.root.maf.blocks.cols.end.create_index()
        self.h5.root.maf.scores.cols.block_id.create_index()
        self.h5.root.maf.pads.cols.block_id.create_index()
        self.h5.root.maf.pads.cols.species.create_index()
        
    def load(self, path):
        """
        Open an HDF5 formatted MAF data file.
        :param path: path to mfa.h5 file.
        """
        self.h5 = tbl.open_file(path, 'r')
    
    def query_interval(self, ref_species, ref_chrom, ref_start, ref_end, select_species = []):
        """
        This is a generator that searches MAF blocks in the desired reference
        species' genomic interval and yields lists of MAFRecord() instances 
        for each MAF block.
        """
        blocks = self.h5.root.maf.blocks.row
        seqs = self.h5.root.seqs
        for ref_hit in self.h5.root.maf.blocks.where("(species == ref_species) & (chrom == ref_chrom) & (start < ref_end) & (end > ref_start)"):
            rows = []
            ref_id = ref_hit['block_id']
            
            for spc_hit in self.h5.root.maf.blocks.where("(block_id == ref_id)"):
                if select_species and not (spc_hit['species'] in select_species):
                    continue
                rows.append( MAFRecord(spc_hit, seqs[spc_hit['seq_id']]) )
                
            yield rows
    
if __name__ == "__main__":
    MAF = MAFBlockDB()
    MAF.create_from_gz(sys.argv[1])
    #MAF.load("/scratch/data/maf/mm10_60way/maf/chrY.maf.h5")
    MAF.get_maf_blocks_for_interval('mm10', 'chrY', 122550, 136258)
