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
import os, re
import logging
from time import time

class MAFBlockScores(tbl.IsDescription):
    """
    Keeps the alignment score assigned to a MAF block. Currently unused but may be 
    useful to break ties between overlapping MAF blocks
    """
    block_id = tbl.UInt32Col()
    score    = tbl.Float32Col()


class MAFBlockCoords(tbl.IsDescription):
    """
    As there are no foreign keys in HDF5 we use the block_id to hold together
    different records for the same MAF block. This holds the coordinates and is 
    thus the main table for PyMAF.
    """
    block_id = tbl.UInt32Col()

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

    pre_code  = tbl.StringCol(1)
    pre_num   = tbl.UInt32Col()
    post_code = tbl.StringCol(1)
    post_num  = tbl.UInt32Col()
    

class MAFRecord(object):
    """
    Thin wrapper around each species record in a MAF block.
    """
    def __init__(self, species, coords, seq, coords_attrs = ['chrom','species','start','end']):
        self.chrom = coords['chrom']
        self.species = species
        self.start = coords['start']
        self.end = coords['end']
        
        # would do the same but more cryptic.
        #for a in coords_attrs:
            #setattr(self,a,coords.__getitem__(a))
            
        if coords['minus']:
            self.sense = '-'
        else:
            self.sense = '+'

        self.seq = seq
        
    def __str__(self):
        return "MAFRecord({self.species}  {self.chrom}:{self.start}-{self.end}{self.sense}) = '{self.seq}')".format(self=self)


class CursorCache(object):
    def __init__(self, group, cls, h5):
        self.group = group
        self.cache = {}
        self._table_cache = {}
        self.cls = cls
        self.h5 = h5

    def __getitem__(self, species):
        if not species in self.cache:
            # create new table for self species upon first encounter
            #spc_group = self.h5.create_group(self.group, 'maf', 'MAF records', filters=filters)
            new_tbl = self.h5.create_table(self.group, species, self.cls, "data for species {0}".format(species), expectedrows=1000000)
            self._table_cache[species] = new_tbl
            self.cache[species] = new_tbl#.row
            
        return self.cache[species]

    def flush(self):
        for t in self._table_cache.values():
            t.flush()

class SpeciesListReconstruct(object):
    def __init__(self):
        from collections import defaultdict
        self.before = defaultdict(set)

    def record(self, l):
        for i,s in enumerate(l):
            self.before[s] |= set(l[:i])
    
    def best_guess(self):
        
        def find_first():
            for s,bef in self.before.items():
                if len(bef) == 0:
                    return s

        def remove(s):
            del self.before[s]

            for spc, before in self.before.items():
                before.discard(s)
                self.before[spc] = before
        
        res = []
        while self.before:
            next = find_first()
            res.append(next)
            remove(next)
        
        return res
        
class MAFBlockDB(object):
    """
    The interface between maf.gz (raw MAF), maf.h5 (HDF5) files on the one
    hand, and higher level on the other.
    """
    
    def __init__(self, fname="", in_memory=False):
        self.in_memory = in_memory
        self.h5 = None
        self.caches = set()

        if fname:
            self.load(fname)
    
        import multiprocessing
        tbl.set_blosc_max_threads(multiprocessing.cpu_count()) # TODO: figure out why BLOSC is still single threaded?

    def flush(self, species_names = []):
        self.logger.debug("flushing rows to disk...")
        t0 = time()
        
        if self.caches:
            self.logger.debug("flushing species caches")
            for cache in self.caches:
                # first flush tables that are completely cached and may not even exist in the file, yet
                cache.flush()

        self.h5.root.scores.flush()
        self.h5.root.seqs.flush()
        self.h5.root.species_names.flush()

        if not species_names:
            species_names = self.species_names

        for species in species_names:
            t = getattr(self.h5.root.coords, species, None)
            if t:
                #self.logger.debug("flushing /coords/{0}".format(species))
                t.flush()
            else:
                self.logger.warning("no coords table for '{0}'".format(species) )
            
            t = getattr(self.h5.root.pads, species, None)
            if t: 
                #self.logger.debug("flushing /pads/{0}".format(species))
                t.flush()
            else:
                self.logger.warning("no pads table for '{0}'".format(species) )
                
        t1 = time()
        self.logger.debug("flush completed in {0:.1f}ms".format( (t1-t0)*1000.) )
 
 
    def build_indices(self, species_names = []):
        self.flush()

        if not species_names:
            species_names = self.species_names

        t0 = time()
        for species in species_names:
            t = getattr(self.h5.root.coords,species, None)
            if t:
                self.logger.debug("creating indices for '{0}'".format(species))
                cols = t.cols
                cols.chrom.create_csindex()
                cols.block_id.create_csindex()
                cols.start.create_csindex()
                cols.end.create_csindex()
            else:
                self.logger.warning("no coords table for '{0}'".format(species) )
            
        t1 = time()
        self.logger.debug("indexing completed in {0:.1f}ms".format((t1-t0)*1000.) )


    def create_from_gz(self, path, regular_flushes=10000, blosc_max_threads=8, complevel=0):
        """
        Converts a raw, gzip compressed MAF file to HDF5 format.
        
        :param path: path to maf.gz file
        :param regular_flushes: number of MAF blocks to process before pushing
            data to disk by PyTables.table.flush()
        :param blosc_max_threads: number of threads for parallel 
            (de-)compression by BLOSC (currently has no effect?)
        """

        tbl.set_blosc_max_threads(blosc_max_threads) # TODO: figure out why BLOSC is still single threaded?
        
        dirname, basename = os.path.split(path)
        name_parts = basename.split('.')[:-1]
        h5path = os.path.join(dirname, ".".join( name_parts + ['hdf5'] ))
        
        self.logger = logging.getLogger("MAFBlockDB({0})".format(h5path))
        self.logger.info("creating hdf5 table '{0}' from '{1}'".format(h5path, path) )
        
        filters = tbl.Filters(complevel=complevel, complib='blosc')
        self.h5 = tbl.open_file(h5path, mode = "w", title = "MAF {0}".format(path), filters=filters)
        
        seqs            = self.h5.create_vlarray(self.h5.root, 'seqs', atom = tbl.VLStringAtom(), title = "all sequences in all MAF blocks", filters=filters)
        names_table     = self.h5.create_vlarray(self.h5.root, 'species_names', atom = tbl.VLStringAtom(), title = "all species names mentioned in the MAF file")
        score_tbl       = self.h5.create_table(self.h5.root, 'scores', MAFBlockScores, "a records for each MAF block", filters=filters)
        species_coords  = self.h5.create_group(self.h5.root, 'coords', 's records for each species', filters=filters)
        species_pads    = self.h5.create_group(self.h5.root, 'pads', 'i records for each species', filters=filters)
        
        coords_lookup = CursorCache(species_coords, MAFBlockCoords, self.h5)
        pads_lookup = CursorCache(species_pads, MAFBlockPads, self.h5)
        
        self.caches = [coords_lookup, pads_lookup]
        scores = score_tbl.row
        
        maf_block_id = 0
        n_seqs = 0

        T0 = time()
        T_last = T0

        all_species = SpeciesListReconstruct()
        covered_species = []

        from gzip import GzipFile
        for maf_line in GzipFile(path):
            if not maf_line.strip():
                continue
            
            if maf_line.startswith('a'):
                maf_block_id += 1 # a new block!
                scores['block_id'] = maf_block_id
                scores['score'] = float(maf_line.split('=')[1])
                scores.append()

                all_species.record(covered_species)
                covered_species = []

                if maf_block_id and (maf_block_id % regular_flushes == 0):
                    self.flush(species_names = coords_lookup.cache.keys() )

                    T = time()
                    dT = T - T_last
                    T_last = T
                    kbps = (regular_flushes / 1000.) / dT
                    
                    self.logger.debug("processed {0:.0f}k MAF blocks ({1:.1f}k blocks per second)".format(maf_block_id/1000., kbps) )
                
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

                covered_species.append(species)

                coords = coords_lookup[species]
                coords.append( [(maf_block_id, chrom, start, end, (strand == '-'), n_seqs)] )
                                
                #coords['block_id'] = maf_block_id
                #coords['chrom'] = chrom
                #coords['start'] = start
                #coords['end'] = end
                #coords['minus'] = (strand == '-')
                
                # store the alignment row in the VLStringArray. 
                # the link is through keeping the n_seqs value
                seqs.append(seq.encode('ascii'))
                #coords['seq_id'] = n_seqs
                n_seqs += 1
                #coords.append()

            elif maf_line.startswith('i'):
                parts = re.split(r'\s+',maf_line)
                loc,pre_code,pre_num,post_code,post_num = parts[1:6]
                species,chrom = loc.split('.',1)
                
                pads = pads_lookup[species]
                #pads['block_id'] = maf_block_id
                #pads['pre_code'] = pre_code
                #pads['pre_num'] = pre_num
                #pads['post_code'] = post_code
                #pads['post_num'] = post_num
                pads.append( [(maf_block_id, pre_code, pre_num, post_code, post_num)] )
            elif maf_line.startswith('e'):
                continue # currently ignored
            else:
                print "ignoring unknown MAF line '{0}'".format(maf_line.strip())
                

        self.logger.info("done processing {0} MAF blocks in {1:.1f}sec.".format(maf_block_id, (time() - T0)) )

        self.species_names = all_species.best_guess()
        self.logger.debug("storing {0} species names ({1})".format(len(self.species_names), ",".join(self.species_names[:100])))

        for s in self.species_names:
            names_table.append(s.encode('ascii'))
        
        self.build_indices()
        self.caches = []

        
    def load(self, path):
        """
        Open an HDF5 formatted MAF data file.
        :param path: path to mfa.h5 file.
        """

        self.logger = logging.getLogger("MAFBlockDB({0})".format(path))
        if self.in_memory:
            self.h5 = tbl.open_file(path, 'r',driver="H5FD_CORE")
        else:
            self.h5 = tbl.open_file(path, 'r')

        self.species_names = self.h5.root.species_names[:]
        self.logger.info("opened '{0}', with data from {1} species".format(path, len(self.species_names)))
    
    def query_interval(self, ref_species, ref_chrom, ref_start, ref_end, select_species = []):
        """
        This is a generator that searches MAF blocks in the desired reference
        species' genomic interval and yields lists of MAFRecord() instances 
        for each MAF block.
        """
        t0 = time()
        if not select_species:
            select_species = self.species_names
            
        t = getattr(self.h5.root.coords,ref_species)
        t3 = t0
        n_blocks = 0
        n_rows = 0
        for ref_hit in t.where("(chrom == ref_chrom) & (start < ref_end) & (end > ref_start)"):
            n_blocks += 1
            rows = []
            ref_id = ref_hit['block_id']
            t1 = time()
            self.logger.debug("ref_hit found in {0:.1f}ms.".format( (t1-t3)*1000.) )

            for species in select_species:
                for spc_hit in getattr(self.h5.root.coords,species).where("(block_id == {0})".format(ref_id)):
                    rows.append( MAFRecord(species, spc_hit, self.h5.root.seqs[spc_hit['seq_id']]) )
                    n_rows += 1
                
            t2 = time()
            yield rows
            t3 = time()
            
            self.logger.debug("collected rows for block in {0:.1f}ms. process after yield in {1:.2f}ms".format( (t2-t1)*1000., (t3-t2)*1000. ) )
            
        t4 = time()
        elapsed = (t4 - t0)*1000
        self.logger.debug("query_interval({ref_species} {ref_chrom}:{ref_start}-{ref_end}) found {n_rows} rows in {n_blocks} blocks. Took {elapsed:.1f}ms.".format( **locals() ) )

    def close(self):
        if self.h5:
            self.h5.close()
