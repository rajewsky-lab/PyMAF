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
import newick
import numpy as np

def find_tree(path, levels=3):
    from glob import glob
    base = os.path.dirname(path)
    for i in range(levels):
        search_dir = os.path.join(base,"../"*i)
        for fname in glob(os.path.join(search_dir,"*.nh")):
            if not "commonNames" in fname:
                return fname

def species_list_from_tree(fname):
    species = []
    from newick import load
    tree = load(file(fname))
    for node in tree[0].walk():
        if node.name:
            species.append(node.name)

    return species

#print species_list_from_tree( find_tree('/scratch/data/maf/mm10_60way/maf/chr1.maf.gz') )
#import sys
#sys.exit(0)

class MAFCoords(tbl.IsDescription):
    """
    As there are no foreign keys in HDF5 we use the block_id to hold together
    different records for the same MAF block. This holds the coordinates and is 
    thus the main table for PyMAF.
    """
    linkv_rownum = tbl.UInt32Col(pos=1)
    
    chrom = tbl.StringCol(32, pos=2)
    start = tbl.UInt32Col(pos=3)
    end   = tbl.UInt32Col(pos=4)
    minus = tbl.BoolCol(pos=5)

    
class MAFPads(tbl.IsDescription):
    """
    Wraps the padding information on how a block connects to neighboring blocks in
    each species (large gaps, etc.). currently this data is not used.
    """
    pre_code  = tbl.StringCol(1, pos=1)
    pre_num   = tbl.UInt32Col(pos=2)
    post_code = tbl.StringCol(1, pos=3)
    post_num  = tbl.UInt32Col(pos=4)
    

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
        
        self.h5.root.scores.flush()
        self.h5.root.seqs.flush()
        self.h5.root.coord_rows_vector.flush()
        self.h5.root.seq_rows_vector.flush()
        
        if not species_names:
            species_names = self.species_names

        for species in species_names:
            t = getattr(self.h5.root.coords, species, None)
            t.flush()
            #else:
                #self.logger.warning("no coords table for '{0}'".format(species) )
            
            #t = getattr(self.h5.root.pads, species, None)
            #if t: 
                ##self.logger.debug("flushing /pads/{0}".format(species))
                #t.flush()
            #else:
                #self.logger.warning("no pads table for '{0}'".format(species) )
                
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
                cols.start.create_csindex()
                cols.end.create_csindex()
            else:
                self.logger.warning("no coords table for '{0}'".format(species) )
            
        t1 = time()
        self.logger.debug("indexing completed in {0:.1f}ms".format((t1-t0)*1000.) )


    def create_from_gz(self, path, tree_path="", regular_flushes=10000, blosc_max_threads=8, complevel=0):
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
        names_table = self.h5.create_vlarray(self.h5.root, 'species_names', atom = tbl.VLStringAtom(), title = "all species names mentioned in the MAF file in the correct order")

        if not tree_path:
            tree_path = find_tree(path)

        self.species_names = species_list_from_tree(tree_path)
        n_species = len(self.species_names)
        self.species_index = {}
        for i,species in enumerate(self.species_names):
            self.species_index[species] = i
            names_table.append(species.encode('ascii'))

        names_table.flush()
        self.logger.debug("stored {0} species names".format(n_species))
        
        seqs   = self.h5.create_vlarray(self.h5.root, 'seqs', atom = tbl.VLStringAtom(), title = "all sequences in all MAF blocks", filters=filters)
        scores = self.h5.create_earray(self.h5.root, 'scores', atom = tbl.Float32Atom(), shape=(0,), title = "alignment scores for each block", filters=filters, expectedrows=1000000)
        self.h5.create_group(self.h5.root, 'coords' )
        
        coords_tables = {}
        coords_curs = {}
        for species in self.species_names:
            table = self.h5.create_table(self.h5.root.coords, species, MAFCoords, "coordinates for species '{0}'".format(species), expectedrows=1000000)
            coords_tables[species] = table
            coords_curs[species] = table.row
            table.flush()

        #pads_tab
        #pads        = self.h5.create_table(self.h5.root, 'pads', MAFPads, "block padding info for each species", expectedrows=1000000)
        
        coord_rows_vector = self.h5.create_earray(self.h5.root, 'coord_rows_vector', atom = tbl.Int32Col(), shape=(0,n_species), title = "vector with row numbers of coords record for each species", filters=filters, expectedrows=1000000)
        seq_rows_vector   = self.h5.create_earray(self.h5.root, 'seq_rows_vector', atom = tbl.Int32Col(), shape=(0,n_species), title = "vector with row numbers of sequences for each species", filters=filters, expectedrows=1000000)
        #pad_rows_vector   = self.h5.create_earray(self.h5.root, 'pad_rows_vector', atom = tbl.Int64Col(shape=n_species), title = "vector with row numbers of pads record for each species", filters=filters, expectedrows=1000000)        
        T0 = time()
        T_last = T0

        curr_coords = {}
        curr_seqs = {}
        
        from gzip import GzipFile
        for maf_line in GzipFile(path):
            if not maf_line.strip():
                continue
            
            if maf_line.startswith('a'):
                scores.append([(float(maf_line.split('=')[1]))])
                if curr_coords:
                    cvec = np.array([curr_coords.get(species,-1) for species in self.species_names])
                    svec = np.array([curr_seqs.get(species,-1) for species in self.species_names])
                    coord_rows_vector.append([cvec])
                    seq_rows_vector.append([svec])
                    
                    curr_coords = {}
                    curr_seqs = {}

                if scores.nrows and (scores.nrows % regular_flushes == 0):
                    self.flush()

                    T = time()
                    dT = T - T_last
                    T_last = T
                    kbps = (regular_flushes / 1000.) / dT
                    
                    self.logger.debug("processed {0:.0f}k MAF blocks ({1:.1f}k blocks per second)".format(scores.nrows/1000., kbps) )
                
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

                curr_coords[species] = coords_tables[species].nrows
                n = coords_tables[species].nrows
                coords_tables[species].append( [(coord_rows_vector.nrows, chrom, start, end, (strand == '-'))] )
                assert coords_tables[species].nrows == n + 1 # need flush??
                
                                
                # store the alignment row in the VLStringArray. 
                # the link is through keeping the seqs current row number
                curr_seqs[species] = seqs.nrows
                seqs.append(seq.encode('ascii'))

            elif maf_line.startswith('i'):
                continue
                #parts = re.split(r'\s+',maf_line)
                #loc,pre_code,pre_num,post_code,post_num = parts[1:6]
                #species,chrom = loc.split('.',1)
                #pads = pads_lookup[species]
                #pads.append( [(maf_block_id, pre_code, pre_num, post_code, post_num)] )
            elif maf_line.startswith('e'):
                continue # currently ignored
            else:
                print "ignoring unknown MAF line '{0}'".format(maf_line.strip())
                
        if curr_coords:
            cvec = np.array([curr_coords.get(species,-1) for species in self.species_names])
            svec = np.array([curr_seqs.get(species,-1) for species in self.species_names])
            coord_rows_vector.append([cvec])
            seq_rows_vector.append([svec])

        self.logger.info("done processing {0} MAF blocks in {1:.1f}sec.".format(scores.nrows, (time() - T0)) )
        self.build_indices()
        
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
        self.species_index = {}
        for i,s in enumerate(self.species_names):
            self.species_index[s] = i

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

        select_indices = [self.species_index[s] for s in select_species]
        t = getattr(self.h5.root.coords, ref_species)
        t3 = t0
        n_blocks = 0
        n_rows = 0
        for ref_hit in t.where("(chrom == ref_chrom) & (start < ref_end) & (end > ref_start)"):
            n_blocks += 1
            rows = []
            vec_row = ref_hit['linkv_rownum']
            t1 = time()
            self.logger.debug("ref_hit found in {0:.1f}ms.".format( (t1-t3)*1000.) )
            seq_rows = self.h5.root.seq_rows_vector[vec_row]
            coord_rows = self.h5.root.coord_rows_vector[vec_row]
            print seq_rows
            print coord_rows
            for species in select_species:
                i = self.species_index[species]
                if seq_rows[i] == -1:
                    # missing coverage in this species
                    continue
                
                coords = getattr(self.h5.root.coords,species)[coord_rows[i]]
                seq = self.h5.root.seqs[seq_rows[i]]
                
                rows.append( MAFRecord(species, coords, seq) )
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
