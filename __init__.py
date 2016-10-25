"""
Part of the PyMAF package. This file holds useful high-level representations of
multiple species alignments and the top-level interface to query the MAF-block
database and obtain alignments for arbitrary regions
"""

__author__ = "Marvin Jens"
__copyright__ = "Copyright 2016, Massachusetts Institute of Technology"
__credits__ = ["Marvin Jens"]
__license__ = "MIT"
__version__ = "0.7"
__email__ = "mjens@mit.edu"
__status__ = "beta"


import os,re,sys
from byo import rev_comp,complement
from byo.track import Track, Accessor
from byo.io import fasta_chunks
import logging

from PyMAF.maf_hdf5 import MAFBlockDB
from byo.io.genome_accessor import GenomeCache
import urllib2

class RemoteCache(object):
    def __init__(self, urlbase):
        self.urlbase = urlbase
        self.logger = logging.getLogger("RemoteCache('{0}')".format(urlbase))

    def __getitem__(self, genome):
        self.logger.debug('getitem {0}'.format(genome))
            
        class Proxy(object):
            no_data = False
            def get_oriented(this, chrom, start, end, strand):
                url = "{base}/{genome}/{chrom}:{start}-{end}{strand}".format(
                    base = self.urlbase, 
                    genome=genome, 
                    chrom=chrom, 
                    start=start, 
                    end=end, 
                    strand=strand
                )
                self.logger.debug('requesting "{0}"'.format(url))
                res = urllib2.urlopen(url)
                return res.read()

        return Proxy()
        
class Alignment(object):
    def __init__(self, mfa, ref=None, acc=None):
        self.logger = logging.getLogger('pymaf.Alignment')
        self.nogaps = {}
        self.by_species = {}
        self.species = []
        self.headers = []
        self.spacers = []
        self.species_index = {}
        self.acc = acc

        for i,(fa_id,seq) in enumerate(fasta_chunks(mfa.split('\n'))):
            self.headers.append(fa_id)
            species = fa_id.split()[0]
            self.species.append(species)
            self.nogaps[species] = seq.replace('-','')
            self.by_species[species] = list(seq)
            self.species_index[species] = i

        self.ref = ref
        if self.ref == None and self.species:
            # assume first species is ref:
            self.ref = self.species[0]

        self.spc_to_mfa = {}
        self.mfa_to_spc = {}
        for spc in self.species:
            self._init_mappings(spc)
        
        self.n_cols = len(self.by_species.get(self.ref,""))

    def lift_over(self, spc_A, pos_A, spc_B):
        """
        finds coordinate in spc_B that corresponds to coordinate pos_A in 
        species spc_A. Coordinates (for now) have to be relative to the start 
        position of the aligned sequences.
        """
        mfa_pos = self.spc_to_mfa[spc_A][pos_A]
        return self.mfa_to_spc[spc_B][mfa_pos]

    def _init_mappings(self,species):
        #self.logger.debug('_init_mappings()')
        x_spc = 0
        x_mfa = 0
        
        self.mfa_to_spc[species] = {}
        self.spc_to_mfa[species] = {}

        row = self.by_species[species]
        while x_mfa < len(row):
            self.mfa_to_spc[species][x_mfa] = x_spc
            if row[x_mfa] != '-':
                self.spc_to_mfa[species][x_spc] = x_mfa
                x_spc += 1
            x_mfa += 1

    def find_closest_inframe_ATG(self,species,spc_start):
        seq = self.nogaps[species]
        L = len(seq)
        frame = spc_start % 3
        
        all_codons = []
        i = 0
        while i+frame < L-3:
            j = i+frame
            if seq[j:j+3] == 'ATG':
                all_codons.append((abs(j-spc_start),j))
            i += 3

        if not all_codons:
            return -1
        else:
            min_dist,start = sorted(all_codons,reverse=True)[0]
            return start
            
        
    def highlight_ORF(self,ref_start,species=None):
        if species == None:
            species = self.ref

        spc_start = self.mfa_to_spc[species][self.spc_to_mfa[self.ref][ref_start]]
        codon = self.nogaps[species][spc_start:spc_start+3]
        if not codon == 'ATG':
            spc_start = self.find_closest_inframe_ATG(species,spc_start)
            
        if spc_start < 0:
            # no ORF found in this species
            return
        
        orf_i = 0
        L = len(self.nogaps[species])
        
        while codon.upper() not in ['TAA','TGA','TAG'] and orf_i + spc_start < L-3:
            #print species,orf_i,spc_start,orf_i+spc_start,L,codon
            p1 = self.spc_to_mfa[species][spc_start + orf_i]
            p2 = self.spc_to_mfa[species][spc_start + orf_i + 1]
            p3 = self.spc_to_mfa[species][spc_start + orf_i + 2]
            
            self.by_species[species][p2] = self.by_species[species][p2].lower()
            self.by_species[species][p3] = self.by_species[species][p3].lower()
        
            orf_i += 3
            codon = self.nogaps[species][spc_start + orf_i:spc_start + orf_i + 3]

    def MUSCLE(self,n_iter=16):
        """
        Crazy little wrapper that builds a full multiple species alignment with
        MUSCLE (Edgar, R.C. Nucleic Acids Res 32(5), 1792-97).
        """

        from subprocess import Popen,PIPE
        # keep correct order for later
        if not self.species:
            return ""
        
        # pipe through MUSCLE
        muscle = Popen(['muscle','-quiet','-maxiters',str(n_iter),'-in','/dev/stdin'],stdin=PIPE,stdout=PIPE)
        mfa = "\n".join(['>{spc}\n{seq}'.format(spc=spc, seq=self.nogaps[spc]) for spc in self.species])
        self.logger.info('passing {size:.2f}kb of sequence ({n} species) through MUSCLE (maxiters={maxiters})'.format(n=len(self.species), size=len(mfa)/1000., maxiters=n_iter) )
        stdout,stderr = muscle.communicate(mfa)

        # sanitize output and re-order
        unordered = {}
        for fa_id,seq in fasta_chunks(stdout.split('\n')):
            unordered[fa_id] = seq
        
        return Alignment("\n".join([">%s\n%s" % (species,unordered[species]) for species in self.species]), ref=self.ref, acc=self.acc)
        

    def __str__(self):
        buf = []
        for species,header in zip(self.species,self.headers):
            buf.append(">{header}".format(**locals()))
            row = list(self.by_species[species])
            for i,spc in enumerate(self.spacers):
                row.insert(spc+i,'|')
            buf.append("".join(row))
        return "\n".join(buf)

    def __add__(self,other):
        """
        overloaded "+" operator. Returns a new, concatenated alignment of 
        species that are covered in both input Alignments.
        """
        if other.ref != self.ref:
            raise ValueError("can't concatenate Alignment() instances with different reference species! '%s' != '%s'" % (self.ref,other.ref))

        # start with a new, empty alignment
        new = Alignment("",ref = self.ref, acc = self.acc)

        # only keep species that are in both blocks for now
        new.species = [s for s in self.species if s in other.species]
        new.headers = ["%s__plus__%s" % (self.headers[self.species_index[s]],other.headers[other.species_index[s]]) for s in new.species]
        
        for i,s in enumerate(new.species):
            new.nogaps[s] = self.nogaps[s] + other.nogaps[s]
            new.by_species[s] = self.by_species[s] + other.by_species[s]
            new._init_mappings(s)
            new.species_index[s] = i

        new.n_cols = self.n_cols + other.n_cols
        new.spacers.append(self.n_cols)

        return new

    def refine_region(self, species, start, end):
        """
        Talks to the upstream MAFBlockMultiGenomeAccessor object to find the blocks within
        the start and end coordinates in the requested species. This works by assuming that 
        the relevant block(s) can be found between the original, reference species start and 
        end coordinate determined blocks and is carried out by bisecting.
        EXPERIMENTAL!
        """
        self.logger.debug(self.headers[0])
        M = re.match(r"(?P<species>\S+) (?P<chrom>\S+)\:(?P<start>\d+)\-(?P<end>\d+)(?P<strand>[\+,\-])",self.headers[0])
        d = M.groupdict()
        self.logger.debug(d)
        
        M_spc = re.match(r"(?P<species>\S+) (?P<chrom>\S+)\:(?P<start>\d+)\-(?P<end>\d+)(?P<strand>[\+,\-])",self.headers[self.species_index[species]])
        s = M_spc.groupdict()
        return self.acc.find_inside_region_by_species(d['chrom'],int(d['start']), int(d['end']), d['strand'], species, s['chrom'], start, end)
        

class MAFCoverageCollector(object):
    """
    A helper class to aggregate the information from multiple MAF blocks
    spanning a reference region. The main task is to keep track of the
    start and end coordinates of the orthologous sequences in other
    species.
    """
    def __init__(self, ref, ref_start, ref_end, genome_provider, excess_threshold=2, lack_threshold=0.1, min_len=1, species=[]):
        from collections import defaultdict
        self.logger = logging.getLogger('pymaf.MAFCoverageCollector')
        self.excess_threshold = excess_threshold
        self.lack_threshold = lack_threshold
        self.min_len = min_len
        self.ref = ref
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.genome_provider = genome_provider
        self.left_col = None
        self.right_col = None
        self.species_min_max = {}
        self.species_intervals = defaultdict(list)
        self.species_strands = defaultdict(set)
        self.species_chroms = defaultdict(set)
        self.species_left_adjust = defaultdict(int)
        self.species_right_adjust = defaultdict(int)
        self.species = []
        self.select_species = set(species)
        self.logger.debug("select_species = {0}".format(species))
        self.ref_start_thisblock = None
        self.ref_end_thisblock = None
        
    def add_record(self,row):
        """
        this is used to feed data. Changes in MAF block are detected by change
        in ref_start_thisblock. This requires that the reference species is
        always passed first! Currently not checked, but silently assumed!!!
        """
        #self.logger.debug("adding '{0}', self.ref={3} self.ref_start={1} self.ref_end={2}".format(species, self.ref_start, self.ref_end, self.ref))
        if row.species == self.ref:
            # we are looking at the reference species!
            self.ref_start_thisblock = row.start
            self.ref_end_thisblock = row.end

            # no ref_start given, use bounds of this block!
            if self.ref_start == None:
                self.ref_start = row.start
            if self.ref_end == None:
                self.ref_end = row.end

            if row.start <= self.ref_start:
                # how many columns do we have to skip until we reach 
                # the desired start position (gaps!)?
                delta = self.ref_start - row.start
                j = 0
                while delta > 0:
                    if row.seq[j] != '-':
                        delta -= 1
                    j += 1

                self.left_col = j
                #print "ref",self.left_col,seq[:self.left_col+1],seq[self.left_col:]

            if row.end >= self.ref_end:
                # same for the end, what column corresponds to ref_end?
                delta = row.end - self.ref_end
                j = len(row.seq)-1
                while delta > 0:
                    if row.seq[j] != '-':
                        delta -= 1
                    j -= 1
                self.right_col = j+1
        else:
            # for other species, how many genomic positions do we need to skip 
            # until we are at left_col?
            if self.left_col != None and self.ref_start_thisblock <= self.ref_start:
                # count non-gap positions before left_col/ref_start
                self.species_left_adjust[row.species] = self.left_col - row.seq[:self.left_col].count('-')
                
            if self.right_col != None and self.ref_end_thisblock >= self.ref_end:
                # count non-gap positions after the ref_end/right_col
                self.species_right_adjust[row.species] = (len(row.seq)-self.right_col) - row.seq[self.right_col:].count('-') 
                
        # record the encountered species, their contigs and strands
        strand = row.sense
        if not row.species in self.species_strands:
            self.species.append(row.species)

        self.species_strands[row.species].add(strand)
        self.species_chroms[row.species].add(row.chrom)
        
        # build a "maximal cover", the largest chunk of orthologous sequence
        # according to the MAF blocks
        if not row.species in self.species_min_max:
            self.species_min_max[row.species] = (row.start,row.end)

        S,E = self.species_min_max[row.species]
        self.species_min_max[row.species] = (min(row.start,S),max(row.end,E))
        self.species_intervals[row.species].append( (row.start, row.end) )

    def get_sequences(self,sense):
        """
        This is a generator, yielding (species,chrom,start,end,strand,seq)
        for as many species as it can. Actual genome access is deferred to
        a GenomeProviderCache instance passed to the constructor.
        """
        ref_len = self.ref_end - self.ref_start

        for species in self.species:
            if not species in self.select_species and self.select_species:
                #self.logger.debug("skipping species '{0}' because it was not selected".format(species))
                continue

            if len(self.species_chroms[species]) > 1:
                chrom_list = ','.join(sorted(self.species_chroms[species]))
                self.logger.warning('sequence is split across different contigs/chroms {chrom_list} in species {species}'.format(**locals()))
                continue
            else:
                chrom = list(self.species_chroms[species])[0]
                
            if len(self.species_strands[species]) > 1:
                self.logger.warning('sequence is split across different strands in species {species}'.format(**locals()))
                continue
            else:
                strand = list(self.species_strands[species])[0]

            start,end = self.species_min_max[species]
            #print "getting",species,chrom,strand
            genome = self.genome_provider[species]
            if species == self.ref:
                start,end = self.ref_start,self.ref_end
            else:
                if strand == '+':
                    start += self.species_left_adjust[species]
                    end -= self.species_right_adjust[species]
                else:
                    # Dr. minus strikes again...
                    end -= self.species_left_adjust[species]
                    start += self.species_right_adjust[species]

            if end - start > ref_len* self.excess_threshold:
                self.logger.warning('homologous sequence exceeds {0}x times the reference sequence for {1}. skipping!'.format(self.excess_threshold, species))
                continue

            if end - start < ref_len* self.lack_threshold:
                self.logger.warning('skipping homologous sequence shorter than {0}x times the reference sequence for {1}. skipping!'.format(self.lack_threshold, species))
                continue

            if end - start < self.min_len:
                self.logger.warning('skipping homologous sequence of length {0} < min_len of {1} for {2}!'.format(end-start, self.min_len, species))
                continue

            #print "retrieving for species='{species}' chrom='{chrom}' sense='{strand}'".format(**locals())
            seq = genome.get_oriented(chrom.split('.')[0],start,end,strand).upper()
            if genome.no_data:
                self.logger.warning("skipping {species} due to missing genome".format(**locals()))
                continue
            else:
                if sense == '-':
                    # Dr. minus again. This time my own convention makes it 
                    # even more weird. But it makes sense if you think about
                    # it: get_data always preserves the order "left-to-right"
                    # in chromosome logic. If you expect rev_comp, use
                    # get_oriented() instead.
                    seq = complement(seq)

                yield species,chrom,start,end,strand,seq


class MAFBlockMultiGenomeAccessor(Accessor):
    """
    This class uses an indexed HDF5 version of the raw MAF blocks to find
    the MAF blocks overlapping the start and end of an arbitrary genomic
    span in O(log(N)) time. It then parses these MAF blocks, extracts the
    orthologous coordinates for the aligned species and retrieves the
    orthologous genomic sequences (with the help of MAFCoverageCollector).
    """

    def __init__(self,maf_path,chrom,sense,sense_specific=False,empty="",genome_provider=None,excess_threshold=5.,lack_threshold=.0,min_len=0, aln_class=Alignment, muscle_iterations=0, in_memory=False, **kwargs):
        super(MAFBlockMultiGenomeAccessor,self).__init__(maf_path,chrom,sense,sense_specific=False,**kwargs)
        self.logger = logging.getLogger("pymaf.MAFBlockMultiGenomeAccessor")

        self.maf_path = maf_path
        self.empty = empty
        self.reference = self.system
        self.aln_class = aln_class
        self.genome_provider = genome_provider
        self.excess_threshold = excess_threshold
        self.lack_threshold = lack_threshold
        self.min_len = min_len
        self.muscle_iterations = muscle_iterations
        h5path = os.path.join(self.maf_path,chrom+".maf.h5")
        self.maf_file = MAFBlockDB(h5path, in_memory=in_memory)
        

    def find_inside_region_by_species(self, ref_chrom, ref_start, ref_end, ref_sense, species, spc_chrom, spc_start, spc_end, select_species = []):
        # TODO implement binary search for relevant maf block such that a MUSCLE alignment can be requested for a smaller region.
        self.logger.debug("collecting candidate blocks inside {ref_start}-{ref_end}".format(ref_start=ref_start, ref_end=ref_end) )
        
        start_i = self.data[ref_start]
        end_i = self.data[ref_end]
        
        block_starts = set()
        for comb in self.index[start_i:end_i+1]:
            block_starts |= set(comb.split(','))

        block_starts = sorted([int(bs,16) for bs in block_starts])
        n_rows = len(block_starts)
        self.logger.debug('find_inside_region_by_species({species} {spc_start}-{spc_end}) binary-searching {n_rows} blocks'.format(**locals()) )        

        # these are the only two we need to consider
        select_species = [self.reference, species]
        
        def search_block(block_starts, pos):
            
            def recurse(block_starts, pos):
                #print "recursion:", block_starts, pos
                if not block_starts:
                    return None

                if len(block_starts) == 1:
                    blk = self.load_block(block_starts[0], species=select_species)
                    s,e = blk.species_min_max.get(species, (-1,-1))
                    #print "checked block", block_starts, s,e, "?", pos
                    if (s == e) or blk.species_chroms[species] != set([spc_chrom]):
                        #print "ignore", s,e, blk.species_chroms[species], spc_chrom
                        return "ignore",block_starts[0],blk
                    elif s <= pos <= e:
                        #print "hit"
                        return "hit",block_starts[0],blk
                    elif e < pos:
                        #print "left"
                        return "left",block_starts[0],blk
                    elif s < pos:
                        #print "right"
                        return "right",block_starts[0],blk
                else:
                    #print "splitting interval"
                    mid = len(block_starts)/2
                    code, blk_start, blk = recurse([block_starts[mid]], pos)
                    #print "result=", code, blk_start
                    if code == "hit":
                        return code, blk_start, blk
                    elif code == "left":
                        return recurse(block_starts[mid:], pos)
                    elif code == "right":
                        return recurse(block_starts[:mid], pos)
                    elif code == "ignore":
                        block_starts.pop(mid)
                        return recurse(block_starts, pos)
        
            res = recurse(list(block_starts), pos)
            #print "search received", res
            if not res:
                return None
            
            code, blk_start, blk = res
            if code == 'hit' : 
                return blk_start, blk
            else:
                return None
            
        #print "-----searching start",block_starts
        start_search = search_block(block_starts, spc_start)
        #print "-----searching end",block_starts
        end_search   = search_block(block_starts, spc_end)

        #print start_search
        #print end_search

        if not (start_search and end_search):
            return 
        
        start_row, start_cov = start_search
        end_row, end_cov = end_search
        
        start_mfa = "\n".join([">{species} {chrom}:{start}-{end}{strand}\n{seq}".format(**locals()) for species,chrom,start,end,strand,seq in start_cov.get_sequences(ref_sense)])
        start_aln = self.aln_class(start_mfa, ref=self.reference, acc=self)

        end_mfa = "\n".join([">{species} {chrom}:{start}-{end}{strand}\n{seq}".format(**locals()) for species,chrom,start,end,strand,seq in end_cov.get_sequences(ref_sense)])
        end_aln = self.aln_class(end_mfa, ref=self.reference, acc=self)

        start_ofs_spc = start_cov.species_min_max[species][0]
        start_ofs_ref = start_cov.species_min_max[self.reference][0]
        
        end_ofs_spc = end_cov.species_min_max[species][0]
        end_ofs_ref = end_cov.species_min_max[self.reference][0]
        
        ref_start = start_aln.lift_over(species, spc_start - start_ofs_spc, self.reference) + start_ofs_ref
        ref_end   =   end_aln.lift_over(species, spc_end   - end_ofs_spc,   self.reference) + end_ofs_ref

        return ref_start, ref_end

       
    def get_data(self, chrom, ref_start, ref_end, ref_sense, select_species = []):
        coverage = MAFCoverageCollector(
            self.reference,
            ref_start,
            ref_end,
            self.genome_provider,
            excess_threshold=self.excess_threshold,
            lack_threshold=self.lack_threshold,
            min_len=self.min_len,
            species=select_species
        )
        
        for block in self.maf_file.query_interval(self.reference, self.chrom, ref_start, ref_end, select_species):
            for row in block:
                #print row
                coverage.add_record(row)

        return list(coverage.get_sequences(ref_sense))

    def get_oriented(self,chrom,start,end,sense, select_species=[]):
        res = self.get_data(chrom,start,end,sense, select_species=select_species)
        if sense == '-':
            # get_data returned already complement(seq). Only need to reverse.
            res = [(species,chrom,start,end,strand,seq[::-1]) for species,chrom,start,end,strand,seq in res]

        mfa = "\n".join([">{species} {chrom}:{start}-{end}{strand}\n{seq}".format(**locals()) for species,chrom,start,end,strand,seq in res])
        aln = self.aln_class(mfa, ref=self.reference, acc=self)
        print "HERE:",aln
        if self.muscle_iterations:
            return aln.MUSCLE(n_iter=self.muscle_iterations)
        else:
            return aln
    
    def get_dummy(self,chrom,start,end,sense):
        return []
    
    def flush(self):
        self.maf_file.close()

def get_track(genome_path, maf_path, reference_species, excess_threshold=2, lack_threshold=.0, min_len=1, new_chrom_flush=False, in_memory=False, muscle_iterations = 0):
    print ">>>> GET TRACK"
    if genome_path.startswith('http://'):
        genome_provider = RemoteCache(genome_path)
        print ">>>.REMOTE"
    else:
        genome_provider = GenomeCache(genome_path)

    MAF_track = Track(
        maf_path,
        MAFBlockMultiGenomeAccessor,
        auto_flush=new_chrom_flush,
        sense_specific=False,
        genome_provider=genome_provider,
        excess_threshold=excess_threshold,
        lack_threshold=lack_threshold,
        min_len=min_len,
        system=reference_species,
        in_memory = in_memory,
        muscle_iterations = muscle_iterations,
    )
    return MAF_track
    
