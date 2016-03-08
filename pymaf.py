import os,re,sys
from byo import rev_comp,complement
from byo.track import Track
from byo.io.genome_accessor import GenomeAccessor
from byo.io.array_accessor import ArrayAccessor
from byo.io import fasta_chunks
import numpy as np
import logging

class GenomeProvider(object):
    def __init__(self,path):
        self.cached = {}
        self.path = path
        self.logger = logging.getLogger("pymaf.GenomeProvider")

    def __getitem__(self,name):
        if not name in self.cached:
            logging.debug("{0} not in genome cache".format(name))
            self.cached[name] = Track(self.path,GenomeAccessor,system=name,split_chrom='.')
        else:
            logging.debug("genome cache hit! {0}".format(name))

        return self.cached[name]

def run_through_MUSCLE(mfa,min_len=0):
    """
    Crazy little wrapper that builds a multiple species alignment with
    MUSCLE (Edgar, R.C. Nucleic Acids Res 32(5), 1792-97).
    """

    from subprocess import Popen,PIPE
    # keep correct order for later
    ordered_species = [fa_id for fa_id,seq in fasta_chunks(mfa.split('\n')) if len(seq) > min_len]
    if not ordered_species:
        return ""
    # pipe through MUSCLE
    muscle = Popen(['muscle','-quiet','-in','/dev/stdin'],stdin=PIPE,stdout=PIPE)
    stdout,stderr = muscle.communicate(mfa)

    # sanitize output and re-order
    unordered = {}
    for fa_id,seq in fasta_chunks(stdout.split('\n')):
        unordered[fa_id] = seq

    return "\n".join([">%s\n%s" % (species,unordered[species]) for species in ordered_species])
    
    
class Alignment(object):
    def __init__(self,mfa,ref=None):
        self.nogaps = {}
        self.by_species = {}
        self.species = []
        self.headers = []
        self.spacers = []
        self.species_index = {}

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

    def _init_mappings(self,species):
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
        new = Alignment("",ref = self.ref)

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


class MAFCoverageCollector(object):
    """
    A helper class to aggregate the information from multiple MAF blocks
    spanning a reference region. The main task is to keep track of the
    start and end coordinates of the orthologous sequences in other
    species.
    """
    def __init__(self,ref,ref_start,ref_end,genome_provider):
        from collections import defaultdict
        self.logger = logging.getLogger('pymaf.MAFCoverageCollector')
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
        
        self.ref_start_thisblock = None
        self.ref_end_thisblock = None
        
    def add_MAF_line(self,maf_line):
        """
        this is used to feed data. Changes in MAF block are detected by change
        in ref_start_thisblock. This requires that the reference species is
        always passed first! Currently not checked, but silently assumed!!!
        """
        #print maf_line.rstrip()
        if not maf_line.startswith('s'):
            return

        parts = re.split(r'\s+',maf_line)
        loc,start,size,strand,total,seq = parts[1:7]
        
        species,chrom = loc.split('.',1)
        start = int(start)
        size = int(size)
        total = int(total)
        if strand == '+':
            end = start + size
        else:
            start,end = total - (start + size),total - start
            # if you think this is sick, go tell the evil master MAF and his
            # sidekick Dr. minus, the inventor of the minus strand, to their 
            # faces. ;)
            
        if species == self.ref:
            # we are looking at the reference species!
            self.ref_start_thisblock = start
            self.ref_end_thisblock = end

            if start <= self.ref_start:
                # how many columns do we have to skip until we reach 
                # the desired start position (gaps!)?
                delta = self.ref_start - start
                j = 0
                while delta > 0:
                    if seq[j] != '-':
                        delta -= 1
                    j += 1

                self.left_col = j
                #print "ref",self.left_col,seq[:self.left_col+1],seq[self.left_col:]

            if end >= self.ref_end:
                # same for the end, what column corresponds to ref_end?
                delta = end - self.ref_end
                j = len(seq)-1
                while delta > 0:
                    if seq[j] != '-':
                        delta -= 1
                    j -= 1
                self.right_col = j+1
        else:
            # for other species, how many genomic positions do we need to skip 
            # until we are at left_col?
            if self.left_col != None and self.ref_start_thisblock <= self.ref_start:
                # count non-gap positions before left_col/ref_start
                self.species_left_adjust[species] = self.left_col - seq[:self.left_col].count('-')
                
            if self.right_col != None and self.ref_end_thisblock >= self.ref_end:
                # count non-gap positions after the ref_end/right_col
                self.species_right_adjust[species] = (len(seq)-self.right_col) - seq[self.right_col:].count('-') 
                
        # record the encountered species, their contigs and strands
        if not species in self.species_strands:
            self.species.append(species)

        self.species_strands[species].add(strand)
        self.species_chroms[species].add(chrom)
        
        # build a "maximal cover", the largest chunk of orthologous sequence
        # according to the MAF blocks
        if not species in self.species_min_max:
            self.species_min_max[species] = (start,end)

        S,E = self.species_min_max[species]
        self.species_min_max[species] = (min(start,S),max(end,E))
        self.species_intervals[species].append( (start,end) )


    def get_sequences(self,sense):
        """
        This is a generator, yielding (species,chrom,start,end,strand,seq)
        for as many species as it can. Actual genome access is deferred to
        a GenomeProvider instance passed to the constructor.
        """
        ref_len = self.ref_end - self.ref_start

        for species in self.species:
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

            if end - start > ref_len*5:
                self.logger.warning('homologous sequence exceeds five times the reference sequence for {0}. skipping!'.format(species))
                continue

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

class MAFBlockMultiGenomeAccessor(ArrayAccessor):
    """
    This class uses an mmap'ed, sparse lookup table and an index to find
    the MAF blocks overlapping the start and end of an arbitrary genomic
    span in constant time. It then parses these MAF blocks, extracts the
    orthologous coordinates for the aligned species and retrieves the
    orthologous genomic sequences (with the help of MAFCoverageCollector).
    """

    def __init__(self,maf_path,chrom,sense,sense_specific=False,dtype=np.uint32,empty="",genome_provider=None,muscle=False,aln_class=Alignment,**kwargs):
        super(MAFBlockMultiGenomeAccessor,self).__init__(maf_path,chrom,sense,dtype=dtype,sense_specific=False,ext=".comb_bin",**kwargs)
        self.logger = logging.getLogger("pymaf.MAFBlockMultiGenomeAccessor")
        self.logger.debug("mmap'ing '%s' lookup sparse-files and indices for chromosome %s" % (str(dtype),chrom))

        self.maf_path = maf_path
        self.empty = empty
        self.reference = self.system
        self.aln_class = aln_class
        self.muscle = muscle
        self.genome_provider = genome_provider        

        index_file = os.path.join(self.maf_path,chrom+".comb_idx")
        if self.load_index(index_file,empty):
            self.maf_file = file(os.path.join(self.maf_path,chrom+".maf"))
        else:
            self.logger.warning("could not find MAF index file '{0}'".format(index_file) )
            self.maf_file = None
           

    def load_index(self,fname,empty):
        self.logger.debug("loading index from '%s'" % fname)
        try:
            self.index = [empty] + [line.rstrip() for line in file(fname)]
        except IOError:
            self.logger.warning("Could not access '%s'. Switching to dummy mode (only empty)" % fname)
            self.index = [empty]
            self.get_data = self.get_dummy
            return False

        self.logger.debug("loaded %d feature combinations" % len(self.index))
        return True

    def get_data(self,chrom,start,end,sense):
        maf_starts = set()
        comb_code_indices = [self.data[start],self.data[end]]
        #comb_code_indices = self.data[start:end]
        #print comb_code_indices
        comb_codes = [self.index[i] for i in comb_code_indices if i]
        #print comb_codes
        for comb in set(comb_codes):
            maf_starts |= set(comb.split(','))
        #print maf_starts
        coverage = MAFCoverageCollector(self.reference,start,end,self.genome_provider)
        for m_start in sorted([int(m,16) for m in maf_starts]):
            self.maf_file.seek(m_start)
            for line in self.maf_file:
                #print start,end,line.rstrip()
                try:
                    coverage.add_MAF_line(line)
                except:
                    import traceback
                    self.logger.error("unhandled exception while parsing MAF-line '{0}'".format(line.rstrip()))
                    exc = traceback.format_exc()
                    logging.error(exc)

                if not line.strip():
                    break
                
        return list(coverage.get_sequences(sense))

    def get_oriented(self,chrom,start,end,sense):
        res = self.get_data(chrom,start,end,sense)
        if sense == '-':
            # get_data returned already complement(seq). Only need to reverse.
            res = [(species,chrom,start,end,strand,seq[::-1]) for species,chrom,start,end,strand,seq in res]

        mfa = "\n".join([">{species} {chrom}:{start}-{end}{strand}\n{seq}".format(**locals()) for species,chrom,start,end,strand,seq in res])
        if self.muscle:
            mfa = run_through_MUSCLE(mfa)
            
        aln = self.aln_class(mfa)
        return aln
    
    def get_dummy(self,chrom,start,end,sense):
        return []



def process_ucsc(src,system=None,**kwargs):
    from byo.gene_model import transcripts_from_UCSC
    from byo.protein import find_ORF

    for tx in transcripts_from_UCSC(src,system=system):
        alignments = []
        for exon in tx.exons:
            aln = MAF_track.get_oriented(exon.chrom,exon.start,exon.end,exon.sense)
            alignments.append(aln)

        aln = alignments[0]
        if len(alignments) > 1:
            for a in alignments[1:]:
                aln = aln + a
        
        n_cols = aln.n_cols
        #aln = aln + aln

        # double for circRNA
        aln.headers[0] += " from {tx.exon_count} exons of {tx.name} gene_name={tx.gene_id} n_cols = {n_cols}".format(**locals())
        #cds_start,cds_ens = tx.map_block_to_spliced(tx.CDS.start,tx.CDS.end)
        
        cds_start,cds_end = tx.map_block_to_spliced(tx.CDS.start,tx.CDS.end)
        ref_seq = aln.nogaps[aln.ref]
        if ref_seq[cds_start+1:cds_start+3].upper() == 'TG':
            # proper start-codon found?
            orf_start = cds_start
        else:
            # de novo search for longest ORF
            aa,orf_start,orf_end = find_ORF(ref_seq,len_thresh=10)

        if orf_start >= 0:
            for species in aln.species:
                aln.highlight_ORF(orf_start,species)
                
        yield aln, tx.name

def process_bed6(src,**kwargs):
    from byo.io.bed import bed_importer
    logging.info('process_bed6()')
    for bed in bed_importer(src):
        logging.info("...retrieving sequences for {bed.name} {bed.chrom}:{bed.start}-{bed.end}:{bed.strand}".format(bed=bed))
        aln = MAF_track.get_oriented(bed.chrom,bed.start,bed.end,bed.strand)
        yield aln, bed.name
        
def not_implemented(*argc,**kwargs):
    logging.error("INPUT FORMAT NOT IMPLEMENTED!")
    sys.exit(1)

if __name__ == '__main__':
    from optparse import *
    usage = """
    usage: %prog [options] <input_file.bed|gff|ucsc>
    """

    parser = OptionParser(usage=usage)
    parser.add_option("-S","--system",dest="system",type=str,default=None,help="model system/reference species (hg19|dm6|...)")
    parser.add_option("-M","--maf-path",dest="maf_path",type=str,default="",help="path to indexed MAF files (default='./')")
    parser.add_option("-G","--genome-path",dest="genome_path",type=str,default="./",help="path to genomes (default='<maf-path>/genomes')")
    parser.add_option("-o","--output-path",dest="output_path",type=str,default="./",help="path to write output files to. if you pass '-', it prints on stdout instead (default='./')")
    parser.add_option("","--muscle",dest="muscle",default=False,action="store_true",help="activate re-alignment through MUSCLE. Warning, this can take a long time for large sequences! (default=Off)")
    parser.add_option("","--debug",dest="debug",default=False,action="store_true",help="activate extensive debug output (default=Off)")
    parser.add_option("-i","--input-format",dest="input_format",default="bed6",choices=["bed6","gff","bed12","ucsc"],help='which format does the input have? ["bed6","gff","bed12","ucsc"] default is bed6')
    options,args = parser.parse_args()

    # prepare logging system
    FORMAT = '%(asctime)-20s\t%(levelname)s\t%(name)s\t%(message)s'
    if options.debug:
        logging.basicConfig(level=logging.DEBUG,format=FORMAT)
    else:
        logging.basicConfig(level=logging.INFO,format=FORMAT)

    if not options.genome_path:
        genome_path = os.path.join(options.maf_path,'genomes')
    else:
        genome_path = options.genome_path

    genome_provider = GenomeProvider(genome_path)
    MAF_track = Track(
        options.maf_path,
        MAFBlockMultiGenomeAccessor,
        sense_specific=False,
        genome_provider=genome_provider,
        muscle=options.muscle,
        system=options.system
    )

    if not args:
        src = sys.stdin
    else:
        src = file(args[0])
    
    handler = { 
        'bed6' : process_bed6,
        'bed12' : not_implemented,
        'gff' : not_implemented,
        'ucsc' : process_ucsc,
    }[options.input_format]
    
    for aln,name in handler(src, muscle=options.muscle, system=options.system):
        mfa = str(aln)

        if options.output_path == '-' or not options.output_path:
            print mfa
        else:
            out_file = os.path.join(options.output_path, "{0}.fa".format(name))
            logging.info("storing {name} in {out_file}".format(name=name, out_file=out_file))
            file(out_file,'w').write(mfa)

    