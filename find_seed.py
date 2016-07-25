import os,sys
import pymaf
import byo.gene_model
from byo.io import fasta_chunks
from glob import glob

#unal_path = "/scratch/transcript_test/*.UTR3.fa"
#let7a = "uGAGGUAGUagguuguauaguu"
#seed = "ACTACCTC"
#minimal = "TACCTC"

#for fname in glob(unal_path):
    #species_hit = set()
    #species_min = set()
    #maxlen = 0
    #for fa_id,seq in fasta_chunks(file(fname)):
        #seq = seq.upper()
        #maxlen = max(maxlen,len(seq))
        #if seed in seq:
            #species_hit.add(fa_id.split()[0])
        #if minimal in seq:
            #species_min.add(fa_id.split()[0])
    
    #if len(species_hit) < 5:
        #continue
    
    #if not 'hg19' in species_min:
        #flags = "CANDIDATE"
        #print fname, flags, sorted(species_hit), maxlen

import logging
FORMAT = '%(asctime)-20s\t%(levelname)s\t%(name)-25s\t%(message)s'
logging.basicConfig(level=logging.DEBUG,format=FORMAT)

from PyMAF.pymaf import *
genome_path = "/home/mjens/data/genomes"
maf_path = "/home/mjens/data/maf/mm10_60way/maf"

genome_provider = GenomeProvider(genome_path)

MAF_track = Track(
    maf_path,
    MAFBlockMultiGenomeAccessor,
    sense_specific=False,
    genome_provider=genome_provider,
    system='mm10'
)


from byo.gene_model import transcripts_from_UCSC
logger = logging.getLogger('find_seed')
from byo.systems import hg19, mm10
from byo import rev_comp

let7 = "UGAGGUAGUAGGUUGUAUAGUU".replace('U','T')
seed8 = rev_comp(let7[1:9])
seed6 = rev_comp(let7[1:7])

print "seeds", seed8, seed6
def get_seed_alignment(chain, species = []):
    alignments = []
    for exon in chain.exons:
        logger.debug('retrieving exon {exon.chrom}:{exon.start}-{exon.end}:{exon.sense} length={l}'.format(exon=exon, l = exon.end - exon.start))
        aln = MAF_track.get_oriented(exon.chrom,exon.start,exon.end,exon.sense, species = species)
        alignments.append(aln)
        logger.debug('got {0} columns of alignment of {1} species'.format(aln.n_cols,len(aln.species)))
        
    # concatenate alignments, if more than one exon
    aln = alignments[0]
    if len(alignments) > 1:
        for a in alignments[1:]:
            aln = aln + a
    
    if not aln.n_cols:
        logger.error("Received empty Alignment({chain.name} {chain.chrom}:{chain.start}-{chain.end}:{chain.sense}). Region not covered by MAF? Skipping.".format(chain=chain))
        return None

    #aln.headers[0] += " from {chain.exon_count} exons of {chain.name} gene_name={chain.gene_id} n_cols = {aln.n_cols}".format(**locals())
    return aln

def test_seed(aln, seed, spc):
    if spc in aln.nogaps:
        if seed in aln.nogaps[spc].upper():
            return True
    return False
    

for tx in transcripts_from_UCSC(sys.stdin,system=mm10):
    if not tx.UTR3:
        continue

    chain = tx.UTR3
    #print "looking for seeds"
    seq = chain.spliced_sequence.upper()
    #print "#",seq
    for M in re.finditer(seed8, seq):
        start, end = chain.map_block_from_spliced(*M.span())
        site = chain.cut(start-2, end + 2)
        print "found seed", site, tx.name
        aln = get_seed_alignment(site, species=['mm10','hg19'])
        
        print aln
        if test_seed(aln, seed6, 'hg19'):
            print "human conserved :("
        else:
            print "human let-7 seed match lost!"
            context = chain.cut(start-6, end + 6)
            aln_detail = get_seed_alignment(context, species=[])
            print aln_detail

