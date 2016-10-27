#!/usr/bin/env python
import os,sys,re
import byo.gene_model
from glob import glob
import logging
import numpy as np

FORMAT = '%(asctime)-20s\t%(levelname)s\t%(name)-25s\t%(message)s'
logging.basicConfig(level=logging.DEBUG,format=FORMAT)
logger = logging.getLogger('find_seed')

from byo.gene_model import transcripts_from_UCSC
from byo.systems import hg19, mm10
from byo.io import fasta_chunks
from byo import rev_comp

let7 = "UGAGGUAGUAGGUUGUAUAGUU".replace('U','T')
seed8 = rev_comp(let7[1:9])
seed6 = rev_comp(let7[1:7])

import PyMAF as pm
#genome_path = "/scratch/data/genomes"
genome_path = "http://localhost:8000"
maf_path = "/scratch/data/maf/mm10_60way/maf"
MAF_track = pm.get_track(genome_path, maf_path, 'mm10', new_chrom_flush=True, in_memory=False, muscle_iterations = 8)

print "seeds", seed8, seed6
def get_seed_alignment(chain, species = []):
    alignments = []
    for exon in chain.exons:
        logger.debug('retrieving exon {exon.chrom}:{exon.start}-{exon.end}:{exon.sense} length={l}'.format(exon=exon, l = exon.end - exon.start))
        aln = MAF_track.get_oriented(exon.chrom,exon.start,exon.end,exon.sense, select_species = species)
        alignments.append(aln)
        #print aln, type(aln)
        logger.debug('got {0} columns of alignment of {1} species'.format(aln.n_cols,len(aln.species)))
        
    # concatenate alignments, if more than one exon
    aln = alignments[0]
    if len(alignments) > 1:
        for a in alignments[1:]:
            aln = aln + a
    
    #aln.headers[0] += " from {chain.exon_count} exons of {chain.name} gene_name={chain.gene_id} n_cols = {aln.n_cols}".format(**locals())
    return aln

def test_seed(aln, seed, spc):
    if spc in aln.nogaps:
        if seed in aln.nogaps[spc].upper():
            return True
    return False
    

for tx in transcripts_from_UCSC(sys.stdin, system=mm10):
    if not tx.UTR3:
        continue

    chain = tx.UTR3
    logger.debug("looking for seeds in '{0}'".format(tx.UTR3.name))
    seq = chain.spliced_sequence.upper()
    #print "#",seq
    for M in re.finditer(seed8, seq):
        start, end = chain.map_block_from_spliced(*M.span())
        site = chain.cut(start-5, end + 5)
        #print "found seed", site, tx.name, tx.key
        aln = get_seed_alignment(site, species=['mm10','hg19'])

        if not aln:
            logger.warning("Received empty Alignment({chain.name} {chain.chrom}:{chain.start}-{chain.end}:{chain.sense}). Region not covered by MAF. Skipping.".format(chain=chain))
            continue
        
        if not test_seed(aln, seed6, 'hg19'):
            context = chain.cut(start-6, end + 6)
            aln_detail = get_seed_alignment(context, species=[])
            
            present_mask = np.array([test_seed(aln_detail, seed6, spc) for spc in aln_detail.species])
            n_present = present_mask.sum()
            if n_present > 3:
                print "seed match at {site.key} interesting!".format(site=site)
                print "present in {0} species, but not human!".format(n_present)
                print aln_detail
                print present_mask

