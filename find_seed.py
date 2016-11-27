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
from byo import systems
from byo.io import fasta_chunks
from byo import rev_comp

let7 = "UGAGGUAGUAGGUUGUAUAGUU".replace('U','T')
seed6 = rev_comp(let7[1:7])
seed7 = rev_comp(let7[1:8])
seed8 = rev_comp(let7[1:9])


import PyMAF as pm
#genome_path = "/scratch/data/genomes"
genome_path = "http://localhost:8000"
maf_path = "/home/mjens/maf/mm10_60way"
MAF_track = pm.get_track(genome_path, maf_path, 'mm10', new_chrom_flush=True, in_memory=False, muscle_iterations = 0, lack_threshold=.1)
flank = 20

#print "# seeds", seed8, seed7, seed6
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
    
# do not scan the same region twice
scanned = set()

for tx in transcripts_from_UCSC(sys.stdin, system=byo.systems.mm10):
    if not tx.UTR3:
        continue

    chain = tx.UTR3
    logger.debug("looking for seeds in '{0}'".format(tx.UTR3.name))
    seq = chain.spliced_sequence.upper()
    #print "#",seq
    for M in re.finditer(seed6, seq):
        start, end = chain.map_block_from_spliced(*M.span())
        site = chain.cut(start-flank, end + flank)
        if site.key in scanned:
            continue
        
        scanned.add(site.key)
        
        #print "found seed", site, tx.name, tx.key
        aln = get_seed_alignment(site, species=['mm10','hg19'])
        #print aln
        if not aln:
            logger.warning("Received empty Alignment({chain.name} {chain.chrom}:{chain.start}-{chain.end}:{chain.sense}). Region not covered by MAF. Skipping.".format(chain=chain))
            continue
        
        if not test_seed(aln, seed6, 'hg19'):
            #print "# seed absent in human? getting more detail"
            context = chain.cut(start-flank, end + flank)
            aln_detail = get_seed_alignment(context, species=[]).MUSCLE()
            
            m6 = np.array([test_seed(aln_detail, seed6, spc) for spc in aln_detail.species])
            m7 = np.array([test_seed(aln_detail, seed7, spc) for spc in aln_detail.species])
            m8 = np.array([test_seed(aln_detail, seed8, spc) for spc in aln_detail.species])
            n_present = m6.sum()
            if n_present > 3:
                print ">seed_lost {chain.name} | {site.key_str} | n_present={n_present}".format(**locals())
                print "m6={0}".format(m6)
                print "m7={0}".format(m7)
                print "m8={0}".format(m8)
                print aln_detail

        #else:
            #print "# seed present in human"

