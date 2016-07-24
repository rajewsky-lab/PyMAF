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
maf_path = "/home/mjens/data/maf/hg19_100way/maf"

genome_provider = GenomeProvider(genome_path)

MAF_track = Track(
    maf_path,
    MAFBlockMultiGenomeAccessor,
    sense_specific=False,
    genome_provider=genome_provider,
    system='hg19'
)

#aln = MAF_track.get_oriented('chr19',8027068,8027687,'-', species=['hg19','mm10'])
aln = MAF_track.get_oriented('chrY',14971391,14971660,'+', species=['hg19','canFam3','mm10'])
print aln

r_start, r_end = aln.refine_region('canFam3', 35672500,35672529)
refine = MAF_track.get_oriented('chrY',r_start,r_end,'+')
print refine
print refine.MUSCLE()



    

