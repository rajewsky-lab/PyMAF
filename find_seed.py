import pymaf
import byo.gene_model
from byo.io import fasta_chunks
from glob import glob

unal_path = "/scratch/transcript_test/*.UTR3.fa"
let7a = "uGAGGUAGUagguuguauaguu"
seed = "ACTACCTC"
minimal = "TACCTC"

for fname in glob(unal_path):
    species_hit = set()
    species_min = set()
    maxlen = 0
    for fa_id,seq in fasta_chunks(file(fname)):
        seq = seq.upper()
        maxlen = max(maxlen,len(seq))
        if seed in seq:
            species_hit.add(fa_id.split()[0])
        if minimal in seq:
            species_min.add(fa_id.split()[0])
    
    if len(species_hit) < 5:
        continue
    
    if not 'hg19' in species_min:
        flags = "CANDIDATE"
        print fname, flags, sorted(species_hit), maxlen

        
    
    

