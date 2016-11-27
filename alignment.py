"""
Part of the PyMAF package. The Alignment() class wraps multiple-species FASTA 
data and provides useful high-level functionality such as 'lift-over'-style
coordinate mapping between species, re-doing the alignment with MUSCLE, etc..
"""

__author__ = "Marvin Jens"
__copyright__ = "Copyright 2016, Massachusetts Institute of Technology"
__credits__ = ["Marvin Jens"]
__license__ = "MIT"
__version__ = "0.7"
__email__ = "mjens@mit.edu"
__status__ = "beta"

import logging
import numpy as np
from byo.io import fasta_chunks

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

    def lift_over(self, spc_A, pos_A, spc_B):
        """
        finds coordinate in spc_B that corresponds to coordinate pos_A in 
        species spc_A. Coordinates (for now) have to be relative to the start 
        position of the aligned sequences.
        """
        mfa_pos = self.spc_to_mfa[spc_A][pos_A]
        return self.mfa_to_spc[spc_B][mfa_pos]

    def __len__(self):
        """
        returns n_cols. Also makes things like 'if not aln:' work.
        """
        return self.n_cols
    
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

    def __getitem__(self, spc):
        """
        Return the alignment row corresponding to the requested species
        """
        #TODO: make this return columns if spc is numerical
        return self.by_species[spc]


    def columns(self, start=None, end=None, select_species=[]):
        """
        Generator. Instead of rows, iterates over columns (each yielded as a list)
        """
        if select_species:
            aln = self.subset(select_species)
        else:
            aln = self

        if start == None:
            start = 0
        if end == None:
            end = aln.n_cols
            
        for j in range(start, end):
            col = [seq[j] for seq in aln]
            yield col


    def variability(self, start=None, end=None, select_species=[]):
        """
        Returns the "fuzziness" of each column on a scale between 0 and 1
        by scaling the Shannon entropy of A,C,G,T,- frequencies to the 
        maximum occuring at uniform frequency. H/H_max
        """
        from collections import Counter
        
        def entropy(seq):
            """Shannon entropy of 'ACGT-' frequencies"""
            cnt = Counter(seq)
            counts = np.array(cnt.values(), dtype=float)
            
            p = counts / counts.sum()
            H = - (p * np.log2(p)).sum()
            return H

        max_H = - np.log2(0.2)
        return [entropy(c)/max_H for c in self.columns(start, end, select_species)]

        
    def __iter__(self):
        """
        Iterate over alignment rows.
        """
        for s in self.species:
            yield self[s]


    def subset(self, select_species):
        # start with a new, empty alignment
        new = Alignment("",ref = self.ref, acc = self.acc)

        # only keep those selected species that are actually there
        new.species = [s for s in select_species if s in self.species]
        new.headers = [self.headers[self.species_index[s]] for s in new.species]

        # remove gap only ('-') columns from the alignment
        drop = []
        for i, col in enumerate(zip(*[self.by_species[s] for s in new.species])):
            if set(col) == set('-'):
                drop.append(i)

        def drop_gaps(s):
            l = list(s)
            for d in drop[::-1]:
                l.pop(d)
            return "".join(l)

        for i,s in enumerate(new.species):
            new.nogaps[s] = self.nogaps[s]
            new.by_species[s] = drop_gaps(self.by_species[s])
            new.species_index[s] = i
            new._init_mappings(s)

        new.n_cols = self.n_cols
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


    def MUSCLE(self,n_iter=16):
        """
        Crazy little wrapper that builds a full multiple species alignment with
        MUSCLE (Edgar, R.C. Nucleic Acids Res 32(5), 1792-97).
        """

        from subprocess import Popen,PIPE
        # keep correct order for later
        if not self.species:
            return Alignment("")
        
        # pipe through MUSCLE
        muscle = Popen(['muscle','-maxiters',str(n_iter),'-in','/dev/stdin'],stdin=PIPE,stdout=PIPE,stderr=PIPE)
        mfa = "\n".join(['>{spc}\n{seq}'.format(spc=spc, seq=self.nogaps[spc]) for spc in self.species])
        self.logger.info('passing {size:.2f}kb of sequence ({n} species) through MUSCLE (maxiters={maxiters})'.format(n=len(self.species), size=len(mfa)/1000., maxiters=n_iter) )

        stdout,stderr = muscle.communicate(mfa)
        if muscle.returncode != 0:
            self.logger.error("MUSCLE died with error {0}. stderr='{1}'".format(muscle.returncode, stderr))
            import random
            uniqname = 'breaks_MUSCLE_{0:030x}.fa'.format(random.randrange(16**30))
            file(uniqname,"w").write(mfa)
            self.logger.error("saved offending MFA as '{0}'".format(os.path.abspath(uniqname)) )
            
            return Alignment("")
        else:
            # sanitize output and re-order
            unordered = {}
            for fa_id,seq in fasta_chunks(stdout.split('\n')):
                unordered[fa_id] = seq
            
            return Alignment("\n".join([">%s\n%s" % (species,unordered[species]) for species in self.species]), ref=self.ref, acc=self.acc)
        
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

