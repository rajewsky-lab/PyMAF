import logging
import sys
FORMAT = '%(asctime)-20s\t%(levelname)s\t%(name)-25s\t%(message)s'
logging.basicConfig(level=logging.DEBUG,format=FORMAT)

from PyMAF import get_track
from PyMAF.maf_hdf5 import MAFBlockDB
db = MAFBlockDB('/scratch/data/maf/mm10_60way/maf/chrM.maf.hdf5')
for x in db.query_interval('mm10', 'chrM', 0, 9):
    print ">>>block"
    for y in x:
        print str(y)
    
#for x in db.query_interval('mm10', 'chrM', 65, 100):
    #print ">>>block"
    #for y in x:
        #print str(y)

sys.exit(1)
maf = get_track('/scratch/data/genomes', '/scratch/data/maf/mm10_60way/maf', 'mm10')
print maf.get_oriented('chrM', 15, 50, '+')

print maf.get_oriented('chrM', 15, 50, '+')