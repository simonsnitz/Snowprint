import sys
import os
import os.path
import shutil

d = os.getcwd()
os.chdir("..")


#works on Linux, but not windows. Uggh!
sys.path.insert(1, '/home/simon/git/tuSeek')

from acc2homologs import homologs2residueFrequency as a2r

#ramr
#acc = "WP_000113609"

#camr
acc = "BAA03510"

#seq = a2h.accID2sequence(acc)
#print(seq)

homologs = a2r(acc,70)

print(homologs)
