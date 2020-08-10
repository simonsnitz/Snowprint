import sys

sys.path.insert(1, '/home/simon/git/tuSeek')

import acc2homologs as a2h

#ramr
acc = "WP_000113609.1"

#camr
#acc = "BAA03510"

#seq = a2h.accID2sequence(acc)
#print(seq)

homologs = a2h.acc2homolog_list(acc,70)

print(homologs)
