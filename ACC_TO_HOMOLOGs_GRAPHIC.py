import sys
import acc2homologs as a2h





acc = str(sys.argv[1])

homologDict = a2h.acc2homolog_list(acc, 50)

print(homologDict)
