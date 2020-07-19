import sys
import acc2homologs as a2h
import acc2operon as a2o
from display import test_operon2graphic as o2g



acc = str(sys.argv[1])
max_regulators = int(sys.argv[2])

homologDict = a2h.acc2homolog_list(acc, 50)

accessions = [ i["accession"] for i in homologDict if i["identity"] <= 59]
identities = [ i["identity"] for i in homologDict ]

print(len(accessions))
operons_with_regulators = a2o.enzyme_acc2regulator(accessions, max_regulators)

o2g.create_html(operons_with_regulators, "display/html_pages/"+acc)

