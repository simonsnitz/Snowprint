import accID2operon as a2o
from getIntergenic import getIntergenicSeq

import pickle



with open('homolog_list.pkl', mode="rb") as f:
    homologList = pickle.load(f)

for i in range(0,len(homologList)):
    protein = homologList[i]
    data = a2o.acc2operon(protein["accession"])
    intergenic = getIntergenicSeq(data["operon"], data["regIndex"], data["genome"])
    homologList[i]["intergenic"] = intergenic
    print('got intergenic region')
    print(intergenic)
    #with open(f'homolog_list2.pkl', mode="wb") as f:
    #    pickle.dump(homologList,f)



