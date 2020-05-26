import accID2operon as a2o
from getIntergenic import operon2Intergenic

import pickle



with open('ramr_homolog_list.pkl', mode="rb") as f:
    homologList = pickle.load(f)

print(len(homologList))

for i in range(0,len(homologList)):
    protein = homologList[i]
    data = a2o.acc2operon(protein["accession"])
    try:
        intergenic = operon2Intergenic(data["operon"], data["regIndex"], data["genome"])
        homologList[i]["intergenic"] = intergenic["intergenicSeq"]
        homologList[i]["regType"] = intergenic["regType"]
    except:
        print("no data for intergenic region")
    print('got intergenic region')
    print(intergenic)
    with open(f'homolog_ramrIntergenic.pkl', mode="wb") as f:
        pickle.dump(homologList,f)



