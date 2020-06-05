import pickle

with open("MarR_ACCs_100K.pkl", mode="rb") as f:
    accs = pickle.load(f)

accList = ""

for i in accs:
    i = i.split(',')
    for j in i:
        if j != "":
            accList += j+'\n'

with open("MarR_ACCs_List.txt", mode="w+") as outfile:
    outfile.write(accList)
