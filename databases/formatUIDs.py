import requests
import json
import time
import pickle

numberOfLines = 0
for line in open("tetrUIDs.txt").readlines( ): numberOfLines += 1

print(numberOfLines)

with open("tetrUIDs.txt", mode="r") as f:
    accList = []
    counter = 0
    temp = ""
    for i in f.readlines():
        i = i.replace('\n',",")
        temp += i
        counter += 1
        if counter % 300 == 0:
            accList.append(temp)
            temp = ""
            counter = 0


with open("formattedTetRs.pkl", mode="wb") as f:
    pickle.dump(accList, f)
