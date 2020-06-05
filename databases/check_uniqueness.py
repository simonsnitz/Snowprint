import sys

filename = sys.argv[1]

keys = set()

with open(filename, mode="r") as f:
    no_duplicates = True
    for line in f:

        if line[0] == ">":
            accession = line[1:15]

            if accession in keys:
                print('got a duplicate record')
                no_duplicates = False
            else:
                keys.add(accession)

    if no_duplicates == True:
        print("No duplicates found")
