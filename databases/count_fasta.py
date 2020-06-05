import sys

inputFile = sys.argv[1]

with open(inputFile, mode="r") as f:
    counter = 0
    for i in f:
        if i[0] == ">":
            counter += 1

    print(counter)
