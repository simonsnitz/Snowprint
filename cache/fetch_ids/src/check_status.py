with open("finished_clusters.txt", "r") as f:
    c = 0
    for line in f.readlines():
        c += 1

    print(str(c)+" out of 9577 regulators finished")
    percent = round(100*(c/9577),3)
    print(str(percent)+"% complete")
