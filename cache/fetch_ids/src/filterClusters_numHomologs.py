

    # for testing purposes
def create_truncated_file(num):
    with open("50perc_clusters_trimmed.clstr", "r") as f:

        with open("truncated.clstr", "w+") as out:

            c = 0
            output = ""
            for line in f.readlines():
                if c < num:
                    output += line
                c += 1

            out.write(output)



    # Filter clusters by the number of homolog sequences it has

def create_filtered_accIDs(min_num_seqs=20):
    #with open("truncated.clstr", "r") as f:
    with open("cache/filtered_clusters.clstr", "r") as f:

        with open("cache/final_accID_list.txt","w+") as out:

            acc = ""
            c = 0
            clusters = []

            for line in f.readlines():
            
                if line[0] == ">":
                    clusters.append([acc,c])
                    acc = ""
                    c = 0
                else:
                    c += 1
                    if line[-2] == "*":
                        acc = line.split(">")[1].split("...")[0]

            for i in clusters:
                if i[1] > min_num_seqs:
                    out.write(i[0]+"\n")


#create_truncated_file(100000)
#create_filtered_accIDs()
