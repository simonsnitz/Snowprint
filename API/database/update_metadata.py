
def update_metadata(Cluster, data, database, cluster_identity, family):

    db_EMBLs = [i.rep_EMBL for i in Cluster.query.all()]

    unique_EMBLs = [{"rep_EMBL": i["EMBL"], "rep_genome": i["genome"], "reg_index": i["regIndex"]}
                        for i in data if i["EMBL"] not in db_EMBLs]

    print(str(len(unique_EMBLs))+" unique EMBLs found")

    if len(unique_EMBLs) != 0:
        for entry in unique_EMBLs:
            newCluster = Cluster(family = family, cluster_percent_identity = cluster_identity,
                rep_EMBL = entry["rep_EMBL"], rep_genome = entry["rep_genome"], reg_index = entry["reg_index"])
            
            database.session.add(newCluster)
            
        database.session.commit()
        print("added "+str(len(unique_EMBLs))+" new entries")
