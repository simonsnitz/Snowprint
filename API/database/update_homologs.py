

def update_homologs(Cluster, data, database):

    """ clusters_without_homologs = Cluster.query.filter_by(homologs = None)

    print("found "+str(clusters_without_homologs.count())+" clusters without homolog data")

    c = 0
    for cluster in clusters_without_homologs:
        EMBL = cluster.rep_EMBL
        entry = [i for i in data if i["EMBL"] == EMBL]
        if len(entry) == 1 and "homologs" in entry[0]:
            cluster.homologs = entry[0]["homologs"]
            c += 1

    #if c > 0:
        #database.session.commit()
    print("homologs for "+str(c)+" entries added") """

    entries_with_homologs = [i for i in data if "homologs" in i]

    c = 0
    for entry in entries_with_homologs:
        EMBL = entry["EMBL"]
        db_cluster = Cluster.query.filter_by(rep_EMBL = EMBL).first()
        if not db_cluster.homologs:
            db_cluster.homologs = entry["homologs"]
            c += 1

    
    database.session.commit()
    print("homologs for "+str(c)+" entries added")