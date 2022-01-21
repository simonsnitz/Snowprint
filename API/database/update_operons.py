

def update_operons(Cluster, data, database):

    """ clusters_without_operons = Cluster.query.filter_by(rep_operon = None)

    print("found "+str(clusters_without_operons.count())+" clusters without a representative operon")

    c = 0
    for cluster in clusters_without_operons:
        EMBL = cluster.rep_EMBL
        entry = [i for i in data if i["EMBL"] == EMBL]
        if len(entry) == 1:
            cluster.rep_operon = entry[0]["operon"]
            c += 1

    if c > 0:
        database.session.commit()
    print("operons for "+str(c)+" entries added") """


    entries_with_operons = [i for i in data if "operon" in i]

    c = 0
    for entry in entries_with_operons:
        EMBL = entry["EMBL"]
        db_cluster = Cluster.query.filter_by(rep_EMBL = EMBL).first()
        if not db_cluster.rep_operon:
            db_cluster.rep_operon = entry["operon"]
            c += 1

    
    database.session.commit()
    print("operons for "+str(c)+" entries added")