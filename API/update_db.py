from api import db, app
from api.models import Cluster

from database.update_metadata import update_metadata
from database.update_homologs import update_homologs
from database.update_operons import update_operons
#from database.update_dna_motifs import update_dna_motifs


import pickle


db_file = "database/data/filtered_TetRs.pkl"
cluster_identity = 50
family = "TetR"

with open(db_file, mode="rb") as f:
    data = pickle.load(f)


#update_metadata(Cluster, data, db, cluster_identity, family)
#update_operons(Cluster, data, db)
#update_homologs(Cluster, data, db)