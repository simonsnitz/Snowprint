from api import db, app
from api.models import Cluster

# import update functions here

import pickle


db_file = "database/data/filtered_TetRs.pkl"

with open(db_file, mode="rb") as f:
    data = pickle.load(f)

print(data[0])
