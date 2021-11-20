from flask import render_template, request, flash, redirect, url_for, jsonify
import time
import json
from api import db, app
from api.models import Cluster
from api.render_operon_graphic import createGraphic
import requests


@app.route('/', methods = ['GET'])
def index():
    return render_template('index.html', title="home page")

@app.route('/cluster/<cluster_id>', methods = ['GET'])
def cluster(cluster_id):
    
    cluster = Cluster.query.filter_by(id=str(cluster_id)).first_or_404()

    operon = createGraphic(cluster.rep_operon, cluster.reg_index)

    # clusterData = {"id":cluster.id, "EMBL":cluster.rep_EMBL, "family":cluster.family, 
    #     "genome":cluster.rep_genome, "reg_index":cluster.reg_index}

    return render_template('cluster.html', title="cluster", cluster=cluster, operon=operon)