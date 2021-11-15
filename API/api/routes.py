from flask import render_template, request, flash, redirect, url_for, jsonify
import time
import json
from api import app
#from api import db, app
#from api.models import Sensor, Ligand, Reference
import requests


@app.route('/', methods = ['GET'])
def index():
    #sensors = Sensor.query.all()
    #sensorData = {s.id: {"alias":s.alias, "uniprotID":s.uniprotID, "family":s.family, "accession":s.accession} for s in sensors}

    return render_template('index.html', title="home page")