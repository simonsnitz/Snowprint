from flask import render_template, request, flash, redirect, url_for, jsonify
import json
from api import db, app
from api.models import Alignment, Regulator, Association, Operator



@app.route('/', methods = ['GET'])
def index():
    return render_template('index.html', title="home page")


@app.route('/GroovIO/<accession_id>', methods = ['GET'])
def GroovIO(accession_id):

        # Extract Alignment record associated with accession ID argument
    record = Alignment.query.filter_by(query_id=accession_id).first()

    if record == None:    
        output = 'WARNING: No alignment found for '+str(accession_id)
        return render_template('prediction.html', title="GroovIO", output=output)

        # Extract accession IDs for homologs within the alignment
    else:
        homologs = json.loads(record.homologs)
        accessions = [i['accession'] for i in homologs]

        # Pull out regulators associated with the alignment
    regulators = [Regulator.query.filter_by(prot_id=acc).first() for acc in accessions]
        # Filter out empty Regulator records
    regulators = [reg for reg in regulators if reg != None]
        # Pull out association objects linked with each regulator
    assoc_list = [Association.query.filter_by(regulator_id=reg.id).first() for reg in regulators]
        # Filter out empty association records
    assoc_list = [assoc for assoc in assoc_list if assoc != None]
        # Pull out the operator object
    operator_id = regulators[0].operator_id
    operator = Operator.query.filter_by(id=operator_id).first()


    # Compile data for the output

        # Operator motif
    motif = operator.motif
        # Consensus score
    score = operator.consensus_score
        # Number of aligned sequences
    num_seqs = operator.number_seqs

        # All relevant data on each homolog used to make the consensus motif
    aligned_homologs = json.loads(operator.aligned_seqs)

    for h in aligned_homologs:
        acc = h["accession"]
        association = Regulator.query.filter_by(prot_id=acc).first().operons[0]

        h["regulated_seq"] = association.regulated_seq
        h["operon"] = json.loads(association.operon.operon)
    

        # Final output
    output = {"accession":accession_id,"motif":motif, "score":score, \
        "num_seqs":num_seqs, "aligned_homologs": aligned_homologs}

    output = json.dumps(output)


    return render_template('prediction.html', title="GroovIO", output=output)