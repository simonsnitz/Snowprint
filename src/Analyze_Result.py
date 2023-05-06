import json
from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker
import os
from pprint import pprint

from Bio.pairwise2 import align
from src.definitions.define_operators import findBestPalindrome, complement

def pull_operator(acc: str):

        # Point to the sqlite database
    engine = create_engine('sqlite:///cache/Snowprint.db')

        # Connect to db, create session, and access the Alignment table
    conn = engine.connect()
    Session = sessionmaker(bind=engine)
    s = Session()
    meta_data = MetaData(bind=engine)
    MetaData.reflect(meta_data)
    Alignment = meta_data.tables["alignment"]
    Regulator = meta_data.tables["regulator"]
    Association = meta_data.tables["association"]
    Operon = meta_data.tables["operon"]
    Operator = meta_data.tables["operator"]

    '''
        Pull relevant data on regulator
    '''

            # Extract Alignment record associated with accession ID argument
    record = s.query(Alignment).filter_by(query_id=acc).first()


        # TODO: have do deal with this edge case.
    if record == None:    
        print('WARNING: No alignment found for '+str(acc))
        return None

        # Extract accession IDs for homologs within the alignment
    else:
        homologs = json.loads(record.homologs)
        accessions = [i['accession'] for i in homologs]


        # Pull out regulators associated with the alignment
    regulators = [s.query(Regulator).filter_by(prot_id=acc).first() for acc in accessions]
        # Filter out empty Regulator records
    regulators = [reg for reg in regulators if reg != None]


        # Pull out association objects linked with each regulator
    assoc_list = [s.query(Association).filter_by(regulator_id=reg.id).first() for reg in regulators]
        # Filter out empty association records
    assoc_list = [assoc for assoc in assoc_list if assoc != None]


        # Pull the operator object
    operator_id = regulators[0].operator_id
    operator = s.query(Operator).filter_by(id=operator_id).first()


        # Display the regulated seq for the target regulator
    #print("Regulated sequence: \n"+str(assoc_list[0].regulated_seq))


        # Display the operator motif
    #print("Operator motif: \n"+str(operator.motif))
    #print("Consensus score: \n"+str(operator.consensus_score))
    #print("number of aligned sequences: \n"+str(operator.number_seqs))

    consensus = "".join(i["base"] for i in json.loads(operator.motif))
    # print("Consensus motif: "+str(consensus))

    #WT_operator = str(json.loads(operator.aligned_seqs)[0])
    #print("WT operator: "+WT_operator)

    try:
        all_seqs = json.loads(operator.aligned_seqs)
        
        entry = {
                "accession": str(acc), 
                "score": operator.consensus_score,
                "sequencesAligned": operator.number_seqs,
                "organism": all_seqs[0]["organism"],
                "intergenic": operator.intergenic,
                "data": all_seqs,
                "consensus": consensus
                }

        return entry
    except:
        return None





    # function to add to JSON
def write_frontend_json(new_data, filename="./frontend/public/data.json"):
    # If data.json doesn't already exist, create it.
    if not os.path.exists(filename):
        print('no existing file')
        with open(filename, 'w+') as file:
            file.write("[]")
    with open(filename,'r+') as file:
        # First we load existing data into a dict.
        file_data = json.load(file)
        # See if entry already exists
        acc = new_data["accession"]
        all_accs = [i["accession"] for i in file_data]
        if acc in all_accs:
            print("Entry already exists in data.json")
            return
        else:
            # Join new_data with file_data inside emp_details
            file_data.append(new_data)
            # Sets file's current position at offset.
            file.seek(0)
            # convert back to json.
            json.dump(file_data, file, indent = 4)
            print("Appended new entry to data.json")



#############################################################################################



    # function to collect perfomance metrics on a model
def assess_model(entry, known_operator: str):


    known_operator = known_operator.upper()



        # What is the "inverted repeat character" of the known operator?
    IR = findBestPalindrome(
        intergenic=known_operator, 
        shortest=4, 
        longest=15, 
        winScore=1, 
        lossScore=-1
        )[0]


    IR_score = IR["score"]/2





        # Is the known operator within the intergenic sequence?
    intergenic = entry["intergenic"]
    comp_known_operator = complement(known_operator)

        # catch an edge case
    if intergenic == None:
        region_align_score = 0
    else:
            # alignment with known operator
        op_align, int_align, score, startPos, endPos = \
            align.localms(known_operator, intergenic, 1, 0, -100, 0)[0]
            # alignment with complement of known operator
        Cop_align, Cint_align, Cscore, CstartPos, CendPos = \
            align.localms(comp_known_operator, intergenic, 1, 0, -100, 0)[0]


        region_align_score = (max(score,Cscore)/len(known_operator))*100


    

        # How well does the predicted operator match the known operator?
            #updated this.
    predicted_operator = (entry["consensus"]).upper()
    # predicted_operator = (entry["data"][0]["predicted_operator"]).upper()

    upstr_align, op_align, score, startPos, endPos = \
        align.localms(known_operator, predicted_operator, 1, 0, -100, 0)[0]

    min_len = min(len(predicted_operator), len(known_operator))
    score = (score/min_len)*100

    operator_align_score = round(score,2)



    performance_metrics = {
        "IR score": IR_score,
        "Region align score": region_align_score,
        "Operator align score": operator_align_score
    }

    return performance_metrics








if __name__ == "__main__":
    lmrr = "WP_011834386.1"
    ebrr = "WP_003976902"
    mexr = "WP_003114897.1"
    ladr = "WP_003721913.1"
    vcer = "WP_001264144.1"
    mtrr = "WP_003693763.1"
    acrr = "WP_000101737"
    mepr = "WP_000397416.1"
    mdtr = "WP_003242592.1"
    mexz = "WP_003088626.1"

    operator_analysis(lmrr)