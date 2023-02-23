from Bio.pairwise2 import align
import json
from pprint import pprint

from src.definitions.define_operators import findBestPalindrome

from sqlalchemy import create_engine, MetaData, insert, update
from sqlalchemy.orm import sessionmaker



def findOperatorInIntergenic(intergenic, operator, ext_length=0):

    operator_length = len(operator)
    
        # This function is needed to extract the ENTIRE aligned region.
        # Otherwise, mismatched ends will not be extracted.
    
    def extractOperator(intergenic, op_align, ext_length):
        begin = 0
        for i in op_align:
            # Find starting position of operator within intergenic region
            if i == '-':
                begin += 1
            else:
                break
        end = (begin + operator_length)
            # Can change indexing of output sequence to include more or less of alignment
        try:
            upstream = intergenic[begin-ext_length:begin].lower()
            mid = intergenic[begin:end]
            downstream = intergenic[end:end+ext_length].lower()
            operator = upstream+mid+downstream
        except:
            operator = intergenic[begin:end]
        return operator


    try:
        upstr_align, op_align, score, startPos, endPos = \
            align.localms(intergenic, operator, 2, -0.5, -100, 0)[0]
    except:
        print("WARNING: Regulated sequence alignment failed")
        return None
    
        # Heavily penalizing gap opening avoids errors with downstream data processing, but may miss out on interesting biological features
        # Returns the aligned operator sequence if a similarity threshold is met. Score threshold (7) should be tuned.
        # Set score cutoff to be 10% of max. Arbitrary, but seems reasonable.
    max_score = 2*operator_length
    score_cutoff = max_score*0.1

    if score > score_cutoff:
        operator = extractOperator(upstr_align, op_align, ext_length)
        return {"operator":operator, "score":score}
    else:
        print('WARNING: Alignment score is not above cutoff')
        return None





def getConsensus(metrics):
    
        # Filter by identity and by alignment score. Some alignments from ~80% homologs have crap scores
    allOperators = [ i["predicted_operator"] for i in metrics 
             if i["align_score"] != 0 ]

    num_seqs = len(allOperators)
	    # Initialize list of dictionaries
    baep = [{base:1}
            for base in allOperators[0] 
        ]

	    # Populate dataset with base representation from all input operators
    for operator in allOperators[1:]:
        if len(operator) == len(allOperators[0]):
            for pos in range(0, len(operator)):
                base = operator[pos]

                try:
                    baep[pos][base] +=1
                except:
                    baep[pos].update({base:1})


    max_values = [max(baep[pos].values()) for pos in range(0,len(baep))]
    max_score = max(max_values)
	    # Convert base conservation scores as a percent of max
    max_values_percent = [ round(i/max_score,2) for i in max_values ] 

    def get_key(my_dict,val):
        for key, value in my_dict.items():
            if val == value:
                return key

        return "key doesn't exist"

	    # Create a list of most conserved bases at each position
    consensusSeq = [ get_key(baep[pos], max_values[pos])
        for pos in range(0, len(baep))
    ]

	    # Dictionary containing the base and it's score at each position
    consensus_data = [{"base":consensusSeq[i] , "score":max_values_percent[i]} 
            for i in range(0,len(max_values))
        ]

    return {"motif_data":consensus_data, "num_seqs":num_seqs}




def get_consensus_score(operator, consensus_data, ext_length):
    
    max_score = 0
    consensus_score = 0

    for i in range(0,len(operator)):
        if operator[i].isupper():
            max_score += 1
            consensus_score += consensus_data["motif_data"][i+ext_length]['score']**2
    
    score = round((consensus_score/max_score)*100, 3)

    return score




def fetch_operator_data(homolog_metadata, acc, ext_length=5):

    regulated_seqs = [h["regulated_seq"] for h in homolog_metadata]

        # Iterates the palindrome locater function through multiple relevant scoring parameters
    operators = []

    test_params = [{"w":2,"l":-2}, {"w":2,"l":-3}, {"w":2,"l":-4}]

    for i in test_params:
        ops = [findBestPalindrome(intergenic=regulated_seqs[0], \
        shortest=5, longest=15, winScore=i["w"], lossScore=i["l"])][0]
        for operator in ops:
            operators.append(operator)



        # Output data to be returned 
    operator_data = { 
        "accession": str(acc),
        "aligned_seq": "None",
        "num_seqs": "None",
        "consensus_score": 0,
        "motif": "None",
        "aligned_seqs": "None",
        "intergenic": regulated_seqs[0],
    }
    

    for operator in operators:

        metrics = []
        for h in homolog_metadata:
            i = h["regulated_seq"]
            homolog = {}
            op = findOperatorInIntergenic(i, operator["seq"], ext_length=ext_length)
            if op != None:
                homolog["predicted_operator"] =  op["operator"]
                homolog["align_score"] =  op["score"]
                homolog["identity"] = h["identity"]
                homolog["coverage"] = h["coverage"]
                homolog["accession"] = h["accession"]
                #homolog["regulated_seq"] = i
                homolog["organism"] = h["organism"]

                metrics.append(homolog)

            # Create the consensus 'motif' and calculate its quality score
        consensus = getConsensus(metrics)
        consensus_score = get_consensus_score(operator["seq"], consensus, ext_length)
            # Pull out the predicted operator from the original query protein
        operator["seq"] = findOperatorInIntergenic(regulated_seqs[0], \
            operator["seq"], ext_length=5)["operator"]

        # Warning: Only the CONSENSUS SCORE is used to identify the 
            # best operator. This should also incorporate the
            # NUMBER OF ALIGNMENTS as a metric to make this decision.

        if consensus_score > operator_data["consensus_score"]:
            operator_data["consensus_score"] = consensus_score
            operator_data["aligned_seq"] = operator["seq"]
            operator_data["num_seqs"] = consensus["num_seqs"]
            operator_data["motif"] = consensus["motif_data"]
            operator_data["aligned_seqs"] = metrics

    return operator_data





def create_operators(acc: str):

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
    Operator = meta_data.tables["operator"]



                # Extract Alignment record associated with accession ID argument
    
    record = s.query(Alignment).filter_by(query_id=acc).first()

    if record == None:    
        print('FATAL: No alignment found for '+str(acc))
        return

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


        # Create a dictionary with relevant metadata on all proteinhomologs
    homolog_metadata = [{
        "regulated_seq": assoc.regulated_seq,
        "accession": s.query(Regulator).filter_by(id=assoc.regulator_id).first().prot_id,
        "organism": s.query(Regulator).filter_by(id=assoc.regulator_id).first().organism
    } for assoc in assoc_list]
    for i in homolog_metadata:
        homolog = [h for h in homologs if h["accession"] == i["accession"]]
        if len(homolog) != 0:
            i["identity"] = homolog[0]["identity"]
            i["coverage"] = homolog[0]["coverage"]
    
        # Remove homologs that don't have a regulated seqeunce associated with them
    homolog_metadata = [h for h in homolog_metadata if h["regulated_seq"] != None]
    
    
    has_operator = regulators[0].operator_id


    if has_operator == None:

        operator_data = fetch_operator_data(homolog_metadata, acc)
        
            # Create a new operon record
        new_operator = (
            insert(Operator).values(
            motif= json.dumps(operator_data["motif"]),
            number_seqs= operator_data["num_seqs"],
            consensus_score= operator_data["consensus_score"],
            aligned_seqs= json.dumps(operator_data["aligned_seqs"]),
            intergenic= operator_data["intergenic"]
                    )
            )

        conn.execute(new_operator)

        operator = s.query(Operator).filter_by(motif=json.dumps(operator_data["motif"])).first()

        for reg in regulators:
                # Add data to the association record
            link_reg_to_operator = (
                update(Regulator).
                where(Regulator.c.id == reg.id).
                values(operator_id = operator.id)
                )
            conn.execute(link_reg_to_operator)

        print('SUCCESS: Added an operator entry for '+str(operator_data["accession"]))
        
        conn.close()

    else:
        print("NOTE: An operator already exists for "+str(regulators[0].prot_id))
        conn.close()
