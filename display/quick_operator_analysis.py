import json
from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker


def operator_analysis(acc: str):

        # Point to the sqlite database
    engine = create_engine('sqlite:///API/GroovIO.db')

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

    if record == None:    
        print('WARNING: No alignment found for '+str(acc))
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


        # Pull the operator object
    operator_id = regulators[0].operator_id
    operator = s.query(Operator).filter_by(id=operator_id).first()


        # Display the regulated seq for the target regulator
    #print("Regulated sequence: \n"+str(assoc_list[0].regulated_seq))


        # Display the operator motif
    #print("Operator motif: \n"+str(operator.motif))
    print("Consensus score: \n"+str(operator.consensus_score))
    print("number of aligned sequences: \n"+str(operator.number_seqs))

    consensus = "".join(i["base"] for i in json.loads(operator.motif))
    print("Consensus motif: "+str(consensus))

    WT_operator = str(json.loads(operator.aligned_seqs)[0])
    print("WT operator: "+WT_operator)


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