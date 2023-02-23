import requests
import xml.etree.ElementTree as ET
import json
from pathlib import Path

from sqlalchemy import create_engine, MetaData, insert
from sqlalchemy.orm import sessionmaker




headers = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.163 Safari/537.36'} 

    # file locations:
tmp = Path("./cache/tmp")
metadata_tmp = tmp / "regulator_metadata.xml"


    # efetch with a list of IDs
def batch_acc2MetaData(prot_acc_list: list):
    
    PROTacc = "".join(i+"," for i in prot_acc_list)[:-1]

    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id='+PROTacc+'&rettype=ipg&retmode=xml')
    if response.ok:
        data = response.content
        
        with open(metadata_tmp, mode='wb') as f:
            f.write(data)
            print('NOTE: Metadata cached')

        
        tree = ET.parse(metadata_tmp)
        root = tree.getroot()

        metadata = []

        for i in root:
            try:
                prot = i[0].attrib['accver']
                attribs = i[1][0][0][0].attrib
                accver = attribs['accver']
                start = attribs['start']
                stop = attribs['stop']
                strand = attribs['strand']
                organism = attribs['org']
                org_id = attribs['taxid']
                data = {'protein_acc':prot,'genome_acc':accver, 'start':start, 'stop':stop, \
                    'strand':strand, 'organism': organism, 'org_id': org_id}
                metadata.append(data)
            except:
                pos = prot_acc_list[len(metadata)]
                print("WARNING: No data for "+str(pos))

        return metadata

    else:
            # I'll likely encounter a 'ProteinList KeyError' at some point and will need to deal with it.
        print('WARNING: eFetch API request failed')




def create_regulators(acc: str):

        # Point to sqlite database
    engine = create_engine('sqlite:///cache/Snowprint.db')

        # Connect to db, create session, and access the Alignment table
    conn = engine.connect()
    Session = sessionmaker(bind=engine)
    session = Session()
    meta_data = MetaData(bind=engine)
    MetaData.reflect(meta_data)
    Alignment = meta_data.tables["alignment"]

        # Extract Alignment record associated with accession ID argument
    record = session.query(Alignment).filter_by(query_id=acc).first()

    if record == None:    
        print('NOTE: No alignment found for '+str(acc))

        # Extract accession IDs for homologs within the alignment
    else:
        homologs = json.loads(record.homologs)
        accessions = [i['accession'] for i in homologs]


            # Fetch data for all regulators associated with the alignment
        reg_metadata = batch_acc2MetaData(accessions)


            # Access the Regulator table
        Regulator = meta_data.tables["regulator"]


            # For each regulator, check if record already exists in the database
                # If not, create a new record and fill in fetched data
        for reg in reg_metadata:
            
            prot_acc = reg["protein_acc"]
            genome_acc = reg["genome_acc"]
            organism = reg["organism"]
            organism_id = reg["org_id"]
            start = reg["start"]
            stop = reg["stop"]
            strand = reg["strand"]

            regulator = session.query(Regulator).filter_by(prot_id= prot_acc).first()

            if regulator == None:    

                    # Create a new record with fetched metadata
                new_row = (
                    insert(Regulator).values(
                        prot_id= prot_acc,
                        genome_id= genome_acc,
                        organism= organism,
                        organism_id = organism_id,
                        start_pos= start,
                        stop_pos= stop,
                        strand= strand
                        )
                )
                    # Add the new record and commit it to the DB
                conn.execute(new_row)
                print('UPDATE: Added a regulator entry for '+str(prot_acc))

            else:
                continue

    conn.close()
