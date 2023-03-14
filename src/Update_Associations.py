import requests
import json
from pprint import pprint
from sqlalchemy import create_engine, MetaData, update
from sqlalchemy.orm import sessionmaker


headers = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.163 Safari/537.36'} 



def operon2Intergenic(operon, regIndex, genome_id):


    if operon[regIndex]["direction"] == "+":
        queryGenes = list(reversed(operon[0:regIndex]))
        index = regIndex
        if len(queryGenes) == 0:
            # print("WARNING: Tiny operon with too few genes. This entry will be omitted.")
            return
        for i in queryGenes:
            if i["direction"] == "-":
                startPos = i["stop"]
                stopPos = operon[index]["start"]
                regType = 1
                break
            else:
                start = operon[regIndex-1]["stop"]
                stop = operon[regIndex]["start"]
                testLength = int(stop) - int(start)
                    # Setting this to 100bp is somewhat arbitrary. 
                    # Most intergenic regions >= 100bp. May need to tweak.
                if testLength > 100:
                    startPos = start
                    stopPos = stop
                    regType = 2
                    break
                else:
                    if index == 1:
                        # print('WARNING: Reached end of operon. This entry will be omitted')
                        return None
                    index -= 1

    elif operon[regIndex]["direction"] == "-":
        queryGenes = operon[regIndex+1:]
        index = regIndex
        if len(queryGenes) == 0:
            # print("WARNING: Tiny operon with too few genes. This entry will be omitted.")
            return
        for i in queryGenes:
            if i["direction"] == "+":
                stopPos = i["start"]
                startPos = operon[index]["stop"]
                regType = 1
                break
            else:
                    # Counterintunitive use of "stop"/"start" ...
                    # Start < stop always true, regardless of direction
                start = operon[regIndex]["stop"]
                stop = operon[regIndex+1]["start"]
                testLength = int(stop) - int(start)
                if testLength > 100:
                    startPos = start
                    stopPos = stop
                    regType = 2
                    break
                else:
                    if index == len(operon)-2:
                        # print('WARNING: Reached end of operon. This entry will be omitted')
                        return None
                    else:
                        index += 1
  
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="+str(genome_id)+"&seq_start="+str(startPos)+"&seq_stop="+str(stopPos)+"&strand=1&rettype=fasta"
    response = requests.get(URL)

    if response.ok:
        intergenic = response.text
    else:
        print('FATAL: Bad eFetch request')

         # 800bp cutoff for an inter-operon region. 
         # A region too long makes analysis fuzzy and less accurate.
    output  = ""
    for i in intergenic.split('\n')[1:]:
        output += i
    if len(output) <= 800:
        return {"regulated_seq": output, "reg_type": regType}
    else:
        # print('WARNING: Intergenic region is over 800bp')
        return None




def update_associations(acc: str):

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


            # Extract Alignment record associated with accession ID argument
    record = s.query(Alignment).filter_by(query_id=acc).first()

    if record == None:    
        # print('WARNING: No alignment found for '+str(acc))
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


        # Add regulated seq data for association records
    for assoc in assoc_list:

            # Check if extra data already present in association record
        if assoc.reg_type == None and assoc.regulated_seq == None:

            operon_rec = s.query(Operon).filter_by(id=assoc.operon_id).first()

                # Fetch regulated seq and reg_type data for the associated operon
            operon = json.loads(operon_rec.operon)
            data = operon2Intergenic(operon, assoc.reg_index, operon_rec.genome_id)  
            
            if data != None:
                    # Add data to the association record
                add_regulated_seq = (
                    update(Association).
                    where(Association.c.id == assoc.id).
                    values(reg_type = data["reg_type"], \
                        regulated_seq = data["regulated_seq"])
                )
                conn.execute(add_regulated_seq)
                # print('UPDATE: updated association entry for '+str(assoc.regulator_id))

        else:
            continue

    conn.close()

