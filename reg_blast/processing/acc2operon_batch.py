from urllib.request import OpenerDirector
import requests
import os
import re
import time
from Bio import Entrez
import pickle
import xml.etree.ElementTree as ET
import math

''' available functions: 

batch_acc2MetaData([protein_acc1, protein_acc2, ...])
batch_NC2genome([genome_acc1, genome_acc2, ...])
parseGenome(genome, seq_start, seq_stop)
fasta2MetaData(allGenes[index])
getOperon(allGenes, index, seq_start, direction)

'''    

Entrez.email = os.getenv("EMAIL", "doelsnitz@utexas.edu")
headers = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.163 Safari/537.36'} 




def batch_acc2MetaData(prot_acc_list: list):
    
    PROTacc = "".join(i+"," for i in prot_acc_list)[:-1]
    #PROTacc = "WP_000113609,WP_014859138,NP_266817.1"

    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id='+PROTacc+'&rettype=ipg&retmode=xml')
    if response.ok:
        data = response.content

        with open('../cache/tmp/metadata.xml', mode='wb') as f:
            f.write(data)
            print('metadata cached')
        
        
        tree = ET.parse('../cache/tmp/metadata.xml')
        root = tree.getroot()

        items = []

        for i in root:
            prot = i[0].attrib['accver']
            metadata = i[1][0][0][0].attrib
            accver = metadata['accver']
            start = metadata['start']
            stop = metadata['stop']
            strand = metadata['strand']
            data = {'protein_acc':prot,'genome_acc':accver, 'start':start, 'stop':stop, 'strand':strand}
            items.append(data)
        
        return items

    else:
        print('efetch API request failed')
        # I'll likely encounter a 'ProteinList KeyError' at some point and will need to deal with it.




def batch_NC2genome(genome_acc_list: list):
    # this has failed in the past due to error: 'requests.exceptions.NetworkError'. Failed to establist new connection. Find fix for this?.
    
    # added sleep for 1/4 sec because I was getting 429 HTML errors that returned "efetch query unsuccessful..."
    genome_accs = "".join(i+"," for i in genome_acc_list)[:-1]
    
    time.sleep(0.25) 
    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='+genome_accs+'&rettype=fasta_cds_aa')

    if response.ok:
        data = response.text
        with open('../cache/tmp/genomes.txt', mode='w+') as f:
            f.write(data)
            print('genomes cached')
    else:
        print(response.status_code)
        print('efetch query unsuccessful. Genome could not be found. re-trying ...')
        
        success = False
        counter = 1
        while(success == False):
            time.sleep(counter) 
            response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='+genome_accs+'&rettype=fasta_cds_aa')

            if response.ok:
                data = response.text
                with open('../cache/tmp/genomes.txt', mode='w+') as f:
                    f.write(data)
                success = True
                print('genomes cached')
            else:
                print('attempt number '+str(counter+1)+' failed')
                counter += 1
        

    with open('../cache/tmp/genomes.txt', mode='r') as f:
        genomes = f.read().split("\n\n")[:-1]
        genomes = [genome.split("\n") for genome in genomes]

    return genomes




def parseGenome(genome, start, stop):
    re1 = re.compile(start)
    re2 = re.compile(stop)
    geneIndex = 0
    allGenes = []
    for i in genome:
        if i[0] == '>':
            if re1.search(i):
                if re2.search(i):
                    regIndex = geneIndex
            geneIndex += 1
            allGenes.append(i)
    try:
        return allGenes, regIndex
    except:
        print('regulator not found in genome')




def fasta2MetaData(fasta):
        # this is a jank way of extracting information (indexing through a string)
    metaData = {}
    regulator = fasta.split(' [')
    
    for i in regulator:
        if i[:10] == 'locus_tag=':
            metaData['alias'] = i[10:-1]
        elif i[:8] == 'protein=':
            metaData['description'] = i[8:-1].replace("'", "")
        elif i[:11] == 'protein_id=':
            metaData['link'] = i[11:-1]
        elif i[:9] == 'location=':
            if i[9:20] == 'complement(':
                metaData['direction'] = '-'
                location = i[20:-2]
                location = location.split('..')
                metaData['start'] = int(re.sub("\D", "", location[0]))
                metaData['stop'] = int(re.sub("\D", "", location[1]))
            else:
                metaData['direction'] = '+'
                location = i[9:-1]
                location = location.split('..')
                metaData['start'] = int(re.sub("\D", "", location[0]))
                metaData['stop'] = int(re.sub("\D", "", location[1]))
    
    if 'link' not in metaData.keys():
        metaData['link'] = ""   #eventually should be empty. Jinja checks if it should set a link or not
    
    return metaData




def getOperon(allGenes, index, seq_start, strand):
    '''
    Rules for inclusion/exclusion of genes from operon:
        - always take immediately adjacent genes
        - if query gene is in same direction as regulator, include it.
        - if query gene is expressed divergently from regulator, 
                grab all adjacent genes that are expressed divergently (change strand direction for next genes)
        - if query gene is expressed divergently from a co-transcribed neighbor of the regulaor, 
                grab that gene. (it may be another regulator. Important to know).
        - if query gene direction converges with regulator, exclude it.
    '''

    def getGene(geneStrand, direction, nextGene, geneList, index):
        
        while geneStrand == nextGene['direction']:
            if direction == '+':
                nextIndex = index+1
            elif direction == '-':
                nextIndex = index-1
                
            try:
                nextGene = fasta2MetaData(allGenes[nextIndex])

                if abs(seq_start - nextGene['start']) > 8000:       #added this. break if too far away
                    break
                elif geneStrand == '-' and nextGene['direction'] == '+' and direction == '+':
                    geneList.append(nextGene)
                elif geneStrand == '+' and nextGene['direction'] == '-' and direction == '-':
                    geneList.append(nextGene)
                elif geneStrand == nextGene['direction']:
                    geneList.append(nextGene)
                index = nextIndex
            except:
                break

    geneStrand = strand
    
    #attempt to get downstream genes, if there are any genes downstream
    try:
        indexDOWN = index-1
        downGene = fasta2MetaData(allGenes[indexDOWN])
        #if seq_start > downGene['start']:
        if strand == '+' and downGene['direction'] == '-':
            geneStrand = downGene['direction']
    
        downgenes = [downGene]
        getGene(geneStrand,'-',downGene, downgenes, indexDOWN)
    
        geneArray = list(reversed(downgenes))
    except:
        geneArray = []

    geneArray.append(fasta2MetaData(allGenes[index]))
    regulatorIndex = (len(geneArray)-1)

    geneStrand = strand
    
    #attempt to get upstream genes, if there are any genes upstream
    try:
        indexUP = index+1
        upGene = fasta2MetaData(allGenes[indexUP])
        #if seq_start > upGene['start']:
        if strand == '-' and upGene['direction'] == '+':
            geneStrand = upGene['direction']
        
        geneArray.append(upGene)

        getGene(geneStrand, '+', upGene, geneArray, indexUP)
    except:
        return geneArray, regulatorIndex

    return geneArray, regulatorIndex




def accession2UniprotID(alias):

    payload = {'from': 'P_REFSEQ_AC',
                'to': 'ACC',
                'format': 'list',
                'query': alias,
                }

    response = requests.get('https://www.uniprot.org/uploadlists/', params=payload, headers=headers)
    if response.ok:
            #may have multiple UniprotIDs. in that case, take the first one.
        try:
            uniprotID = response.text.split('\n')[0]
        except:
            uniprotID = response.text[:-1]
        return uniprotID
    else:
        response.raise_for_status()




def map_retrieve(ids2map, source_fmt='ACC+ID',
    target_fmt='ACC', output_fmt='list'):

    URL_ENDPOINT = 'https://www.uniprot.org/uploadlists/'

    #if hasattr(ids2map, 'pop'):
    ids2map = ' '.join(ids2map)
    payload = {'from': source_fmt,
    'to': target_fmt,
    'format': output_fmt,
    'query': ids2map,
    }
    response = requests.get(URL_ENDPOINT, params=payload)
    if response.ok:
        return response.text.split('\n')[:-1]
    else:
        response.raise_for_status()




def batch_acc2operon(accessions: list):
    
    startTime = time.time()
    metaData = batch_acc2MetaData(accessions) # modify this to input a list
    if (i != "EMPTY" for i in metaData):
        time.sleep(0.25)
        genome_accs = [i["genome_acc"] for i in metaData]
        genomes = batch_NC2genome(genome_accs)
        for i in range(0,len(metaData)):
            allGenes, index = parseGenome(genomes[i], metaData[i]["start"], metaData[i]["stop"])
            reg = fasta2MetaData(allGenes[index])
            operon, regIndex = getOperon(allGenes, index, reg['start'], reg['direction'])
            data = {"operon": operon, "regIndex": regIndex, "genome": metaData[i]["genome_acc"] }
                # data = {"operon": operon, "regIndex": regIndex, "genome": metaData[i]["genome_acc"], "organism":metaData["org"] }
            metaData[i]["data"] = data

    endTime = time.time()
    print("fetching operon batch took "+str(endTime-startTime)+" seconds")
            
    return metaData




def append_operons(input_file, batch_size):
    with open(input_file, mode="rb") as f:
        db = pickle.load(f)

        num_batches = int(math.ceil(len(db)/batch_size))

        #print(db[0])

        accessions = []
        for i in range(0,num_batches):
            accessions.append( [entry["EMBL"] for entry in db[i*batch_size : i*batch_size+batch_size]] )

        # print(accessions[0])
        # print(num_batches)

        print("batch size: "+str(batch_size))
        metaData = batch_acc2operon(accessions[0])




if __name__ == "__main__":
    
    accessions =  ["WP_000113609","WP_014859138","NP_266817.1"]
    
    genome_accs = ['NC_003197.2','NC_018220.1','NC_002662.1']

    #batch_acc2operon(accessions)

    append_operons('../cache/all_the_regulators/metadata/filtered_TetRs.pkl', 100)



    