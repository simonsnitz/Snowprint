import requests
import os
from pathlib import Path
from os.path import exists
import re
import time
import pickle
import xml.etree.ElementTree as ET
import math

   

headers = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.163 Safari/537.36'} 


    # file locations:
p = Path("./cache/tmp")
metadata_tmp = p / "metadata.xml"
genomes_tmp = p / "genomes.txt"




    # just efetch with a list of IDs
def batch_acc2MetaData(prot_acc_list: list):
    
    PROTacc = "".join(i+"," for i in prot_acc_list)[:-1]
    #PROTacc = "WP_000113609,WP_014859138,NP_266817.1"

    startTime = time.time()

    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id='+PROTacc+'&rettype=ipg&retmode=xml')
    if response.ok:
        data = response.content

        #with open('../cache/tmp/metadata.xml', mode='wb') as f:
        with open(metadata_tmp, mode='wb') as f:
            f.write(data)
            print('metadata cached')
        
        
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
                data = {'protein_acc':prot,'genome_acc':accver, 'start':start, 'stop':stop, 'strand':strand}
                metadata.append(data)
            except:
                pos = prot_acc_list[len(metadata)]
                print("no data for "+str(pos))
                    # 'U00096.3' is used as a placeholder to avoid a downstream error
                    # this is jank, but shouldn't bias the result (interoperon won't align)
                # data = {'protein_acc':'None','genome_acc':'U00096.3', 'start':'None', 'stop':'None', 'strand':'None'}

        

        endTime = time.time()
        print("efetching operon batch took "+str(endTime-startTime)+" seconds")
        return metadata

    else:
        print('efetch API request failed')
        # I'll likely encounter a 'ProteinList KeyError' at some point and will need to deal with it.



    # Returns a list of full genome sequences. BIG file.
def batch_NC2genome(acc, genome_acc_list: list):
    # this has failed in the past due to error: 'requests.exceptions.NetworkError'. Failed to establist new connection. Find fix for this?.
    startTime = time.time()
    # added sleep for 1/4 sec because I was getting 429 HTML errors that returned "efetch query unsuccessful..."
    genome_accs = "".join(i+"," for i in genome_acc_list)[:-1]
    
    # time.sleep(0.25) 
    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='+genome_accs+'&rettype=fasta_cds_aa')

    if response.ok:
        data = response.text
        with open(genomes_tmp, mode='w+') as f:
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
                with open(genomes_tmp, mode='w+') as f:
                    f.write(data)
                success = True
                print('genomes cached')
            else:
                print('attempt number '+str(counter+1)+' failed')
                counter += 1
        
        # cache genomes text file
    with open(genomes_tmp, mode='r') as f:
        genomes = f.read().split("\n\n")[:-1]
        genomes = [genome.split("\n") for genome in genomes]

        # cache genomes pickle object
    genomes_obj = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', f"cache/tmp/genomes/{acc}.pkl"))
    with open(genomes_obj, mode="wb") as out:
        pickle.dump(genomes, out)

    endTime = time.time()
    print("efetching genome batch took "+str(endTime-startTime)+" seconds")

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
        return "None", "None"




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




    # This is the main function imported in getIntergenic_batch.

def batch_acc2operon(accessions: list):
    
    acc = accessions[0]

    startTime = time.time()

    metaData = batch_acc2MetaData(accessions)

    genomes_obj = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', f"cache/tmp/genomes/{acc}.pkl"))
    if exists(genomes_obj):
        with open(genomes_obj, mode="rb") as f:
            genomes = pickle.load(f)
            print('cached genomes found')
    
    else:
        print("no cache found. fetching genomes")
        genome_accs = [i["genome_acc"] for i in metaData]
        genomes = batch_NC2genome(acc, genome_accs)


    assert len(metaData) == len(genomes), "different number of metadata/genome entries"

    for i in range(0,len(metaData)):
            # Return metadata on all genes. Find where the target regulator is.
        allGenes, index = parseGenome(genomes[i], metaData[i]["start"], metaData[i]["stop"])
        if allGenes != "None":
                # TODO:
                # Necessary? How is reg['start'] and reg['direction'] different from metadata[i]['start']
                    # and metadata[i]['strand']?
            reg = fasta2MetaData(allGenes[index])
            operon, regIndex = getOperon(allGenes, index, reg['start'], reg['direction'])
            data = {"operon": operon, "regIndex": regIndex, "genome": metaData[i]["genome_acc"] }
                # data = {"operon": operon, "regIndex": regIndex, "genome": metaData[i]["genome_acc"], "organism":metaData["org"] }
            metaData[i]["data"] = data
        else:
            metaData[i]["data"] = "None"
            # how will I process this downstream?


        # will this filtering step be compatible with downstream code? TBD
    metaData = [i for i in metaData if i["data"] != "None"]
    endTime = time.time()
    print("fetching operon batch took "+str(endTime-startTime)+" seconds")
            
    
    # TODO: (get organism ID, from genome_id?)
    # TODO: (get full operon sequence)
    # TODO: (renumber the start/stop positions relative to the operon start? \
    #   to be able to use with the full operon sequence.)


    return metaData



if __name__ == "__main__":
    
    accessions =  ["WP_000113609","WP_014859138","NP_266817.1"]


    