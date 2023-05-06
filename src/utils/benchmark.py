import pandas as pd
import re
from src.Analyze_Result import pull_operator, assess_model
from snowprint import predict_operator
from Bio.Emboss.Applications import WaterCommandline

from io import StringIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


benchmark_file = "../../cache/tmp/Snowprint_benchmarking_04202023.xlsx"



    ###                     IMPORTANT!                          ###

    ###     YOU MUST HAVE BLAST INSTALLED TO RUN THIS SCRIPT    ###

    ### REFER TO THE README FILE FOR INSTALLATION INSTRUCTIONS  ###



def needle_align_score(query_seq, target_seq):
    needle_cline = WaterCommandline(asequence="asis:" + query_seq,
                                     bsequence="asis:" + target_seq,
                                     aformat="simple",
                                     gapopen=10,
                                     gapextend=0.5,
                                     outfile='stdout'
                                     )
    out_data, err = needle_cline()
    out_split = out_data.split("\n")
    p = re.compile("\((.*)\)")
    score =  p.search(out_split[25]).group(1).replace("%", "")
    smallest_seq = min(len(query_seq), len(target_seq))
    bases_compared = int(out_split[23].split(" ")[2])
    coverage = round((( bases_compared/ smallest_seq)*100),2)
    
    out = {"score": score, "coverage": coverage, "bases_compared": bases_compared}
    return out





def blast_align(sequence1, sequence2):

    # Create two sequence files
    seq1 = SeqRecord(Seq(sequence1),
                    id="seq1")
    seq2 = SeqRecord(Seq(sequence2),
                    id="seq2")
    SeqIO.write(seq1, "seq1.fasta", "fasta")
    SeqIO.write(seq2, "seq2.fasta", "fasta")

    # Run BLAST and parse the output as XML
    output = NcbiblastnCommandline(query="seq1.fasta", subject="seq2.fasta", outfmt=5, word_size=4)()[0]
    blast_result_record = NCBIXML.read(StringIO(output))

    # Print some information on the result
    data = []
    for alignment in blast_result_record.alignments:
        for hsp in alignment.hsps:

            identity = round((hsp.identities/hsp.align_length)*100, 2)
            coverage = round((hsp.align_length/blast_result_record.query_length)*100, 2)

            entry = {"length": alignment.length,
                    "e value": hsp.expect,
                    "identity": identity,
                    "coverage": coverage}
            data.append(entry)

    try:
        best_e = min([i["e value"] for i in data])
        best_align = [i for i in data if i["e value"] == best_e][0]
        return best_align

    except:
        entry = {"length": None,
        "e value": None,
        "identity": None,
        "coverage": None}
        return entry



def benchmark():
    with open(benchmark_file, "rb+") as f:

        df = pd.read_excel(f)

        IDs = df.loc[:,"Protein ID"].values
        Operators = df.loc[:,"Known operator"].values

        for i in range(0, len(IDs)):

            if pd.isna(df.loc[i,"Predicted operator"]):

                print("starting: "+str(IDs[i]))
                
                predict_operator(IDs[i])
                operator_data = pull_operator(IDs[i])
                if operator_data != None:

                    metrics = assess_model(operator_data, Operators[i])


                    df.loc[i,"Predicted operator"] = operator_data["consensus"]                
                    df.loc[i,"Inverted repeat score"] = metrics["IR score"]
                    df.loc[i,"Region align score"] = metrics["Region align score"]
                    df.loc[i,"Operator align score"] = metrics["Operator align score"]
                    df.loc[i,"Number aligned seqs"] = operator_data["sequencesAligned"]
                    df.loc[i,"Consensus Score"] = operator_data["score"]

                    op_align = blast_align(Operators[i], operator_data["consensus"])

                    df.loc[i,"Op-align Evalue"] = op_align["e value"]
                    df.loc[i,"Op-align identity"] = op_align["identity"]
                    df.loc[i,"Op-align coverage"] = op_align["coverage"]

                    df.to_excel(benchmark_file)
                    print("Updated Snowprint_benchmarking.xlsx")
                    
                else:
                    print("Snowprint failed for "+str(IDs[i]))




if __name__ == "__main__":

    benchmark()