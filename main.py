from src.Create_Alignment import create_alignment
from src.Create_Regulators import create_regulators
from src.Create_Operons import create_operons
from src.Update_Associations import update_associations
from src.Create_Operators import create_operators
from src.Analyze_Result import pull_operator, write_frontend_json, assess_model

from alive_progress import alive_bar
import pandas as pd
import time
import sys
from pprint import pprint




    # Perform the Snowprint workflow
def predict_operator(acc):

    with alive_bar(5, dual_line=True, title='Snowprint') as bar:
      
        bar.text = f'-> Collecting homologs for {acc}'
        create_alignment(acc)
        bar()

        bar.text = f'-> Collecting metadata for {acc} homologs'
        create_regulators(acc)
        bar()

        bar.text = f'-> Collecting genetic context for {acc} homologs'
        create_operons(acc)
        bar()

        bar.text = f'-> Extracting operators for {acc} homologs'
        update_associations(acc)
        bar()

        bar.text = f'-> Aligning operators and generating metrics for {acc}'
        create_operators(acc)
        bar()





    # Update data.json in the react app's public folder
# def update_frontend_data(acc):
#     operator = pull_operator(acc)
#     write_frontend_json(operator)






def benchmark():

    with open("Snowprint_benchmarking.xlsx", "rb+") as f:

        xl = pd.ExcelFile(f)
        families = xl.sheet_names

        for family in families:

            df = pd.read_excel(f, sheet_name=family)

            IDs = df.loc[:,"Protein ID"].values
            Operators = df.loc[:,"Known operator"].values

            for i in range(0, len(IDs)):

                start = time.time()
                
                predict_operator(IDs[i])
                operator_data = pull_operator(IDs[i])
                metrics = assess_model(operator_data, Operators[i])

                end = time.time()
                elapsed_time = end-start

                df.loc[i,"Predicted operator"] = operator_data["data"][0]["predicted_operator"]             

                df.loc[i,"Time to complete"] = elapsed_time    

                df.loc[i,"Inverted repeat score"] = metrics["IR score"]
                
                df.loc[i,"Region align score"] = metrics["Region align score"]

                df.loc[i,"Operator align score"] = metrics["Operator align score"]

                df.loc[i,"Number aligned seqs"] = operator_data["sequencesAligned"]

                df.loc[i,"Consensus Score"] = operator_data["score"]

                df.to_excel("Snowprint_benchmarking.xlsx", sheet_name=family)
                print("Updated Snowprint_benchmarking.xlsx")


benchmark()


    # RamR (40%)
# r = "WP_000113609.1"
# known_operator = "TATAATGAGTGAGTAAGCACTCATTATAA"
# predict_operator(r)
# analyze_model(r, known_operator)





#     # Open regulator ID file and predict operators for each
# with open("cache/ids/reg_ids.txt", "r") as f:
#     regs = [r.strip() for r in f.readlines() if r[0] != ">"]

# for r in regs:
#     print(r)
#     predict_operator(r)
#     operator_analysis(r)



    # perform snowprint analysis. update data.json for frontend display
def snowprint(acc):
    predict_operator(acc)
    data = pull_operator(acc)
    write_frontend_json(data)


# acc = sys.argv[1]
# snowprint(acc)