from src.Create_Alignment import create_alignment
from src.Create_Regulators import create_regulators
from src.Create_Operons import create_operons
from src.Update_Associations import update_associations
from src.Create_Operators import create_operators
from src.Analyze_Result import pull_operator, write_frontend_json, assess_model

import pandas as pd
import time
import sys
from pprint import pprint


    # Perform the Snowprint workflow
def predict_operator(acc):
    startTime = time.time()

    create_alignment(acc)
    create_regulators(acc)
    create_operons(acc)
    update_associations(acc)
    create_operators(acc)

    endTime = time.time()
    print("total processing time: "+str(endTime-startTime)+" seconds")



    # Update data.json in the react app's public folder
# def update_frontend_data(acc):
#     operator = pull_operator(acc)
#     write_frontend_json(operator)






def benchmark():

    with open("Snowprint_metrics.xlsx", "rb+") as f:

        sheet_name = "LacI"

        df = pd.read_excel(f, sheet_name=sheet_name)

        IDs = df.loc[:,"Protein ID"].values
        Operators = df.loc[:,"Known operator"].values

        for i in range(0, len(IDs)):
            
            predict_operator(IDs[i])
            operator_data = pull_operator(IDs[i])
            metrics = assess_model(operator_data, Operators[i])

            df.loc[i,"Inverted repeat score"] = metrics["IR score"]
            
            df.loc[i,"Region align score"] = metrics["Region align score"]

            df.loc[i,"Operator align score"] = metrics["Operator align score"]

            df.loc[i,"Number aligned seqs"] = operator_data["sequencesAligned"]

            df.loc[i,"Consensus Score"] = operator_data["score"]

            df.to_excel("Snowprint_metrics.xlsx", sheet_name=sheet_name)
            print("Updated Snowprint_metrics.xlsx")


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