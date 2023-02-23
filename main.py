from src.Create_Alignment import create_alignment
from src.Create_Regulators import create_regulators
from src.Create_Operons import create_operons
from src.Update_Associations import update_associations
from src.Create_Operators import create_operators
from src.Analyze_Result import pull_operator, write_frontend_json, assess_model

import time
import sys


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


    # for benchmarking purposes only
def analyze_model(r, known_operator):
    operator = pull_operator(r)
    print("accession: "+str(r))
    print("predicted operator: "+str(operator["data"][0]["predicted_operator"]))
    print("known operator: "+str(known_operator))
    assess_model(operator, known_operator)







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




def get_snowprint(acc):
    predict_operator(acc)
    data = pull_operator(acc)
    write_frontend_json(data)


acc = sys.argv[1]
get_snowprint(acc)