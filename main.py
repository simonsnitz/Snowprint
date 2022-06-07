from src.Create_Alignment import create_alignment
from src.Create_Regulators import create_regulators
from src.Create_Operons import create_operons
from src.Update_Associations import update_associations
from src.Create_Operators import create_operators
from display.quick_operator_analysis import operator_analysis

import time




def predict_operator(acc):

    startTime = time.time()

    create_alignment(acc)
    create_regulators(acc)
    create_operons(acc)
    update_associations(acc)
    create_operators(acc)

    endTime = time.time()

    print("total processing time: "+str(endTime-startTime)+" seconds")




    # Open regulator ID file and predict operators for each
with open("cache/ids/reg_ids.txt", "r") as f:
    regs = [r.strip() for r in f.readlines() if r[0] != ">"]

for r in regs:
    print(r)
    predict_operator(r)
    operator_analysis(r)
