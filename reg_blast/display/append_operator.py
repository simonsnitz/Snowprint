from pathlib import Path
import json


def create_operator_html(operators):
    
    
    operators = [ "["+json.dumps(i)+"]" for i in operators ]

    graphic_js_variable = "var operator_graphic = ["+str(operators)[1:-1]+"]"



    JS_string = """

            // Operator object //


            """+graphic_js_variable+"""


            //Deals with inputting list of consensus motif data with different percent identity cutoffs

            for (i in range(0, operator_graphic.length)){
                var operator_data = JSON.parse(operator_graphic[i]);
                for (j in range(0, operator_data.length)){
                    renderOperator(operator_data[j]);
                }
            }

        """

        # TODO: CHANGE THIS
    p = Path("/home/simonsnitz/projects/tuSeek/reg_blast/display/Webpage")
    fp = p / "operator_data.js"


       # Overwrite javascript file.
    with open(str(fp), mode="w") as f:
        f.write(JS_string)



