import pickle
from pprint import pprint
import re
import json
import io


def create_operator_html(operators, regulator_name):
    
    #graphic = [ createGraphic(i["operon"], i["regIndex"]) for i in operons ]
    
    operators = [ json.dumps(i) for i in operators ]

    graphic_js_variable = "var graphic = ["+str(operators)[1:-1]+"]"


    newHTML = io.open(regulator_name + ".html", "w")

	#make a copy of template file
    for line in io.open("display/operator_template.html", 'r'):
        newHTML.write(line)

    append2HTML = """

            // new section appended below //

            """+graphic_js_variable+"""


            //Deals with inputting list of consensus motif data with different percent identity cutoffs

            for (i in range(0, graphic.length)){
                var data = JSON.parse(graphic[i]);
                for (j in range(0, data.length)){
                    renderOperator(data[j]);
                }
            }

    </script>

    </html>
    """

    newHTML.write(append2HTML)

    newHTML.close()

if __name__ == "__main__":

    inputFile = ""
    regulator_name = "alkx"

    create_operator_html(inputFile, regulator_name)


