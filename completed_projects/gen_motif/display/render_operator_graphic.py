import pickle
from pprint import pprint
import re
import json
import io


def create_operator_html(operators, html_file):
    
    #graphic = [ createGraphic(i["operon"], i["regIndex"]) for i in operons ]
    
    operators = [ json.dumps(i) for i in operators ]

    graphic_js_variable = "var operator_graphic = ["+str(operators)[1:-1]+"]"


    #HTML_FILE = io.open( html_file, "w")

    with open(html_file, "a") as HTML_FILE:


        append2HTML = """

            // operator variable appended below //

            """+graphic_js_variable+"""


            //Deals with inputting list of consensus motif data with different percent identity cutoffs

            for (i in range(0, operator_graphic.length)){
                var operator_data = JSON.parse(operator_graphic[i]);
                for (j in range(0, operator_data.length)){
                    renderOperator(operator_data[j]);
                }
            }

        </script>

        </html>
        """

        HTML_FILE.write(append2HTML)

if __name__ == "__main__":

    inputFile = ""
    regulator_name = "alkx"

    create_operator_html(inputFile, regulator_name)


