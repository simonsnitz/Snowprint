import pickle
from pprint import pprint
import re
import json
import io


def create_operator_html(operator, regulator_name):
    
    #graphic = [ createGraphic(i["operon"], i["regIndex"]) for i in operons ]
    
    operator = json.dumps(operator)

    graphic_js_variable = "var graphic = '["+str(operator)[1:-1]+"]'"


    newHTML = io.open(regulator_name + ".html", "w")

	#make a copy of template file
    for line in io.open("display/operator_display_template.html", 'r'):
        newHTML.write(line)

    append2HTML = """

            // new section appended below //

            """+graphic_js_variable+"""

            var data = JSON.parse(graphic)
                renderOperator(data)


    </script>

    </html>
    """

    newHTML.write(append2HTML)

    newHTML.close()

if __name__ == "__main__":

    inputFile = ""
    regulator_name = "alkx"

    create_operator_html(inputFile, regulator_name)


