"""
  This is a prototype subsection for the 'test_operon2graphic.py' script.
  Everything in here was pasted into the final code of that file.
  
  This script takes in the name of a JSON file(?) and outputs an HTML page that displays a list of operons.
"""

import sys
import io

    #creates HTML file
newHTML = io.open(sys.argv[1] + ".html", "w")

    #copies code from template HTML file
for line in io.open("template_displayOperon.html", 'r'):
    newHTML.write(line)

    #appends new section with 'operon-specific data'?
append2HTML = """

        // new section appended below //

        var len = graphic.length
        for (i in range(0, len)){
            var data = JSON.parse(graphic[i])
            renderOperon(data)
        }


</script>

</html>
"""

newHTML.write(append2HTML)

newHTML.close()
