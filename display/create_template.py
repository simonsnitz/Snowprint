import sys
import io

newHTML = io.open(sys.argv[1] + ".html", "w")

for line in io.open("template_displayOperon.html", 'r'):
    newHTML.write(line)

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
