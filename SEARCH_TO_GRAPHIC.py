import search2accessionList as s2accs
import acc2operon
from display import test_operon2graphic as o2g
import sys



search_term = sys.argv[1]

max_regulators = sys.argv[2]



application = "create_graphic"

accList4graphic = s2accs.searchTerms2ACCs(search_term, application)

operons_with_regulators = acc2operon.enzyme_acc2regulator(accList4graphic, max_regulators)

o2g.create_html(operons_with_regulators, "display/html_pages/"+str(search_term).replace(" ","_"))

