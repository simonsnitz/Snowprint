import search2accessionList as s2accs
import acc2operon
from display import test_operon2graphic as o2g
import sys

"""
    This function takes in a search term and # of max regulators to limit search.
    Returns an HTML page with a list of operons containing at least 1 regulator and a gene that fits the search term
"""

search_term = sys.argv[1]

max_regulators = sys.argv[2]



application = "create_graphic"

accList4graphic = s2accs.searchTerms2ACCs(search_term, application)

operons_with_regulators = acc2operon.enzyme_acc2regulator(accList4graphic, max_regulators)

o2g.create_html(operons_with_regulators, "display/html_pages/"+str(search_term).replace(" ","_"))

