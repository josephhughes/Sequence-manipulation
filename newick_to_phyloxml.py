import sys
import re
from StringIO import StringIO

from ete2 import Phyloxml, phyloxml

#Creates empty phyloxml document
project = Phyloxml()

# Loads newick tree
phylo = phyloxml.PhyloxmlTree(newick=sys.argv[1])

# Set basic tree info as a phyloxml phylogeny object
phylo.phyloxml_phylogeny.set_name("test_tree")
if len(phylo.children) <= 2:
    phylo.phyloxml_phylogeny.set_rooted("true")
else:
    phylo.phyloxml_phylogeny.set_rooted("false")
    
# Add the tree to the phyloxml project
project.add_phylogeny(phylo)

# Export phyloxml document
OUTPUT = StringIO()
project.export(OUTPUT)

# Some ad-hoc changes to the phyloxml formatted document to meet the schema definition
text = OUTPUT.getvalue()
text = text.replace("phy:", "")
text = re.sub('branch_length_attr="[^"]+"', "", text)
header = """
 <phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
 xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd"
 xmlns="http://www.phyloxml.org">
"""
text = re.sub('<Phyloxml[^>]+>', header, text)
text = text.replace('Phyloxml', 'phyloxml')

# Save result 
open(sys.argv[1]+".phyloxml", "w").write(text)
