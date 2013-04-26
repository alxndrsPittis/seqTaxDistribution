import os
import sqlite3
from sys import exit


#link = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz"

#os.system("wget "+link)
#os.system("gunzip -f idmapping_selected.tab.gz")
#print "Selecting columns from Uniprot idmapping_selected.tab..."
#os.system("cut -f 1,2,14 idmapping_selected.tab > idmapping_selected.3tab")


CMD = open("commands.tmp", "w")
cmd = """
DROP TABLE IF EXISTS idmapping;
CREATE TABLE idmapping(ac VARCHAR(50) PRIMARY KEY, id VARCHAR(50), taxon INT);
.separator "\t"
.import idmapping_selected.3tab idmapping
"""

CMD.write(cmd)
CMD.close()

#OUT = open("commands.file", "w")
# 
#print >>OUT, 'delete from maps;'
#print >>OUT, '.separator "\\t"'
#print >>OUT, '.import idmapping_selected.3tab maps'
# 
#OUT.close()
print "Updating idmapping database..."
os.system("sqlite3 idmapping.sqlite < commands.tmp")
os.system("rm commands.tmp")
print "idmapping updated."
