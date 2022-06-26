#! python3

import cgitb
cgitb.enable()

# below two lines are compulsory, don't delete them
print("Content-Type: text/html")
print()

# don't change its order
import os
import sys
import cgi
import uuid
import operations


# getting user inputs
form = cgi.FieldStorage()
gene_name = form.getvalue('genename')
species_name = form.getvalue('speciesname')
is_direct = form.getvalue('structure')


# input validation
if not gene_name and species_name:
    print("<h1>Please enter all the parameters!!!</h1>")
    sys.exit(1)

# operate
try:
    protein_2d_list = operations.get_top_two_protein_matches(gene_name,species_name, is_direct)
    alignment = operations.do_pairwise_alignment(protein_2d_list[0], protein_2d_list[1])
    operations.do_protein_superimposing([protein_2d_list[0][0], protein_2d_list[1][0]], is_direct)
except Exception as e:
    print(e)
