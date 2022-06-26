#!/usr/bin/env python3
import cgi;


import os
import sys
import logging
import Bio.PDB
import Bio.PDB.PDBParser
from Bio.PDB import PDBList, MMCIFParser
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO
#import nglview as nv


from Bio.Blast.Applications import NcbiblastpCommandline
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


Entrez.email = "ahmed@mu.edu.tr"


def get_top_two_protein_matches(gene, species, is_direct):
    """
    -expects gene name and species name
    -finds the relevent proteins with direct structure
    -returns the top two proteins [[pid1,pseq1],[pid2, pseq2]]
    """
    try:
        
        query = f"{gene} AND {species}[porgn] AND protein "
        
        if(is_direct == 'true'):
            query += "structure direct[Filter]"
        else:
            query += "structure[Filter]"
            
        handle=Entrez.esearch(db="protein",term=query, retmax=200, sort="relevance" ) #
        esearch_result=Entrez.read(handle)
      
        efetch_results = []
        handle=Entrez.efetch(db="protein",
                              id= esearch_result["IdList"][:2],
                              rettype="gb")

        for r in SeqIO.parse(handle,"genbank"):
            efetch_results.append([r.id, r.seq]) 
            seq_arr = [ str(r.seq[i:i+60]) for i in range(0, len(r.seq), 60) ]
            print(f"""
                    <pre style="">
                    {' - ' * 20}
                    {' - ' * 8}Protein Info{' - ' * 8}
                    {' - ' * 20}
                    Id:           <span style="color:#fa7268">{r.id}</span>
                    Description:  <span style="color:#fa7268">{r.description}</span>
                    Sequence:</pre>
                    
            """)
            for line in seq_arr:
                    print(f"""
                    <pre style="color:red">
                    {line}</pre>
                    """)
                
           
        return efetch_results
    except Exception as e:
        print("Error occurred while finding top proteins match")
        print(e)
        sys.exit(1)
        
        
def do_pairwise_alignment(protein1, protein2):
    try:
        # Create two sequence files
        seq1 = SeqRecord(protein1[1],
                           id=protein1[0])
        seq2 = SeqRecord(protein2[1],
                           id=protein2[0])

        seq1_file = f"{protein1[0]}.fasta"
        seq2_file = f"{protein2[0]}.fasta"
        SeqIO.write(seq1, seq1_file, "fasta")
        SeqIO.write(seq2, seq2_file, "fasta")

        # Run BLAST and parse the output as XML
        output = NcbiblastpCommandline(query=seq1_file, subject=seq2_file, outfmt=5)()[0]
        blast_result_record = NCBIXML.read(StringIO(output))
        
        hsps_arr = []
        
        # Print some information on the result
        for alignment in blast_result_record.alignments:
            ihsp = 0
            for hsp in alignment.hsps:
                hsps_arr.append(hsp)
   
                print(f"""
                    <pre>
                    {' - ' * 20}
                    {' - ' * 6}Blast Pairwise Alignment - {ihsp+1}{' - ' * 6}
                    {' - ' * 20}
                    Subject:  <span style="color:#fa7268">{protein2[0]}</span>
                    Query:    <span style="color:#fa7268">{protein1[0]}</span>
                    Length:   <span style="color:#fa7268">{alignment.length}</span>
                    E value:  <span style="color:#fa7268">{hsp.expect}</span>
                    Identity: <span style="color:#fa7268">{hsp.identities}</span>
                    {' - ' * 20}
                    </pre>
                """)

                chunk_size =  60
                query_arr = [ hsp.query[i:i+chunk_size] for i in range(0, len(hsp.query), chunk_size) ]
                subject_arr = [ hsp.sbjct[i:i+chunk_size] for i in range(0, len(hsp.sbjct), chunk_size) ]
                match_arr = [ hsp.match[i:i+chunk_size] for i in range(0, len(hsp.match), chunk_size) ]
                ihsp+=1
                for i in range(len(subject_arr)):
                    print(f"""
                    <pre>
                    <span style="color:red">{subject_arr[i]}</span>
                    <span style="color:green">{match_arr[i]}</span>
                    <span style="color:blue">{query_arr[i]}</span>
                    </pre>
                    """)
                
                
        return hsps_arr
        
    except Exception as e:
        print("Error occurred while aligning the proteins")
        print(e)
        #sys.exit(1)
        

def retrieve_pdb_file(pdb_ref, file_type):
    
    pdbl = PDBList(verbose=False)
    file_path = pdbl.retrieve_pdb_file(pdb_ref, pdir=f"./", file_format=file_type, overwrite=True)
    
    return file_path

    
def store_structure_as_image(pdb_file):
    try:
        import __main__
        __main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
        import pymol
        pymol.finish_launching()
        
        pdb_name =pdb_file.split('.')[0]
        pymol.cmd.feedback("disable","all","everything")
        pymol.cmd.load(pdb_file, pdb_name)
        
        
        
        #pymol.cmd.disable("all")
        pymol.cmd.enable(pdb_name)
        
        #pymol.cmd.hide('all')
        #pymol.cmd.show('cartoon')
        #pymol.cmd.set('ray_opaque_background', 0)
        
        pymol.cmd.png("%s.png"%(pdb_name))
        
        pymol.cmd.quit()
        
        return 0
    except Exception as e:
        print(" | Error: Unable to find pymol module | ")
        return 1

   
def do_protein_superimposing_v2(ids, is_direct=False):
    try:
        #print(ids)
        handle = Entrez.elink(dbfrom="protein", id=ids, linkname="protein_structure", dbto='structure')
        record = Entrez.read(handle)
        handle.close()

        if(len(record[0]["LinkSetDb"][0]["Link"]) <= 1):
            return

        struct_1_id = record[0]["LinkSetDb"][0]["Link"][0]['Id']
        struct_2_id = record[0]["LinkSetDb"][0]["Link"][1]['Id']

        #print(f"{struct_1_id} --- {struct_2_id}")

        handle = Entrez.esummary(db="structure", id=f"{struct_1_id},{struct_2_id}")
        record = Entrez.read(handle)
        struct_1_pdb_id, struct_2_pdb_id = record[0]['PdbAcc'], record[1]['PdbAcc']

        #print(f"{struct_1_pdb_id} --- {struct_2_pdb_id}")

        print(f"""
            <pre>
            {' - ' * 20}
            {' - ' * 6}Structure Superimposition{' - ' * 6}
            {' - ' * 20}

            PDB Structure 1:  <span style="color:#fa7268">{struct_1_pdb_id}\t|\t{record[0]['Id']}\t|\t{record[0]['PdbDescr']}</span>
            PDB Structure 2:  <span style="color:#fa7268">{struct_2_pdb_id}\t|\t{record[1]['Id']}\t|\t{record[1]['PdbDescr']}</span>

            </pre>
        """)

        struct_1_path = retrieve_pdb_file(struct_1_pdb_id, file_type='pdb')
        struct_2_path = retrieve_pdb_file(struct_2_pdb_id, file_type='pdb')

        #print(f"{struct_1_path} --- {struct_2_path}")
        parser = Bio.PDB.PDBParser(QUIET=True)
        structure1 = parser.get_structure(struct_1_pdb_id,struct_1_path)
        structure2 = parser.get_structure(struct_2_pdb_id,struct_2_path)


        # SUPERIMPOSE
        import __main__
        __main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
        import pymol
        pymol.finish_launching()
        pymol.cmd.feedback("disable","all","everything")


        prot1 = os.path.basename(struct_1_path).split('.')[0]
        prot2 = os.path.basename(struct_2_path).split('.')[0]


        pymol.cmd.load(struct_1_path, prot1)
        pymol.cmd.load(struct_2_path, prot2)

        pymol.cmd.enable(prot1)
        pymol.cmd.enable(prot2)


        align_obj = f"{prot1}-{prot2}"
        align1=pymol.cmd.align(prot2, prot1, object=align_obj)
        pymol.cmd.disable(align_obj)
        # PART 2 ----

        #PDB-to-PNG
        out_png_file = f"{prot1}-{prot2}-super.png"
        pymol.cmd.png(out_png_file)
        pymol.cmd.quit()

        # SHOW PNG
        print(f"<img src='./{out_png_file}'/>")
    except Exception as e:
        print(str(e))

    
# DEPRECATED
def do_protein_superimposing(ids, is_direct):
    try:
        #print(ids)
        handle = Entrez.elink(dbfrom="protein", id=ids, linkname="protein_structure", dbto='structure')
        record = Entrez.read(handle)
        handle.close()

        if(len(record[0]["LinkSetDb"][0]["Link"]) <= 1):
            return 
            
        struct_1_id = record[0]["LinkSetDb"][0]["Link"][0]['Id']
        struct_2_id = record[0]["LinkSetDb"][0]["Link"][1]['Id']
        
        
         
        handle = Entrez.esummary(db="structure", id=f"{struct_1_id},{struct_2_id}")
        record = Entrez.read(handle)
        struct_1_pdb_id, struct_2_pdb_id = record[0]['PdbAcc'], record[1]['PdbAcc'] 
        
        

        print(f"""
            <pre>
            {' - ' * 20}
            {' - ' * 6}Structure Superimposition{' - ' * 6}
            {' - ' * 20}
            
            PDB Structure 1:  <span style="color:#fa7268">{struct_1_pdb_id}\t|\t{record[0]['Id']}\t|\t{record[0]['PdbDescr']}</span>
            PDB Structure 2:  <span style="color:#fa7268">{struct_2_pdb_id}\t|\t{record[1]['Id']}\t|\t{record[1]['PdbDescr']}</span>
                        
            </pre>
        """)
        
        struct_1_path = retrieve_pdb_file(struct_1_pdb_id, file_type='pdb')
        struct_2_path = retrieve_pdb_file(struct_2_pdb_id, file_type='pdb')


        parser = Bio.PDB.PDBParser(QUIET=True)
        structure1 = parser.get_structure(struct_1_pdb_id,struct_1_path)
        structure2 = parser.get_structure(struct_2_pdb_id,struct_2_path)
        
        
        ppb=Bio.PDB.PPBuilder()
        query = ppb.build_peptides(structure1)[0]
        target = ppb.build_peptides(structure2)[0]
        
        _max = min(len(query), len(target)) - 1
        query_atoms = [ r['CA'] for r in query ][:_max]
        target_atoms = [ r['CA'] for r in target ][:_max]
        
        superimposer = Bio.PDB.Superimposer()
        superimposer.set_atoms(query_atoms, target_atoms)

        superimposer.apply(structure2.get_atoms())

        # Write modified structures to one file
        #outfile_path = "X-modified.pdb"
        #outfile=open(outfile_path, "w") 
        #io=Bio.PDB.PDBIO() 
        #io.set_structure(structure2) 
        #io.save(outfile) 
        #outfile.close()
        
        pdb_io = Bio.PDB.PDBIO()
        outfile_path = f'{struct_1_pdb_id}_{struct_2_pdb_id}_proteins_struc.pdb'
        pdb_file = open(outfile_path, 'w')
        for struct in [structure1, structure2]:
                pdb_io.set_structure(struct[0])
                pdb_io.save(pdb_file)
        pdb_file.close()
        
        
        
        store_status = store_structure_as_image(outfile_path)
                    
        if(store_status == 0): # success
            print("<img src='./%s.png'/>"%(outfile_path.split('.')[0]))
        else:
            print(f'<p style="color:red;font-weight:bold"> [!] To view the 3D superimposed structures; open {outfile_path} in Pymol</>')
            #print('<p style="color:blue;font-weight:600">Static Image of 7P6D and 7P6E protein Structures</p>')
            if(is_direct == 'true'):
                print("<img src='./superstructure.png' height='350px' width='850px'/>")
            else:
                print("<img src='./6XLI_6PXR_proteins_struc.png' height='350px' width='850px'/>")
        #
        #os.remove(outfile_path)
        #os.remove(struct_1_path)
        #os.remove(struct_2_path)
    except Exception as e:
        print(str(e))
    
        
