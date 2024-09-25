'''
brew install pymol

conda create --name proteinfolding
conda activate proteinfolding

conda install -c conda-forge python ipython
pip install reqests biopython pytest


git init



'''


uniprot_id = "P0DOX5"













# import pymol
# from pymol import cmd

# # Initialize PyMOL
# pymol.finish_launching()

# # Load the PDB structures for different species
# species_files = ["species1.pdb", "species2.pdb", "species3.pdb"]
# for idx, species_file in enumerate(species_files):
#     cmd.load(species_file, f"species_{idx}")

# # Align structures to the first one (species_0)
# for idx in range(1, len(species_files)):
#     cmd.align(f"species_{idx}", "species_0")

# # Calculate RMSD and display differences
# cmd.color("blue", "species_0")
# for idx in range(1, len(species_files)):
#     cmd.color("red", f"species_{idx}")
#     rmsd = cmd.rms_cur(f"species_{idx}", "species_0")
#     print(f"RMSD between species_{idx} and species_0: {rmsd}")

# # Save the aligned structures for visualization
# cmd.save("aligned_structures.pse")



# from Bio import ExPASy, SwissProt, PDB
# import requests

# def fetch_pdb_from_uniprot(uniprot_id):
#     uniprot_url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
#     response = requests.get(uniprot_url)
#     pdb_ids = []
    
#     if response.status_code == 200:
#         from xml.etree import ElementTree as ET
#         tree = ET.fromstring(response.content)
#         for ref in tree.findall(".//dbReference[@type='PDB']"):
#             pdb_id = ref.get('id')
#             pdb_ids.append(pdb_id)
    
#     return pdb_ids


# uniprot_id = "P0DOX5"
# fetch_pdb_from_uniprot(uniprot_id)






## source: https://stackoverflow.com/questions/71600329/get-pdb-id-chain-id-from-uniprot-id

import urllib.parse
import urllib.request

url = 'https://www.uniprot.org/uploadlists/'

params = {
'from': 'ACC+ID',
'to': 'PDB_ID',
'format': 'tab',
'query': uniprot_id
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with open('UniProt_PDB_IDs.txt', 'a') as f:
   with urllib.request.urlopen(req) as q:
      response = q.read()
      f.write(response.decode('utf-8'))








# https://labstructbioinf.github.io/localpdb/

# !pip install localpdb

# Setup the database and sync protein structures in the mmCIF format:

# !localpdb_setup -db_path ~/projects/python/proteinfolding/my_local_pdb --fetch_cif









# https://www.blopig.com/blog/2024/02/working-with-pdb-structures-in-pandas/


# https://www.uniprot.org/help/retrieve_3d
# https://www.uniprot.org/uniprotkb?query=%28database%3Apdb%29+AND+%28reviewed%3Atrue%29+AND+P0DOX5



import requests

# Define the UniProt accession number (in this case P0DOX5)
accession = "P0DOX5"

# Make a request to UniProt's API to get the PDB ID
url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
response = requests.get(url)



######################################################
# Check if the request was successful




# Extract the PDB IDs from the UniProt cross-references
if response.status_code == 200:
    data = response.json()
    print(f"PDB IDs obtained for {accession}")
else:
    print(f"Failed to retrieve data for {accession}")

pdb_entries = [ref['id'] for ref in data['uniProtKBCrossReferences'] if ref['database'] == 'PDB']

print(f"PDB IDs: {pdb_entries}")

'''
PDB IDs: ['1N0X', '3PGF', '4R2G', '5O4E', '5VJ6', '5VU0', '5VZX', '5VZY', '5W5L', '5WAV', '5XJE', '5XJF', '5XMH', '5Y56', '5YC5', '6APD', '6ARP', '6ARU', '6B70', '6B7Z', '6BF7', '6BF9', '6BFT', '6BGT', '6BKB', '6BKC', '6BZ4', '6DKJ', '6EAQ', '6FCZ', '6FGO', '6G1E', '6IFJ', '6IQG', '6IQH', '6KA7', '6MB3', '6MSY', '6MU3', '6MUB', '6N2X', '6N32', '6N35', '6OGE', '6OKQ', '6P6D', '6UBI', '6UGW', '6UGX', '6UGY', '6UOE', '6V8Z', '6VSL', '6VSZ', '6X3I', '6YSC', '6YT7', '6YTB', '7CZT', '7CZU', '7CZV', '7T17', '7URU', '7X13', '7XXL', '8A49', '8A64', '8DAO', '8DBZ', '8DV1', '8DV2', '8ECQ', '8ECV', '8ECZ', '8ED1', '8EDF', '8GHR']
'''


for pdb_id in pdb_entries:
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(pdb_url)
    
    if response.status_code == 200:
        with open(f"{pdb_id}.pdb", 'wb') as file:
            file.write(response.content)
        print(f"PDB file {pdb_id}.pdb downloaded successfully.")
    else:
        print(f"Failed to download PDB file for {pdb_id}")
'''

PDB file 1N0X.pdb downloaded successfully.
PDB file 3PGF.pdb downloaded successfully.
PDB file 4R2G.pdb downloaded successfully.
PDB file 5O4E.pdb downloaded successfully.
PDB file 5VJ6.pdb downloaded successfully.
PDB file 5VU0.pdb downloaded successfully.
PDB file 5VZX.pdb downloaded successfully.
PDB file 5VZY.pdb downloaded successfully.
PDB file 5W5L.pdb downloaded successfully.
PDB file 5WAV.pdb downloaded successfully.
PDB file 5XJE.pdb downloaded successfully.
PDB file 5XJF.pdb downloaded successfully.
PDB file 5XMH.pdb downloaded successfully.
PDB file 5Y56.pdb downloaded successfully.
PDB file 5YC5.pdb downloaded successfully.
PDB file 6APD.pdb downloaded successfully.
PDB file 6ARP.pdb downloaded successfully.
PDB file 6ARU.pdb downloaded successfully.
PDB file 6B70.pdb downloaded successfully.
PDB file 6B7Z.pdb downloaded successfully.
PDB file 6BF7.pdb downloaded successfully.
PDB file 6BF9.pdb downloaded successfully.
PDB file 6BFT.pdb downloaded successfully.
PDB file 6BGT.pdb downloaded successfully.
PDB file 6BKB.pdb downloaded successfully.
PDB file 6BKC.pdb downloaded successfully.
PDB file 6BZ4.pdb downloaded successfully.
PDB file 6DKJ.pdb downloaded successfully.
PDB file 6EAQ.pdb downloaded successfully.
PDB file 6FCZ.pdb downloaded successfully.
PDB file 6FGO.pdb downloaded successfully.
PDB file 6G1E.pdb downloaded successfully.
PDB file 6IFJ.pdb downloaded successfully.
PDB file 6IQG.pdb downloaded successfully.
Failed to download PDB file for 6IQH
PDB file 6KA7.pdb downloaded successfully.
PDB file 6MB3.pdb downloaded successfully.
PDB file 6MSY.pdb downloaded successfully.
PDB file 6MU3.pdb downloaded successfully.
PDB file 6MUB.pdb downloaded successfully.
PDB file 6N2X.pdb downloaded successfully.
PDB file 6N32.pdb downloaded successfully.
PDB file 6N35.pdb downloaded successfully.
PDB file 6OGE.pdb downloaded successfully.
PDB file 6OKQ.pdb downloaded successfully.
PDB file 6P6D.pdb downloaded successfully.
PDB file 6UBI.pdb downloaded successfully.
PDB file 6UGW.pdb downloaded successfully.
PDB file 6UGX.pdb downloaded successfully.
PDB file 6UGY.pdb downloaded successfully.
PDB file 6UOE.pdb downloaded successfully.
PDB file 6V8Z.pdb downloaded successfully.
PDB file 6VSL.pdb downloaded successfully.
PDB file 6VSZ.pdb downloaded successfully.
PDB file 6X3I.pdb downloaded successfully.
PDB file 6YSC.pdb downloaded successfully.
PDB file 6YT7.pdb downloaded successfully.
PDB file 6YTB.pdb downloaded successfully.
PDB file 7CZT.pdb downloaded successfully.
PDB file 7CZU.pdb downloaded successfully.
PDB file 7CZV.pdb downloaded successfully.
PDB file 7T17.pdb downloaded successfully.
PDB file 7URU.pdb downloaded successfully.
PDB file 7X13.pdb downloaded successfully.
PDB file 7XXL.pdb downloaded successfully.
PDB file 8A49.pdb downloaded successfully.
Failed to download PDB file for 8A64
PDB file 8DAO.pdb downloaded successfully.
PDB file 8DBZ.pdb downloaded successfully.
PDB file 8DV1.pdb downloaded successfully.
PDB file 8DV2.pdb downloaded successfully.
PDB file 8ECQ.pdb downloaded successfully.
PDB file 8ECV.pdb downloaded successfully.
PDB file 8ECZ.pdb downloaded successfully.
PDB file 8ED1.pdb downloaded successfully.
PDB file 8EDF.pdb downloaded successfully.
PDB file 8GHR.pdb downloaded successfully.
'''

# filter for successfully donwloaded
orig_pdb_entries = pdb_entries[:]
pdb_entries = [x for x in pdb_entries if os.path.exists(f"{pdb_directory}{pdb}.pdb")]

import os
from Bio.PDB import PDBParser

def analyze_pdb(pdb_file):
    # Check if the PDB file exists
    if not os.path.exists(pdb_file):
        print(f"File not found: {pdb_file}")
        return

    # Create a PDB parser
    parser = PDBParser(QUIET=True)
    
    # Parse the structure
    structure = parser.get_structure("protein", pdb_file)
    
    for model in structure:
        for chain in model:
            residues = [residue for residue in chain.get_residues()]
            if residues:
                start_res = residues[0].get_id()[1]  # Get start residue
                end_res = residues[-1].get_id()[1]  # Get end residue
                print(f"Chain {chain.id} in {pdb_file} spans residues {start_res}-{end_res}")


for pdb in pdb_entries:
    analyze_pdb(f"data/{pdb}.pdb")
'''
Chain L in data/1N0X.pdb spans residues 1-963
Chain H in data/1N0X.pdb spans residues 1-1388
Chain M in data/1N0X.pdb spans residues 1-976
Chain K in data/1N0X.pdb spans residues 1-1083
Chain P in data/1N0X.pdb spans residues 1-732
Chain R in data/1N0X.pdb spans residues 1-917
Chain A in data/3PGF.pdb spans residues 1-559
Chain H in data/3PGF.pdb spans residues 2-399
Chain L in data/3PGF.pdb spans residues 1-395
Chain B in data/3PGF.pdb spans residues 1-2
Chain E in data/4R2G.pdb spans residues 92-519
Chain F in data/4R2G.pdb spans residues 1-175
Chain P in data/4R2G.pdb spans residues 5-210
Chain Q in data/4R2G.pdb spans residues 2-301
Chain B in data/4R2G.pdb spans residues 1-175
Chain C in data/4R2G.pdb spans residues 5-210
Chain D in data/4R2G.pdb spans residues 2-301
Chain H in data/4R2G.pdb spans residues 1-176
Chain I in data/4R2G.pdb spans residues 5-210
Chain J in data/4R2G.pdb spans residues 2-301
Chain L in data/4R2G.pdb spans residues 3-175
Chain M in data/4R2G.pdb spans residues 5-210
Chain N in data/4R2G.pdb spans residues 1-301
Chain O in data/4R2G.pdb spans residues 92-518
Chain A in data/4R2G.pdb spans residues 92-521
Chain K in data/4R2G.pdb spans residues 92-516
Chain G in data/4R2G.pdb spans residues 1-10
Chain R in data/4R2G.pdb spans residues 1-10
Chain S in data/4R2G.pdb spans residues 1-10
Chain T in data/4R2G.pdb spans residues 1-2
Chain U in data/4R2G.pdb spans residues 1-2
Chain V in data/4R2G.pdb spans residues 1-10
Chain A in data/5O4E.pdb spans residues 236-659
Chain B in data/5O4E.pdb spans residues 236-652
Chain C in data/5O4E.pdb spans residues 240-640
Chain D in data/5O4E.pdb spans residues 237-665
Chain E in data/5O4E.pdb spans residues 13-315
Chain F in data/5O4E.pdb spans residues 13-311
Chain G in data/5O4E.pdb spans residues 1-8
Chain H in data/5O4E.pdb spans residues 1-8
Chain I in data/5O4E.pdb spans residues 1-2
Chain J in data/5O4E.pdb spans residues 1-6
Chain A in data/5VJ6.pdb spans residues 518-664
Chain B in data/5VJ6.pdb spans residues 518-664
Chain C in data/5VJ6.pdb spans residues 518-664
Chain D in data/5VJ6.pdb spans residues 33-505
Chain E in data/5VJ6.pdb spans residues 33-505
Chain F in data/5VJ6.pdb spans residues 33-505
Chain H in data/5VJ6.pdb spans residues 3-233
Chain M in data/5VJ6.pdb spans residues 1-213
Chain N in data/5VJ6.pdb spans residues 1-213
Chain O in data/5VJ6.pdb spans residues 1-213
Chain P in data/5VJ6.pdb spans residues 1-213
Chain Q in data/5VJ6.pdb spans residues 1-213
Chain R in data/5VJ6.pdb spans residues 1-213
Chain L in data/5VJ6.pdb spans residues 2-208
Chain A in data/5VU0.pdb spans residues 228-1213
Chain B in data/5VU0.pdb spans residues 229-1188
Chain C in data/5VU0.pdb spans residues 3-353
Chain D in data/5VU0.pdb spans residues 1-8
Chain E in data/5VU0.pdb spans residues 1-8
Chain H in data/5VZX.pdb spans residues 1-429
Chain L in data/5VZX.pdb spans residues 1-442
Chain E in data/5VZX.pdb spans residues 1-426
Chain I in data/5VZX.pdb spans residues 1-437
Chain H in data/5VZY.pdb spans residues 1-355
Chain L in data/5VZY.pdb spans residues 1-379
Chain A in data/5VZY.pdb spans residues 13-108
Chain A in data/5W5L.pdb spans residues 236-705
Chain B in data/5W5L.pdb spans residues 235-755
Chain C in data/5W5L.pdb spans residues 1-9
Chain D in data/5W5L.pdb spans residues 1-9
Chain A in data/5WAV.pdb spans residues 237-681
Chain B in data/5WAV.pdb spans residues 239-650
Chain C in data/5WAV.pdb spans residues 1-8
Chain D in data/5WAV.pdb spans residues 1-8
Chain A in data/5XJE.pdb spans residues 232-1137
Chain B in data/5XJE.pdb spans residues 230-1159
Chain C in data/5XJE.pdb spans residues 10-1113
Chain D in data/5XJE.pdb spans residues 1-9
Chain E in data/5XJE.pdb spans residues 1-9
Chain F in data/5XJE.pdb spans residues 1-2
Chain G in data/5XJE.pdb spans residues 1-3
Chain A in data/5XJF.pdb spans residues 230-1139
Chain B in data/5XJF.pdb spans residues 230-1155
Chain C in data/5XJF.pdb spans residues 10-1115
Chain D in data/5XJF.pdb spans residues 1-9
Chain E in data/5XJF.pdb spans residues 1-9
Chain F in data/5XJF.pdb spans residues 1-5
Chain A in data/5XMH.pdb spans residues 237-619
Chain B in data/5XMH.pdb spans residues 238-609
Chain L in data/5XMH.pdb spans residues 1-318
Chain H in data/5XMH.pdb spans residues 1-313
Chain D in data/5XMH.pdb spans residues 1-302
Chain C in data/5XMH.pdb spans residues 1-303
Chain E in data/5XMH.pdb spans residues 1-8
Chain F in data/5XMH.pdb spans residues 1-6
Chain A in data/5Y56.pdb spans residues 236-659
Chain B in data/5Y56.pdb spans residues 236-636
Chain C in data/5Y56.pdb spans residues 1-9
Chain D in data/5Y56.pdb spans residues 1-7
Chain A in data/5YC5.pdb spans residues 234-1118
Chain B in data/5YC5.pdb spans residues 232-1110
Chain C in data/5YC5.pdb spans residues 4-207
Chain D in data/5YC5.pdb spans residues 1-9
Chain E in data/5YC5.pdb spans residues 1-8
Chain A in data/6APD.pdb spans residues 27-515
Chain B in data/6APD.pdb spans residues 27-517
Chain C in data/6APD.pdb spans residues 27-517
Chain D in data/6APD.pdb spans residues 1-215
Chain E in data/6APD.pdb spans residues 1-212
Chain F in data/6APD.pdb spans residues 1-216
Chain G in data/6APD.pdb spans residues 1-214
Chain H in data/6APD.pdb spans residues 1-215
Chain I in data/6APD.pdb spans residues 1-211
Chain J in data/6APD.pdb spans residues 1-214
Chain K in data/6APD.pdb spans residues 1-214
Chain L in data/6APD.pdb spans residues 2-209
Chain M in data/6APD.pdb spans residues 2-209
Chain N in data/6APD.pdb spans residues 1-214
Chain O in data/6APD.pdb spans residues 2-209
Chain C in data/6ARP.pdb spans residues 1-661
Chain D in data/6ARP.pdb spans residues 1-590
Chain A in data/6ARP.pdb spans residues 1-642
Chain B in data/6ARP.pdb spans residues 1-563
Chain A in data/6ARU.pdb spans residues 4-710
Chain B in data/6ARU.pdb spans residues 1-210
Chain C in data/6ARU.pdb spans residues 1-301
Chain D in data/6ARU.pdb spans residues 1-8
Chain E in data/6ARU.pdb spans residues 1-2
Chain A in data/6B70.pdb spans residues 46-1011
Chain B in data/6B70.pdb spans residues 46-1011
Chain C in data/6B70.pdb spans residues 4-221
Chain D in data/6B70.pdb spans residues 2-212
Chain E in data/6B70.pdb spans residues 4-221
Chain F in data/6B70.pdb spans residues 2-212
Chain a in data/6B70.pdb spans residues 1-21
Chain c in data/6B70.pdb spans residues 1-15
Chain A in data/6B7Z.pdb spans residues 46-1011
Chain B in data/6B7Z.pdb spans residues 46-1011
Chain C in data/6B7Z.pdb spans residues 4-221
Chain D in data/6B7Z.pdb spans residues 2-212
Chain E in data/6B7Z.pdb spans residues 4-221
Chain F in data/6B7Z.pdb spans residues 2-212
Chain A in data/6BF7.pdb spans residues 46-1011
Chain B in data/6BF7.pdb spans residues 46-1011
Chain C in data/6BF7.pdb spans residues 4-221
Chain D in data/6BF7.pdb spans residues 2-212
Chain E in data/6BF7.pdb spans residues 4-221
Chain F in data/6BF7.pdb spans residues 2-212
Chain A in data/6BF9.pdb spans residues 46-1011
Chain B in data/6BF9.pdb spans residues 46-1011
Chain C in data/6BF9.pdb spans residues 4-221
Chain D in data/6BF9.pdb spans residues 2-212
Chain E in data/6BF9.pdb spans residues 4-221
Chain F in data/6BF9.pdb spans residues 2-212
Chain B in data/6BFT.pdb spans residues 1-420
Chain A in data/6BFT.pdb spans residues 1-324
Chain L in data/6BFT.pdb spans residues 1-436
Chain H in data/6BFT.pdb spans residues 1-435
Chain G in data/6BFT.pdb spans residues 13-208
Chain C in data/6BFT.pdb spans residues 14-208
Chain C in data/6BGT.pdb spans residues 1-828
Chain A in data/6BGT.pdb spans residues 1-301
Chain B in data/6BGT.pdb spans residues 1-218
Chain H in data/6BKB.pdb spans residues 2-304
Chain E in data/6BKB.pdb spans residues 421-705
Chain L in data/6BKB.pdb spans residues 2-311
Chain A in data/6BKB.pdb spans residues 1-2
Chain H in data/6BKC.pdb spans residues 3-312
Chain L in data/6BKC.pdb spans residues 1-331
Chain E in data/6BKC.pdb spans residues 421-804
Chain A in data/6BZ4.pdb spans residues 237-624
Chain B in data/6BZ4.pdb spans residues 237-634
Chain C in data/6BZ4.pdb spans residues 1-7
Chain D in data/6BZ4.pdb spans residues 1-6
Chain H in data/6DKJ.pdb spans residues 1-541
Chain L in data/6DKJ.pdb spans residues 1-502
Chain C in data/6DKJ.pdb spans residues 29-284
Chain A in data/6DKJ.pdb spans residues 1-619
Chain B in data/6DKJ.pdb spans residues 1-477
Chain D in data/6DKJ.pdb spans residues 29-380
Chain A in data/6EAQ.pdb spans residues 235-746
Chain B in data/6EAQ.pdb spans residues 233-731
Chain C in data/6EAQ.pdb spans residues 5-372
Chain D in data/6EAQ.pdb spans residues 1-7
Chain E in data/6EAQ.pdb spans residues 1-8
Chain F in data/6EAQ.pdb spans residues 1-3
Chain G in data/6EAQ.pdb spans residues 1-3
Chain A in data/6FCZ.pdb spans residues 90-375
Chain B in data/6FCZ.pdb spans residues 92-343
Chain C in data/6FCZ.pdb spans residues 89-376
Chain H in data/6FCZ.pdb spans residues 232-502
Chain K in data/6FCZ.pdb spans residues 235-502
Chain A in data/6FGO.pdb spans residues 237-1176
Chain E in data/6FGO.pdb spans residues 3-218
Chain B in data/6FGO.pdb spans residues 237-1128
Chain F in data/6FGO.pdb spans residues 1-226
Chain C in data/6FGO.pdb spans residues 237-1123
Chain G in data/6FGO.pdb spans residues 4-201
Chain D in data/6FGO.pdb spans residues 237-1105
Chain H in data/6FGO.pdb spans residues 5-602
Chain I in data/6FGO.pdb spans residues 1-8
Chain J in data/6FGO.pdb spans residues 1-7
Chain K in data/6FGO.pdb spans residues 1-8
Chain L in data/6FGO.pdb spans residues 1-5
Chain A in data/6G1E.pdb spans residues 237-1275
Chain B in data/6G1E.pdb spans residues 237-1244
Chain C in data/6G1E.pdb spans residues 1-9
Chain D in data/6G1E.pdb spans residues 1-9
Chain A in data/6IFJ.pdb spans residues 236-771
Chain B in data/6IFJ.pdb spans residues 235-763
Chain C in data/6IFJ.pdb spans residues 1-110
Chain D in data/6IFJ.pdb spans residues 1-106
Chain E in data/6IFJ.pdb spans residues 1-7
Chain F in data/6IFJ.pdb spans residues 1-7
Chain A in data/6IQG.pdb spans residues 237-629
Chain B in data/6IQG.pdb spans residues 237-622
Chain C in data/6IQG.pdb spans residues 603-615
Chain D in data/6IQG.pdb spans residues 602-616
Chain E in data/6IQG.pdb spans residues 1-8
Chain F in data/6IQG.pdb spans residues 1-8
File not found: data/6IQH.pdb
Chain A in data/6KA7.pdb spans residues 7-301
Chain B in data/6KA7.pdb spans residues 17-259
Chain C in data/6KA7.pdb spans residues 238-601
Chain D in data/6KA7.pdb spans residues 238-604
Chain E in data/6KA7.pdb spans residues 1-8
Chain F in data/6KA7.pdb spans residues 1-8
Chain E in data/6MB3.pdb spans residues 103-192
Chain H in data/6MB3.pdb spans residues 2-113
Chain L in data/6MB3.pdb spans residues 2-107
Chain A in data/6MB3.pdb spans residues 2-113
Chain K in data/6MB3.pdb spans residues 2-107
Chain B in data/6MB3.pdb spans residues 2-113
Chain M in data/6MB3.pdb spans residues 2-107
Chain C in data/6MB3.pdb spans residues 2-113
Chain N in data/6MB3.pdb spans residues 2-107
Chain D in data/6MB3.pdb spans residues 2-113
Chain O in data/6MB3.pdb spans residues 2-107
Chain F in data/6MB3.pdb spans residues 2-113
Chain P in data/6MB3.pdb spans residues 2-107
Chain G in data/6MB3.pdb spans residues 2-113
Chain Q in data/6MB3.pdb spans residues 2-107
Chain I in data/6MB3.pdb spans residues 2-113
Chain R in data/6MB3.pdb spans residues 2-107
Chain J in data/6MB3.pdb spans residues 2-113
Chain S in data/6MB3.pdb spans residues 2-107
Chain L in data/6MSY.pdb spans residues 1-472
Chain H in data/6MSY.pdb spans residues 1-489
Chain A in data/6MSY.pdb spans residues 1-4
Chain L in data/6MU3.pdb spans residues 2-331
Chain H in data/6MU3.pdb spans residues 1-417
Chain K in data/6MU3.pdb spans residues 2-318
Chain M in data/6MU3.pdb spans residues 1-416
Chain A in data/6MU3.pdb spans residues 1-4
Chain B in data/6MU3.pdb spans residues 1-4
Chain L in data/6MUB.pdb spans residues 1-302
Chain H in data/6MUB.pdb spans residues 1-402
Chain K in data/6MUB.pdb spans residues 1-213
Chain M in data/6MUB.pdb spans residues 1-402
Chain A in data/6MUB.pdb spans residues 1-5
Chain B in data/6MUB.pdb spans residues 1-5
Chain L in data/6N2X.pdb spans residues 1-213
Chain H in data/6N2X.pdb spans residues 1-228
Chain K in data/6N2X.pdb spans residues 1-213
Chain M in data/6N2X.pdb spans residues 1-228
Chain A in data/6N2X.pdb spans residues 1-11
Chain B in data/6N2X.pdb spans residues 1-11
Chain L in data/6N32.pdb spans residues 2-465
Chain H in data/6N32.pdb spans residues 1-346
Chain M in data/6N32.pdb spans residues 1-499
Chain K in data/6N32.pdb spans residues 1-367
Chain L in data/6N35.pdb spans residues 2-517
Chain H in data/6N35.pdb spans residues 1-739
Chain K in data/6N35.pdb spans residues 1-688
Chain M in data/6N35.pdb spans residues 1-517
Chain A in data/6N35.pdb spans residues 1-2
Chain A in data/6OGE.pdb spans residues 23-709
Chain B in data/6OGE.pdb spans residues 1-214
Chain C in data/6OGE.pdb spans residues 1-216
Chain D in data/6OGE.pdb spans residues 1-214
Chain E in data/6OGE.pdb spans residues 1-220
Chain F in data/6OGE.pdb spans residues 1-4
Chain A in data/6OKQ.pdb spans residues 1-213
Chain B in data/6OKQ.pdb spans residues 2-210
Chain C in data/6OKQ.pdb spans residues 2-210
Chain D in data/6OKQ.pdb spans residues 1-213
Chain E in data/6OKQ.pdb spans residues 2-210
Chain F in data/6OKQ.pdb spans residues 1-213
Chain A in data/6P6D.pdb spans residues 238-755
Chain B in data/6P6D.pdb spans residues 238-740
Chain C in data/6P6D.pdb spans residues 1-8
Chain D in data/6P6D.pdb spans residues 1-8
Chain A in data/6UBI.pdb spans residues 1-539
Chain B in data/6UBI.pdb spans residues 1-506
Chain C in data/6UBI.pdb spans residues 512-614
Chain D in data/6UBI.pdb spans residues 1-487
Chain E in data/6UBI.pdb spans residues 1-515
Chain F in data/6UBI.pdb spans residues 512-602
Chain A in data/6UGW.pdb spans residues 240-752
Chain B in data/6UGW.pdb spans residues 1-8
Chain A in data/6UGX.pdb spans residues 240-1289
Chain B in data/6UGX.pdb spans residues 240-750
Chain C in data/6UGX.pdb spans residues 1-8
Chain D in data/6UGX.pdb spans residues 1-7
Chain A in data/6UGY.pdb spans residues 241-740
Chain B in data/6UGY.pdb spans residues 1-8
Chain H in data/6UOE.pdb spans residues 1-611
Chain L in data/6UOE.pdb spans residues 1-1102
Chain P in data/6UOE.pdb spans residues 68-133
Chain A in data/6V8Z.pdb spans residues 32-609
Chain B in data/6V8Z.pdb spans residues 506-703
Chain D in data/6V8Z.pdb spans residues 1-211
Chain E in data/6V8Z.pdb spans residues 5-212
Chain C in data/6V8Z.pdb spans residues 1-210
Chain F in data/6V8Z.pdb spans residues 1-213
Chain G in data/6V8Z.pdb spans residues 32-609
Chain H in data/6V8Z.pdb spans residues 506-703
Chain J in data/6V8Z.pdb spans residues 1-211
Chain K in data/6V8Z.pdb spans residues 5-212
Chain I in data/6V8Z.pdb spans residues 1-210
Chain L in data/6V8Z.pdb spans residues 1-213
Chain M in data/6V8Z.pdb spans residues 32-609
Chain N in data/6V8Z.pdb spans residues 506-703
Chain P in data/6V8Z.pdb spans residues 1-211
Chain Q in data/6V8Z.pdb spans residues 5-212
Chain O in data/6V8Z.pdb spans residues 1-210
Chain R in data/6V8Z.pdb spans residues 1-213
Chain S in data/6V8Z.pdb spans residues 1-2
Chain T in data/6V8Z.pdb spans residues 1-3
Chain U in data/6V8Z.pdb spans residues 1-2
Chain V in data/6V8Z.pdb spans residues 1-2
Chain W in data/6V8Z.pdb spans residues 1-2
Chain X in data/6V8Z.pdb spans residues 1-2
Chain Y in data/6V8Z.pdb spans residues 1-2
Chain Z in data/6V8Z.pdb spans residues 1-8
Chain a in data/6V8Z.pdb spans residues 1-4
Chain b in data/6V8Z.pdb spans residues 1-2
Chain c in data/6V8Z.pdb spans residues 1-2
Chain d in data/6V8Z.pdb spans residues 1-2
Chain e in data/6V8Z.pdb spans residues 1-2
Chain f in data/6V8Z.pdb spans residues 1-3
Chain g in data/6V8Z.pdb spans residues 1-2
Chain h in data/6V8Z.pdb spans residues 1-2
Chain i in data/6V8Z.pdb spans residues 1-2
Chain j in data/6V8Z.pdb spans residues 1-2
Chain k in data/6V8Z.pdb spans residues 1-2
Chain l in data/6V8Z.pdb spans residues 1-8
Chain m in data/6V8Z.pdb spans residues 1-4
Chain n in data/6V8Z.pdb spans residues 1-2
Chain o in data/6V8Z.pdb spans residues 1-2
Chain p in data/6V8Z.pdb spans residues 1-2
Chain q in data/6V8Z.pdb spans residues 1-2
Chain r in data/6V8Z.pdb spans residues 1-3
Chain s in data/6V8Z.pdb spans residues 1-2
Chain t in data/6V8Z.pdb spans residues 1-2
Chain u in data/6V8Z.pdb spans residues 1-2
Chain v in data/6V8Z.pdb spans residues 1-2
Chain w in data/6V8Z.pdb spans residues 1-2
Chain x in data/6V8Z.pdb spans residues 1-8
Chain y in data/6V8Z.pdb spans residues 1-4
Chain z in data/6V8Z.pdb spans residues 1-2
Chain 0 in data/6V8Z.pdb spans residues 1-2
Chain 1 in data/6V8Z.pdb spans residues 1-2
Chain A in data/6VSL.pdb spans residues 235-719
Chain B in data/6VSL.pdb spans residues 237-670
Chain C in data/6VSL.pdb spans residues 1-9
Chain D in data/6VSL.pdb spans residues 1-9
Chain A in data/6VSZ.pdb spans residues 236-609
Chain B in data/6VSZ.pdb spans residues 236-611
Chain C in data/6VSZ.pdb spans residues 1-8
Chain D in data/6VSZ.pdb spans residues 1-8
Chain A in data/6X3I.pdb spans residues 236-636
Chain B in data/6X3I.pdb spans residues 1-2
Chain A in data/6YSC.pdb spans residues 237-648
Chain B in data/6YSC.pdb spans residues 237-633
Chain C in data/6YSC.pdb spans residues 1-8
Chain D in data/6YSC.pdb spans residues 1-8
Chain A in data/6YT7.pdb spans residues 237-685
Chain B in data/6YT7.pdb spans residues 237-720
Chain C in data/6YT7.pdb spans residues 1-8
Chain D in data/6YT7.pdb spans residues 1-8
Chain A in data/6YTB.pdb spans residues 237-693
Chain B in data/6YTB.pdb spans residues 237-740
Chain C in data/6YTB.pdb spans residues 1-8
Chain D in data/6YTB.pdb spans residues 1-8
Chain A in data/7CZT.pdb spans residues 27-1409
Chain B in data/7CZT.pdb spans residues 27-1411
Chain C in data/7CZT.pdb spans residues 27-1408
Chain H in data/7CZT.pdb spans residues 1-221
Chain L in data/7CZT.pdb spans residues 1-219
Chain I in data/7CZT.pdb spans residues 1-221
Chain M in data/7CZT.pdb spans residues 1-219
Chain D in data/7CZT.pdb spans residues 1-2
Chain E in data/7CZT.pdb spans residues 1-2
Chain F in data/7CZT.pdb spans residues 1-2
Chain G in data/7CZT.pdb spans residues 1-2
Chain J in data/7CZT.pdb spans residues 1-2
Chain K in data/7CZT.pdb spans residues 1-2
Chain N in data/7CZT.pdb spans residues 1-2
Chain O in data/7CZT.pdb spans residues 1-2
Chain P in data/7CZT.pdb spans residues 1-2
Chain Q in data/7CZT.pdb spans residues 1-2
Chain R in data/7CZT.pdb spans residues 1-2
Chain S in data/7CZT.pdb spans residues 1-2
Chain T in data/7CZT.pdb spans residues 1-2
Chain U in data/7CZT.pdb spans residues 1-2
Chain V in data/7CZT.pdb spans residues 1-2
Chain W in data/7CZT.pdb spans residues 1-2
Chain X in data/7CZT.pdb spans residues 1-2
Chain Y in data/7CZT.pdb spans residues 1-2
Chain Z in data/7CZT.pdb spans residues 1-2
Chain a in data/7CZT.pdb spans residues 1-2
Chain b in data/7CZT.pdb spans residues 1-2
Chain A in data/7CZU.pdb spans residues 27-1409
Chain B in data/7CZU.pdb spans residues 27-1411
Chain C in data/7CZU.pdb spans residues 27-1408
Chain J in data/7CZU.pdb spans residues 1-229
Chain N in data/7CZU.pdb spans residues 1-214
Chain H in data/7CZU.pdb spans residues 1-229
Chain K in data/7CZU.pdb spans residues 1-214
Chain D in data/7CZU.pdb spans residues 1-2
Chain E in data/7CZU.pdb spans residues 1-2
Chain F in data/7CZU.pdb spans residues 1-2
Chain G in data/7CZU.pdb spans residues 1-2
Chain I in data/7CZU.pdb spans residues 1-2
Chain L in data/7CZU.pdb spans residues 1-2
Chain M in data/7CZU.pdb spans residues 1-2
Chain O in data/7CZU.pdb spans residues 1-2
Chain P in data/7CZU.pdb spans residues 1-2
Chain Q in data/7CZU.pdb spans residues 1-2
Chain R in data/7CZU.pdb spans residues 1-2
Chain S in data/7CZU.pdb spans residues 1-2
Chain T in data/7CZU.pdb spans residues 1-2
Chain U in data/7CZU.pdb spans residues 1-2
Chain V in data/7CZU.pdb spans residues 1-2
Chain W in data/7CZU.pdb spans residues 1-2
Chain X in data/7CZU.pdb spans residues 1-2
Chain Y in data/7CZU.pdb spans residues 1-2
Chain Z in data/7CZU.pdb spans residues 1-2
Chain a in data/7CZU.pdb spans residues 1-2
Chain b in data/7CZU.pdb spans residues 1-2
Chain A in data/7CZV.pdb spans residues 27-1411
Chain B in data/7CZV.pdb spans residues 27-1410
Chain C in data/7CZV.pdb spans residues 27-1410
Chain H in data/7CZV.pdb spans residues 1-229
Chain K in data/7CZV.pdb spans residues 1-214
Chain I in data/7CZV.pdb spans residues 1-229
Chain M in data/7CZV.pdb spans residues 1-214
Chain J in data/7CZV.pdb spans residues 1-229
Chain N in data/7CZV.pdb spans residues 1-214
Chain D in data/7CZV.pdb spans residues 1-2
Chain E in data/7CZV.pdb spans residues 1-2
Chain F in data/7CZV.pdb spans residues 1-2
Chain G in data/7CZV.pdb spans residues 1-2
Chain L in data/7CZV.pdb spans residues 1-2
Chain O in data/7CZV.pdb spans residues 1-2
Chain P in data/7CZV.pdb spans residues 1-2
Chain Q in data/7CZV.pdb spans residues 1-2
Chain R in data/7CZV.pdb spans residues 1-2
Chain S in data/7CZV.pdb spans residues 1-2
Chain T in data/7CZV.pdb spans residues 1-2
Chain U in data/7CZV.pdb spans residues 1-2
Chain V in data/7CZV.pdb spans residues 1-2
Chain W in data/7CZV.pdb spans residues 1-2
Chain X in data/7CZV.pdb spans residues 1-2
Chain Y in data/7CZV.pdb spans residues 1-2
Chain Z in data/7CZV.pdb spans residues 1-2
Chain a in data/7CZV.pdb spans residues 1-2
Chain b in data/7CZV.pdb spans residues 1-2
Chain c in data/7CZV.pdb spans residues 1-2
Chain H in data/7T17.pdb spans residues 1-105
Chain I in data/7T17.pdb spans residues 1-106
Chain J in data/7T17.pdb spans residues 1-123
Chain K in data/7T17.pdb spans residues 1-110
Chain F in data/7T17.pdb spans residues 1-123
Chain G in data/7T17.pdb spans residues 1-110
Chain A in data/7T17.pdb spans residues 1-392
Chain C in data/7T17.pdb spans residues 1-392
Chain E in data/7T17.pdb spans residues 1-392
Chain A in data/7URU.pdb spans residues 233-505
Chain B in data/7URU.pdb spans residues 232-505
Chain C in data/7URU.pdb spans residues 5-1104
Chain D in data/7URU.pdb spans residues 1-8
Chain E in data/7URU.pdb spans residues 1-8
Chain A in data/7X13.pdb spans residues 236-459
Chain B in data/7X13.pdb spans residues 236-458
Chain C in data/7X13.pdb spans residues 236-458
Chain D in data/7X13.pdb spans residues 236-457
Chain E in data/7X13.pdb spans residues 236-457
Chain F in data/7X13.pdb spans residues 236-457
Chain G in data/7X13.pdb spans residues 236-456
Chain H in data/7X13.pdb spans residues 236-458
Chain I in data/7X13.pdb spans residues 236-458
Chain J in data/7X13.pdb spans residues 236-458
Chain K in data/7X13.pdb spans residues 236-456
Chain L in data/7X13.pdb spans residues 236-456
Chain C in data/7XXL.pdb spans residues 3-211
Chain A in data/7XXL.pdb spans residues 1-229
Chain B in data/7XXL.pdb spans residues 333-601
Chain A in data/8A49.pdb spans residues 238-506
Chain C in data/8A49.pdb spans residues 102-1131
Chain D in data/8A49.pdb spans residues 103-1117
Chain B in data/8A49.pdb spans residues 238-514
Chain E in data/8A49.pdb spans residues 1-8
File not found: data/8A64.pdb
Chain C in data/8DAO.pdb spans residues 1-126
Chain D in data/8DAO.pdb spans residues 1-109
Chain E in data/8DAO.pdb spans residues 114-216
Chain F in data/8DAO.pdb spans residues 108-213
Chain I in data/8DAO.pdb spans residues 812-823
Chain A in data/8DAO.pdb spans residues 1-126
Chain B in data/8DAO.pdb spans residues 1-109
Chain G in data/8DAO.pdb spans residues 114-216
Chain H in data/8DAO.pdb spans residues 108-213
Chain J in data/8DAO.pdb spans residues 812-823
Chain A in data/8DBZ.pdb spans residues 22-368
Chain D in data/8DBZ.pdb spans residues 1-125
Chain E in data/8DBZ.pdb spans residues 1-109
Chain B in data/8DBZ.pdb spans residues 1-127
Chain C in data/8DBZ.pdb spans residues 1-112
Chain F in data/8DBZ.pdb spans residues 1-103
Chain G in data/8DBZ.pdb spans residues 1-107
Chain H in data/8DBZ.pdb spans residues 1-103
Chain I in data/8DBZ.pdb spans residues 1-107
Chain A in data/8DV1.pdb spans residues 330-1301
Chain D in data/8DV1.pdb spans residues 21-1006
Chain A in data/8DV2.pdb spans residues 330-1301
Chain D in data/8DV2.pdb spans residues 21-1106
Chain L in data/8ECQ.pdb spans residues 3-463
Chain H in data/8ECQ.pdb spans residues 1-467
Chain L in data/8ECV.pdb spans residues 3-411
Chain H in data/8ECV.pdb spans residues 1-443
Chain A in data/8ECV.pdb spans residues 3-430
Chain B in data/8ECV.pdb spans residues 1-440
Chain L in data/8ECZ.pdb spans residues 3-212
Chain H in data/8ECZ.pdb spans residues 1-302
Chain A in data/8ECZ.pdb spans residues 3-211
Chain B in data/8ECZ.pdb spans residues 2-303
Chain L in data/8ED1.pdb spans residues 3-404
Chain H in data/8ED1.pdb spans residues 1-403
Chain L in data/8EDF.pdb spans residues 3-210
Chain H in data/8EDF.pdb spans residues 2-267
Chain A in data/8EDF.pdb spans residues 334-601
Chain A in data/8GHR.pdb spans residues 187-1403
Chain B in data/8GHR.pdb spans residues 187-1403
Chain C in data/8GHR.pdb spans residues 2-124
Chain D in data/8GHR.pdb spans residues 2-124
Chain E in data/8GHR.pdb spans residues 1-2
Chain F in data/8GHR.pdb spans residues 1-2
Chain G in data/8GHR.pdb spans residues 1-2
Chain H in data/8GHR.pdb spans residues 1-2
'''


# pip install pymolPy3
# pip install ipymol
## no these are wrong

# brew --prefix pymol     # this shows root
# /opt/homebrew/opt/pymol/libexec/lib/python3.12/site-packages/
# is the address to add to PYTHONPATH

from pymol import cmd, CmdException

# # List of PDB entries without the '.pdb' extension
# pdb_entries = ['1N0X', '3PGF', '4R2G', '5O4E', '5VJ6']  # replace with your list

# # Define the path to the directory where the PDB files are stored
# pdb_directory = "data/"

# # Step 1: Load PDB files and determine the longest structure
# longest_structure = None
# max_residues = 0

# for pdb in pdb_entries:
#     pdb_file = f"{pdb_directory}{pdb}.pdb"

#     cmd.load(pdb_file)

#     # Find the length of the structure by checking the number of residues
#     cmd.select("all_residues", f"{pdb} and polymer")
#     num_residues = cmd.count_atoms("all_residues")

#     if num_residues > max_residues:
#         max_residues = num_residues
#         longest_structure = pdb      # '6V8Z'



from pymol import cmd

pdb_directory = "data/"
longest_structure = None
longest_chain = None
max_continuous_residues = 0

# Function to find the longest continuous chain of residues
def get_longest_continuous_chain(pdb):
    chains = cmd.get_chains(pdb)
    max_chain = None
    max_length = 0

    for chain in chains:
        # Select all polymer residues in this chain
        cmd.select("current_chain", f"{pdb} and chain {chain} and polymer")
        resi_list = cmd.get_model("current_chain").get_residues()
        
        # Find continuous stretches of residues
        current_length = 0
        previous_resi = None
        longest_stretch = 0
        
        for resi in resi_list:
            resi_number = int(resi.resi)
            if previous_resi is None or resi_number == previous_resi + 1:
                current_length += 1
            else:
                current_length = 1
            
            # Update the longest stretch found so far
            if current_length > longest_stretch:
                longest_stretch = current_length
            
            previous_resi = resi_number
        
        # Update max continuous chain if this one is longer
        if longest_stretch > max_length:
            max_length = longest_stretch
            max_chain = chain

    return max_chain, max_length

# Iterate over PDB entries to find the structure with the longest continuous chain
for pdb in pdb_entries:
    pdb_file = f"{pdb_directory}{pdb}.pdb"
    
    # Load the PDB file into PyMOL
    cmd.load(pdb_file)

    # Find the longest continuous chain in this structure
    chain, length = get_longest_continuous_chain(pdb)

    # Check if this structure has the longest continuous chain found so far
    if length > max_continuous_residues:
        max_continuous_residues = length
        longest_structure = pdb
        longest_chain = chain

# Output the result
print(f"Structure with the longest continuous chain: {longest_structure}, Chain: {longest_chain}, Length: {max_continuous_residues} residues")

# Step 2: Align all structures to the longest structure
for pdb in pdb_entries:
    if pdb != longest_structure:
        try:
            cmd.align(pdb, longest_structure)
        except CmdException:
            print(f"Aligning was not possible for {pdb}")
'''
 Selector-Error: Invalid selection name "6IQH".
6IQH<--
Aligning was not possible for 6IQH
 Selector-Error: Invalid selection name "8A64".
8A64<--
Aligning was not possible for 8A64
'''

# Step 3: Color each structure differently
colors = ['red', 'green', 'blue', 'yellow', 'magenta', 'cyan', 'orange', 'purple']
for idx, pdb in enumerate(pdb_entries):
    try:
        cmd.color(colors[idx % len(colors)], pdb)  # Cycle through colors if there are many structures
    except CmdException:
        print(f"Coloring was not possible for {pdb}")
'''
Coloring was not possible for 6IQH
Coloring was not possible for 8A64
'''

# Step 4: Save the session so you can visually inspect it later in PyMOL
cmd.save('aligned_structures.pse')


colors = [
    'green', 'blue', 'yellow', 'magenta', 'cyan', 'orange', 'purple',
    'pink', 'lime', 'teal', 'salmon', 'violet', 'forest', 'deepblue', 'firebrick'
]

for idx, pdb in enumerate(pdb_entries):
    try:
        # Color each PDB structure differently
        cmd.color(colors[idx % len(colors)], pdb)
        
        # Get the longest chain and color it red
        # longest_chain = get_longest_chain(pdb)
        longest_chain = longest_structure
        
        if longest_chain:
            cmd.color('red', f'{pdb} and chain {longest_chain}')
        
        # Make all other chains in this structure transparent
        cmd.set('transparency', 0.6, f'{pdb} and not chain {longest_chain}')
    
    except CmdException:
        print(f"Coloring or transparency adjustment was not possible for {pdb}")
'''
Coloring or transparency adjustment was not possible for 6IQH
Coloring or transparency adjustment was not possible for 8A64
'''
        
# Step 4: Save the session so you can visually inspect it later in PyMOL
cmd.save('aligned_structures_1.pse')

'''
delete all
reinitialize
load ~/projects/python/proteinfolding/aligned_structures_1.pse
'''

'''
hide all
6V8Z is longest
'''













# Define the UniProt accession number (in this case P0DOX5)
accession = "P0DOX5"

# Make a request to UniProt's API to get the PDB ID
url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
response = requests.get(url)



######################################################
# Check if the request was successful

















import os
import requests
from pathlib import Path
from Bio.PDB import PDBParser
import numpy as np
import matplotlib.pyplot as plt

class ProteinAnalyzer:
    def __init__(self, uniprot_accession, out_dir="/Users/josephus/Downloads/pdb_results/"):
        """
        Initialize the ProteinAnalyzer with a UniProt accession number.

        Parameters:
            uniprot_accession (str): The UniProt accession number of the protein.
        """
        self.uniprot_accession = uniprot_accession
        self.pdb_ids = []
        self.out_dir = out_dir
        self.pdb_files = []
        self.structure = None
        self.residue_b_factors = []
        self.residue_ids = []
        self.high_flex_residues = []
        self.threshold = None

        # ensure existence of output directory
        os.makedirs(self.out_dir, exist_ok=True)

    def fetch_pdb_ids(self):
        """
        Fetch PDB IDs associated with the UniProt accession number.
        """
        # Uniprot's API
        response = requests.get(f'https://rest.uniprot.org/uniprotkb/{self.uniprot_accession}.json')

        if response.status_code == 200:
            data = response.json()
            self.pdb_ids = [ref['id'] for ref in data['uniProtKBCrossReferences'] if ref['database'] == 'PDB']
            print(f"PDB IDs associated with UniProt accession {self.uniprot_accession}: {self.pdb_ids}")
        else:
            print(f"Failed to retrieve data for {self.uniprot_accession}")
            self.pdb_ids = []
        
    def download_pdb_files(self):
        """
        Download PDB files for the fetched PDB IDs.
        """
        if not self.pdb_ids:
            print("No PDB IDs to download.")
            return
         
        for pdb_id in self.pdb_ids:
            pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            pdb_response = requests.get(pdb_url)
            
            if pdb_response.status_code == 200:
                # Save the PDB file
                filename = Path(f"{self.out_dir}/{pdb_id}.pdb")
                with open(filename, 'wb') as file:
                    file.write(pdb_response.content)
                self.pdb_files.append(filename)
                print(f"Downloaded PDB file for {pdb_id}")
            else:
                print(f"Failed to download PDB file for {pdb_id}")

    def parse_pdb_file(self, pdb_file):
        """
        Parse a PDB file and store its structure.

        Parameters:
            pdb_file (Path): The path to the PDB file.
        """
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure(pdb_file.stem, pdb_file)

    def analyze_flexibility(self):
        """
        Analyze the parsed PDB structure to identify regions of high flexibility based on B-factors.
        """
        if not self.structure:
            print("No structure to analyze.")
            return
        self.residue_b_factors = []
        self.residue_ids = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    # Skip hetero atoms and water molecules
                    if residue.id[0] != ' ':
                        continue
                    b_factors = [atom.get_bfactor() for atom in residue]
                    avg_b_factor = np.mean(b_factors)
                    self.residue_b_factors.append(avg_b_factor)
                    # Create a unique residue identifier
                    res_id = f"{chain.id}_{residue.id[1]}_{residue.get_resname()}"
                    self.residue_ids.append(res_id)
        # Identify residues with high B-factors
        mean_b = np.mean(self.residue_b_factors)
        std_b = np.std(self.residue_b_factors)
        self.threshold = mean_b + std_b
        self.high_flex_residues = [
            (res_id, b_factor) for res_id, b_factor in zip(self.residue_ids, self.residue_b_factors) if b_factor > self.threshold
        ]
        print("Residues with high flexibility:")
        for res_id, b_factor in self.high_flex_residues:
            print(f"{res_id}: B-factor = {b_factor:.2f}")

   def plot_b_factors(self, save_path=None):
       """
       Plot the average B-factors along the sequence and indicate the high flexibility threshold.
       If save_path is provided, save the plot to the specified file. Otherwise, the plot is shown.
       """
       if not self.residue_b_factors:
           print("No B-factors to plot. Please run analyze_flexibility() first.")
           return
       plt.figure(figsize=(10, 4))
       plt.plot(self.residue_b_factors, label='Average B-factor')
       if self.threshold is not None:
           plt.axhline(y=self.threshold, color='r', linestyle='--', label='High flexibility threshold')
       plt.xlabel('Residue Index')
       plt.ylabel('Average B-factor')
       plt.title(f'B-factors along the sequence of {self.structure.id}')
       plt.legend()
       if save_path:
           plt.savefig(save_path)
           print(f"Plot saved to {save_path}")
           plt.close()
       else:
           plt.show()



# Create an instance of the ProteinAnalyzer class
analyzer = ProteinAnalyzer("P0DOX5")

# Step 1: Fetch PDB IDs associated with the UniProt accession
analyzer.fetch_pdb_ids()

# Step 2: Download the PDB files for the fetched PDB IDs
analyzer.download_pdb_files()

# Step 3: Analyze flexibility for each downloaded PDB file
for pdb_file in analyzer.pdb_files:
    analyzer.parse_pdb_file(pdb_file)
    analyzer.analyze_flexibility()
    # Step 4: Plot the B-factors
    analyzer.plot_b_factors()

