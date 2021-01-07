### PUBMED SEARCH INFORMATION:
chemical = ''
date_from = "1900/01/01" 
date_to = "2020/03/01" 
###

### USER INFORMATION:
user = '' # Windows user name
entrez_email = '' # Email used at Entrez account
ndex_cred = {'user': '', 'password': ''} # Ndex credentials

#%% CREATION OF FOLDERS
""" If they still don't exist, it creates folders where all the obtained files will be saved:
    - folder 'pmids' at Desktop
    - folder 'chemicalname' at Desktop\\pmids
    - folders 'full_pmids_lists' and 'json_files' at Desktop\\pmids\\chemicalname
    - one folder for each query at Desktop\\pmids\\chemicalname"""
import os

queries = []
for n in list(range(1,8)):
    queries.append('pmids\\' + chemical + '\\query' + str(n))
root_path = 'C:\\Users\\' + user + '\\Desktop\\'
folders = ['pmids', ('pmids\\' + chemical), ('pmids\\' + chemical + '\\full_pmids_lists'), ('pmids\\' + chemical + '\\json_files')] + queries

for folder in folders: 
    if not os.path.exists(os.path.join(root_path, folder)):
        os.makedirs(os.path.join(root_path, folder))
        
chemical_folder = root_path + 'pmids\\' + chemical + '\\'

#%% PMIDs LISTS RETRIVAL 
"""Pubmed search for each query"""
from Bio import Entrez

def search(query):
    Entrez.email = entrez_email
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax='200000',
                            retmode='xml', 
                            term=query)
    results = Entrez.read(handle)
    return results

date = ' AND (' + date_from + '[PDat] : ' + date_to + '[PDat])'

terms = []
terms.append(chemical + '[Title] AND ("pharmacokinetics"[MeSH Terms] OR "pharmacokinetics"[Subheading] OR "absorption"[MeSH Terms] OR "distribution"[Title] OR "excretion"[All Fields])' + date)
terms.append(chemical + '[Title] AND ("Mutation"[MeSH] OR "Cytogenetic Analysis"[MeSH] OR "Mutagens"[MeSH] OR "Oncogenes"[MeSH] OR "Genetic Processes"[MeSH] OR "genomic instability"[MeSH] OR chromosom* OR clastogen* OR "genetic toxicology" OR "strand break" OR "unscheduled DNA synthesis" OR "DNA damage" OR "DNA adducts" OR "SCE" OR "chromatid" OR micronucle* OR mutagen* OR "DNA repair" OR "UDS" OR "DNA fragmentation" OR "DNA cleavage")' + date)
terms.append(chemical + '[Title] AND ("rna"[MeSH] OR "epigenesis, genetic"[MesH] OR rna OR "rna, messenger"[MeSH] OR "rna" OR "messenger rna" OR mrna OR "histones"[MeSH] OR histones OR epigenetic OR miRNA OR methylation)' + date)
terms.append(chemical + '[Title] AND ("reactive oxygen species"[MeSH Terms] OR "reactive oxygen species"[All Fields] OR "oxygen radicals"[All Fields] OR "oxidative stress"[MeSH Terms] OR "oxidative"[All Fields] OR "oxidative stress"[All Fields] OR "free radicals"[All Fields])' + date)
terms.append(chemical + '[Title] AND (inflamm* OR immun* OR chemokine OR cytokine OR leukocyte OR white blood cell)' + date)
terms.append(chemical + '[Title] AND ("Hormones, Hormone Substitutes, and Hormone Antagonists"[MeSH] OR "Endocrine Disruptors"[MeSH] OR "Thyroid Hormones"[MeSH] OR "Estrogens"[MeSH] OR "Progesterone"[MeSH] OR "Receptors, Estrogen"[MeSH] OR "Receptors, Androgen"[MeSH] OR "Receptors, Progesterone"[MeSH] OR "Receptors, Thyroid Hormone"[MeSH] OR "Receptors, Aryl Hydrocarbon"[MeSH] OR "Peroxisome Proliferator-Activated Receptors"[MeSH] OR "constitutive androstane receptor"[Supplementary Concept] OR "farnesoid X-activated receptor"[Supplementary Concept] OR "liver X receptor"[Supplementary Concept] OR "Retinoid X Receptors"[MeSH])' + date)
terms.append(chemical + '[Title] AND ("Cell Transformation, Neoplastic"[MeSH] OR "Cell Proliferation"[MeSH] OR apoptosis OR "necrosis"[MeSH] OR "DNA Replication"[MeSH] OR "Cell Cycle"[MeSH] OR brdu OR thymidine OR angiogenesis)' + date)

number = []
for i in range(0,7):
    query = search(terms[i])
    l = []
    for q in query["IdList"]:
        l.append(str(q))
    number.append(len(l))
    with open(root_path + 'pmids\\' + chemical + '\\full_pmids_lists\\pubmed_result_query' + str(i+1) + '.txt', "w") as f:
        for x in l:
            f.write(x + "\n")
            
# GENERATE TABLE
""""Saves an excel table with the search terms and the number of pmids for each query"""
import pandas as pd

qn = []
for n in list(range(1,8)):
    qn.append('query' + str(n))

data = {'Query': qn, 'Search Terms': terms, '# pmids': number}

df = pd.DataFrame(data)
df.to_excel(chemical_folder + "queries " + chemical + " PubMed.xlsx", index = False)

#%% GENERATION OF RAW INDRA STATEMENTS
from indra import literature
from indra.sources import reach
from indra.tools import assemble_corpus as ac

for i in range(0,7):
    query = qn[i]
    print('Processing %s' % query)
    # Retrieve article text from pmids
    paper_contents = {} 
    with open(chemical_folder + 'full_pmids_lists\\pubmed_result_' + query + '.txt') as art_pmid:
            for pmid in art_pmid:
                print(pmid.rstrip())
                try:
                    content, content_type = literature.get_full_text(pmid.rstrip(), 'pmid')
                    paper_contents[pmid] = (content, content_type)
                except AttributeError:
                    content = literature.pubmed_client.get_abstract(pmid.rstrip(), prepend_title=True)
                    content_type = 'abstract'
                    paper_contents[pmid] = (content, content_type)
                    print('Attribute Error')
                    
    # Save the text of all retrieved articles in files (text)
    # Generate statements (REACH)
    read_offline = True
    literature_stmts = []
    
    for pmid, (content, content_type) in paper_contents.items():
        rp = None
        literature_stmts = []
        print('Reading %s' % pmid)
        if content_type != None:
            if content_type == 'abstract':
                print('abstract')
                rp = reach.api.process_text(content, citation=pmid, offline=read_offline)
            elif content_type == 'pmc_oa_xml':
                print('pmc_oa_xml')
                rp = reach.api.process_nxml_str(content, citation=pmid, offline=read_offline) 
            elif content_type == 'elsevier_xml':
                print('elsevier_xml')
                txt = literature.elsevier_client.extract_text(content) 
                if txt:
                    print('text')
                    rp = reach.api.process_text(txt, citation=pmid, offline=read_offline)
            if rp is not None:
                print((rp))
                literature_stmts = rp.statements
                fr = (chemical_folder + query + '\\' + str(pmid.rstrip()) + '.pickle')
                ac.dump_statements(literature_stmts, fr, protocol=4)
                
#%% GENERATE NETWORK
from itertools import chain

from indra.statements import statements
from indra.statements import RegulateAmount, RegulateActivity
from indra.assemblers.cx import CxAssembler 
from indra.databases import ndex_client


for i in list(range(0,7)):
    query = qn[i]
    stmts = []
    # Load statements from processed pmids
    folder = 'C:\\Users\\' + user + '\\Desktop\\pmids\\' + chemical + '\\' + query + '\\'
    for file in os.listdir(folder):
        if file.endswith(".pickle"):
            fr = open(os.path.join(folder, file), 'rb')
            stmts.append(ac.load_statements(os.path.join(folder, file), as_dict=False))
    stmts = list(chain.from_iterable(stmts))


    # Process statements and preassembly  
    stmts = ac.filter_no_hypothesis(stmts)  # Optional: filter out hypothetical statements
    stmts = ac.map_grounding(stmts)         # Map grounding
    stmts = ac.filter_grounded_only(stmts)  # Optional: filter out ungrounded agents
    stmts = ac.map_sequence(stmts)          # Map sequence
    
    stmts_filtered = [s for s in stmts if (not isinstance(s, (RegulateAmount, RegulateActivity)) or s.obj is not None)]
    preassembled_stmts = ac.run_preassembly(stmts_filtered, return_toplevel=False) # Run preassembly 
    
    #Save preassembled statements in a json file
    f = statements.stmts_to_json_file(preassembled_stmts, 'C:\\Users\\' + user + '\\Desktop\\pmids\\' + chemical + '\\json_files\\preassembled_' + chemical + '_' + query + '.json')
    
    # Assemble statements and generate network
    cxa = CxAssembler(preassembled_stmts, network_name = chemical + '_' + query)
    cx_str = cxa.make_model()
    network_id = ndex_client.create_network(cx_str, ndex_cred) 
    
    print(query + ': ' + network_id)
