import os
import glob
import argparse
import pandas as pd
from Bio import Entrez
from indra import literature
from indra.sources import reach
from indra.tools import assemble_corpus as ac
from indra.statements import statements
from indra.statements import RegulateAmount, RegulateActivity
from indra.assemblers.cx import CxAssembler
from indra.databases import ndex_client


### USER INFORMATION
entrez_email = ''  # Email used at Entrez account
ndex_cred = {'user': '', 'password': ''}  # Ndex credentials

### PUBMED SEARCH INFORMATION:
chemical = ''
date_from = "1900/01/01" 
date_to = "2020/03/01" 
###

# Get the home folder on any platform
HOME = os.path.expanduser('~')
root_path = os.path.join(HOME, 'Desktop')
chemical_folder = os.path.join(root_path, 'pmids', chemical)


date = ' AND (' + date_from + '[PDat] : ' + date_to + '[PDat])'

query_patterns = [
    '[Title] AND ("pharmacokinetics"[MeSH Terms] OR "pharmacokinetics"[Subheading] OR "absorption"[MeSH Terms] OR "distribution"[Title] OR "excretion"[All Fields])',
    '[Title] AND ("Mutation"[MeSH] OR "Cytogenetic Analysis"[MeSH] OR "Mutagens"[MeSH] OR "Oncogenes"[MeSH] OR "Genetic Processes"[MeSH] OR "genomic instability"[MeSH] OR chromosom* OR clastogen* OR "genetic toxicology" OR "strand break" OR "unscheduled DNA synthesis" OR "DNA damage" OR "DNA adducts" OR "SCE" OR "chromatid" OR micronucle* OR mutagen* OR "DNA repair" OR "UDS" OR "DNA fragmentation" OR "DNA cleavage")',
    '[Title] AND ("rna"[MeSH] OR "epigenesis, genetic"[MesH] OR rna OR "rna, messenger"[MeSH] OR "rna" OR "messenger rna" OR mrna OR "histones"[MeSH] OR histones OR epigenetic OR miRNA OR methylation)',
    '[Title] AND ("reactive oxygen species"[MeSH Terms] OR "reactive oxygen species"[All Fields] OR "oxygen radicals"[All Fields] OR "oxidative stress"[MeSH Terms] OR "oxidative"[All Fields] OR "oxidative stress"[All Fields] OR "free radicals"[All Fields])',
    '[Title] AND (inflamm* OR immun* OR chemokine OR cytokine OR leukocyte OR white blood cell)',
    '[Title] AND ("Hormones, Hormone Substitutes, and Hormone Antagonists"[MeSH] OR "Endocrine Disruptors"[MeSH] OR "Thyroid Hormones"[MeSH] OR "Estrogens"[MeSH] OR "Progesterone"[MeSH] OR "Receptors, Estrogen"[MeSH] OR "Receptors, Androgen"[MeSH] OR "Receptors, Progesterone"[MeSH] OR "Receptors, Thyroid Hormone"[MeSH] OR "Receptors, Aryl Hydrocarbon"[MeSH] OR "Peroxisome Proliferator-Activated Receptors"[MeSH] OR "constitutive androstane receptor"[Supplementary Concept] OR "farnesoid X-activated receptor"[Supplementary Concept] OR "liver X receptor"[Supplementary Concept] OR "Retinoid X Receptors"[MeSH])',
    '[Title] AND ("Cell Transformation, Neoplastic"[MeSH] OR "Cell Proliferation"[MeSH] OR apoptosis OR "necrosis"[MeSH] OR "DNA Replication"[MeSH] OR "Cell Cycle"[MeSH] OR brdu OR thymidine OR angiogenesis)',
]
terms = ['%s%s%s' % (chemical, qpattern, date) for qpattern in query_patterns]


#%% CREATION OF FOLDERS
def make_folders(chemical, terms):
    """Create all folders that we will put files into
        - folder 'pmids' at Desktop
        - folder 'chemicalname' at Desktop/pmids
        - folders 'full_pmids_lists' and 'json_files' at Desktop/pmids/chemicalname
        - one folder for each query at Desktop/pmids/chemicalname"""
    queries = []
    for n in range(len(terms)):
        path = os.path.join('pmids', chemical, 'query' + str(n + 1))
        queries.append(path)

    folders = ['pmids',
               os.path.join('pmids', chemical),
               os.path.join('pmids', chemical, 'full_pmids_lists'),
               os.path.join('pmids', chemical, 'json_files')] + queries

    for folder in folders:
        if not os.path.exists(os.path.join(root_path, folder)):
            os.makedirs(os.path.join(root_path, folder))


def search(query):
    """Run a specific PubMed query with BioPython to get PMIDs"""
    Entrez.email = entrez_email
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax='200000',
                            retmode='xml', 
                            term=query)
    results = Entrez.read(handle)
    return results


def search_pmids(terms):
    """Return PMIDs for a list of search terms."""
    number = []
    for i, term in enumerate(terms):
        query = search(term)
        l = [str(q) for q in query["IdList"]]
        number.append(len(l))
        fname = 'pubmed_result_query' + str(i+1) + '.txt'
        path = os.path.join(root_path, 'pmids', chemical, 'full_pmids_lists', fname)
        with open(path, 'w') as fh:
            fh.write('\n'.join(l))
    return number


def get_content(pmid):
    """Retrieve text content and content type for a given PMID."""
    pmid = pmid.strip()
    print(pmid)
    try:
        content, content_type = literature.get_full_text(pmid, 'pmid')
    except AttributeError:
        content = literature.pubmed_client.get_abstract(pmid,
                                                        prepend_title=True)
        content_type = 'abstract'
        print('Attribute Error')
    return content, content_type


def read_content(pmid, content, content_type, query):
    """Process text content for a given PMID and query and dump statements
    into a pickle file."""
    read_offline = True
    rp = None
    literature_stmts = []
    print('Reading %s' % pmid)
    if content_type is not None:
        if content_type == 'abstract':
            print('abstract')
            rp = reach.api.process_text(content, citation=pmid,
                                        offline=read_offline)
        elif content_type == 'pmc_oa_xml':
            print('pmc_oa_xml')
            rp = reach.api.process_nxml_str(content, citation=pmid,
                                            offline=read_offline)
        elif content_type == 'elsevier_xml':
            print('elsevier_xml')
            txt = literature.elsevier_client.extract_text(content)
            if txt:
                print('text')
                rp = reach.api.process_text(txt, citation=pmid,
                                            offline=read_offline)
        if rp is not None:
            print(rp)
            fr = os.path.join(chemical_folder, query, pmid + '.pickle')
            ac.dump_statements(rp.statements, fr, protocol=4)


def load_query_statements(query):
    stmts = []
    # Load statements from processed pmids
    folder = os.path.join(root_path, 'pmids', chemical, query)
    for file in glob.glob(os.path.join(folder, '*.pickle')):
        stmts += ac.load_statements(os.path.join(folder, file),
                                    as_dict=False)
    return stmts


def assemble_statements(stmts):
    """Process statements and run preassembly."""
    stmts = ac.filter_no_hypothesis(stmts)  # Filter out hypothetical statements
    stmts = ac.map_grounding(stmts)  # Map grounding
    stmts = ac.filter_grounded_only(stmts)  # Filter out ungrounded agents
    stmts = ac.map_sequence(stmts)  # Map sequence

    stmts_filtered = [s for s in stmts if
                      (not isinstance(s, (RegulateAmount, RegulateActivity))
                       or s.obj is not None)]
    preassembled_stmts = ac.run_preassembly(stmts_filtered,
                                            return_toplevel=False)  # Run preassembly
    return preassembled_stmts


if __name__ == '__main__':
    # Prepare folders
    make_folders(chemical, terms)

    # GENERATE TABLE
    # Saves an excel table with the search terms and the number
    # of pmids for each query
    number = search_pmids(terms)
    qn = ['query%d' % (n+1) for n in range(len(terms))]
    data = {'Query': qn, 'Search Terms': terms, '# pmids': number}
    df = pd.DataFrame(data)
    fname = 'queries' + chemical + ' PubMed.xlsx'
    df.to_excel(os.path.join(chemical_folder, fname), index=False)

    # GENERATE RAW INDRA STATEMENTS
    for i, query in enumerate(qn):
        print('Processing %s' % query)
        # Retrieve article text from pmids
        paper_contents = {}
        fname = os.path.join(chemical_folder, 'full_pmids_lists',
                             'pubmed_result_' + query + '.txt')
        with open(fname, 'r') as fh:
            for line in fh.readlines():
                pmid = line.strip()
                paper_contents[pmid] = get_content(pmid)

        # Save the text of all retrieved articles in files (text)
        # Generate statements (REACH)
        for pmid, (content, content_type) in paper_contents.items():
            read_content(pmid, content, content_type, query)

    # GENERATE NETWORK
    for i, query in enumerate(qn):
        stmts = load_query_statements(query)
        preassembled_stmts = assemble_statements(stmts)
        # Save preassembled statements in a json file
        fname = os.path.join(root_path, 'pmids', chemical, 'json_files',
                             'preassembled_' + chemical + '_' + query + '.json')
        statements.stmts_to_json_file(preassembled_stmts, fname)

        # Assemble statements and generate network
        cxa = CxAssembler(preassembled_stmts, network_name=chemical + '_' + query)
        cx_str = cxa.make_model()
        network_id = ndex_client.create_network(cx_str, ndex_cred)

        print(query + ': ' + network_id)