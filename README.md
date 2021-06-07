# Automated network assembly of mechanistic literature to support carcinogen hazard identification with INDRA

## Setup and usage
Running the `INDRA_pipeline.py` script requires  Python 3.6 or above.
To install dependencies, run
```
pip install -r requirements.txt
```
Then edit the script to set a `chemical`, as well as credentials for NDEx
in case the networks are to be uploaded there. Then run the script as
```
python INDRA_pipeline.py
```

## Workflow and results
Given a chemical, `INDRA_pipeline.py` script runs the following operations:
1. Search PubMed with a set of search terms constructed for the chemical,
   and save the PMIDs into a file.
2. For each PMID, use INDRA to retrieve full text content via PubMed Central or
   an abstract via PubMed, depending on availability.
3. Call Reach via INDRA on all retrieved text content to produce "raw" INDRA
   Statements, then save these into pickle files.
4. Run an INDRA assembly pipeline on the "raw" statements and save the
   assembled statements into a JSON file.
5. Using INDRA, generate a CX network from the assembled statements and
   upload it to NDEx.

The intermediate and final results are written into files in the user's
home folder.