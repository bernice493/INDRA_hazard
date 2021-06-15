# Automated network assembly of mechanistic literature to support carcinogen hazard identification with INDRA

## Setup and usage
Running the `INDRA_pipeline.py` script requires  Python 3.6 or above.
To install dependencies, run
```
pip install -r requirements.txt
```

To set up Reach for reading locally, follow instructions [here](https://indra.readthedocs.io/en/latest/modules/sources/reach/index.html#indra-using-a-reach-jar-through-a-python-java-bridge-indra-sources-reach-reader).


Then run the script as
```
python INDRA_pipeline.py [chemical]
```
As seen above, the script takes the name of a chemical as its first argument
and takes a number of other optional arguments as follows

```
usage: INDRA_pipeline.py [-h] [--date_from DATE_FROM] [--date_to DATE_TO] [--entrez_email ENTREZ_EMAIL] [--ndex_user NDEX_USER] [--ndex_password NDEX_PASSWORD]
                         chemical

positional arguments:
  chemical              The name of the chemical to query for

optional arguments:
  -h, --help            show this help message and exit
  --date_from DATE_FROM
                        The starting date from which to searchfor PMIDs
  --date_to DATE_TO     The end date up to which to searchfor PMIDs
  --entrez_email ENTREZ_EMAIL
                        The email to use when querying PubMed via BioPython.
  --ndex_user NDEX_USER
                        The user name to use if uploading networks to NDEx
  --ndex_password NDEX_PASSWORD
                        The password to use if uploading networks to NDEx
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