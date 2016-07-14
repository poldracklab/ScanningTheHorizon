"""
Sample size analysis for "Scanning the horizon:..." by Poldrack et al.
Copright Russell Poldrack, 2016

Load the spreadsheet of sample sizes manually annotated by Joe Wexler
based on the neurosynth automated estimates

"""

import numpy,pandas
from get_pmid_data_from_pmid import get_pmid_data_from_pmid


def load_worksheet(fname='estimated_n_format_real.csv',verbose=True):
    """load data and clean up"""
    d=pandas.read_csv(fname)
    # drop entries with no PMID
    d=d[numpy.isfinite(d.PMID)]
    return d

wsdata=load_worksheet()
print('found %d records and %d unique PMIDS in worksheet'%(wsdata.shape[0],len(set(wsdata.PMID))))
print('found %d single-group studies'%numpy.sum(wsdata.kind=='study'))
def get_pubmed_records(d):
    """ get records for all unique PMIDS in set"""
    pmidrecs={}
    problem_pmids=[]
    for pmid in list(set(d.PMID)):
        try:
            pmidrecs[pmid],r=get_pmid_data_from_pmid(pmid)
        except Exception as e:
            print('problem with',pmid,e)
            problem_pmids.append(pmid)
    return pmidrecs,problem_pmids

pmidrecs=get_pubmed_records(wsdata)
