"""
Sample size analysis for "Scanning the horizon:..." by Poldrack et al.
Copright Russell Poldrack, 2016

Load the spreadsheet of sample sizes manually annotated by Joe Wexler
based on the neurosynth automated estimates

"""
import os
import numpy,pandas
from get_pmid_data_from_pmid import get_pmid_data_from_pmid
from filter_fMRI_terms import filter_fMRI_terms

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

if os.path.exists('fullset_pmid_records'):
    pmidrecs=pickle.load(open('fullset_pmid_records','rb'))
else:
    pmidrecs=get_pubmed_records(wsdata)
    pmidrecs_fmri=filter_fMRI_terms(pmidrecs)
    pickle.dump(pmidrecs,open('fullset_pmid_records','wb'))

# get single-study and group sizes

def get_study_and_group_sizes(wsdata):
    """ split records into single-group studies ('study') and
    multi-group studies ('group')
    """
    
    study_sizes={}
    group_sizes={}

    for pmid,kind,n_actual in zip(wsdata.PMID,wsdata.kind,wsdata.n_actual):
        if kind=='study':
            # make sure there are no duplicates
            assert pmid not in study_sizes
            study_sizes[pmid]=n_actual
        elif kind=='group':
            if not pmid in study_sizes:
                study_sizes[pmid]=[]
            study_sizes[pmid].append(n_actual)
    return study_sizes,group_sizes

study_sizes,group_sizes=get_study_and_group_sizes(wsdata)
