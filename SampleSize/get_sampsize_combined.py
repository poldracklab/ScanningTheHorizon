"""
Sample size analysis for "Scanning the horizon:..." by Poldrack et al.
Copright Russell Poldrack, 2016

Load the spreadsheet of sample sizes manually annotated by Joe Wexler
based on the neurosynth automated estimates

"""
import os,pickle
import numpy,pandas
from get_pmid_data_from_pmid import get_pmid_data_from_pmid
from filter_fMRI_terms import filter_fMRI_terms

def load_worksheet(fname='estimated_n_format_real.csv',verbose=True):
    """load data and clean up"""
    d=pandas.read_csv(fname)
    # drop entries with no PMID
    d=d[numpy.isfinite(d.PMID)]
    d=d[d.n_actual!='n/a']

    return d

wsdata=load_worksheet()
print('found %d records and %d unique PMIDS in worksheet'%(wsdata.shape[0],len(set(wsdata.PMID))))
print('found %d single-group studies'%numpy.sum(wsdata.kind=='study'))

def get_pubmed_records(d,pmidrecs={}):
    """ get records for all unique PMIDS in set"""
    problem_pmids=[]
    for pmid in list(set(d.PMID)):
        if not pmid in pmidrecs:
            try:
                pmidrecs[pmid],r=get_pmid_data_from_pmid(pmid,verbose=True)
            except Exception as e:
                print('problem with',pmid,e)
                problem_pmids.append(pmid)
    return pmidrecs,problem_pmids

if os.path.exists('fullset_pmid_records.pkl'):
    pmidrecs=pickle.load(open('fullset_pmid_records.pkl','rb'))
else:
    pmidrecs,problem_pmids=get_pubmed_records(wsdata)
    pmidrecs=filter_fMRI_terms(pmidrecs)
    pickle.dump(pmidrecs,open('fullset_pmid_records.pkl','wb'))

pubyears=[]
for i in range(wsdata.shape[0]):
    try:
        pubyears.append(pmidrecs[wsdata.PMID.iloc[i]]['PubDate'])
    except KeyError:
        pubyears.append('n/a')
wsdata['pubdate']=numpy.array(pubyears)
wsdata.to_csv('fullset_with_pubyears.csv')

# get single-study and group sizes

def get_study_and_group_sizes(wsdata,verbose=False):
    """ split records into single-group studies ('study') and
    multi-group studies ('group')
    """

    study_sizes={}
    group_sizes={}

    for pmid,kind,n_actual in zip(wsdata.PMID,wsdata.kind,wsdata.n_actual):
        if verbose:
            print(pmid,kind,n_actual)
        if kind=='study':
            # make sure there are no duplicates
            assert pmid not in study_sizes
            try:
                study_sizes[pmid]=int(n_actual)
            except ValueError:
                print('skipping',pmid,n_actual)
                continue
        elif kind=='group':
            if not pmid in group_sizes:
                group_sizes[pmid]=[]
            try:
                group_sizes[pmid].append(int(n_actual))
            except ValueError:
                print('skipping',pmid,n_actual)
                continue
    return study_sizes,group_sizes

study_sizes,group_sizes=get_study_and_group_sizes(wsdata)

bad_groups=[]
for k in group_sizes:
    if len(group_sizes[k])<2:
        print('only one n for',k)
        bad_groups.append(k)
for b in list(set(bad_groups)):
    del group_sizes[b]

dates=numpy.arange(2011,2016)

study_data=[]
for pmid in study_sizes:
    if pmid in study_sizes and pmid in pmidrecs:
        if pmidrecs[pmid]['PubDate'] in dates:
            study_data.append([pmid,pmidrecs[pmid]['PubDate'],study_sizes[pmid]])

group_data=[]
for pmid in group_sizes:
    if pmid in group_sizes and pmid in pmidrecs:
        if pmidrecs[pmid]['PubDate'] in dates:
            for n in group_sizes[pmid]:
                group_data.append([pmid,pmidrecs[pmid]['PubDate'],n])

numpy.savetxt('neurosynth_study_data.txt',numpy.array(study_data))
numpy.savetxt('neurosynth_group_data.txt',numpy.array(group_data))
