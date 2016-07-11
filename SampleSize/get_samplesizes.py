"""
analysis of sample sizes from David et al. PLOS one data
and neurosynth automatic estimates

Copyright Russell Poldrack, 2016
"""

import os
import numpy
import xlrd
import requests
from Bio import Entrez
import pickle
import time

# should change this to your own email address
Entrez.email='slacker@harvard.edu'

# first load David et al data from excel spreadsheet provided by Sean David
def load_david_worksheet(fname='fMRI MA significance bias database 03-24-13.xlsx',
    verbose=False):
    workbook = xlrd.open_workbook(fname)
    sheet=workbook.sheet_by_name('Original Sheet')

    studies={}
    all_ids=[]

    # load identifier and N value from spreadsheet cells

    for i in range(1,sheet.nrows):
        idval=sheet.cell(i,5).value

        if idval == xlrd.empty_cell.value:
            continue

        # check for cases with no DOI or PMID
        if str(idval).find('No')==0:
            print('skipping',idval)
            continue

        # determine whether it's a doi or pmid and treat accordingly
        try:
            id=int(idval)
            is_pmid=True
        except ValueError:
            is_pmid=False

        if not is_pmid:
            try:
                id=idval.split('doi:')[1].replace(' ','') #.replace("doi:",'')
            except:
                print('bad identifier',idval)
                continue

        if not id in studies.keys():
            studies[id]=[]

        # get sample size value
        n=sheet.cell(i,6).value
        if n==xlrd.empty_cell.value:
            if verbose:
                print('no value for',id)
            continue

        try:
            nval=int(n)
        except ValueError:
            print('skipping bad value: %s'%n)

        all_ids.append(id)
        studies[id].append(sheet.cell(i,6).value)

    unique_ids=list(studies.keys())
    print('found %d  ids'%len(all_ids))
    print('found %d unique ids'%len(unique_ids))


    # exclude studies with multiple N values
    good_unique_ids=unique_ids
    for id in unique_ids:
        if len(set(studies[id]))>1:
            if verbose:
                print('excluding:',(id,studies[id]))
            good_unique_ids.remove(id)
    return good_unique_ids

david_unique_ids=load_david_worksheet() #verbose=True)


# we need pubmed IDs in order to be able to obtain info about the studies
# and match them to the automated estimates, but not all of the entries
# in the spreadsheet have a PMID associated with them. for those, we
# use the DOI to grab the PMID via the CrossRef API

# first we need to split the ids into dois and pmids
def split_pmid_and_doi(ids):
    pmid=[]
    doi=[]
    for id in ids:
        # if the unique ID is castable as an integer then assume it's a PMID,
        # otherwise assume it's a DOI
        try:
            pmid.append(int(id))
        except:
            doi.append(id)
    return pmid,doi

david_pmids,david_dois=split_pmid_and_doi(david_unique_ids)


# get DOI records using crossref api


def get_doi_records_from_crossref(dois,outfile='doi_records.pkl'):
    doi_records={}
    bad_doi=[]
    for id in dois:
        doi_s=id.replace('doi:','').replace(' ','').split('/')
        assert len(doi_s)>1

        url='http://api.crossref.org/works/%s/%s'%(doi_s[0],'/'.join(doi_s[1:]))
        print(url)
        resp = requests.get(url)
        try:
            doi_records[id]=resp.json()
            print('success:',id)
        except:
            print('Bad resp:',id,resp)
            bad_doi.append(id)

    pickle.dump(doi_records,open(outfile,'wb'))
    return doi_records,bad_doi

if os.path.exists('doi_records.pkl'):
    doi_records=pickle.load(open('doi_records.pkl','rb'))
else:
    doi_records,bad_dois=get_doi_records_from_crossref(dois)


# get pmids for DOIs

def get_pmids_for_dois(doi_records,verbose=False):
    print('Grabbing information from pubmed')
    print('This will take a while because we have to throttle our request rate')
    doi_pmid_cvt={}
    doi_pmids=[]
    bad_cvt=[]
    for id in doi_records.keys():

        time.sleep(0.5) # slow down requests so that we don't get locked out
        # first try searching using the DOI
        handle = Entrez.esearch(db="pubmed", retmax=10, term=id)
        record = Entrez.read(handle)
        handle.close()
        # if DOI search fails, then try searching using title
        if len(record['IdList'])!=1:
            if verbose:
                print('%d matches for doi, trying title'%len(record['IdList']))
            handle = Entrez.esearch(db="pubmed", retmax=10, term=doi_records[id]['message']['title'][0])
            record = Entrez.read(handle)

        if len(record['IdList'])==1:
            doi_pmid_cvt[id]=record['IdList'][0]
            doi_pmids.append(record['IdList'][0])
            if verbose:
                print(record['IdList'])
        else:
            print('no/bad PMID for %s (%d records)'%(id,len(record['IdList'])))
            if verbose:
                print(record['IdList'])
            bad_cvt.append(id)

        if 0:
            if verbose:
                print('problem searching for %s'%id)
            bad_cvt.append(id)
    pickle.dump(doi_pmid_cvt,open( 'doi_pmid_cvt.pkl','wb'))
    pickle.dump(doi_pmids,open( 'doi_pmids.pkl','wb'))
    return doi_pmids,doi_pmid_cvt,bad_cvt

try:
    open('foosdf')
    doi_pmid_cvt=pickle.load(open( 'doi_pmid_cvt.pkl','rb'))
    doi_pmids=pickle.load(open( 'doi_pmids.pkl','rb'))
except FileNotFoundError:
    doi_pmids,doi_pmid_cvt,bad_cvt=get_pmids_for_dois(doi_records,verbose=True)

asdf


# get paper information from PMID using Entrez tools
def get_pmid_data_from_pmid(pmid, delay=0.34,verbose=True):

    pmid=int(pmid)

    handle = Entrez.efetch("pubmed", id="%s"%pmid, retmode="xml")
    time.sleep(delay)
    records=Entrez.parse(handle)
    rec=[i for i in records]
    data={}
    data['Journal']=rec[0]['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
    data['Title']=rec[0]['MedlineCitation']['Article']['ArticleTitle']
    data['Abstract']=rec[0]['MedlineCitation']['Article']['Abstract']['AbstractText']
    data['PubDate']=rec[0]['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
    data['DateCreated']=rec[0]['MedlineCitation']['DateCreated']['Year']
    data['MeshTerms']=[]
    if len(rec[0]['MedlineCitation']['MeshHeadingList'])>0:
        for i in range(len(rec[0]['MedlineCitation']['MeshHeadingList'])):
            data['MeshTerms'].append(str(rec[0]['MedlineCitation']['MeshHeadingList'][i]['DescriptorName']))
    print((pmid,data))
    return data,rec


pmids=david_pmids + doi_pmids

if os.path.exists('pmid_data.pkl'):
    pmid_data=pickle.load(open('pmid_data.pkl','rb'))

else:
    pmid_records={}
    pmid_data={}
    problem_pmid=[]
    for id in pmids:
        d,r=get_pmid_data_from_pmid(id)
        if not d is None:
            pmid_data[id]=d
        else:
            problem_pmid.append(id)
    pickle.dump(pmid_data,open('pmid_data.pkl','wb'))


# check for specific mention of fMRI-related terms
# to filter out non-fMRI papers

def filter_fMRI_terms(pmids,fmri_terms=['fMRI','functional MRI','functional magnetic resonance']):
    good_pmids={}
    for pmid in pmids.keys():
        goodkey=0
        if 'MeshTerms' in pmids[pmid].keys():
            if 'Magnetic Resonance Imaging' in pmids[pmid]['MeshTerms']:
                goodkey=1
        for t in fmri_terms:
            if pmids[pmid]['Abstract'][0].find(t)>-1:
                goodkey=1
            if pmids[pmid]['Title'].find(t)>-1:
                goodkey=1
        if goodkey:
            good_pmids[pmid]=pmids[pmid]
    return good_pmids

print('original David list: %d items'%len(pmid_data))
pmid_data=filter_fMRI_terms(pmid_data)
print('filtered David list: %d items'%len(pmid_data))

 # get N for each

pmid_doi_cvt={}
for k in doi_pmid_cvt.keys():
    pmid_doi_cvt[int(doi_pmid_cvt[k])]=k

sampsize={}
pubyear={}
for id in pmid_data.keys():
    iid=int(id)
    if not iid in pmid_data.keys():
        pmid_data[iid]=pmid_data[id]
        del pmid_data[id]

    # the date can be stored in one of two fields in the pubmed record
    if 'PubDate' in pmid_data[iid]:
        try:
            pubyear[iid]=int(pmid_data[iid]['PubDate']['Year'])
        except:
            pubyear[iid]=int(pmid_data[iid]['PubDate']['MedlineDate'].split(' ')[0])
    if iid in pmid_doi_cvt:
        k=pmid_doi_cvt[iid]
    else:
        k=iid
    if k in studies:
        try:
            sampsize[iid]=int(studies[k][0])
        except:
            print(('bad sample size:',studies[k][0]))

    else:
        print(k)

sampsizes=[]
pubyears=[]
alldata=[]
for k in pubyear.keys():
    try:
        if k in sampsize:
            pubyears.append(pubyear[k])
            sampsizes.append(sampsize[k])
            alldata.append([k,pubyear[k],sampsize[k]])
    except:
        pass
numpy.savetxt('david_sampsizedata.txt',alldata)



# load neurosynth data
f=open('estimated_n.txt')
header=f.readline()
lines=[i.strip().split('\t') for i in f.readlines()]

# get all N estimates for each pmid
# previously we had just taken the first but that was problematic

try:
    ns_n_estimates
except NameError:
    ns_n_estimates={}
    for l in lines:
        print(l[:2])
        pmid=int(l[0])
        if not pmid in ns_n_estimates.keys():
            ns_n_estimates[pmid]=[]
        ns_n_estimates[pmid].append(int(l[1]))


try:
    tal_data=pickle.load(open('tal_data.pkl','rb'))
except FileNotFoundError:
    tal_data={}
    tal_problem_pmid=[]
    for l in lines:
        pmid=int(l[0])
        # only include abstracts with a single N estimate
        # there does not appear to be a good strategy to accurately
        # assess sample size when there are multiple N's reported

        if len(ns_n_estimates[pmid])>1:
            continue
        if not pmid in tal_data:
            tal_data[pmid],_=get_pmid_data_from_pmid(pmid)
        print(tal_data[pmid],l[1:])

        try:
            tal_data[pmid]['nest']=ns_n_estimates[pmid]
        except:
            del tal_data[pmid]
            continue

        tal_data[pmid]['pubdate']=None
        try:
            tal_data[pmid]['pubdate']=int(tal_data[pmid]['PubDate']['Year'])
        except:
            tal_data[pmid]['pubdate']=int(tal_data[pmid]['PubDate']['MedlineDate'].split(' ')[0])

    pickle.dump(tal_data,open('tal_data.pkl','wb'))

match_pmid=[]
tal_est=[]
david_est=[]
for l in lines:
    pmid=int(l[0])

    if pmid in pmid_data and pmid in tal_data:
        match_pmid.append(pmid)
        david_est.append(sampsize[pmid])
        tal_est.append(tal_data[pmid]['nest'])

print('original neurosynth list: %d items'%len(tal_data))
tal_data=filter_fMRI_terms(tal_data)
print('filtered neurosynth list: %d items'%len(tal_data))

tal_data_mtx=[]
for pmid in tal_data.keys():
    if 'nest' in tal_data[pmid]:
        tal_data_mtx.append([int(tal_data[pmid]['pubdate']),tal_data[pmid]['nest'][0]])
tal_data_mtx=numpy.array(tal_data_mtx)

numpy.savetxt('neurosynth_sampsizedata.txt',tal_data_mtx)
