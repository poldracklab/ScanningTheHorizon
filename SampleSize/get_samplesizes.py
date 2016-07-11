"""
analysis of sample sizes from David et al. PLOS one data
and neurosynth automatic estimates

Copyright Russell Poldrack, 2016
"""

import numpy
import xlrd
import requests
from Bio import Entrez
import pickle
import time

Entrez.email='poldrack@stanford.edu'

# first load David et al data from excel spreadsheet provided by Sean David

workbook = xlrd.open_workbook('fMRI MA significance bias database 03-24-13.xlsx')
sheet=workbook.sheet_by_name('Original Sheet')


id=[]
studies={}
ids=[]

for i in range(1,sheet.nrows):
    try:
        id=sheet.cell(i,5).value.replace("doi:",'').replace(' ','')
    except:
        id=int(sheet.cell(i,5).value)

    if id == xlrd.empty_cell.value:
        continue

    n=sheet.cell(i,6)
    try:
        int(n.value)
    except:
        print('bad value: %s'%n.value)

    if not n.value == xlrd.empty_cell.value:
        if not id in studies:
            studies[id]=[]
        ids.append(id)

        studies[id].append(sheet.cell(i,6).value)

unique_ids=list(set(ids))
print('found %d  ids'%len(ids))
print('found %d unique ids'%len(unique_ids))

good_unique_ids=list(set(ids))

for id in unique_ids:
    if len(set(studies[id]))>1:
        print((id,studies[id]))
        good_unique_ids.remove(id)

# we need pubmed IDs in order to be able to obtain info about the studies
# and match them to the automated estimates, but not all of the entries
# in the spreadsheet have a PMID associated with them. for those, we
# use the DOI to grab the PMID

has_pmid=[]
has_doi=[]
for id in good_unique_ids:
    try:
        has_pmid.append(int(id.value))
    except:
        has_doi.append(id)

# get DOI records using crossref api

import pickle

try:
    doi_records=pickle.load(open('doi_records.pkl','rb'))
except:
    doi_records={}
    bad_doi=[]
    for id in has_doi:
        try:
            doi_s=id.replace('doi:','').replace(' ','').split('/')
            assert len(doi_s)>1
        except:
            print('found PMID',id)
            has_pmid.append(id)
            continue
        #time.sleep(0.5)
        url='http://api.crossref.org/works/%s/%s'%(doi_s[0],'/'.join(doi_s[1:]))
        print(url)
        resp = requests.get(url)
        try:
            doi_records[id]=resp.json()
            print('success:',id)
        except:
            print('Bad resp:',id,resp)
            bad_doi.append(id)

    pickle.dump(doi_records,open('doi_records.pkl','wb'))



# get pmids for DOIs

try:
    doi_pmid_cvt=pickle.load(open( 'doi_pmid_cvt.pkl','rb'))
    doi_pmids=pickle.load(open( 'doi_pmids.pkl','rb'))
except:
    print('Grabbing information from pubmed')
    print('This will take a while because we have to throttle our request rate')
    doi_pmid_cvt={}
    doi_pmids=[]
    bad_cvt=[]
    for id in doi_records.keys():
        try:
            time.sleep(0.5) # slow down requests so that we don't get locked out
            handle = Entrez.esearch(db="pubmed", retmax=10, term=id)
            record = Entrez.read(handle)
            handle.close()
            if len(record['IdList'])<1:
                print('no match, trying title')
                handle = Entrez.esearch(db="pubmed", retmax=10, term=doi_records[id]['message']['title'][0])
                record = Entrez.read(handle)
            if len(record['IdList'])==1:
                doi_pmid_cvt[id]=record['IdList'][0]
                doi_pmids.append(record['IdList'][0])
                print(record['IdList'])
            else:
                print('no/bad PMID for %s'%id)
                print(record['IdList'])

        except:
            print('problem searching for %s'%id)
    pickle.dump(doi_pmid_cvt,open( 'doi_pmid_cvt.pkl','wb'))
    pickle.dump(doi_pmids,open( 'doi_pmids.pkl','wb'))

# get paper information from PMID using Entrez tools
def get_pmid_data_from_pmid(pmid, delay=0.34,verbose=True):
    pmid=int(pmid)
    try:
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

    except:
        print (pmid,'None')
        return None,None



has_pmid=has_pmid + doi_pmids
try:
    pmid_data=pickle.load(open('pmid_data.pkl','rb'))
except:

    pmid_records={}
    pmid_data={}
    problem_pmid=[]
    for id in has_pmid:
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
except:
    ns_n_estimates={}
    for l in lines:
        print(l[:2])
        pmid=int(l[0])
        if not pmid in ns_n_estimates.keys():
            ns_n_estimates[pmid]=[]
        ns_n_estimates[pmid].append(int(l[1]))


try:
    tal_data=pickle.load(open('tal_data.pkl','rb'))
except:
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
