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
    """load the David et al. data from excel workbook
    and return ID and sample size"""

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
        if len(studies[id])>1:
            if verbose:
                print('excluding:',(id,studies[id]))
            good_unique_ids.remove(id)
            del studies[id]
    return good_unique_ids,studies

david_unique_ids,studies=load_david_worksheet() #verbose=True)


# we need pubmed IDs in order to be able to obtain info about the studies
# and match them to the automated estimates, but not all of the entries
# in the spreadsheet have a PMID associated with them. for those, we
# use the DOI to grab the PMID via the CrossRef API

# first we need to split the ids into dois and pmids
def split_pmid_and_doi(ids):
    """split the list into integers (which are assumed to be PMIDS)
    and strings (which are assumed to be DOIs)"""
    pmid=[]
    doi=[]
    for id in ids:
        # if the unique ID is castable as an integer then assume it's a PMID,
        # otherwise assume it's a DOI
        try:
            pmid.append(int(id))
        except ValueError:
            doi.append(id)
    return pmid,doi

david_pmids,david_dois=split_pmid_and_doi(david_unique_ids)


# get DOI records using crossref api


def get_doi_records_from_crossref(dois,outfile='doi_records.pkl'):
    """get information about a DOI"""
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
        except Exception as e:
            # TBD: fix this exception to be more specific
            print(e)
            print('Bad resp:',id,resp)
            bad_doi.append(id)

    pickle.dump(doi_records,open(outfile,'wb'))
    return doi_records,bad_doi

if os.path.exists('doi_records.pkl'):
    doi_records=pickle.load(open('doi_records.pkl','rb'))
else:
    doi_records,bad_dois=get_doi_records_from_crossref(david_dois)


# get pmid for each DOI

def get_pmids_for_dois(doi_records,verbose=False):
    """ use entrez API to identify PMID from the DOI"""
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
    doi_pmid_cvt=pickle.load(open( 'doi_pmid_cvt.pkl','rb'))
    doi_pmids=pickle.load(open( 'doi_pmids.pkl','rb'))
except FileNotFoundError:
    doi_pmids,doi_pmid_cvt,bad_cvt=get_pmids_for_dois(doi_records,verbose=True)



# get paper information from PMID using Entrez tools
def get_pmid_data_from_pmid(pmid, delay=0.34,verbose=False):
    """ use entrez API tool to get full record for a PMID"""

    pmid=int(pmid)
    handle = Entrez.efetch("pubmed", id="%s"%pmid, retmode="xml")
    time.sleep(delay)
    records=Entrez.parse(handle)
    rec=[i for i in records]
    if len(rec)>1:
        print('Warning: record of length >1',pmid,rec)
    if verbose:
        print(rec)
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
    if verbose:
        print(pmid,data)
    return data,rec


pmids=david_pmids + doi_pmids
def get_all_pmid_data(pmids):
    """loop through a set of PMIDS and grab the record data for each"""

    pmid_data={}
    problem_pmid=[]
    for id in pmids:
        try:
            d,r=get_pmid_data_from_pmid(id)
        except KeyError:
            problem_pmid.append(id)
        else:
            pmid_data[id]=d
    pickle.dump(pmid_data,open('david_pmid_data.pkl','wb'))
    return pmid_data,problem_pmid

if os.path.exists('david_pmid_data.pkl'):
    pmid_data=pickle.load(open('david_pmid_data.pkl','rb'))
else:
    pmid_data,problem_pmid=get_all_pmid_data(pmids)


# check for specific mention of fMRI-related terms
# to filter out non-fMRI papers

def filter_fMRI_terms(pmids,fmri_terms=['fMRI','functional MRI',
                        'functional magnetic resonance']):
    """ return pmids that include fMRI-related terms
    in MESH keywords or abstract or title"""

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

def get_pubdate(pmidrec):
        try:
            pubdate=int(pmidrec['PubDate']['Year'])
        except KeyError:
            pubdate=int(pmidrec['PubDate']['MedlineDate'].split(' ')[0])
        return pubdate

 # get N for each
def get_n_and_date_for_david_pmids(doi_pmid_cvt,pmid_data):
    """get sample size and date for each PMID"""

    # create the reverse dictionary
    pmid_doi_cvt={}
    for k in doi_pmid_cvt.keys():
        pmid_doi_cvt[int(doi_pmid_cvt[k])]=k

    sampsize={}
    pubyear={}
    for id in pmid_data.keys():
        # convert string ids to integer
        iid=int(id)
        if not iid in pmid_data.keys():
            pmid_data[iid]=pmid_data[id]
            del pmid_data[id]

        # the date can be stored in one of two fields in the pubmed record
        if 'PubDate' in pmid_data[iid]:
            pubyear[iid]=get_pubdate(pmid_data[iid])

        # WTF?
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
            print(k,'not found in studies')
    return pubyear,sampsize

pubyear,sampsize=get_n_and_date_for_david_pmids(doi_pmid_cvt,pmid_data)

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

def load_neurosynth_data(infile='estimated_n.txt'):
    f=open(infile)
    header=f.readline()
    lines=[i.strip().split('\t') for i in f.readlines()]
    return lines

lines=load_neurosynth_data()

# get all N estimates for each pmid
# previously we had just taken the first but that was problematic


def get_neurosynth_n_estimates(lines,verbose=False):
    """extract pmid and estimated n from lines of text file"""
    ns_n_estimates={}
    for l in lines:
        if verbose:
            print(l[:2])
        pmid=int(l[0])
        if not pmid in ns_n_estimates.keys():
            ns_n_estimates[pmid]=[]
        ns_n_estimates[pmid].append(int(l[1]))
    return ns_n_estimates





def get_neurosynth_data_from_pmids(lines,verbose=False):
    ns_n_estimates=get_neurosynth_n_estimates(lines)

    tal_data={}
    tal_problem_pmid=[]
    for l in lines:
        pmid=int(l[0])
        # only include abstracts with a single N estimate
        # there does not appear to be a good strategy to accurately
        # assess sample size when there are multiple N's reported

        if len(ns_n_estimates[pmid])>1:
            continue
        try:
            assert pmid in ns_n_estimates
        except AssertionError:
            if verbose:
                print('missing n estimate:',pmid)
            continue

        if not pmid in tal_data:
            try:
                tal_data[pmid],_=get_pmid_data_from_pmid(pmid)
            except KeyError:
                continue
        else:
            raise Exception('duplicate pmids - check this out')

        tal_data[pmid]['nest']=ns_n_estimates[pmid]

        if verbose:
            print(tal_data[pmid],l[1:])

        tal_data[pmid]['pubdate']=get_pubdate(tal_data[pmid])

    pickle.dump(tal_data,open('tal_data.pkl','wb'))
    return tal_data

if os.path.exists('tal_data.pkl'):
    tal_data=pickle.load(open('tal_data.pkl','rb'))
else:
    tal_data=get_neurosynth_data_from_pmids(lines)


print('original neurosynth list: %d items'%len(tal_data))
tal_data=filter_fMRI_terms(tal_data)
print('filtered neurosynth list: %d items'%len(tal_data))

tal_data_mtx=[]
for pmid in tal_data.keys():
    if 'nest' in tal_data[pmid]:
        tal_data_mtx.append([int(tal_data[pmid]['pubdate']),tal_data[pmid]['nest'][0]])
tal_data_mtx=numpy.array(tal_data_mtx)

numpy.savetxt('neurosynth_sampsizedata.txt',tal_data_mtx)
