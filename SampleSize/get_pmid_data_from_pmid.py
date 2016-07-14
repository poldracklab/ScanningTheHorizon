"""
get paper information from PMID using Entrez tools
"""

import xlrd
import time
from Bio import Entrez

def get_pubdate(pubdaterec):
        try:
            pubdate=int(pubdaterec['Year'])
        except KeyError:
            pubdate=int(pubdaterec['MedlineDate'].split(' ')[0])
        return pubdate

def get_pmid_data_from_pmid(pmid, delay=0.34,verbose=False,email='slacker@harvard.edu'):
    """ use entrez API tool to get full record for a PMID"""
    Entrez.email=email
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
    pubdaterec=rec[0]['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
    data['PubDate']=get_pubdate(pubdaterec)
    data['DateCreated']=rec[0]['MedlineCitation']['DateCreated']['Year']
    data['MeshTerms']=[]
    if 'MeshHeadingList' in rec[0]['MedlineCitation']:
        if len(rec[0]['MedlineCitation']['MeshHeadingList'])>0:
            for i in range(len(rec[0]['MedlineCitation']['MeshHeadingList'])):
                data['MeshTerms'].append(str(rec[0]['MedlineCitation']['MeshHeadingList'][i]['DescriptorName']))
                

    if verbose:
        print(pmid,data)
    return data,rec
