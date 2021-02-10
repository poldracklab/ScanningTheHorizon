# get all fMRI abstracts from pubmed
# and estimate sample size using Tal's code from NRN analysis

from pathlib import Path
import pandas as pd
from autocv import pubmed
import json
from get_ns_sample_sizes import estimate_n
from Bio import Entrez


def get_pubmed_pmids(query, email, retmax=1000):
    Entrez.email = email
    print(f'using {email} for Entrez service')
    print('searching for', query)
    handle = Entrez.esearch(db="pubmed", retmax=retmax, term=query)
    record = Entrez.read(handle)
    handle.close()
    pmids = [int(i) for i in record['IdList']]
    print('found %d matches' % len(pmids))
    return(pmids)


def divide_list_into_chunks(list_to_split, chunk_size):   
    # looping till length l 
    for i in range(0, len(list_to_split), chunk_size):  
        yield list_to_split[i:(i + chunk_size)] 


def get_pubmed_data(pmids, email, chunk_size=5000):
    Entrez.email = email
    pubmed_records = []
    # load full records in batches of 5000 - since pubmed seems to
    # only return 10000 regardless of larger retmax
    for chunk in divide_list_into_chunks(pmids, chunk_size):
        handle = Entrez.efetch(db="pubmed", id=",".join(['%d' % i for i in chunk]),
                               retmax = chunk_size+1, retmode="xml")
        chunk_records = Entrez.read(handle)['PubmedArticle']
        pubmed_records += chunk_records
        print(len(chunk), len(chunk_records), len(pubmed_records))
    return(pubmed_records)



def get_pubmed_records(basedir='./', force_reload=False,
                       pickle_filename='fmri_pubmed_records.json'):
    basedir = Path(basedir)
    records_pickle_file = basedir / pickle_filename
    if force_reload or not records_pickle_file.exists():
        query = '("fMRI" OR "functional MRI" OR "functional magnetic resonance imaging") AND brain'
        pubmed_ids = get_pubmed_pmids(query, 'poldrack@stanford.edu', 100000)
        pubmed_records = get_pubmed_data(pubmed_ids, 'poldrack@stanford.edu', 5000) 
        with open(records_pickle_file, 'w') as f:
            json.dump(pubmed_records, f)
    # reload from json so that it's a proper dict
    with open(records_pickle_file, 'r') as f:
        pubmed_records = json.load(f)
    return(pubmed_records)


def get_pubmed_year(record):
    year = None
    if 'Year' in record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
        year = int(record['MedlineCitation']['Article'][
            'Journal']['JournalIssue']['PubDate']['Year'])
    elif 'MedlineDate' in record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
        date_parts = record['MedlineCitation']['Article'][
            'Journal']['JournalIssue']['PubDate']['MedlineDate'].split(' ')
        for part in date_parts:
            try:
                year = int(part)
            except ValueError:
                pass
    return(year)


def parse_pubmed_records(pubmed_records):
    pubs = {}
    skipped = []
    no_estimate = []
    for i in pubmed_records:
        pub={}
        pmid = pubmed.get_pubmed_pmid(i)
        pub['year'] = get_pubmed_year(i)
        pub['abstract'] = pubmed.get_pubmed_abstract(i)
        if pub['abstract'] is None:
            skipped.append(pmid)
            continue
        pub['n'] = estimate_n(pub['abstract'])
        if len(pub['n']) == 0:
            no_estimate.append(pmid)
            continue
        pubs[pmid] = pub
    return(pubs)


def get_total_n_data_frame(parsed_records):
    pmid = set(parsed_records.keys())
    year = [parsed_records[key]['year'] for key in pmid]
    group_n = [parsed_records[key]['n'] for key in pmid]
    total_n = [sum([int(i) for i in n]) for n in group_n]
    df = pd.DataFrame({'year': year, 'total_n': total_n, 'group_n': group_n}, index=pmid)
    return(df)

def get_group_n_data_frame(parsed_records):
    pmids = set(parsed_records.keys())
    group_n_df = pd.DataFrame(columns=['pmid', 'year', 'group', 'n'])
    for pmid in pmids:
        year = parsed_records[pmid]['year']
        for idx, group_n in enumerate(parsed_records[pmid]['n']):
            pmid_rec = {
                'year': year,
                'pmid': pmid,
                'n': int(group_n),
                'group': idx
            }
            df_idx = f'{pmid}_{idx}'
            group_n_df.loc[df_idx, :] = pmid_rec
    return(group_n_df)


if __name__ == "__main__":
    pubmed_records = get_pubmed_records()
    parsed_records = parse_pubmed_records(pubmed_records)
    df = get_total_n_data_frame(parsed_records)
    df.to_csv('total_n_data.csv')
    group_df = get_group_n_data_frame(parsed_records)
    group_df.to_csv('group_n_data.csv')