# check for specific mention of fMRI-related terms
# to filter out non-fMRI papers


def filter_fMRI_terms(pmids,fmri_terms=['fMRI','functional MRI',
                        'functional magnetic resonance'],usemesh=False):
    """ return pmids that include fMRI-related terms
    in MESH keywords or abstract or title"""

    good_pmids={}
    for pmid in pmids.keys():
        goodkey=0
        if usemesh:
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
