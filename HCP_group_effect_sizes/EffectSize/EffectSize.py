import pandas as pd
import os
import sys
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import csv

Contrast = sys.argv[1]
ConDir = sys.argv[2]
FslDir = os.environ.get('FSLDIR')
OutDir = sys.argv[3]

Masks={
"tfMRI_MOTOR":[["C",6],
                ["C",25],
                ["S",5],
                ["S",16]],
'tfMRI_WM':[["C",3]],
'tfMRI_LANGUAGE':[["C",6]],
'tfMRI_EMOTION':[["S",9],
                ["S",19]],
'tfMRI_GAMBLING':[["S",10],
                ["S",20]]
}

# read cope + varcope
CopeFile = os.path.join(ConDir,'group.gfeat','cope1.feat','stats','cope1.nii.gz')
Cope = nib.load(CopeFile).get_data()
VarCopeFile = os.path.join(ConDir,'group.gfeat','cope1.feat','stats','varcope1.nii.gz')
VarCope = nib.load(VarCopeFile).get_data()

# read masks
HarvardOxfordCortFile = os.path.join(FslDir,'data','atlases','HarvardOxford','HarvardOxford-cort-maxprob-thr25-2mm.nii.gz')
HarvardOxfordSubFile = os.path.join(FslDir,'data','atlases','HarvardOxford','HarvardOxford-sub-maxprob-thr25-2mm.nii.gz')
HarvardOxford = {
    "C":nib.load(HarvardOxfordCortFile).get_data(),
    "S":nib.load(HarvardOxfordSubFile).get_data()
}

# Compute Cohens D
D = Cope/VarCope/np.sqrt(186)

# extract coordinates of masks
mask_ind = pd.DataFrame()
for msk,no in Masks[Contrast]:
    indices = list(np.where(HarvardOxford[msk]==no+1))
    ind_df = pd.DataFrame(indices).transpose()
    ind_df.columns=['x','y','z']
    mask_ind = mask_ind.append(ind_df)

mask_ind = mask_ind.drop_duplicates()
if Contrast == "tfMRI_LANGUAGE":
    mask_ind = mask_ind[mask_ind.x<45]

# export values
#PERCENTBOLD = np.mean(Cope[mask_ind.x,mask_ind.y,mask_ind.z])/10000
COHENSD = np.nanmedian(D[mask_ind.x,mask_ind.y,mask_ind.z])
COHENS25 = np.nanpercentile(D[mask_ind.x,mask_ind.y,mask_ind.z],25)
COHENS75 = np.nanpercentile(D[mask_ind.x,mask_ind.y,mask_ind.z],75)
SIZE = len(mask_ind)

with open(os.path.join(OutDir,"es.csv"),'a') as output_file:
    dict_writer = csv.writer(output_file)
    dict_writer.writerow([SIZE,COHENS25,COHENSD,COHENS75])
