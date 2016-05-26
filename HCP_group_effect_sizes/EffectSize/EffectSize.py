import pandas as pd
import os
import sys
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import csv
import shutil

Contrast = sys.argv[1]
ConDir = sys.argv[2]
OutDir = sys.argv[3]
FslDir = os.environ.get('FSLDIR')

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
HarvardOxfordCortFileAll = os.path.join(FslDir,'data','atlases','HarvardOxford','HarvardOxford-cort-prob-2mm.nii.gz')
HarvardOxfordSubFileAll = os.path.join(FslDir,'data','atlases','HarvardOxford','HarvardOxford-sub-prob-2mm.nii.gz')
HarvardOxford = {
    "C":nib.load(HarvardOxfordCortFileAll).get_data(),
    "S":nib.load(HarvardOxfordSubFileAll).get_data()
}
HarvardOxfordFile = {
    "C":HarvardOxfordCortFileAll,
    "S":HarvardOxfordSubFileAll
}


# Compute Cohens D
D = Cope/VarCope/np.sqrt(186)

# extract coordinates of masks
for msk,no in Masks[Contrast]:

    print("Contrast : ",Contrast)
    print("Masktype : ",msk)
    print("Atlas number : ",no)
    indd = np.where(HarvardOxford[msk][:,:,:,no]>0)
    mask_ind = pd.DataFrame()
    mask_ind['x'] = indd[0]
    mask_ind['y'] = indd[1]
    mask_ind['z'] = indd[2]
    #print(mask_ind)

    # export values
    COHENSD = np.nanmedian(D[mask_ind.x,mask_ind.y,mask_ind.z])
    COHENS10 = np.nanpercentile(D[mask_ind.x,mask_ind.y,mask_ind.z],10)
    COHENS90 = np.nanpercentile(D[mask_ind.x,mask_ind.y,mask_ind.z],90)
    SIZE = len(mask_ind)

    # compute percent bold

    tmpquer = ConDir+"group.gfeat/cope1.feat/featquery"
    featroi_cmd = "fslroi %s %s %s 1" %(HarvardOxfordFile[msk],tmpquer,no)
    os.popen(featroi_cmd).read()
    featquery_cmd = "featquery 1 %s/group.gfeat/cope1.feat 2 stats/pe1 stats/cope1 featquery -p -s -w -b %s" %(ConDir,tmpquer)
    os.popen(featquery_cmd).read()

    featquery_file = ConDir+"group.gfeat/cope1.feat/featquery/report.txt"
    for line in open(featquery_file):
        list = line.split(" ")

    # devide by 2 (except motor) because the contrastmatrices are [1 -1]
    # motor is just an average effect
    if Contrast == "tfMRI_MOTOR":
        PERCENT = float(list[6])
        PERCENT10 = float(list[4])
        PERCENT90 = float(list[7])
    else:
        PERCENT = float(list[6])/2
        PERCENT10 = float(list[4])/2
        PERCENT90 = float(list[7])/2

    shutil.rmtree(ConDir+"group.gfeat/cope1.feat/featquery/")

    with open(os.path.join(OutDir,"es.csv"),'a') as output_file:
        dict_writer = csv.writer(output_file)
        dict_writer.writerow([Contrast,msk,no,SIZE,COHENS10,COHENSD,COHENS90,PERCENT10,PERCENT,PERCENT90])
