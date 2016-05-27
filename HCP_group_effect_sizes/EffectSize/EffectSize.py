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
    "MOTOR":[["Precentral Gyrus",1,os.path.join(OutDir,"CreateMasks","masks","MOTOR_1_CJ.nii.gz")],
            ["Supplementary motor cortex",2,os.path.join(OutDir,"CreateMasks","masks","MOTOR_2_CJ.nii.gz")],
            ["Left putamen",3,os.path.join(OutDir,"CreateMasks","masks","MOTOR_3_CJ.nii.gz")],
            ["Right putamen",4,os.path.join(OutDir,"CreateMasks","masks","MOTOR_4_CJ.nii.gz")]],
    'WM':[["Middle frontal gyrus",1,os.path.join(OutDir,"CreateMasks","masks","WM_1_CJ.nii.gz")]],
    'EMOTION':[["Left amygdala",1,os.path.join(OutDir,"CreateMasks","masks","EMOTION_1_CJ.nii.gz")],
            ["Right amygdala",2,os.path.join(OutDir,"CreateMasks","masks","EMOTION_2_CJ.nii.gz")]],
    'GAMBLING':[["Left accumbens",1,os.path.join(OutDir,"CreateMasks","masks","GAMBLING_1_CJ.nii.gz")],
            ["Right accumbens",2,os.path.join(OutDir,"CreateMasks","masks","GAMBLING_2_CJ.nii.gz")]]
}


# read cope + varcope
CopeFile = os.path.join(ConDir,'group.gfeat','cope1.feat','stats','cope1.nii.gz')
Cope = nib.load(CopeFile).get_data()
VarCopeFile = os.path.join(ConDir,'group.gfeat','cope1.feat','stats','varcope1.nii.gz')
VarCope = nib.load(VarCopeFile).get_data()

# Compute Cohens D
D = Cope/VarCope/np.sqrt(186)

# extract coordinates of masks
for a in range(len(Masks[Contrast])):

    mask = Masks[Contrast][a][2]
    indd = np.where(nib.load(mask).get_data()>0)
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

    featquery_cmd = "featquery 1 %s/group.gfeat/cope1.feat 2 stats/pe1 stats/cope1 featquery -p -s -w -b %s" %(ConDir,mask)
    os.popen(featquery_cmd).read()

    featquery_file = ConDir+"group.gfeat/cope1.feat/featquery/report.txt"
    for line in open(featquery_file):
        list = line.split(" ")

    # devide by 2 (except motor) because the contrastmatrices are [1 -1]
    # motor is just an average effect
    if Contrast == "MOTOR":
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
        dict_writer.writerow([Contrast,Masks[Contrast][a][0],SIZE,COHENS10,COHENSD,COHENS90,PERCENT10,PERCENT,PERCENT90])
