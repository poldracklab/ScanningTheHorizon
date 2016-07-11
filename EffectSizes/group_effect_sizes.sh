#!/bin/bash

#SBATCH --job-name=EffectSizes
#SBATCH --output=error/out.EffectSizes
#SBATCH --error=error/err.EffectSizes
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --partition=russpold
#SBATCH --qos=russpold

# Directory of Connectome in a box
ConnectomeInABoxDir=$HCPDIR
# Experiment to be analyzed
Paradigm=("MOTOR" "WM" "EMOTION" "GAMBLING")
# Which contrast from the experiment
Contrast=(7 11 3 6)
# Working Directory
WorkDir=/home/jdurnez/effect_sizes/HCP_group_effect_sizes/

# File with subjects to be included and disks on which they appear
SubjectsFile=$WorkDir/SubjectSelection/IDs_all_cons_and_unrelated.txt

for exp in {0..3} ; do

  ConDir=$WorkDir/03_GroupAnalyses/tfMRI_${Paradigm[exp]}/

  # Extract Effectsizes
  python $WorkDir/04_EffectSize/EffectSize.py ${Paradigm[exp]} $ConDir $WorkDir/

done
