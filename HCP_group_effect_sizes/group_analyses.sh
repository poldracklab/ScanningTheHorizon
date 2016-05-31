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
Paradigm=("tfMRI_MOTOR" "tfMRI_WM" "tfMRI_EMOTION" "tfMRI_GAMBLING")
# Which contrast from the experiment
Contrast=(7 11 3 6)
# Working Directory
WorkDir=/home/jdurnez/effect_sizes/HCP_group_effect_sizes/

# File with subjects to be included and disks on which they appear
SubjectsFile=$WorkDir/SubjectSelection/IDs_all_cons_and_unrelated.txt

for exp in {0..3} ; do

  ConDir=$WorkDir/GroupAnalyses/${Paradigm[exp]}/
  mkdir $ConDir

  # Make design.fsf and design.grp
  python $WorkDir/GroupAnalyses/make_design.py $SubjectsFile $ConnectomeInABoxDir ${Paradigm[exp]} ${Contrast[exp]} $ConDir

  # Analyze group stats
  feat $ConDir/design

done
