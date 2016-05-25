# Directory of Connectome in a box
ConnectomeInABoxDir=$HCPDIR
# Experiment to be analyzed
Paradigm=("tfMRI_MOTOR" "tfMRI_WM" "tfMRI_LANGUAGE" "tfMRI_EMOTION" "tfMRI_GAMBLING")
# Which contrast from the experiment
Contrast=(7 11 2 3 6)
# Working Directory
WorkDir=`pwd`
# File with subjects to be included and disks on which they appear
SubjectsFile=$WorkDir/SubjectSelection/IDs_all_cons_and_unrelated.txt

for exp in {0..4} ; do

  ConDir=$WorkDir/GroupAnalyses/${Paradigm[exp]}/
  mkdir $ConDir

  # Make design.fsf and design.grp
  python $WorkDir/GroupAnalyses/make_design.py $SubjectsFile $ConnectomeInABoxDir ${Paradigm[exp]} ${Contrast[exp]} $ConDir

  # Analyze group stats
  feat $ConDir/design

  # Extract Effectsizes
  python $WorkDir/EffectSize/EffectSize.py ${Paradigm[exp]} $ConDir $WorkDir/EffectSize/

done
