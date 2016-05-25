#!/bin/bash

#######################################################
##                                                   ##
## This script loops over all subjects and saves the ##
## subject ID's of subjects who have results for all ##
## contrasts, to avoid different effect sizes due to ##
## the use of different subjects.                    ##
##                                                   ##
## The script uses the program s3cmd, which is       ##
## available here:                                   ##
## https://github.com/pcorliss/s3cmd-modification    ##
##                                                   ##
#######################################################

# Make an iterable array of all ID's

SubjString=$(s3cmd ls s3://hcp-openaccess/HCP_900/ | # list folder content
  sed -e 's!                       DIR   s3://hcp-openaccess/HCP_900/!!' | #strip folder
  head -n 898 | #strip last line
  sed -e 's!/!!' #remove last dash after each subject
  )
SubjList=($SubjString)

# Make an iterable array of paradigms and contrasts per paradigm

ParsList=('SOCIAL' 'EMOTION' 'RELATIONAL' 'LANGUAGE' 'GAMBLING' 'WM' 'MOTOR');
ConsList=(6 6 6 6 6 30 26)

# Loop over all subjects and check which ones have all contrasts

for sub in {0..897} ; do
  echo ${SubjList[sub]}

  cnt=0 #reset counter

  for par in {0..6}; do

    # Count how many copes this subject has for given paradigm
    NO=$(( $(s3cmd ls s3://hcp-openaccess/HCP_900/${SubjList[sub]}/MNINonLinear/Results/tfMRI_${ParsList[par]}/tfMRI_${ParsList[par]}_hp200_s4_level2vol.feat/ | wc | awk '{print $1}') - 1 ))

    # If number of copes is equal to expected number of copes --> increment variable
    if [ $NO == ${ConsList[par]} ]; then
      cnt=$((cnt+1))
    fi

  done

  # If subjects has all copes for all 7 paradigms --> write subject id to file
  if [ $cnt == "7" ]; then
    echo ${SubjList[sub]} >> IDs_all_cons.txt
  fi

done
