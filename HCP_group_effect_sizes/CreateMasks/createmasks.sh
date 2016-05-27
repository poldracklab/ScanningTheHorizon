#NEUROSYNTH
fslmaths masks/neurosynth/emotion_pAgF_z_FDR_0.01.nii.gz -bin masks/EMOTION_NS.nii.gz
fslmaths masks/neurosynth/gambling_pAgF_z_FDR_0.01.nii.gz -bin masks/GAMBLING_NS.nii.gz
fslmaths masks/neurosynth/motor_pAgF_z_FDR_0.01.nii.gz -bin masks/MOTOR_NS.nii.gz
fslmaths masks/neurosynth/WM_pAgF_z_FDR_0.01.nii.gz -bin masks/WM_NS.nii.gz

#ANATOMICAL
fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz masks/MOTOR_1_HO.nii.gz 6 1
fslmaths masks/MOTOR_1_HO -thr 0.25 -bin masks/MOTOR_1_HO.nii.gz
fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz masks/MOTOR_2_HO.nii.gz 25 1
fslmaths masks/MOTOR_2_HO -thr 0.25 -bin masks/MOTOR_2_HO.nii.gz
fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz masks/MOTOR_3_HO.nii.gz 5 1
fslmaths masks/MOTOR_3_HO -thr 0.25 -bin masks/MOTOR_3_HO.nii.gz
fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz masks/MOTOR_4_HO.nii.gz 16 1
fslmaths masks/MOTOR_4_HO -thr 0.25 -bin masks/MOTOR_4_HO.nii.gz
fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz masks/WM_HO.nii.gz 3 1
fslmaths masks/WM_HO -thr 0.25 -bin masks/WM_HO.nii.gz
fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz masks/EMOTION_1_HO.nii.gz 9 1
fslmaths masks/EMOTION_1_HO -thr 0.25 -bin masks/EMOTION_1_HO.nii.gz
fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz masks/EMOTION_2_HO.nii.gz 19 1
fslmaths masks/EMOTION_2_HO -thr 0.25 -bin masks/EMOTION_1_HO.nii.gz
fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz masks/GAMBLING_1_HO.nii.gz 10 1
fslmaths masks/GAMBLING_2_HO -thr 0.25 -bin masks/GAMBLING_1_HO.nii.gz
fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz masks/GAMBLING_2_HO.nii.gz 20 1
fslmaths masks/GAMBLING_2_HO -thr 0.25 -bin masks/GAMBLING_2_HO.nii.gz

#CONJUNCTION
fslmaths masks/EMOTION_1_HO -mul masks/EMOTION_NS masks/EMOTION_1_CJ
fslmaths masks/EMOTION_2_HO -mul masks/EMOTION_NS masks/EMOTION_2_CJ
fslmaths masks/GAMBLING_1_HO -mul masks/GAMBLING_NS masks/GAMBLING_1_CJ
fslmaths masks/GAMBLING_2_HO -mul masks/GAMBLING_NS masks/GAMBLING_2_CJ
fslmaths masks/WM_HO -mul masks/GAMBLING_NS masks/WM_CJ
fslmaths masks/MOTOR_1_HO -mul masks/MOTOR_NS masks/MOTOR_1_CJ
fslmaths masks/MOTOR_2_HO -mul masks/MOTOR_NS masks/MOTOR_2_CJ
fslmaths masks/MOTOR_3_HO -mul masks/MOTOR_NS masks/MOTOR_3_CJ
fslmaths masks/MOTOR_4_HO -mul masks/MOTOR_NS masks/MOTOR_4_CJ
