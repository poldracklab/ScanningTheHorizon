# HCP_group_effect_sizes

In this repository, we explore the effect sizes of common psychological paradigms.  We focus on four experiments administered by the Human Connectome Project: an emotion task, gambling task, working memory task and motor task.  The pipeline for analysis:

1. **SubjectSelection**: We analyzed the data from 186 independent subjects.  In this folder, the code is shared on how we've selected these subjects so that (1) all subjects have results for all tasks and (2) there are no genetically related subjects in the analysis.
2. **GroupAnalyses**: The first level analyses were shared by the Human Connectome Project, we've performed second level analyses with FSL's flame: a linear mixed effects regression at each voxel, using generalized least squares with a local estimate of random effects variance. (description from [cobidas report, Nichols 2016)](http://biorxiv.org/content/biorxiv/early/2016/05/20/054262.full.pdf#page=40).  The specific contrasts that have been tested are:
  - MOTOR: average
  - EMOTION: faces - houses
  - GAMBLING: reward - punishment
  - WORKING MEMORY: 2-back - 0-back
3. **CreateMasks**:  The masks used for the analyses are the intersections of anatomical and functional masks for each contrast.
  - Functional:  We've created masks using [neurosynth](www.neurosynth.org).  We've used the search terms ["Motor","Emotion","Gambling","Working memory"] and have used meta-analyses with FDR-control at 0.01 and Forward Inference.
  - Anatomical: We've used Harvard-Oxford probabilistic atlas at p>0.25.  We've selected the following anatomical regions for the contrasts:

 | Contrast | Mask (intersected with neurosynth meta-analyse map) |
 | -------- | ---- |
 | Motor | Precentral gyrus |
 || Supplementary motor cortex |
 || Left Putamen |
 || Right Putamen |
 | Working memory | Middle frontal gyrus |
 | Emotion | Left amygdala |
 || Right amygdala |
 | Gambling | Left accumbens |
 || Right accumbens |

4. **EffectSize**: This script computes Cohen's D over the whole brain and takes the median (and the 10 and 90 percentile) within the masks created before.  Featquery is called to compute %BOLD change within the masks (and outputs as well the 10 and 90 percentile).
