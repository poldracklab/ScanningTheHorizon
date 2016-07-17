This folder contains code to process the sample size estimates for use in
generating Figure 1. This must be run to generate the files used in the
Figure 1 notebook.

To run all analyses, use <p>make all</p>

### Individual analyses:

#### David et al. data

Manually annotated sample size data were obtained from Sean David, Pooja Loftus, and
Isabella Chu, which were originally created for this paper:

[David SP, Ware JJ, Chu IM, Loftus PD, Fusar-Poli P, Radua J, et al. (2013) Potential Reporting Bias in fMRI Studies of the Brain. PLoS ONE 8(7): e70104. doi:10.1371/journal.pone.0070104](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0070104)


The original data are in the file titled "fMRI MA significance bias database 03-24-13.xlsx"

These data are processed using get_sampsize_david.py which generates the text
file "david_sampsizedata.txt" which is used for the Figure 1 notebook. In particular,
because pubmed IDs were not available for all entries, we used the CrossRef API to
obtain Pubmed IDs for all papers (excluding those for whom a unique Pubmed ID could
  not be obtained).  The publication date for each paper was then obtained from
  the Pubmed record using the Entrez API.

#### Neurosynth data

Automatically annotated data were obtained from neurosynth using the
code in get_ns_sample_sizes.py based on the data in abstracts.txt,
which generated the file estimated_n.txt. These
were manually annotated by Joe Wexler.  This is a description of Joe's
annotation process:

>1) I wrote a script to change the color of all the numbers and number words to make it easier for me to find them visually.
2) I looked through these numbers searching for ones that describe the number of subjects/participants in the study or in a particular group in the study.
3) I made a note of how many groups I thought there were and how many participants are in each. You mentioned not to worry about some of the finer distinctions so I tried to group them according to the criteria that was relevant to the goal of the study, and generally to features that were known/decided before the study rather than discovered differences. I think I also tried to make sure that the sum of the groups equaled the total number of subjects in the study, so if there were multiples ways of dividing up the subjects, I choose what seemed like the best one to avoid including the same subject more than once. I also gave Tal's program the benefit of the doubt a little bit, meaning if I wasn't sure how to divide the subjects into groups, I went with the way the program had done it.
4) I recorded the number of participants per group in column C, with each row corresponding to a group. If there were multiple groups, I wrote the word 'group' in column D for each group, and if there was only one, I wrote 'study'. If my number of groups for a study was greater than that already present, I added rows, filling in columns C and D and leaving columns B and E blank. If I thought there fewer groups, I wrote 'n/a' in the columns C and D for the rows I thought were unnecessary.

These data were saved into the file estimated_n_format_2011-2015numbers.csv.  They were
further processed using get_sampsize_neurosynth.py which retreived the paper
metadata from the Entrez API and also filtered papers based on the inclusion
of terms relevant to fMRI (to exclude PET and structural studies). This generated
the files "neurosynth_study_data.txt" (which includes sample sizes for single-group
 studies) and "neurosynth_group_data.txt" (which includes sample sizes for multiple-group
   studies, one group per line).
