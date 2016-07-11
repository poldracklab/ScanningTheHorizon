## SubjectSelection

For the group analyses, we want:
- the same subjects in all analyses to avoid different effect sizes due to different subjects:
  * `subs_all_contrasts.sh` searches for subjects that have data for all experiments
- independent subjects to avoid an overestimation of the variance due to correlations among familymembers
  * `subs_unrelated.R` searches for unrelated subjects in the HCP data (need restricted access)
- `subs_all_contrasts_and_unrelated.R` combines both previous lists.
- `IDs_all_cons_and_unrelated.txt` is a csv file that lists all subjects used in this analyses and the disks on which they appear in the Connectome-in-a-box
