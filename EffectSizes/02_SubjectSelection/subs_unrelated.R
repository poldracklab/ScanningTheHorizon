########################################
## EXTRACT UNRELATED SUBJECTS FOR HCP ##
########################################

# Read in first 6 columns of restricted data, make a new variable 'fam'
file <- "RESTRICTED_jdurnez_11_17_2015_8_48_33.csv"
restr <- read.table(file,header=TRUE,sep=",",)[,1:6]
restr$fam <- NA
.
# Read in a textfile with the subjects ID's per disk on the Connectome-in-a-box.
# This will allow to create a file that tells for each subject on which disk they are.
# if all subjects are in the same folder, you won't need this file.
#
# to obtain the .txt-files do the following command in bash in the Connectome-in-a-box folder
# for i in {1..5}; do ls $HCPDIR/Disk$i\of5/ > disk$i.txt; done

dsk <- subs <- c()
for(d in 1:5){
  a <- read.table(paste("disk",d,".txt",sep=""))$V1
  dsk <- c(dsk,rep(d,length(a)))
  subs <- c(subs,a)
}

# Loop over all subjects and find out which other subjects share the same mother OR father
# Give all subjects that share a mother or a father the same family ID number (f)
f <- 0
for(sub in 1:541){
  ind <- unique(c(which(restr$Mother_ID[sub] == restr$Mother_ID),which(restr$Father_ID[sub] == restr$Father_ID)))
  if(is.na(restr$fam[sub])){
    f <- f+1
    restr$fam[ind] <- f
  }
}

# Select one sibling per family (the first one)
# Write away the subject IDs and the disks on which these subjects are
first_sib <- c()
disk <- c()
for(fam in 1:length(unique(restr$fam))){
  new <- restr$Subject[which(restr$fam == fam)[1]]
  if (sum(new==subs)>0){
    first_sib <- c(first_sib,new)
    disk <- c(disk,dsk[which(new== subs)])
  }
}

#Write away
write.table(first_sib,"IDs_all_unrelated.txt",col.names=FALSE,row.names=FALSE)
write.table(disk,"IDs_all_unrelated_disk.txt",col.names=FALSE,row.names=FALSE)
