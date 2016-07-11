sub_con <- read.table("IDs_all_cons.txt")
names(sub_con) <- "ID"
sub_unr_disk <- read.table("IDs_all_unrelated_disk.txt")$V1
sub_unr <- read.table("IDs_all_unrelated.txt")$V1
sub_unr <- data.frame(sub_unr,sub_unr_disk)
names(sub_unr) <- c("ID","disk")
total <- merge(sub_con,sub_unr,by="ID")

write.table(total,"IDs_all_cons_and_unrelated.txt",row.names=FALSE)
