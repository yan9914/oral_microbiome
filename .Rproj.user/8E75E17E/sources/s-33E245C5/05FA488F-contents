library(dplyr)
CPS2 <- read.table("CPS2_SraRunTable.txt", sep = ",", header = T)
PLCO <- read.table("PLCO_SraRunTable.txt", sep = ",", header = T)
metadata <- bind_rows(CPS2, PLCO)
complete.cases(metadata) %>% sum

write.csv(metadata, "oral_metadata.csv", row.names = FALSE)
