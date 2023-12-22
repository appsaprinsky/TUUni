# setwd("C:\\holger\\SEM\\modelfit\\stanversion")      HB!!!!!!!!!!!!!
setwd("/Users/admin/Desktop/Research Project/CODE/TUUni/FINAL/old run")
getwd()


library(readr)
library(dplyr)
library(writexl)

# Assuming your data frame is named 'rounded_data_frame'
write_xlsx(rounded_data_frame, "output_file.xlsx")


data1 <- read_csv("run1.csv")
data2 <- read_csv("run2.csv")
data3 <- read_csv("run3.csv")
data4 <- read_csv("run4.csv")



data_new<-merge(data1, data2, all = TRUE)
data_new<-merge(data_new, data3, all = TRUE)
data_new<-merge(data_new, data4, all = TRUE)

data_new



df<-data_new[!duplicated(data_new, fromLast = TRUE), ]
df

df %>% filter( ModelTorF == 0.0)

df %>% filter( ModelTorF == 0.3)




write_xlsx(round(df %>% filter( ModelTorF == 0.0) %>% filter( TimePoint == 1), 2), "output/0_1.xlsx")

write_xlsx(round(df %>% filter( ModelTorF == 0.0) %>% filter( TimePoint == 2), 2), "output/0_2.xlsx")

write_xlsx(round(df %>% filter( ModelTorF == 0.0) %>% filter( TimePoint == 3), 2), "output/0_3.xlsx")

write_xlsx(round(df %>% filter( ModelTorF == 0.0) %>% filter( TimePoint == 4), 2), "output/0_4.xlsx")

write_xlsx(round(df %>% filter( ModelTorF == 0.0) %>% filter( TimePoint == 5), 2), "output/0_5.xlsx")

write_xlsx(round(df %>% filter( ModelTorF == 0.0) %>% filter( TimePoint == 10), 2), "output/0_10.xlsx")

write_xlsx(round(df %>% filter( ModelTorF == 0.0) %>% filter( TimePoint == 15), 2), "output/0_15.xlsx")



write_xlsx(round(df %>% filter( ModelTorF == 0.3) %>% filter( TimePoint == 1), 2), "output/03_1.xlsx")

write_xlsx(round(df %>% filter( ModelTorF == 0.3) %>% filter( TimePoint == 2), 2), "output/03_2.xlsx")

write_xlsx(round(df %>% filter( ModelTorF == 0.3) %>% filter( TimePoint == 3), 2), "output/03_3.xlsx")

write_xlsx(round(df %>% filter( ModelTorF == 0.3) %>% filter( TimePoint == 4), 2), "output/03_4.xlsx")

write_xlsx(round(df %>% filter( ModelTorF == 0.3) %>% filter( TimePoint == 5), 2), "output/03_5.xlsx")

write_xlsx(round(df %>% filter( ModelTorF == 0.3) %>% filter( TimePoint == 10), 2), "output/03_10.xlsx")

write_xlsx(round(df %>% filter( ModelTorF == 0.3) %>% filter( TimePoint == 15), 2), "output/03_15.xlsx")




