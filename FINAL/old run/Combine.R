# setwd("C:\\holger\\SEM\\modelfit\\stanversion")      HB!!!!!!!!!!!!!
setwd("/Users/admin/Desktop/Research Project/CODE/TUUni/FINAL/old run")
getwd()


library(readr)
library(dplyr)
data1 <- read_csv("run1.csv")
data2 <- read_csv("run2.csv")
data3 <- read_csv("run3.csv")
data4 <- read_csv("run4.csv")



data_new<-merge(data1, data2, all = TRUE)
data_new<-merge(data_new, data3, all = TRUE)
data_new<-merge(data_new, data4, all = TRUE)

data_new



df<-data_new[!duplicated(data_new, fromLast = TRUE), ]

df %>% filter( ModelTorF == 0.0) %>% filter( TimePoint == 1)

df %>% filter( ModelTorF == 0.0) %>% filter( TimePoint == 2)

df %>% filter( ModelTorF == 0.0) %>% filter( TimePoint == 3)

df %>% filter( ModelTorF == 0.0) %>% filter( TimePoint == 4)

df %>% filter( ModelTorF == 0.0) %>% filter( TimePoint == 5)

df %>% filter( ModelTorF == 0.0) %>% filter( TimePoint == 10)

df %>% filter( ModelTorF == 0.0) %>% filter( TimePoint == 15)



df %>% filter( ModelTorF == 0.3) %>% filter( TimePoint == 1)

df %>% filter( ModelTorF == 0.3) %>% filter( TimePoint == 2)

df %>% filter( ModelTorF == 0.3) %>% filter( TimePoint == 3)

df %>% filter( ModelTorF == 0.3) %>% filter( TimePoint == 4)

df %>% filter( ModelTorF == 0.3) %>% filter( TimePoint == 5)

df %>% filter( ModelTorF == 0.3) %>% filter( TimePoint == 10)

df %>% filter( ModelTorF == 0.3) %>% filter( TimePoint == 15)

