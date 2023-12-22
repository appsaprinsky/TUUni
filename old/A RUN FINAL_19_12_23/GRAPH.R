setwd("/Users/admin/Desktop/Research Project/CODE/TUUni/FINAL")
getwd()

library(readr)
library(dplyr)
library(ggplot2)
data <- read_csv("global_SIMULATE_Info.csv")





#final_column_num 4, 5, 6, 7 for BRMSEA BGammaHat adjBGammaHat      BMc
# final_column_name_str "BRMSEA"  "BGammaHat" "adjBGammaHat"      "BMc"

Create2graphs <- function(data_MAIN, final_column_num, final_column_name_str){

  BRMSEA_data <- data_MAIN[,c(1,2,3,final_column_num)]
  BRMSEA_TM <- BRMSEA_data %>% filter( ModelTorF == 0.3)
  BRMSEA_TM <- BRMSEA_TM[,c(1,2,4)]
  BRMSEA_FM <- BRMSEA_data  %>% filter(ModelTorF == 0)
  BRMSEA_FM <- BRMSEA_FM[,c(1,2,4)]
  

  my_plot1 <-ggplot(BRMSEA_TM, aes_string(x = "PersonSize", y= final_column_name_str, color = "as.factor(TimePoint)")) +
    geom_line() +
    labs(x = "PersonSize", y = final_column_name_str, color = "TimePoint", title = paste(final_column_name_str, "For True Model", sep = " ") ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

  my_plot2 <- ggplot(BRMSEA_FM, aes_string(x = "PersonSize", y = final_column_name_str, color = "as.factor(TimePoint)")) +
    geom_line() +
    labs(x = "PersonSize", y = final_column_name_str, color = "TimePoint", title = paste(final_column_name_str, "For False Model", sep = " ") ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(paste(final_column_name_str, "T.jpg", sep = "_"), plot = my_plot1, width = 8, height = 6, units = "in")
  ggsave(paste(final_column_name_str, "F.jpg", sep = "_"), plot = my_plot2, width = 8, height = 6, units = "in")
}

#final_column_num 4, 5, 6, 7 for BRMSEA BGammaHat adjBGammaHat      BMc
# final_column_name_str "BRMSEA"  "BGammaHat" "adjBGammaHat"      "BMc"

Create2graphs(data, 4, "BRMSEA" )
Create2graphs(data, 5, "BGammaHat" )
Create2graphs(data, 6, "adjBGammaHat" )
Create2graphs(data, 7, "BMc" )




