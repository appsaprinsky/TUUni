setwd("/Users/admin/Desktop/Research Project/CODE/TUUni/FINAL/old run")
getwd()

library(readr)
library(dplyr)
library(ggplot2)
data <- read_csv("final_combined.csv")





#final_column_num 4, 5, 6, 7 for BRMSEA BGammaHat adjBGammaHat      BMc
# final_column_name_str "BRMSEA"  "BGammaHat" "adjBGammaHat"      "BMc"

Create2graphs <- function(data_MAIN, final_column_num, final_column_name_str){

  BRMSEA_data <- data_MAIN[,c(1,2,3,final_column_num)]
  BRMSEA_TM <- BRMSEA_data %>% filter( ModelTorF == 0)
  BRMSEA_TM <- BRMSEA_TM[,c(1,2,4)]
  BRMSEA_FM <- BRMSEA_data  %>% filter(ModelTorF == 0.3)
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

# Create2graphs(data, 4, "BRMSEA" )
# Create2graphs(data, 5, "BGammaHat" )
# Create2graphs(data, 6, "adjBGammaHat" )
# Create2graphs(data, 7, "BMc" )
 

Create3graphs <- function(data_MAIN, final_column_num, final_column_name_str){
  
  BRMSEA_data <- data_MAIN[,c(1,2,3,final_column_num)]
  mac <- c(1, 2, 3, 4, 5, 10)
  BRMSEA_TM <- BRMSEA_data# %>% filter( ModelTorF == 0)
  
  
  
  
  my_plot1 <-ggplot(BRMSEA_TM, aes_string(x = "PersonSize", y= final_column_name_str, color = "as.factor( TimePoint)")) +
    geom_line() +
    facet_wrap(~ ModelTorF, scales = "free_y", nrow = 3, ncol = 2) +
    labs(x = "PersonSize", y = final_column_name_str, color = "TimePoint", title = paste(final_column_name_str, "For True Model", sep = " ") ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(paste(final_column_name_str, "Time.jpg", sep = "_"), plot = my_plot1, width = 8, height = 6, units = "in")
}

Create3graphs(data, 4, "BRMSEA" )
Create3graphs(data, 5, "BGammaHat" )
Create3graphs(data, 6, "adjBGammaHat" )
Create3graphs(data, 7, "BMc" )



Create4graphs <- function(data_MAIN, final_column_num, final_column_name_str){
  
  BRMSEA_data <- data_MAIN[,c(1,2,3,final_column_num)]
  mac <- c(1, 2, 3, 4, 5, 10)
  BRMSEA_TM <- BRMSEA_data# %>% filter( ModelTorF == 0)
  
  
  

  my_plot1 <-ggplot(BRMSEA_TM, aes_string(x = "PersonSize", y= final_column_name_str, color = "as.factor( ModelTorF)")) +
    geom_line() +
    facet_wrap(~ TimePoint, scales = "free_y", nrow = 3, ncol = 2) +
    labs(x = "PersonSize", y = final_column_name_str, color = "ModelTorF", title = paste(final_column_name_str, "For True Model", sep = " ") ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  

  ggsave(paste(final_column_name_str, "0_03.jpg", sep = "_"), plot = my_plot1, width = 8, height = 6, units = "in")
}

Create4graphs(data, 4, "BRMSEA" )
Create4graphs(data, 5, "BGammaHat" )
Create4graphs(data, 6, "adjBGammaHat" )
Create4graphs(data, 7, "BMc" )



