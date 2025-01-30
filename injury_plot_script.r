###################################
########## Plot Creation ##########
###################################
#
#
#
#
#

# Load Packages
library(ggplot2) # Load ggplot2
library(ggdark) #Load ggdark
library(tidyr) # Load tidyr
library(xgboost) # Load XGBoost
library(xgboostExplainer) # Load XGboost Explainer
library(EnvStats) # Load EnvStats
library(splines) # Load Splines
library(grid) # Load grid
library(gridExtra) # Load grid extra
library(scales) # Load scales


# Create injury plot function
create_inj_plot <- function(data, variable, injury, n_seg = 100, 
                             filter_2 = TRUE){
  #'
  #' This function calculates the empirical distribution
  #' of the injured and healthy samples, it then multiplies
  #' the areas under the distribution by the number of healthy and
  #' injured samples to calculate the relative risk of injury at
  #' each segment of the distribution/across the range of the variable
  #' 
  #' @param data The dataset containing the metric data
  #' @param variable The variable to use for the plot
  #' @param injury A binary vector of injury
  #' @param n_seg The number of segments to split the distribution into
  #' @param filter_2 Remove first and last 5% of values - TRUE/FALSE
  #'
  #' @return A data frame with the start of each segment, the proportion
  #' of the distribution falling in each segment and the relative injury risk
  #' 
  #'
  #'
  
  # Create sequence of the range of the variable
  x <- seq(from = min(data[, variable], na.rm=T), # Start from minimum value 
           to = max(data[,variable], na.rm = T), # To maximum value
           length.out = n_seg) # Number of segments is output
  
  # Calculate empirical distribution for injured samples
  p2 <- pemp(unique(x), data[which(injury == 1), variable])
  # Calculate empirical distribution for healthy samples
  p3 <- pemp(unique(x), data[, variable])
  # Join range and distribution values
  t1 <- cbind.data.frame(unique(x), p2, p3)
  
  # Remove duplicated values
  
    mid_val <- round(mean(1:nrow(t1)))
    t1_1 <- t1[1:mid_val,]
    t1_2 <- t1[(mid_val + 1):nrow(t1),]
    t1_1_fil <- t1_1[!(duplicated(t1_1$p2, fromLast = TRUE) | 
                         duplicated(t1_1$p3, fromLast =TRUE)),]
    t1_2_fil <- t1_2[!(duplicated(t1_2$p2) | duplicated(t1_2$p3)),]
    t1_fil <- rbind.data.frame(t1_1_fil, t1_2_fil)
  
  
  # For each row
  for (v in 1:nrow(t1_fil)){
    # For the injured and healthy distributions
    for(l in 2:3){
      # If duplicated values, take mean value
      if(sum((t1[,l] ==t1_fil[v,l])) > 1){
        val_to_take <- round(mean(which(t1[,l] == t1_fil[v,l])))
        t1_fil[v,] <- val_to_take
      }
    }
  }
  
  # Create empty vectors to store area
  area_inj <- area_all <- rep(NA, nrow(t1_fil)) 
  # For each row
  for(j in 1:nrow(t1_fil)){
    # If first row
    if(j ==1){
      # Take area of intial segment
      area_inj[j] <- t1_fil$p2[j]
      area_all[j] <- t1_fil$p3[j]
      # If very last row
    } else if (j == nrow(t1_fil) + 1){
      area_inj[j] <- 1 - t1_fil$p2[j-1]
      area_all[j] <- 1 - t1_fil$p3[j-1]
    } else {
      # Take difference from previous row
      area_inj[j] <- t1_fil$p2[j] -  t1_fil$p2[j-1]
      area_all[j] <- t1_fil$p3[j] -  t1_fil$p3[j-1]
    }
    
  }
  # For last row
  t1_fil[j,] <- c(t1_fil[j-1, 1] + (t1_fil[j-1, 1] - t1_fil[j-2, 1]), 0, 0)
  # Add values to initial data frame 
  t1_fil$area_inj <- area_inj
  t1_fil$area_all <- area_all
  # If filter 2
  if(filter_2){
    # Calculate bottom five percent of values
    rem_1 <- floor(nrow(t1_fil) * 0.05)
    # Calculate top five percent of values
    rem_2 <- ceiling(nrow(t1_fil) * 0.95)
    # Remove extreme values
    t1 <- t1_fil[rem_1:rem_2, ]
  }
  
  # Add names to dataset
  names(t1_fil) <- c("variable", "inj", "all", "area_inj", "area_all")
  # Calculate injury risk
  t1_fil$test_inj <- ((t1_fil$area_inj * sum(injury == 1)) /
                        (t1_fil$area_all * length(injury)))
  # Return created data, injury vector and original data
  return(list(t1 = t1_fil, new_inj = injury, data = data))
}

# Create plot function
plot.function_2 <- function(dat_list, xgb_res, num_plot = 50, var_data, 
                            path = "injury_plots"){
  #'
  #' This function creates an injury risk zone plot for a given metric
  #' 
  #' 
  #' @param dat_list A list containing the metric dataset as its first item, 
  #' and the injury vector as its second item.
  #' @param xgb_res The variable importance to be used to identify the variables
  #' to plot
  #' @param num_plot The number of variables to produce a plot for
  #' @param var_data 
  #' @param path The path to be used for saving the plot
  #' 
  #' @return A list containing the generated plots
  #'
  #'
  #'
  
  # Identify number of plots to create as either number of plots or number
  # of available variables
  num_plot <- min(num_plot, nrow(xgb_res$imp_mat))
  # Extract variables to plot
  plot_vars <- xgb_res$Feature[1:num_plot]
  # Create a list to store risk zones
  zone_store <- vector(mode ="list", length = num_plot)
  # Create a list to store generated plots
  plot_list_1 <- plot_list_2 <- plot_list_3 <- vector(mode = "list", length = num_plot)
  
  # For each variable
  for(i in 1:length(plot_vars)){
    
    # Calculate injury risk
    temp <- create_inj_plot(data = dat_list[[1]], variable = plot_vars[i], 
                            injury = dat_list[[2]],
                            n_seg = 500, filter_2 =TRUE)
    # Store variable name
    n_name <- plot_vars[i]
    
    # Extract results from injury risk function
    t1 <- temp[[1]]
    n_inj <- temp[[2]]
    n_data <- temp[[3]]
    
    # Set weight as the number of overall samples falling in a given area
    t1$weight <- t1$area_all
    # Calculate proportion of injuries
    m_inj <- sum(n_inj==1)/length(n_inj)
    # Remove rows where area is negative
    t1 <- t1[t1$area_inj >= 0 & t1$area_all >= 0,]
    # Drop first row
    t1 <- t1[-1,]
    
    # Fit ggplot on data range and injury risk
    s_1 <- ggplot(t1, aes(x = variable, y =test_inj)) +
      geom_smooth(aes(weight = weight), n = nrow(t1), # Set weights as proportion of samples
                  method="glm",
                  method.args = list(family = "quasipoisson"), # Set distribution as quasipoisson
                  formula = y ~ ns(x, 2)) # Set formula for smoothing line
    
    # Build plot to extract smoothing line
    b_dat <- ggplot_build(s_1)$data[[1]]
    # Create empty vector to store risk zones
    col_vec <- rep(NA, nrow(b_dat))
    #If confidence region upper bound falls below injury rate
    col_vec[which(b_dat$ymax <= m_inj)] <- "blue"
    # If confidence region lower bound falls above injury rate
    col_vec[which(b_dat$ymin >= m_inj)] <- "red"
    # If smoothing line is below injury rate but confidence upper bound is above
    col_vec[which(b_dat$y <= m_inj & b_dat$ymax >= m_inj)] <- "lightblue"
    # If smoothing line is above injury rate but confidence lower bound is below
    col_vec[which(b_dat$y >= m_inj & b_dat$ymin <= m_inj)] <- "lightred"
    
    # Extract range of values
    x_vals <- b_dat$x
    # Create a data frame of rectables
    rects <- data.frame(xstart = b_dat$x[1:(nrow(b_dat) - 1)],
                        xend = b_dat$x[2:nrow(b_dat)],
                        Zone = col_vec[2:nrow(b_dat)])
    # Set factor levels and add zones which do not appear as levels
    zones <- c("blue", "lightblue", "lightred", "red")
    levels(rects$Zone) <- c(levels(rects$Zone), zones[which(!zones %in% rects$Zone)])
    
    # Create ggplot
    g_1 <- ggplot(t1, # Set data 
                  aes(x = variable, y= test_inj)) + # Set aesthetics
      geom_smooth(aes(weight = weight, linetype = "solid"), # Add smothing line
                  n = nrow(t1), 
                  method = "glm", 
                  method.args = list(family = "quasipoisson"),
                  formula = y ~ ns(x, 2),
                  color = "#5bc0de",
                  fill = "#dedede") +
      theme_classic() + # Set theme classic
      theme(panel.grid.major = element_blank(),  # Turn of grid
            panel.grid.minor = element_blank(),
            legend.background = element_rect(linetype="solid", colour ="lightgrey",
                                             size = 0.01), # Set legend background colors
            legend.title = element_text(size = 7),  # Set legend title size
            legend.text = element_text(size = 7), # Set legend text size
            legend.position = "top", # Set legend position
            legend.justification =  c(0,1),  # Set legend location
            axis.line = element_line(color ="grey50", size = 0.25), # Set axis line colors
            axis.ticks = element_line(color = "darkgrey", size = 0.25), # Set axis tick colors
            axis.title.x =element_text(face = "bold"),# Set axis title font
            axis.text.x = element_text("grey50"), # Set x-axis text color
            axis.text.y = element_text("grey50"), # Set y-axis text color
            axis.line.y.right = element_line(colour = "white"), # Set y-axis on right to white
            plot.title = element_text(face = "bold")) + # Set font for plot title
      geom_hline(linetype = 2,  # Add horizontal line at injury risk level
                 yintercept = m_inj) +
      ggtitle(paste("Injury Rate vs", n_name )) + # Set plot title
      xlab(paste(n_name)) + # Set xlabel
      ylab("Injury Rate") + # Set y-label
      geom_segment(data =rects, aes(y = -0.0005, # Add injury risk zone bar to plot
                                    yend = -0.0005,
                                    x = xstart,
                                    xend = xend,
                                    colour = Zone), 
                   size = 5, inherit.aes = FALSE)+
      scale_color_manual(breaks = c("blue" = "blue", "lightblue" = "lightblue", # Set colors manually
                                    "lightred" = "lightred", "red"= "red",
                                    "black" = "black"),
                         labels = c("blue" = "Low Risk", "lightblue" = "Minor Risk",
                                    "lightred" = "Moderate Risk", "red"= "Severe Risk",
                                    "black" = "Injury Rate"),
                         values =  c("blue" = "#ffdd55", "lightblue" = "#feb24c", 
                                     "lightred" = "#fd8d3c", "red"= "#ff2a2a",
                                     "black" = "black"),
                         drop = FALSE, # Keep all levels
                         aesthetics = c("colour")) +
      scale_linetype_manual(values = c("solid" = "solid"), # Set linetype label and value
                            labels = c("solid" = "Injury Rate Line & \nUncertainty Region")) +
      labs(linetype = "", color ="Risk Zone:") # Set labels
    
    # Add secondary axis to plot for baseline injury rate label
    g_1 <- g_1 + (scale_y_continuous(#breaks = sort(c(as.numeric(ggplot_build(g_1)$layout$panel_params[1]$y.labels))),
                                     labels = scales::percent_format(0.01), 
                                     sec.axis = sec_axis(~., breaks = m_inj,  labels = c("Baseline \nInjury Rate"))))
    # Generate plot
    g <- ggplotGrob(g_1)
    # Move title location to left corner
    g$layout$l[g$layout$name == "title"] <- 1
    # Uncomment line to save plot
    #ggsave(g, file=paste(path,"Injury_plot_", plot_vars[i], ".jpeg", sep = ""), width = 6, height = 6, dpi = 600)
    
    # Store plot in list
    plot_list_1[[i]] <- g
  }
  # Return plot
  return(plot_list_1)
}



# Create vector of variables to use
imp <- data.frame(Feature = c("HSR_ac_28_56", "total_distance_ac_28_56", "HSR_ac_7_21",
                              "total_time_ac_28_56", "decelerations_ac_28_56"))
# Create data list
dat_list <- list(met_dat, inj_dat$injuries_1)

# Create injury plots
inj_plots_1 <- plot.function_2(dat_list = dat_list, xgb_res = imp, num_plot = 5)

# View plots
grid.draw(inj_plots_1[[1]])
grid.draw(inj_plots_1[[2]])
grid.draw(inj_plots_1[[3]])
grid.draw(inj_plots_1[[4]])
grid.draw(inj_plots_1[[5]])


ggsave(inj_plots_1[[1]], file=paste("Injury_plot_", imp$Feature[1], ".jpeg", sep = ""), width = 6, height = 6, dpi = 600)
ggsave(inj_plots_1[[2]], file=paste("Injury_plot_", imp$Feature[2], ".jpeg", sep = ""), width = 6, height = 6, dpi = 600)
ggsave(inj_plots_1[[3]], file=paste("Injury_plot_", imp$Feature[3], ".jpeg", sep = ""), width = 6, height = 6, dpi = 600)
ggsave(inj_plots_1[[4]], file=paste("Injury_plot_", imp$Feature[4], ".jpeg", sep = ""), width = 6, height = 6, dpi = 600)
ggsave(inj_plots_1[[5]], file=paste("Injury_plot_", imp$Feature[5], ".jpeg", sep = ""), width = 6, height = 6, dpi = 600)


