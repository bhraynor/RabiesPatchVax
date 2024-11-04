#############################################################################
#dependencies
#############################################################################
library(dplyr) #tidy coding
library(ggplot2) #pretty figures
library(sf) #simple features
library(viridis) #good color palettes
library(reshape2) #formating


#############################################################################
#Rabies dynamics plots
#############################################################################
#Input: ODE_output, cell.names, init
#         -ODE_output: the matrix that the deSOLVE ODE solver returrns
#         -cell.names: the vector of cell names
#         -init: the vector of initial conditions fed into the ODE solver

#Output:
#   fig1 = list where each list item is a figure of total dynamics per a single spatial cell (sub-division)
#   fig2 = list where each list item is a figure of E/I state dynamics per a single spatial cell (sub-division)
#   fig3 = figure looking at infectious disease states per cell and overall in one plot
#   fig4 = chloropleth map of infected dogs for entire simlation
#   fig5 = chloropleth map of infected dogs for last week of simulation

#----------------------------------------------------
# #Example:
# source("Functions_SimulationVisualization.R")
#----------------------------------------------------

rfun_RabiesDynamics <- function(ODE_output, name_patch){ #work on wrapping in fx

  
  #calculate statistics
  name_comps <- colnames(out%>%select(-time, -sim))
  num_comps <- length(name_comps)
  name_states <- unique(substr(name_comps,1,1))
  num_states <- length(name_states)
  num_patch <- num_comps/num_states
  end_time <- max(ODE_output$time)

  #set color palettes
  Pal1 <- c('Susceptible' = 'gold',
            'Exposed' = 'orange1',
            'Infectious' = 'red3',
            'Vaccinated' = 'dodgerblue3')
  Pal2 <- c('All districts' = 'red4',
            'Single district' = 'red2' )
  
  
  df <- ODE_output
  
  df.dist <- melt(df, id=c('time','sim'))%>%
    mutate(state = substr(variable, 1, 1), 
           district.num = substr(variable, 2, 10))
  
  df.global <- df.dist %>%
    group_by(state,time, sim)%>%
    summarize(
      total =sum(value))%>%
    ungroup()%>%
    tidyr::spread(key=state, value=total)
  
  #ALL STATES
  fig1 <- list()
  fig2 <- list()
  for(i in 1:num_patch){
    sub <- df.dist%>% 
      filter(district.num == i)%>%
      select(-variable, -district.num)%>%
      tidyr::spread(key=state, value=value)
    
    sub <- sub %>%
      mutate(I = ifelse(I==0, NA, I),
             E= ifelse(E==0, NA, E))%>%
      mutate(I = ifelse(is.na(I) & time == 0, 0, I),
             E = ifelse(is.na(E) & time == 0, 0, E))
    

    fig.sub1 <- ggplot()+
      theme_classic()+ #white background w/ out grid
      geom_line(data=sub, aes(x=time, y=S, group=sim, color="Susceptible"), size=1, alpha=0.5)+
      geom_line(data=sub, aes(x=time, y=E, group=sim, color="Exposed"), size=1, alpha=0.5)+
      geom_line(data=sub, aes(x=time, y=I, group=sim, color="Infectious"), size=1, alpha=0.5)+
      geom_line(data=sub, aes(x=time, y=V, group=sim, color="Vaccinated"), size =1, alpha=0.5)+
      scale_color_manual(values = Pal1, name= "Diease state", 
                         breaks=c("Susceptible","Exposed","Infectious","Vaccinated")) +
      theme(text = element_text(size=20),
            axis.text.x = element_text(size=16))+
      labs(y= "Number of dogs", x = "Time (days)")+ #axis labels
      ggtitle(paste0(name_patch[i]))

    fig.sub2 <- ggplot()+
      theme_classic()+ #white background w/ out grid
      #geom_line(data=sub, aes(x=time, y=E, group=sim, color="Exposed"), size=1.5, alpha=0.5)+
      geom_line(data=sub, aes(x=time, y=I, group=sim, color="Infectious"), size=1.5, alpha=0.5)+
      scale_color_manual(values = Pal1, name= "Diease state", 
                         breaks=c("Susceptible","Exposed","Infectious","Vaccinated")) +
      theme(text = element_text(size=20),
            axis.text.x = element_text(size=16))+
      labs(y= "Number of dogs", x = "Time (days)")+ #axis labels
      scale_x_continuous(expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0))+
      ggtitle(paste0(name_patch[i]))
    
    
    fig1[[i]] <- fig.sub1
    fig2[[i]] <- fig.sub2
  }
  
  #Figure 3 - all I lines on one plot
  df.dist_I <- df.dist %>%
    filter(state == "I") %>%
    select(-state, -district.num)%>%
    tidyr::spread(key= variable, value=value)
  df.dist_I$total <- rowSums(df.dist_I%>%select(-time, -sim))
  
  df.I_plot <- df.dist_I %>% 
    mutate(total = ifelse(total==0, NA, total))%>%
    mutate(time = time-3)%>%
    filter(time >= 0)
  
  
  fig3 <- ggplot()+
    theme_classic()+ #white background w/ out grid
    geom_line(data=df.I_plot, 
              aes(x=time, y=total, group=sim, color="Infectious"), 
              size=1, alpha=0.5)+
    scale_color_manual(values = Pal1, name= "Diease state", 
                     breaks=c("Susceptible","Exposed","Infectious","Vaccinated")) +
    #geom_line(data=df.dist_I, aes(x=time, y=total, group="sim", color="All districts"), size=1, alpha=0.8)+
    # scale_color_manual(values = Pal2, name= "Infectious dogs", 
    #                    breaks=c("Single district","All districts")) +
    theme(text = element_text(size=20),
          axis.text.x = element_text(size=16))+
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0))+
    labs(y= "System total", x = "Time (days)") #axis labels
  
  # #Figure 4 - cloropleth maps - total sim cases
  # dist_AQPcity <-st_read("~/RabiesLabPeru/peru_spatial_data/data_modified/01_sf_data/dist_AQPcity.shp")
  # 
  # 
  # df_means <- df%>%
  #   select(-sim)%>%
  #   group_by(time) %>%
  #   summarise_all(funs(mean))
  # 
  # df_dist_means <- melt(df_means, id=c('time'))%>%
  #   mutate(state = substr(variable, 1, 1), 
  #          district.num = substr(variable, 2, 10))
  # 
  # final_I <- df_dist_means%>%
  #   filter(state=="I" & time==end_time)
  # 
  # dist_AQPcity$SimulatedCasesEnd <- final_I$value
  # 
  # fig4 <- ggplot() + 
  #   theme_void()+
  #   ggtitle("End cases")+
  #   geom_sf(data = dist_AQPcity, aes(fill = SimulatedCasesEnd), size=1)+
  #   scale_fill_viridis_c(option = "A", direction = -1)
  
  fig1 <<- fig1
  fig2 <<- fig2
  fig3 <<- fig3
  #fig4 <<- fig4
  }
  