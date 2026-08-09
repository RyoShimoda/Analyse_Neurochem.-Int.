# Analyse for Neurochemistry International

# Library ----
library(car)
library(PerformanceAnalytics)

# Set working directry ----
setwd("Analyse_Neurochem_int")

# My Function & ANOVAKUN----
source("scripts/AUTOGRAPH.txt")
source("scripts/anovakun_489.txt")
source("scripts/MF_IHC.txt")


# Contextual Fear Conditioning ----
  # Create Datasets & Plots ----
    mf_autograph(folder_name = ".",
                 experiment_type = "FC",
                 path = ".\\FC", 
                 plot_number = "_FC")

    mf_FC_annotation(folder_name = ".", 
                          graph = GFC_FC, 
                          y_first_arrow = 25, y_second_arrow = 45, y_third_arrow3 = 65, 
                          legend_posision_y = 65)
    
  # Statistical Analysis ----
    sink(file = "Results/ANOVA_FC.txt", split = T)
    
    anovakun(dataset = sumFC_FC[-1],
             design = "AsB",
             Group = c("SED", "LIE"),
             Time = TimeFC,
             hf = T,
             peta = T)
    sink()
    
# Contextual Fear Extinction ----
  # Create Datasets & Plots ----
    mf_autograph(folder_name = ".", 
                     path = ".\\Ex1", 
                     experiment_type = "Ex1")
    
  # Statistical Analysis ---- 
    AnodataEx1 <- dataEx1 %>% 
      group_by(No, Group) %>% 
      summarise(Freezing = mean(Freezing))
    
    sink(file = "Results/ANOVA_Ex1.txt", split = T)
    
    anovakun(dataset = sumEx1per3[-1], "AsB",
             Group = c("SED", "LIE"),
             Time = Timeper3,
             hf = T, 
             peta = T)
    
    cat("\n-- Shapiro-Wilk test--\n")
    shapiro.test(AnodataEx1$Freezing)
    
    cat("\n-- Levene test --\n")
    leveneTest(AnodataEx1$Freezing, AnodataEx1$Group)
    
    cat("\n-- Welch t test --\n")
    t.test(AnodataEx1$Freezing ~ AnodataEx1$Group)
    sink()
    
# Muscle Sampling ----
  # Create Datasets & Plots ----
    muscle_list <- c("Soleus", "Plantaris", "Adrenal", "Thymus")
    for (i in muscle_list) {
      mf_sampling(folder_name = ".", 
                       dataset = "Sampling_after_extinction_day1.xlsx", 
                       sampname = i)
    }
    
    modSoleus <- gSoleus +
      scale_y_continuous(expand = c(0,0), limits = c(0,60), breaks = c(0,10,20,30,40,50)) +
      theme(axis.line.y  = element_blank()) +
      annotate("segment", x = .3, xend = .3,  y = 0, yend = 50.3, size = 1) +
      annotate("path", x = c(1,1,2,2), y = c(47, 55, 55, 52)) +
      annotate("text", x = 1.5, y = 57, label = "*", size = 6)
    modSoleus
    
    tiff(filename = "Results/Plot/Soleus.tiff", 
         width = 4 * 900, height = 3.5 * 900, res = 900)
    modSoleus
    dev.off()
    
    modPlantaris <- gPlantaris +
      scale_y_continuous(expand = c(0,0), limits = c(0, 150), breaks = c(0,25,50, 75,100,125)) +
      theme(axis.line.y = element_blank()) +
      annotate("segment", x = .3, xend = .3, y = 0, yend = 125.3, size = 1) +
      annotate("path", x = c(1,1,2,2), y = c(105, 125, 125, 115)) +
      annotate("text", x = 1.5, y = 127, label = "*", size = 6)
    modPlantaris
    
    tiff(filename = "Results/Plot/Plantaris.tiff", 
         width = 4 * 900, height = 3.5 * 900, res = 900)
    modPlantaris
    dev.off()
    
  # Statistical Analysis ----
    anomuscle <- Sampdata %>% 
      mutate("Soleus" = Soleus/Body_weight * 100) %>% 
      mutate("Plantaris" = Plantaris/Body_weight * 100) %>% 
      mutate("Adrenal" = Adrenal/Body_weight * 100) %>% 
      mutate("Thymus" = Thymus/Body_weight * 100)
    
    sink(file = "Results/Sampling_ANOVA.txt", split = T)

    cat("\n== Body Weight ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(anomuscle$Body_weight)
    
    cat("\n-- Levene's Test --\n")
    leveneTest(anomuscle$Body_weight, anomuscle$Group)
    
    cat("\n-- unpaired t test --\n")
    t.test(anomuscle$Body_weight ~ anomuscle$Group, var.equal  = TRUE)
    
    cat("\n== Soleus ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(anomuscle$Soleus)
    
    cat("\n-- Levene's Test --\n")
    leveneTest(anomuscle$Soleus, anomuscle$Group)
    
    cat("\n-- unpaired t test --\n")
    t.test(anomuscle$Soleus ~ anomuscle$Group, var.equal  = TRUE)
    
    
    cat("\n== Plantaris ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(anomuscle$Plantaris)
    
    cat("\n-- Levene's Test --\n")
    leveneTest(anomuscle$Plantaris, anomuscle$Group)
    
    cat("\n== unpaired t test ==\n")
    t.test(anomuscle$Plantaris ~ anomuscle$Group, var.equal = TRUE)
    
    
    cat("\n== Adrenal ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(anomuscle$Adrenal)
    
    cat("\n-- Levene's Test --\n")
    leveneTest(anomuscle$Adrenal, anomuscle$Group)
    
    cat("\n ===== unpaired t test ===== \n")
    t.test(anomuscle$Adrenal ~ anomuscle$Group, var.equal = TRUE)
    
    cat("\n== Thymus ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(anomuscle$Thymus)
    
    cat("\n-- Levene's Test --\n")
    leveneTest(anomuscle$Thymus, anomuscle$Group)
    
    cat("\n-- unpaired t test --\n")
    t.test(anomuscle$Thymus ~ anomuscle$Group, var.equal = TRUE)
    
    sink()
    

        
# Immunohistochemistry =========================================================
  # My Function
    # Create Lists----
    RegionList <- c("dDGsp", "dDGip", "dCA3", "dCA2", "dCA1", "vDGsp", "vDGip", "vCA3", "vCA1")
    PFCList <- c("PL", "IL")
    RegionListAll <- c("PL", "IL", "dDGsp", "dDGip", "dCA3", "dCA2", "dCA1", "vDGsp", "vDGip", "vCA3", "vCA1")
    
    # Use My Function----
    mf_IHC(Folder_name = "c-Fos\\Hippocampus",
                Path = ".\\c-Fos\\Hippocampus\\Analyse_G50_R20",
                Region_list = RegionList, Red_name = "c-Fos", Green_name = "NeuN", 
                Result_folder = "Result_Analyse_G50_R20")
    
    mf_IHC(Folder_name = "c-Fos\\PL_IL",
                Path = ".\\c-Fos\\PL_IL\\G40_R20",
                Region_list = PFCList, 
                Red_name = "c-Fos", Green_name = "NeuN", PFC = "", 
                number_of_region = 2, D_V = "", 
                Result_folder = "Result_G40_R20", 
                save_plot_width_all = 3, save_plot_width_i = 6, save_plot_height_i = 6,
                title_size_i = 8, axis_text_size_i = 8, 
                axis_title_size_i = 8, save_histplot_width = 7.5, save_histplot_height = 3,
                legend_position = c(.8, .9))
    
    # Read KatiKati Data ----
    Kachidata_D4 <- read.csv("c-Fos\\Hippocampus\\ForAnalyse.csv", header = F) %>% 
      mutate(Group = ifelse(str_detect(V1, "SED"), "SED", "LIE")) %>% 
      mutate(Region = if_else(str_detect(V1, pattern = "dCA1"), "dCA1",
                              if_else(str_detect(V1, pattern = "dCA2"), "dCA2",
                              if_else(str_detect(V1, pattern = "dCA3"), "dCA3",
                              if_else(str_detect(V1, pattern = "dDGsp"), "dDGsp",
                              if_else(str_detect(V1, pattern = "dDGip"), "dDGip",
                              if_else(str_detect(V1, pattern = "vCA1"), "vCA1",
                              if_else(str_detect(V1, pattern = "vCA3"), "vCA3",
                              if_else(str_detect(V1, pattern = "vDGsp"), "vDGsp",
                              if_else(str_detect(V1, pattern = "vDGip"), "vDGip", "")))))))))) %>% 
      mutate(No = substring(V1 ,1, 6)) %>% 
      rename("cFos" = "V2") %>% 
      group_by(No, Group, Region) %>% 
      summarise(cFos = sum(cFos)) %>% 
      arrange(Region)
    
    # Create Plot Data ----
    GData_Hipp <- GatheringData_Hipp %>% 
      mutate(Region = as.character(Region)) %>% 
      mutate(No = substring(No, 1, 6)) %>% 
      arrange(Region)
    
    plotdata_d4 <- bind_cols(GData_Hipp, Kachidata_D4) %>%
      # Checking rows and columns 
      dplyr::select("No...1", "No...13","Group...2", "Group...14","Region...3", "Region...15","Area_G","cFos") %>%
      dplyr::select("No...1", "Group...2", "Region...3", "Area_G", "cFos") %>% 
      rename("No" = "No...1") %>% 
      rename("Group" = "Group...2") %>%
      rename("Region" = "Region...3") %>% 
      mutate(cFos_st = cFos/Area_G * 1000000 / 0.04) %>%
      mutate(Group = as.factor(Group)) %>%
      mutate(Group = relevel(Group, ref = "SED"))
    
    # Create DG combined Data ----
    DGdata_d4 <- plotdata_d4 %>% 
      filter(str_detect(Region, pattern = "DG")) %>% 
      mutate(DV = if_else(str_detect(Region, pattern = "d"), "Dorsal", "Ventral")) %>% 
      group_by(No, Group, DV) %>% 
      summarise(Area_G = sum(Area_G),
                cFos = sum(cFos)) %>% 
      mutate(Region = if_else(str_detect(DV, pattern = "Dorsal"), "dDG", "vDG")) %>% 
      dplyr::select("No", "Group", "Region", "Area_G","cFos")
    
    # Create Plot Data ----
    plotdata_DGcomb_d4 <- plotdata_d4 %>% 
      filter(str_detect(Region, pattern = "CA")) %>% 
      bind_rows(DGdata_d4) %>% 
      mutate(cFos_st = cFos/Area_G * 1000000 / 0.04) %>% 
      mutate(Group = as.factor(Group)) %>% 
      mutate(Group = relevel(Group, ref = "SED")) %>% 
      dplyr::select("No", "Group", "Region", "cFos_st")
    
    
    # Create DG combined List
    RegionListAll_DGcomb <- c("PL", "IL", "dDG", "dCA3", "dCA2", "dCA1", "vDG", "vCA3", "vCA1")
    RegionList_DGcomb <- c("dDG", "dCA3", "dCA2", "dCA1", "vDG", "vCA3", "vCA1")
    
    
    # Read PL/IL Data ----
    PFC <- GatheringData_PFC %>%
      dplyr::select("No", "Group", "Region", "Count") %>% 
      mutate(cFos_st = Count/0.04) %>% 
      mutate(No = gsub("_1_R", "", No)) %>% 
      dplyr::select(-"Count")
    
    # Create Plot Data ----
    sumHipp_DGcomb <- plotdata_DGcomb_d4 %>% 
      group_by(Group, Region) %>%
      summarise(meancFos = mean(cFos_st),
                secFos = sd(cFos_st)/sqrt(n()-1)) %>%
      mutate(Group = as.factor(Group)) %>%
      mutate(Group = relevel(Group, ref = "SED"))
    
    sumPFC <- PFC %>% 
      group_by(Group, Region) %>% 
      summarise(meancFos = mean(cFos_st),
                secFos = sd(cFos_st)/sqrt(n()-1)) %>% 
      mutate(Group = as.factor(Group)) %>% 
      mutate(Group = relevel(Group, ref = "SED"))
    
    #  Hippocampus Plot ----
    PlotHipp_DGcomb <- ggplot(sumHipp_DGcomb, aes(x = Region, y = meancFos, fill = Group)) +
      geom_bar(stat = 'identity', position = 'dodge', width = .7, colour = "black") + 
      geom_errorbar(aes(ymin = meancFos - secFos,
                        ymax = meancFos + secFos),
                    width = .2, position = position_dodge(.7)) +
      geom_jitter(data = plotdata_DGcomb_d4, aes(x = Region, y = cFos_st, color = Group),
                  size = 1.2, alpha = .7, fill = "black", shape = 21,
                  position = position_jitterdodge(jitter.width = .1, jitter.height = 0)) +
      labs(title = "", y = bquote(paste("c-Fos"^{"+"} ~ "&" ~ "NeuN"^{"+"} ~ "/" ~ "NeuN"^{"+"} ~ " (# /" ~ mm^3 ~")"))) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 60000), breaks = seq(0, 60000, by = 10000)) + 
      scale_x_discrete(limits = RegionList_DGcomb) +
      scale_color_manual(values = c(SED = "black", LIE = "black")) +
      scale_fill_manual(values = c(SED = "grey85", LIE = "skyblue")) +
      theme_classic(base_family = "TNR") +
      theme(
        legend.position = c(.9, 1),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10)) +
      #dCA3
      annotate("path", x = c(1.8,1.8,2.2,2.2), y = c(18000,20000,20000,15000)) +
      annotate("text", x = 2, y = 21000, label = "*", size = 10) +
      #vDG
      annotate("path", x = c(4.8,4.8,5.2,5.2), y = c(20000,22000,22000,17000)) +
      annotate("text", x = 5, y = 23000, label = "*", size = 10)
    
    PlotHipp_DGcomb
    
    ggsave("c-Fos\\Plot_Hippocampus.png", width = 7, height = 3, dpi = 300)
    
    # PFC Plot ----
    PlotPFC <- ggplot(sumPFC, aes(x = Region, y = meancFos, fill = Group)) +
      geom_bar(stat = 'identity', position = 'dodge', width = .7, colour = "black") + 
      geom_errorbar(aes(ymin = meancFos - secFos,
                        ymax = meancFos + secFos),
                    width = .2, position = position_dodge(.7)) +
      geom_jitter(data = PFC, aes(x = Region, y = cFos_st, color = Group),
                  size = 1.2, alpha = .7, fill = "black", shape = 21,
                  position = position_jitterdodge(jitter.width = .1, jitter.height = 0)) +
      labs(title = "", y = bquote(paste("c-Fos"^{"+"} ~ "/ area" ~ "(# /" ~ mm^3 ~")"))) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 80000), breaks = seq(0, 80000, by = 20000)) + 
      scale_x_discrete(limits = PFCList) +
      scale_color_manual(values = c(SED = "black", LIE = "black")) +
      scale_fill_manual(values = c(SED = "grey85", LIE = "skyblue")) +
      theme_classic(base_family = "TNR") +
      theme(
        legend.position = c(.9, 1),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10)) 
    
    PlotPFC
    ggsave("c-Fos\\Plot_PFC.png", width = 3, height = 3, dpi = 300)
    
    
    # Correlation ----
    # Create Datasets ----
    Freezingdata <- dataEx1 %>%
      group_by(No, Group) %>%
      summarise(Freezing = mean(Freezing))
    
    Freezing <-
      Freezingdata[order(Freezingdata$Group),] %>%
      dplyr::select(3)
    
    Hippocampus_cor <- plotdata_DGcomb_d4 %>% 
      spread(key = Region, value = cFos_st)
    
    PFC_cor <- PFC %>% 
      spread(key = Region, value = cFos_st)
      
    cordata <- bind_cols(Hippocampus_cor, PFC_cor, Freezing) %>%
      # Check columns and rows
      dplyr::select("No...1","No...10","No...14","Group...2", "Group...11",RegionListAll_DGcomb, "Freezing") %>% 
      dplyr::select(-c(2:4)) %>%
      rename("No" = "No...1",
             "Group" = "Group...11") %>% 
      arrange(Group)
    
   
    png("c-Fos/Chart_Correlation.png", width = 9*300, height = 9*300, res  = 300)
    chart.Correlation(cordata[3:12])
    dev.off()
     
    # ScatterPlot----
    
    dir.create("c-Fos/ScatterResults", showWarnings = F)

    for(i in RegionListAll_DGcomb){
      
      plotdata <- cordata %>% 
        dplyr::select("Group", i, "Freezing") %>% 
        rename("PlotX" = i)
      
      get_p_value <- function(dataset, x_var, y_var) {
        test <- cor.test(dataset[[x_var]], dataset[[y_var]])
        return(test$p.value)
      }
      
      g <- ggplot(plotdata, aes(x = PlotX, y = Freezing, colour = Group, fill = Group)) +
        geom_point(size = 4, alpha = .7, shape = 21) +
        if (i == "PL" || i == "IL"){
          labs(x = bquote(paste("c-Fos"^{"+"} ~ "/ area" ~ "(# /" ~ mm^3 ~")")),
               y = "Freezing Time (%)")
        }
        else {
          labs(x = bquote(paste("c-Fos"^{"+"} ~ "&" ~ "NeuN"^{"+"} ~ "/" ~ "NeuN"^{"+"} ~ " (# /" ~ mm^3 ~")")), 
               y = "Freezing Time (%)")
        }
      g <- g +
        scale_color_manual(values = c("SED" = "grey85", "LIE" = "skyblue")) +
        scale_fill_manual(values = c("SED" = "grey85", "LIE" = "skyblue")) +
        # theme_classic(base_family = "TNR") +
        theme_classic(base_family = "MEI") +
        theme(plot.title = element_text(size = 14, hjust = .5),
              # legend.position = "none",
              legend.key = element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              axis.text.x = element_text(size = 8, colour = "black"),
              axis.text.y = element_text(size = 10, colour = "black"),
              axis.line = element_line(colour = "black"),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_text(size = 10))
      
      if(get_p_value(plotdata, x_var = "PlotX", y_var = "Freezing") < 0.05){
        modliner <- lm(Freezing ~ PlotX, plotdata)
        predict <- g + geom_abline(intercept = modliner$coefficients[1], slope = modliner$coefficients[2],size = 1)
      ggsave(paste0("c-Fos/ScatterResults/Lined_", i, ".png"), width = 4, height = 3, dpi = 300)
      }
      else{
        lmgroup <- g
      assign(paste0("Plot", i), lmgroup, envir = .GlobalEnv)
      ggsave(paste0("c-Fos/ScatterResults/", i, ".png"), width = 4, height = 3, dpi = 300)
      }
    }
    
    # dCA3 Plot
    mod_PlotdCA3 <- Plot_predictdCA3 +
      scale_color_manual(values = c("SED" = "black", "LIE" = "black")) +
      scale_y_continuous(limits = c(59, 100), breaks = seq(60, 100, by = 10)) +
      scale_x_continuous(limits = c(2500, 15000), breaks = seq(3000,15000, by = 3000)) +
      theme(axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 8))
    mod_PlotdCA3
    ggsave("c-Fos/Mod_dCA3_Ex1.png", width = 4, height = 3, dpi = 300)
    
    # Analyse ==========================================
    sink("c-Fos/Immunohistochemistry_Analyse.txt", split = T)
    
    #PFC
    cat("\n== PL ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(cordata$PL)
    cat("\n-- Levene's Test --\n")
    leveneTest(cordata$PL, cordata$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(cordata$PL ~ cordata$Group, var.equal = TRUE)
    
    cat("\n== IL ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(cordata$IL)
    cat("\n-- Levene's Test --\n")
    leveneTest(cordata$IL, cordata$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(cordata$IL ~ cordata$Group, var.equal = TRUE)
    
    #Hippocampus
    
    cat("\n== dDG ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(cordata$dDG)
    cat("\n-- Levene's Test --\n")
    leveneTest(cordata$dDG, cordata$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(cordata$dDG ~ cordata$Group, var.equal = TRUE)
    
    cat("\n== vDG ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(cordata$vDG)
    cat("\n-- Levene's Test --\n")
    leveneTest(cordata$vDG, cordata$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(cordata$vDG ~ cordata$Group, var.equal = TRUE)
    
    cat("\n== dCA3 ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(cordata$dCA3)
    cat("\n-- Levene's Test --\n")
    leveneTest(cordata$dCA3, cordata$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(cordata$dCA3 ~ cordata$Group, var.equal = TRUE)
    
    cat("\n== dCA2 ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(cordata$dCA2)
    cat("\n-- Levene's Test --\n")
    leveneTest(cordata$dCA2, cordata$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(cordata$dCA2 ~ cordata$Group, var.equal = TRUE)
    
    cat("\n== dCA1 ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(cordata$dCA1)
    cat("\n-- Levene's Test --\n")
    leveneTest(cordata$dCA1, cordata$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(cordata$dCA1 ~ cordata$Group, var.equal = TRUE)
    
    cat("\n== vCA3 ==\n")
    cat("\n-- Shapiro-Wilk Test--\n")
    shapiro.test(cordata$vCA3)
    cat("\n-- Levene's Test --\n")
    leveneTest(cordata$vCA3, cordata$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(cordata$vCA3 ~ cordata$Group, var.equal = TRUE)
    
    cat("\n== vCA1 ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(cordata$vCA1)
    cat("\n-- Levene's Test --\n")
    leveneTest(cordata$vCA1, cordata$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(cordata$vCA1 ~ cordata$Group, var.equal = TRUE)
    
    cat("\n== Correlation ==\n")
    cat("\n-- dCA3 & Freezing Time --\n")
    cor.test(cordata$dCA3, cordata$Freezing)
  
    sink()
    