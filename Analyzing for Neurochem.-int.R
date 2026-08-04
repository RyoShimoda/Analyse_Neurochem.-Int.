# Analyse for Neurochemistry International

# Library ----
library(car)
library(psych)
library(ppcor)
library(PerformanceAnalytics)
library(readr)

# My Function & ANOVAKUN----
source("AUTOGRAPH.txt")
source("anovakun_489.txt")
source("IHCC_cFos.txt")

# Contextual Fear Conditioning ----
  # Create Datasets & Plots ----
    mf_autograph(folder_name = "Sampling_after_extinction_day1", 
                     experiment_type = "FC",
                     path = ".\\Sampling_after_extinction_day1\\FC", 
                     plot_number = "_FC")

    mf_FC_annotation(folder_name = "Sampling_after_extinction_day1_demo", 
                          graph = GFC_FC, 
                          y_first_arrow = 25, y_second_arrow = 45, y_third_arrow3 = 65, 
                          legend_posision_y = 65)
    
  # Statistical Analysis ----
    sink(file = "Sampling_after_extinction_day1/Result/ANOVA_FC.txt", split = T)
    
    anovakun(dataset = sumFC[-1],
             design = "AsB",
             Group = c("SED", "LIE"),
             Time = TimeFC,
             hf = T,
             peta = T)
    sink()
    
# Contextual Fear Extinction ----
  # Create Datasets & Plots ----
    mf_autograph(folder_name = "Sampling_after_extinction_day1_demo", 
                     path = ".\\Sampling_after_extinction_day1_demo\\Ex1", 
                     experiment_type = "Ex1")
    
  # Statistical Analysis ---- 
    AnodataEx1 <- dataEx1 %>% 
      group_by(No, Group) %>% 
      summarise(Freezing = mean(Freezing))
    
    sink(file = "Sampling_after_extinction_day1/Result/ANOVA_Ex1.txt", split = T)
    
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
      mf_sampling(folder_name = "Sampling_after_extinction_day1_demo", 
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
    
    tiff(filename = "Sampling_after_extinction_day1/Result/Plot/Soleus.tiff", 
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
    
    tiff(filename = "Sampling_after_extinction_day1/Result/Plot/Plantaris.tiff", 
         width = 4 * 900, height = 3.5 * 900, res = 900)
    modPlantaris
    dev.off()
    
  # Statistical Analysis ----
    anomuscle <- Sampdata %>% 
      mutate("Soleus" = Soleus/Body_weight * 100) %>% 
      mutate("Plantaris" = Plantaris/Body_weight * 100) %>% 
      mutate("Adrenal" = Adrenal/Body_weight * 100) %>% 
      mutate("Thymus" = Thymus/Body_weight * 100)
    
    sink(file = "Sampling_after_extinction_day1/Result/Sampling_ANOVA.txt", split = T)

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
    #Set Working Directory -----
    setwd("D:\\Soya_lab\\Experimental_Raw_data")
    
    # Create Lists----
    RegionList <- c("dDGsp", "dDGip", "dCA3", "dCA2", "dCA1", "vDGsp", "vDGip", "vCA3", "vCA1")
    PFCList <- c("PL", "IL")
    RegionListAll <- c("PL", "IL", "dDGsp", "dDGip", "dCA3", "dCA2", "dCA1", "vDGsp", "vDGip", "vCA3", "vCA1")
    
    # Use My Function----
    mf_IHCC_cFos(Folder_name = "Immunohistochemistry_demo\\c-Fos\\Hippocampus",
                 Path = "D:\\Soya_lab\\Experimental_Raw_data\\Immunohistochemistry_demo\\c-Fos\\Hippocampus\\Analyse_G50_R20",
                 Region_list = RegionList, Red_name = "c-Fos", Green_name = "NeuN", 
                 Result_folder = "Result_Analyse_G50_R20")
    
    mf_IHCC_cFos(Folder_name = "Immunohistochemistry_demo\\c-Fos\\PL_IL",
                 Path = "D:\\Soya_lab\\Experimental_Raw_data\\Immunohistochemistry_demo\\c-Fos\\PL_IL\\G40_R20",
                 Region_list = PFCList, Red_name = "c-Fos", Green_name = "NeuN", PFC = "", number_of_region = 2, D_V = "", 
                 Result_folder = "Result_G40_R20", save_plot_width_all = 3, save_plot_width_i = 6, save_plot_height_i = 6,
                 title_size_i = 8, axis_text_size_i = 8, axis_title_size_i = 8, save_histplot_width = 7.5, save_histplot_height = 3,
                 legend_position = c(.8, .9))
    
    # Read KatiKati Data ----
    kachidata_d4 <- read.csv("Immunohistochemistry_demo\\c-Fos\\Hippocampus\\ForAnalyse.csv", header = F) %>% 
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
      summarise(cFos = sum(cFos))
    
    Kachidata_D4 <- kachidata_d4[order(kachidata_d4$Region),]
    
    # Read GatheringData ----
    gatheringdata <- read.csv("Immunohistochemistry_demo\\c-Fos\\Hippocampus\\Result_Analyse_G50_R20\\GatheringData.csv")
    Gatheringdata <- gatheringdata[order(gatheringdata$Region),]
    
    # Create Plot Data ----
    plotdata_d4 <- bind_cols(Gatheringdata, Kachidata_D4) %>% 
      dplyr::select(14, 3, 4, 6, 17) %>%
      rename("Group" = "Group...3") %>%
      rename("No" = "No...14") %>% 
      rename("Region" = "Region...4") %>% 
      mutate(cfos = cFos/Area_G * 1000000 / 0.04) %>%
      mutate(Group = as.factor(Group)) %>%
      mutate(Group = relevel(Group, ref = "SED")) %>%
      dplyr::select(1,2,3,6)
    
    individualdata_d4 <- plotdata_d4[order(plotdata_d4$No),] %>% 
      mutate(ID = rep(1:16, each = 9))
    
    # Create Region Data ----
    for (i in RegionList) {
      assign(paste0("Data_", i), plotdata_d4[plotdata_d4$Region == i,], envir = .GlobalEnv )
    }
    
    # Bind Gatharing Data & KatiKatiData ----
    summarisedata_d4 <- bind_cols(Gatheringdata, Kachidata_D4) %>%
      dplyr::select(14, 3, 4, 6, 17) %>%
      rename("Group" = "Group...3") %>%
      rename("No" = "No...14") %>% 
      rename("Region" = "Region...4") %>% 
      mutate(cfos = cFos/Area_G * 1000000 / 0.04) %>%
      group_by(Group, Region) %>%
      summarise(mean_cfos = mean(cfos),
                se_cfos = sd(cfos)/sqrt(n()-1)) %>%
      mutate(Group = as.factor(Group)) %>%
      mutate(Group = relevel(Group, ref = "SED"))
    
    # Create DG combind Data ----
    DGdata_d4 <-bind_cols(Gatheringdata, Kachidata_D4) %>% 
      dplyr::select(14, 3, 4, 6, 17) %>%
      rename("Group" = "Group...3") %>% 
      rename("No" = "No...14") %>% 
      rename("Region" = "Region...4") %>% 
      filter(str_detect(Region, pattern = "DG")) %>% 
      mutate(DV = if_else(str_detect(Region, pattern = "d"), "Dorsal", "Ventral")) %>% 
      group_by(No, Group, DV) %>% 
      summarise(Area_G = sum(Area_G),
                cFos = sum(cFos)) %>% 
      mutate(Region = if_else(str_detect(DV, pattern = "Dorsal"), "dDG", "vDG")) %>% 
      dplyr::select(1,2,6,4,5)
    
    # Create DG Plot Data ----
    plotdata_DGcomb_d4 <- bind_cols(Gatheringdata, Kachidata_D4) %>% 
      dplyr::select(14, 3, 4, 6, 17) %>% 
      rename("Group" = "Group...3") %>%
      rename("No" = "No...14") %>% 
      rename("Region" = "Region...4") %>% 
      filter(str_detect(Region, pattern = "CA")) %>% 
      bind_rows(DGdata_d4) %>% 
      mutate(cfos = cFos/Area_G * 1000000 / 0.04) %>% 
      mutate(Group = as.factor(Group)) %>% 
      mutate(Group = relevel(Group, ref = "SED")) %>% 
      dplyr::select(1,2,3,6)
    
    
    # Create DG combined List
    RegionListAll_DGcomb <- c("PL", "IL", "dDG", "dCA3", "dCA2", "dCA1", "vDG", "vCA3", "vCA1")
    RegionList_DGcomb <- c("dDG", "dCA3", "dCA2", "dCA1", "vDG", "vCA3", "vCA1")
    
    # Calculate mm^2 -> mm^3 ----
    DatadDG <- DGdata_d4 %>% 
      mutate(Count = cFos/Area_G * 1000000 / 0.04) %>% 
      filter(Region == "dDG")
    
    DatavDG <- DGdata_d4 %>% 
      mutate(Count = cFos/Area_G * 1000000 / 0.04) %>% 
      filter(Region == "vDG")
    
    # Read Hippocampus PL/IL Data ----
    Hippocampus_DGcomb <- plotdata_DGcomb_d4 %>% 
      dplyr::select(-1) %>% 
      rename("Count" = "cfos")
    
    PFC <- read.csv("Immunohistochemistry\\c-Fos\\PL_IL\\Result_LG10_R20\\GatheringData.csv") %>%
      dplyr::select(2, 3, 4, 11) %>% 
      mutate(Count = Count/0.04)
    
    # Create Plot Data ----
    PlotAllData_DGcomb <- bind_rows(Hippocampus_DGcomb, PFC) %>% 
      mutate(Group = as.factor(Group)) %>% 
      mutate(Group = relevel(Group, ref = "SED"))
    
    PlotSum_DGcomb <- PlotAllData_DGcomb %>% 
      group_by(Group, Region) %>%
      summarise(meanCount = mean(Count),
                seCount = sd(Count)/sqrt(n()-1)) %>%
      mutate(Group = as.factor(Group)) %>%
      mutate(Group = relevel(Group, ref = "SED"))
    
    # Create Whole Hippocampus Data ----
    wholesum <- dvData_d4 %>% 
      group_by(Group) %>% 
      summarise(mean = mean(cfos), 
                se = sd(cfos)/sqrt(n()-1))
    
    wholejitter <-dvData_d4 %>% 
      group_by(No, Group) %>% 
      summarise(mean = mean(cfos))
    
    # Create Whole PL/IL Data ----
    wholePFCsum <- PFC %>% 
      group_by(Group) %>% 
      summarise(mean = mean(Count), 
                se = sd(Count)/sqrt(n()-1)) %>% 
      mutate(Group = as.factor(Group)) %>% 
      mutate(Group = relevel(Group, ref = "SED"))
    
    wholePFCjitter <- PFC %>% 
      group_by(No, Group) %>% 
      summarise(mean = mean(Count)) %>% 
      mutate(Group = as.factor(Group)) %>% 
      mutate(Group = relevel(Group, ref = "SED"))
    
    # ALL Region Plot ----
    PlotAll_DGcomb <- ggplot(PlotSum_DGcomb, aes(x = Region, y = meanCount, fill = Group)) +
      geom_bar(stat = 'identity', position = 'dodge', width = .7, colour = "black") + 
      geom_errorbar(aes(ymin = meanCount - seCount,
                        ymax = meanCount + seCount),
                    width = .2, position = position_dodge(.7)) +
      geom_jitter(data = PlotAllData_DGcomb, aes(x = Region, y = Count, color = Group),
                  size = 1.2, alpha = .7, fill = "black", shape = 21,
                  position = position_jitterdodge(jitter.width = .1, jitter.height = 0)) +
      labs(title = "", y = expression(paste("c-Fos+ / NeuN+ (cell / ", {mm^3},")"))) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 60000), breaks = c(0, 10000, 20000, 30000, 40000, 50000, 60000)) + 
      scale_x_discrete(limits = RegionListAll_DGcomb) +
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
      annotate("path", x = c(3.8,3.8,4.2,4.2), y = c(18000,20000,20000,15000)) +
      annotate("text", x = 4, y = 21000, label = "*", size = 10) +
      #vDG
      annotate("path", x = c(6.8,6.8,7.2,7.2), y = c(20000,22000,22000,17000)) +
      annotate("text", x = 7, y = 23000, label = "*", size = 10)
    
    PlotAll_DGcomb  
    ggsave("Immunohistochemistry\\c-Fos\\PlotAll_ver2_DGcomb.png", width = 9, height = 3, dpi = 300)
    
    # Correlation ----
    # Create Datasets ----
    Freezingdata <- dataEx1 %>%
      group_by(No, Group) %>%
      summarise(Freezing = mean(Freezing))
    
    Freezing <-
      Freezingdata[order(Freezingdata$Group),] %>%
      dplyr::select(3)
    
    Hippocampus <- plotdata_DGcomb_d4 %>% 
      spread(key = Region, value = cfos)
    
    PFC_cor <- PFC %>% 
      spread(key = Region, value = Count)
      
    cordata <- bind_cols(Hippocampus, PFC_cor, Freezing) %>%
      dplyr::select(2, 13, 12, 6, 3, 4, 5, 9, 7, 8, 15) %>%
      rename("Group" = "Group...2")
    
    chart.Correlation(cordata[2:11])
    
    #ScatterPlot----

    for(i in RegionList_DGcomb){
      
      plotdata <- cordata %>% 
        dplyr::select("Group", i, "Freezing") %>% 
        rename("PlotX" = i)
      
      assign(paste0("Data_",i), plotdata, envir = .GlobalEnv)
      
      g <- ggplot(plotdata, aes(x = PlotX, y = Freezing, colour = Group, fill = Group)) +
        geom_point(size = 4, alpha = .7, shape = 21) +
        labs(x = expression(paste("c-Fos+ / NeuN+ (cell / ", {mm^3},")")), y = "Freezing Time (%)") +
        scale_color_manual(values = c("SED" = "grey85", "LIE" = "skyblue")) +
        scale_fill_manual(values = c("SED" = "grey85", "LIE" = "skyblue")) +
        # theme_classic(base_family = "TNR") +
        theme_classic(base_family = "MEI") +
        theme(plot.title = element_text(size = 14, hjust = .5),
              # legend.position = "none",
              legend.key = element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              axis.text.x = element_text(size = 10, colour = "black"),
              axis.text.y = element_text(size = 10, colour = "black"),
              axis.line = element_line(colour = "black"),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_text(size = 10))
      
      lmgroup <- g
      assign(paste0("Plot", i), lmgroup, envir = .GlobalEnv)
      ggsave(paste0("Immunohistochemistry_demo/c-Fos/ScatterResult/", i, ".png"), width = 3.5, height = 3, dpi = 300)
      
      modliner <- lm(Freezing ~ PlotX, plotdata)
      predict <- g + geom_abline(intercept = modliner$coefficients[1], slope = modliner$coefficients[2],size = 1)
      assign(paste0("Plot_predict", i), predict, envir = .GlobalEnv)
      ggsave(paste0("Immunohistochemistry_demo/c-Fos/ScatterResult/All_", i, ".png"), width = 3.5, height = 3, dpi = 300)
      
    }
    
    mod_PlotdCA3 <- Plot_predictdCA3 +
      scale_color_manual(values = c("SED" = "black", "LIE" = "black")) +
      scale_y_continuous(limits = c(59, 100), breaks = seq(60, 100, by = 10)) +
      scale_x_continuous(limits = c(2500, 15000), breaks = seq(3000,15000, by = 3000)) +
      theme(axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 8))
    mod_PlotdCA3
    ggsave("Immunohistochemistry_demo/c-Fos/ScatterResult/dCA3_Ex1.png", width = 4, height = 3, dpi = 300)
    
    #Analyse ==========================================
    # For correlation plot
    cordataLIE <- cordata %>% 
      filter(Group == "LIE") %>% 
      dplyr::select(-1)
    
    cordataSED <- cordata %>% 
      filter(Group == "SED")%>% 
      dplyr::select(-1)
    
    sink("Immunohistochemistry/c-Fos/Immunohistochemistry_Analyse.txt", split = T)
    
    #PFC
    cat("\n== PL ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(Data_PL$Count)
    cat("\n-- Levene's Test --\n")
    leveneTest(Data_PL$Count, Data_PL$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(Data_PL$Count ~ Data_PL$Group, var.equal = TRUE)
    
    cat("\n== IL ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(Data_IL$Count)
    cat("\n-- Levene's Test --\n")
    leveneTest(Data_IL$Count, Data_IL$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(Data_IL$Count ~ Data_IL$Group, var.equal = TRUE)
    
    #Hippocampus
    
    cat("\n== dDG ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(DatadDG$Count)
    cat("\n-- Levene's Test --\n")
    leveneTest(DatadDG$Count, DatadDG$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(DatadDG$Count ~ DatadDG$Group, var.equal = TRUE)
    
    cat("\n== vDG ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(DatavDG$Count)
    cat("\n-- Levene's Test --\n")
    leveneTest(DatavDG$Count, DatavDG$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(DatavDG$Count ~ DatavDG$Group, var.equal = TRUE)
    
    cat("\n== dCA3 ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(Data_dCA3$cfos)
    cat("\n-- Levene's Test --\n")
    leveneTest(Data_dCA3$cfos, Data_dCA3$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(Data_dCA3$cfos ~ Data_dCA3$Group, var.equal = TRUE)
    
    cat("\n== dCA2 ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(Data_dCA3$cfos)
    cat("\n-- Levene's Test --\n")
    leveneTest(Data_dCA2$cfos, Data_dCA2$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(Data_dCA2$cfos ~ Data_dCA2$Group, var.equal = TRUE)
    
    cat("\n== dCA1 ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(Data_dCA1$cfos)
    cat("\n-- Levene's Test --\n")
    leveneTest(Data_dCA1$cfos, Data_dCA1$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(Data_dCA1$cfos ~ Data_dCA1$Group, var.equal = TRUE)
    
    cat("\n== vCA3 ==\n")
    cat("\n-- Shapiro-Wilk Test--\n")
    shapiro.test(Data_vCA3$cfos)
    cat("\n-- Levene's Test --\n")
    leveneTest(Data_vCA3$cfos, Data_vCA3$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(Data_vCA3$cfos ~ Data_vCA3$Group, var.equal = TRUE)
    
    cat("\n== vCA1 ==\n")
    cat("\n-- Shapiro-Wilk Test --\n")
    shapiro.test(Data_vCA1$cfos)
    cat("\n-- Levene's Test --\n")
    leveneTest(Data_vCA1$cfos, Data_vCA1$Group)
    cat("\n-- unpaired t Test --\n")
    t.test(Data_vCA1$cfos ~ Data_vCA1$Group, var.equal = TRUE)
    
    cat("\n== Correlation ==\n")
    cat("\n-- dCA3 & Freezing Time --\n")
    cor.test(cordata$dCA3, cordata$Freezing)
  
    sink()
    