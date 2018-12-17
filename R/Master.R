
#Set project
project <- "K:/FORSK-Projekt/Projekter/Scientific Projects/110_PhD_Lasse_2015/Projekter/LossOfLifetimeEstimation/"
setwd(project)

#Load required packages
pkgs <- c("survival", "relsurv", "ggplot2", "numDeriv", "mstate", "flexsurv", "xtable", 
          "rootSolve", "quadprog", "matrixcalc", "grid", "gridExtra", "matrixStats", 
          "VGAM", "githubinstall", "rstpm2", "foreign", "readstata13", "cuRe")

for(i in 1:length(pkgs)){
  if(!pkgs[i] %in% rownames(installed.packages())){
    if(pkgs[i] == "rstpm2"){
      gh_install_packages(pkgs[i], ref = "develop")
    }else{
      install.packages(pkgs[i]) 
    }
  }
  library(pkgs[i], character.only = T)
}


##Github commit
#Update commit number
#commit a21daa052c12a3505ab87b27ecc4976d010ddbb6

#Create directories
thesis <- FALSE
if(thesis){
  fig.out <- "C:/Users/sw1y/Dropbox/Apps/ShareLaTeX/Thesis/papers/paperB/Output/Figures/"
  tab.out <- "C:/Users/sw1y/Dropbox/Apps/ShareLaTeX/Thesis/papers/paperB/Output/Tables/"
  data.out <- "GeneratedData/"
  dir.create(fig.out, showWarnings = F, recursive = T)
  dir.create(tab.out, showWarnings = F, recursive = T)   
} else {
  fig.out <- "C:/Users/sw1y/Dropbox/Apps/ShareLaTeX/Loss of lifetime estimation/Output/Figures2/"
  tab.out <- "C:/Users/sw1y/Dropbox/Apps/ShareLaTeX/Loss of lifetime estimation/Output/Tables2/"
  data.out <- "GeneratedData/"
  dir.create(fig.out, showWarnings = F, recursive = T)
  dir.create(tab.out, showWarnings = F, recursive = T) 
}


#Global project settings
mai_par <- par("mai")
mfrow_par <- par("mfrow")
format <- "%.1f"
ayear <- 365.24

#Load Danish lifetime
survexp.dk <- transrate.hmd("../Data/DanishLifeTable/mltper_1x1.txt", 
                            "../Data/DanishLifeTable/fltper_1x1.txt")

#Colour palette for models
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Import functions from the cuRe package
source("R/LoadData2.R", encoding = "utf-8")

