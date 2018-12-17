#Load data and remove CPR-number, birth date and hospital identification
filename <- "GeneratedData/All_data_LYFO.RData"
if(file.exists(filename)){
  load(filename)
}else{
  #Set working directory and load the datasets from the drives
  setwd("../StatCure/ExternalData/")
  LYFO <- read.csv2("F_2017_05_05_lyfodel.csv", header = T, stringsAsFactors = F)
  
  #Format birth date and diagnosis and create an age variable. Remove birth day and hospital
  LYFO$Fødselsdato <- as.Date(LYFO$Fødselsdato, format = "%d-%m-%Y")
  LYFO$diag_date <- as.Date(LYFO$Dato.for.diagnostisk.biopsi, format = "%d-%m-%Y")
  LYFO$age <- as.numeric(LYFO$diag_date - LYFO$Fødselsdato)
  LYFO <- subset(LYFO, select = -c(Fødselsdato, UnitCode))
  
  #Remove decrypted files from folder and save temporary file in project folder
  file.remove("F_2017_05_05_amlformat.csv", "F_2017_05_05_DCCGdel.csv", 
              "F_2017_05_05_DCCG_ekstra.csv", "F_2017_05_05_lyfodel.csv")
  setwd(project)
  save(LYFO, file = filename)
}

#LYFO <- LYFO[LYFO$diag_date > as.Date("2006-01-01"),]
LYFO$Histology <- LYFO$WHO.histologikode
LYFO$diag_date <- as.Date(LYFO$Dato.for.diagnostisk.biopsi, format = "%d-%m-%Y")
LYFO$age_years <- LYFO$age / ayear
LYFO$sex <- factor(LYFO[, "Køn"], levels = c("Mand", "Kvinde"), labels = c("male", "female"))
LYFO$death_date <- LYFO[, "CPR.reg..dødsdato"]
LYFO$death_date[LYFO$death_date == "."] <- LYFO[LYFO$death_date == ".", "Dødsdato"]
LYFO$death_date <- as.Date(LYFO$death_date, format = "%d-%m-%Y")
LYFO$look_up <- as.Date(LYFO[, "Opslagsdato.for.KM.dage..død."], format = "%d-%m-%Y")
LYFO$status <- as.numeric(!is.na(LYFO$death_date))
wh <- is.na(LYFO$death_date)
LYFO$death_date[wh] <- LYFO$look_up[wh]

LYFO$FU <- as.numeric(LYFO$death_date - LYFO$diag_date)
LYFO <- LYFO[!is.na(LYFO$FU), ]
LYFO <- LYFO[LYFO$FU > 0, ]
LYFO <- LYFO[LYFO$age_years >= 18,]

LYFO2 <- LYFO[, c("diag_date", "age", "age_years", 
                  "sex", "FU", "status", "Histology")]
LYFO2$FU_years <- LYFO2$FU / ayear
#LYFO2$gender <- as.numeric(LYFO2$sex) - 1

DLBCL <- LYFO2[LYFO2$Histology == "9680",]
FL <- LYFO2[LYFO2$Histology %in% c("9690", "9691", "9695", "9698"),]
ML <- LYFO2[LYFO2$Histology == "9673",]
#BL <- LYFO2[LYFO2$Histology == "9687",]
#PTCL <- LYFO2[LYFO2$Histology == "9702",]
#HL <- LYFO2[LYFO2$Histology %in% c("9659", "9650", "9663", "9651", "9652", "9653"),]

DLBCL$exp_haz <- general.haz(time = "FU", age = "age", sex = "sex", year = "diag_date", 
                             data = DLBCL, ratetable = survexp.dk)
FL$exp_haz <- general.haz(time = "FU", age = "age", sex = "sex", year = "diag_date", 
                             data = FL, ratetable = survexp.dk)
ML$exp_haz <- general.haz(time = "FU", age = "age", sex = "sex", year = "diag_date", 
                             data = ML, ratetable = survexp.dk)


##Danish cancer registry data
filename <- "GeneratedData/OldCancers.RData"
if(file.exists(filename)){
  load(filename)
}else{
  CR_person <- read.table("ExternalData/t_person.asc",
                          header = T, sep = ";", stringsAsFactors = F)
  CR_tumor <- read.table("ExternalData/t_tumor.asc",
                          header = T, sep = ";", stringsAsFactors = F)
  save(CR_tumor, CR_person, file = filename)
}


# dupps <- which(duplicated(CR_person$v_pnr_enc))
# 
# res <- rep(NA, length(dupps))
# for(i in 1:length(dupps)){
#   a <- CR_person[CR_person$v_pnr_enc == CR_person$v_pnr_enc[dupps[i]],]
#   res[i] <- all(a[1,] == a[2,]) 
# }
# table(res)
# 
# CR_person <- CR_person[!duplicated(CR_person$v_pnr_enc),]
# rownames(CR_person) <- CR_person$v_pnr_enc
# 
# wh <- which(is.na(CR_tumor$C_STATUS) & CR_tumor$k_cprnr_enc %in% CR_person$v_pnr_enc)
# 
# a <- CR_person[CR_tumor$k_cprnr_enc[wh],]


table(CR_tumor$C_DIAGGR)
CR_tumor$disease <- factor(CR_tumor$C_DIAGGR, levels = c(22, 36, 49, 51), 
                           labels = c("Colon cancer", "Breast cancer",
                                      "Bladder cancer", "Melanoma"))
CR_tumor$diag_date <- as.Date(CR_tumor$D_DIAGNOSEDATO)
CR_tumor$D_STATDATO <- as.Date(CR_tumor$D_STATDATO)
CR_tumor$D_FDSDATO <- as.Date(CR_tumor$D_FDSDATO)
CR_tumor$age <- as.numeric(CR_tumor$diag_date - CR_tumor$D_FDSDATO) 
CR_tumor$age_years <- CR_tumor$age / ayear
CR_tumor$sex <- factor(CR_tumor$C_SEX, levels = c(1, 2), labels = c("male", "female"))

#Dont really know what these three types mean, so I'm removing them.
wh <- which(CR_tumor$C_STATUS %in% c(20, 50, 60))
CR_tumor <- CR_tumor[-wh,]

wh <- which(is.na(CR_tumor$D_STATDATO))
CR_tumor$D_STATDATO[wh] <- "2016-12-31"
CR_tumor$C_STATUS[wh] <- 0
table(CR_tumor$C_STATUS)
CR_tumor$status <- ifelse(CR_tumor$C_STATUS == 90, 1, 0)

#Calculate follow-up
CR_tumor$FU <- as.numeric(CR_tumor$D_STATDATO - CR_tumor$diag_date) / ayear
wh <- which(CR_tumor$FU <= 0)
CR_tumor <- CR_tumor[-wh,]

CR_tumor$FU_days <- CR_tumor$FU * ayear
CR_tumor$age_group <- cut(CR_tumor$age_years, breaks = c(0, 60, 70, 80, 200), 
                          labels = c("50-59", "60-69", "70-79", "80+"))


#There are some patients who are more than 108 years of age and still alive today. It seems like they are missing
#information on death since the top oldest people are all censored. I think there might be something wrong here. So
#I am removing these patients.
age_death <- CR_tumor$age_years + CR_tumor$FU
wh <- which(age_death > 108)
sort(age_death[wh])
CR_tumor <- CR_tumor[-wh,]


