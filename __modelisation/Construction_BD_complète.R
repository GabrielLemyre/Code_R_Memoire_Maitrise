# ————————————————————————————————————————————————————————————————————————————————————————————————————
# Hidden Markov Models
# IMPORTATION DES DONNÉES - modélisation
# ————————————————————————————————————————————————————————————————————————————————————————————————————
# written
# Gabriel LEMYRE
# ————————————————————————————————————————————————————————————————————————————————————————————————————
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ————————————————————————————————————————————————————————————————————————————————————————————————————
# First version : 21 décembre, 2020
# Last version : 21 décembre, 2020
# ————————————————————————————————————————————————————————————————————————————————————————————————————

Construction_BD_complete <- function(path.DATA){
  
  # Date pour initialisé la base de données
  oldest.time <- "1919-01-01"
  
  "Matlab2Rdate()"

  # Fonction pour message d'erreur si problème d'importation
  errorImport <- function(error.message,
                          index){
    cat("\n")
    herror()
    cat("Erreur lors de l'importation de",index,":\n")
    hline()
    cat(gsub(": ",":\n    ",as.character(error.message)))
    herror()
    cat("\n")
  }
  # ————————————————————————————————————————————————————————————————————————————————————
  # ////////////////////////////////////////////////////////////////////
  # PRÉPARATION DES DONNÉES - TESTS NORMALITÉ - FAITS STYLISÉS
  # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  # ————————————————————————————————————————————————————————————————————————————————————
  
  # ———————————————————————————————————————————————————————
  # Importation des données ----
  # ———————————————————————————————————————————————————————--
  
  # Table pour recevoir données QUOTIDIENNES
  DATA.daily <- data.frame(Date = seq(as.Date(oldest.time),Sys.Date(),1))
  DATA.daily$Date <- as.Date(DATA.daily$Date, "%Y/%m/%d")
  
  # Table pour recevoir données HEBDOMADAIRES
  DATA.weekly <- data.frame(Date = seq(as.Date(oldest.time),Sys.Date(),1))
  DATA.weekly$Date <- as.Date(DATA.weekly$Date, "%Y/%m/%d")
  
  # Table pour recevoir données MENSUELLES
  DATA.monthly <- data.frame(Date = seq(as.Date(oldest.time),Sys.Date(),1))
  DATA.monthly$Date <- as.Date(DATA.monthly$Date, "%Y/%m/%d")
  
  # Longueur des tables vides pour recevoir les données
  full.length <- length(DATA.daily$Date)
  
  # Provenance des données (Langage de programmation [R ou Matlab])
  data.Origin <- "Matlab"
  
  # import des données
  
  # S&P500 -- DAILY
  # https://finance.yahoo.com/quote/%5EGSPC?p=^GSPC&.tsrc=fin-srch
  tryCatch({
    SP.500 <- read.csv(paste(data.path,"GSPC.csv",sep="/")) # MATLAB
    DATA.daily$SP.500 <- SP.500 %>% fit.to.Dates(keep.name = "Adj.Close", 
                                                 fullset = DATA.daily$Date)
  },
  error = function (err){
    errorImport(err,"SP.500")
  })
  
  
  # S&P/TSX (S&P/TSX Composite index (^GSPTSE)) -- DAILY
  # https://finance.yahoo.com/quote/%5EGSPTSE?p=%5EGSPTSE
  # tryCatch({
  #   SP.TSX <- read.csv(paste(data.path,"GSPTSE.csv",sep="/")) # MATLAB
  #   DATA.daily$SP.TSX <- SP.TSX %>% fit.to.Dates(keep.name = "Adj.Close",
  #                                                fullset = DATA.daily$Date)
  # },
  # error = function (err){
  #   errorImport(err,"SP.TSX")
  # })
  
  
  # VIX (CBOE Volatility Index: VIX) -- DAILY
  # https://fred.stlouisfed.org/series/VIXCLS
  # tryCatch({
  #   VIX <- read.csv(paste(data.path,"VIXCLS.csv",sep="/")) # MATLAB
  #   DATA.daily$VIX <- VIX %>% fit.to.Dates(keep.name = "VIXCLS",
  #                                          date.name = "DATE",
  #                                          fullset = DATA.daily$Date)
  # },
  # error = function (err){
  #   errorImport(err,"VIX")
  # })
  
  
  # TYX (Treasury Yield 30 Years) -- DAILY
  # https://finance.yahoo.com/quote/%5ETYX?p=^TYX&.tsrc=fin-srch
  # tryCatch({
  #   TYX <- read.csv(paste(data.path,"TYX.csv",sep="/")) # MATLAB
  #   DATA.daily$TYX <- TYX %>% fit.to.Dates(keep.name = "Adj.Close",
  #                                          fullset = DATA.daily$Date)
  # },
  # error = function (err){
  #   errorImport(err,"TYX")
  # })
  
  
  # NASDAQCOM -- DAILY
  # https://fred.stlouisfed.org/series/NASDAQCOM
  # tryCatch({
  #   NASDAQCOM <- read.csv(paste(data.path,"NASDAQCOM.csv",sep="/")) # MATLAB
  #   DATA.daily$NASDAQCOM <- NASDAQCOM %>% fit.to.Dates(keep.name = "NASDAQCOM",
  #                                                      date.name = "DATE",
  #                                                      fullset = DATA.daily$Date)
  # },
  # error = function (err){
  #   errorImport(err,"NASDAQCOM")
  # })
  
  
  # CAD-USD taux change -- DAILY
  # https://finance.yahoo.com/quote/CAD%3DX?p=CAD%3DX
  tryCatch({
    CADUSD <- read.csv(paste(data.path,"CADUSD.csv",sep="/")) # MATLAB
    DATA.daily$CADUSD <- CADUSD %>% fit.to.Dates(keep.name = "Adj.Close",
                                                 fullset = DATA.daily$Date)
  },
  error = function (err){
    errorImport(err,"CADUSD")
  })
  
  
  # HRX.TO -- WEEKLY
  # https://finance.yahoo.com/quote/HRX.TO?p=HRX.TO&.tsrc=fin-srch
  # tryCatch({
  #   HRX.TO <- read.csv(paste(data.path,"HRX.TO.csv",sep="/")) # MATLAB
  #   DATA.weekly$HRX.TO <- HRX.TO %>% fit.to.Dates(keep.name = "Adj.Close",
  #                                                 fullset = DATA.weekly$Date,
  #                                                 methode.donnee.manquante = "replicate")
  # },
  # error = function (err){
  #   errorImport(err,"HRX.TO")
  # })
  
  # INDPRO -- MONTHLY
  # https://fred.stlouisfed.org/series/INDPRO
  tryCatch({
    INDPRO <- read.csv(paste(data.path,"INDPRO.csv",sep="/")) # MATLAB
    
    INDPRO.log <- data.frame(Date=INDPRO$DATE[-1],
                             INDPRO.log=(log(INDPRO$INDPRO[-1]/INDPRO$INDPRO[-dim(INDPRO)[1]])))
    
    DATA.monthly$INDPRO <- INDPRO %>% fit.to.Dates(keep.name = "INDPRO",
                                                   date.name = "DATE",
                                                   fullset = DATA.monthly$Date,
                                                   methode.donnee.manquante = "replicate")
    
    DATA.monthly$INDPRO.log <- INDPRO.log %>% fit.to.Dates(keep.name = "INDPRO.log",
                                                           fullset = DATA.monthly$Date,
                                                           methode.donnee.manquante = "replicate")
  },
  error = function (err){
    errorImport(err,"INDPRO")
  })
  
  # ———————————————————————————————————————————————————————
  # Préparation des données ----
  # ———————————————————————————————————————————————————————
  # Binding the data together before transform and trim
  BD.complete <- cbind(if(dim(DATA.daily)[2]>0){DATA.daily},
                       if(dim(DATA.weekly)[2]>0){DATA.weekly[-1]},
                       if(dim(DATA.monthly)[2]>0){DATA.monthly[-1]}
  )
  
  names(BD.complete) <- c(if(dim(DATA.daily)[2]>0){names(DATA.daily)},
                          if(dim(DATA.weekly)[2]>0){names(DATA.weekly)[-1]},
                          if(dim(DATA.monthly)[2]>0){names(DATA.monthly)[-1]})
  
  # Impression des entêtes afin vérifier les données
  if (dim(BD.complete)[2]>0){
    hline()
    print(head(BD.complete))
  }
  
  return(BD.complete)
}