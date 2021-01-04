# ----------------------------------------------------------------------------------------------------
# Data Preparation
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# First version : february 14th, 2020
# Last version : 5 mars 2020
# ----------------------------------------------------------------------------------------------------

# --------------------------------------------------------
# Préparation de l'échantillon
# --------------------------------------------------------
# Fonction manipulant les données
#     L'idée est de prendre l'échantillon et de préparer
#     les données afin d'être utilisées directement dans
#     les modèles
#   Cette fonction retourne l'échantillon avec et sans 
#   transformations, borné aux dates[start.date, end.date]
# --------------------------------------------------------

data.prep <- function(data.raw, data.Origin="R", 
                      data.names=NULL, 
                      start.date=NULL,
                      end.date=NULL,
                      transformations="none", # "troncature"  "logR" "none
                      freq.str=1, # TODO : ajout d'option de compression des données
                      scalingFactor=1,
                      date.col.name="Date",
                      train.test.ratio=1){ # train.test.ratio = 1 => juste un taining set, train.test.ratio=0.6 => 60% train 40% test
  
  # Sélection des colonnes spécifiées par le vecteur 'data.names'
  #   et de la colonne des dates
  if (is.null(data.names)){
    data.names <- names(data.raw)[-1]
  }
  data.select <- data.raw[c(date.col.name,data.names)]
  
  # Retrait des lignes incomplètes (contenant un NA) et enregistrement
  #   dans un data.frame
  data.noNa <- data.frame(data.select[complete.cases(data.select),])
  names(data.noNa) <- names(data.select)
  
  if (data.Origin=="matlab"){ # Si les données viennent de Matlab
    # transformations de la variable date au format approprié
    data.noNa$Date <- as.Date(as.character(Matlab2Rdate(data.noNa$Date)))
  }
  
  hline()
  
  # Date de début effective
  indice.date.debut.USER <- which(data.noNa$Date==start.date)
  if (length(indice.date.debut.USER)==0){
    # Si la première date donnée n'est pas dans l'échantillon, on utilise la première date de l'échantillon sans NAs
    st <- which.min(abs(data.noNa$Date-start.date))
    first.date <- data.noNa$Date[st]
    cat("La date de début donnée n'est pas dans l'échantillon.",
        "\nNous prennons la date la plus près dans l'échantillon. \nDate début effective :",
        weekdays(first.date),as.character(first.date),
        "[",as.numeric(rownames(data.noNa))[st],"]\n# --------------------------------------------------------\n")
  } else {
    st <- indice.date.debut.USER
    first.date <- data.noNa$Date[st]
    cat("Date début :",weekdays(first.date),as.character(first.date),
        "[",as.numeric(rownames(data.noNa))[st],"]\n# --------------------------------------------------------\n")
  }
  
  # Date de fin effective
  indice.date.fin.USER <- which(data.noNa$Date==end.date)
  if (length(indice.date.fin.USER)==0){
    # Si la dernière date donnée n'est pas dans l'échantillon, on utilise la dernière date de l'échantillon sans NAs
    en <- which.min(abs(data.noNa$Date-end.date))
    last.date <- data.noNa$Date[en]
    cat("La date de fin donnée n'est pas dans l'échantillon.",
        "\nNous prennons la date la plus près dans l'échantillon. \nDate fin effective :",
        weekdays(last.date),as.character(last.date),
        "[",as.numeric(rownames(data.noNa))[en],"]\n# --------------------------------------------------------\n")
  } else {
    en <- indice.date.fin.USER
    last.date <- data.noNa$Date[en]
    cat("Date fin :",weekdays(last.date),as.character(last.date),
        "[",as.numeric(rownames(data.noNa))[en],"]\n# --------------------------------------------------------\n")
  }
  
  # Construction de la table propre, sans NA et avec les bonnes dates et de type 'numeric'
  numeric.clean.data <- sapply(data.noNa[st:en,data.names], function(x) {as.numeric(as.character(x))})
  clean.data <- data.frame(numeric.clean.data,row.names=data.noNa$Date[st:en])
  # Donne les bons noms aux colonnes de la table de données
  names(clean.data) <- data.names
  
  # Obtention de la longueur de l'échantillon
  n <- dim(clean.data)[1]
  
  # Compression des données, obtention des valeurs Quotidienne, Hebdomadaire ou mensuelle
  frequency.nb = switch(freq.str, # freq.str contient une lettre correspondant au choix de l'usager
                        "Q" = 1,  # Quotidien
                        "H" = 5,  # Hebdomadaire
                        "M" = 20) # Mensuel
  sequence.index <- rownames(clean.data)[seq(from=frequency.nb, to=n,by=frequency.nb)] # Sequence des indices retenus
  clean.data <- subset(clean.data, rownames(clean.data) %in% sequence.index) # Selection de ces indices dans la BD
  
  # Obtention de la nouvelle longueur de l'échantillon
  n <- dim(clean.data)[1]
  
  # Si pas assez de transformations
  if (length(transformations)==1){
    transformations <- rep(transformations,length(data.names))
  } else {
    if (length(transformations) != length(data.names)){
      stop("Pas assez de transformations pour le nombre de variables.")
    }
  }
  
  DATA <- list()
  for (variable in data.names){
    no.variable <- which(data.names==variable)
    if (transformations[no.variable]=="logR"){
      logR <- log(clean.data[[variable]][-1]/clean.data[[variable]][-n])
      # Multipliction des log-rendements par un facteur 'mult'
      DATA[[variable]] <- logR * scalingFactor[no.variable]
      
    } else if (transformations[no.variable]=="scale") {
      DATA[[variable]] <- clean.data[[variable]][-1] * scalingFactor[no.variable]
      
    } else if (transformations[no.variable]=="diff") {
      diffs <- (clean.data[[variable]][-1]-clean.data[[variable]][-n])/clean.data[[variable]][-n]
      # Multipliction des log-rendements par un facteur 'mult'
      DATA[[variable]] <- diffs * scalingFactor[no.variable]
      
    } else if (transformations[which(data.names==variable)]=="none") {
      DATA[[variable]] <- clean.data[[variable]][-1]
    }
    
  }
  
  
  
  # print(rownames(clean.data)[-1])
  # Sauvegarde de l'échantillon obtenu dans un data.frame
  DATA <- data.frame(DATA, row.names=rownames(clean.data)[-1])
  
  names(DATA) <- data.names
  
  print(head(DATA))
  
  # Retourne les données transformées 'DATA'
  #   et la version avant la 'transformations' sous 'DATA.init'
  nb.train <- ceiling(train.test.ratio*n)
  
  RBIND <- function(datalist) {
    require(plyr)
    temp <- rbind.fill(datalist)
    rownames(temp) <- unlist(lapply(datalist, row.names))
    temp
  }
  
  if (train.test.ratio<1){
    
    train.DATA <- data.frame(DATA[c(1:nb.train),], row.names=rownames(DATA)[c(1:nb.train)])
    names(train.DATA) <- names(data.select)[-1]
    train.DATES <- rownames(train.DATA)
    
    test.DATA <- data.frame(DATA[-c(1:nb.train),], row.names=rownames(DATA)[-c(1:nb.train)])
    names(test.DATA) <- names(data.select)[-1]
    test.DATES <- rownames(test.DATA)
    
    full.DATA <- RBIND(list(train.DATA,test.DATA))
    names(full.DATA) <- names(data.select)[-1]
    full.DATES <- rownames(full.DATA)
    
  } else {
    train.DATA <- as.data.frame(DATA, row.names=rownames(DATA)[c(1:nb.train)])
    names(train.DATA) <- names(data.select)[-1]
    train.DATES <- rownames(train.DATA)
    
    test.DATA <- NULL
    test.DATES <- NULL
    
    full.DATA <- train.DATA
    names(full.DATA) <- names(data.select)[-1]
    full.DATES <- train.DATES
  }
  
  
  # Nombre de données d'entrainement
  nb.train <- if(!is.null(dim(train.DATA)[1])){dim(train.DATA)[1]}else{length(train.DATA)}
  # nb.train.exp.var <- if(!is.null(dim(train.DATA.expl)[1])){dim(train.DATA.expl)[1]}else{length(train.DATA.expl)}
  
  # Nombre de données pour les tests
  nb.test <- if(!is.null(dim(test.DATA)[1])){dim(test.DATA)[1]}else{length(test.DATA)}
  # nb.test.exp.var <- if(!is.null(dim(test.DATA.expl)[1])){dim(test.DATA.expl)[1]}else{length(test.DATA.expl)}
  
  
  # ———————————————————————————————————————————————————————
  # IMPRESSION DES ENTÊTES DES BDs RÉSULTANTES ----
  # ———————————————————————————————————————————————————————
  # ————————————————————————————
  # TRAIN
  # ————————————————————————————
  if (nb.train>0){
    cat("\n")
    hline()
    cat("Entête - Base de donnée pour entraînement :\n")
    hline()
    print(head(train.DATA))
    hline()
  }
  # if (nb.train.exp.var>0){
  #   hline()
  #   cat("Entête - Variable explicative pour entraînement :\n")
  #   hline()
  #   exp.var.w.Dates <- as.data.frame(cbind(train.DATES,as.numeric(as.character(train.DATA.expl))))
  #   names(exp.var.w.Dates) <- c("Date",names(bound.DATA[variable.explicative]))
  #   print(head(exp.var.w.Dates))
  #   hline()
  # }
  
  # ————————————————————————————
  # TEST
  # ————————————————————————————
  if (nb.test>0){
    cat("\n")
    hline()
    cat("Entête - Base de donnée pour tests 'out of sample' :\n")
    hline()
    print(head(test.DATA))
    hline()
  }
  # if (nb.test.exp.var>0){
  #   hline()
  #   cat("Entête - Variable explicative pour 'out of sample' :\n")
  #   hline()
  #   exp.var.w.Dates <- as.data.frame(cbind(test.DATES,as.numeric(as.character(test.DATA.expl))))
  #   names(exp.var.w.Dates) <- c("Date",names(bound.DATA[variable.explicative]))
  #   print(head(exp.var.w.Dates))
  #   hline()
  # }
  
  
  # ———————————————————————————————————————————————————————
  # RETOURNE LES BDs RÉSULTANTES ----
  # ———————————————————————————————————————————————————————
  return(
    list(
      train.DATA=train.DATA,
      test.DATA=test.DATA,
      train.DATES=train.DATES,
      test.DATES=test.DATES,
      full.DATA=full.DATA,
      full.DATES=full.DATES,
      DATA.init=clean.data,
      real.start=first.date,
      real.end=last.date,
      n.IS=nb.train,
      n.OofS=nb.test,
      n.total= max(nb.train,nb.train + nb.test)
    )
  )
}

