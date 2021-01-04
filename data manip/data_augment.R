# ----------------------------------------------------------------------------------------------------
# Data Augmentation
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# First version : 9 mars 2020
# Last version : 9 mars 2020
# ----------------------------------------------------------------------------------------------------

# --------------------------------------------------------
# Augmentation de l'échantillon
# --------------------------------------------------------
# Fonction manipulant les données
#   L'idée est de prendre un échantillon et d'y ajouter
#     les données les plus à jour.
#   Cette fonction retourne l'échantillon sans 
#   transformations
# --------------------------------------------------------

data.augment <- function(Data.Old,
                         index.old,
                         Data.New,
                         index.new,
                         data.Origin.old="Matlab",
                         data.Origin.new="R",
                         old.date.col.name="Date",
                         new.date.col.name="Date") {
  
  # Sélection des colonnes spécifiées par le vecteur 'index'
  #   et de la colonne des dates OLD
  Data.Old.df <- data.frame(Data.Old[c(old.date.col.name,index.old)])
  names(Data.Old.df) <- c("date",index.old)
  
  
  # Sélection des colonnes spécifiées par le vecteur 'index'
  #   et de la colonne des dates NEW
  Data.New.df <- data.frame(Data.New[new.date.col.name],Data.New[index.new])
  names(Data.New.df) <- c("date",index.new)
  print(head(Data.New.df))
  
  if (data.Origin.old=="Matlab"){ # Si les données viennent de Matlab
    # Transformation de la variable date au format approprié
    Data.Old.df$date <- as.Date(as.character(Matlab2Rdate(Data.Old.df$date)))
  } else {
    Data.Old.df$date <- as.Date(as.character(Data.Old.df$date))
  }
  
  if (data.Origin.new=="Matlab"){ # Si les données viennent de Matlab
    # Transformation de la variable date au format approprié
    Data.New.df$date <- as.Date(as.character(Matlab2Rdate(Data.New.df$date)))
  } else {
    Data.New.df$date <- as.Date(as.character(Data.New.df$date))
  }
  
  # Nous nous assurons que les observations soient numériques
  Data.New.df[index.new] <- sapply(Data.New.df[index.new], function(x) {as.numeric(as.character(x))})
  
  # Longueurs de l'ancien et du nouveau data.frame
  len.old <- length(Data.Old.df$date)
  len.new <- length(Data.New.df$date)
  
  # Quelle est la dernière date de l'ancien DF
  last.date <- Data.Old.df$date[len.old]
  
  # Cela correspond à quelle entrée du nouveau DF
  id.missing.date <- which(Data.New.df$date==last.date)+1
  
  # Range où ajouter les nouvelles données
  new.range <- c(id.missing.date:len.new)
  
  # Extraction des nouvelles données
  missing.data <- Data.New.df[new.range,]
  
  # Ajout de lignes vides
  Data.Old.df[new.range,] <- NA
  
  # Ajout des nouvelles données et des dates correspondantes
  Data.Old.df$date[new.range] <- missing.data$date
  Data.Old.df[new.range,index.new] <- c(missing.data[,index.new])
  
  # print(Data.Old.df[new.range-3,])
  
  return(Data=Data.Old.df)
}