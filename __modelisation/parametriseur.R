# ————————————————————————————————————————————————————————————————————————————————————————————————————
# Hidden Markov Models
# PARAMÉTRISE LES MODÈLES
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

# ————————————————————————————————————————————————————————
# ENREGISTRE TOUS LES MODÈLES À ENTRAINER DANS UNE LISTE
# ————————————————————————————————————————————————————————
# Retourne l'objet de modèle général
# ————————————————————————————————————————————————————————

# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# CONSTRUIT LA LISTE DE SPÉCIFICATIONS À ENTRAÎNER AVEC TOUS LEURS PARAMÈTRES ----
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ————————————————————————————————————————————————————————————————————————————————————

Build.liste.modeles <- function(liste,
                                default.test.Transformations = FALSE,
                                default.transformations = "none",        # Pas de transformation
                                default.scalingFactor=1,                 # Pas de scaling
                                default.freq.str = "Q",                  # Quotidien
                                default.type = "HMM",
                                default.Transition.Type = "Homogeneous", # Homogène
                                default.initial.Distribution= NULL,
                                default.tau = 5,
                                admissible.types = c("UNIVARIEE","MELANGE","HMM","HHMM","HHMM.simplifie","DDMS"),
                                default.type.a.rouler = "ALL"){
  
  # PRÉPARATION general.pars
  load.model <- liste$general.pars$load.model
  test.Transformations <- swap.null(liste$general.pars$test.Transformations, default.test.Transformations)
  type.a.rouler <- swap.null(liste$general.pars$type.a.rouler, default.type.a.rouler)
  
  # PRÉPARATION data.pars
  index <- liste$data.pars$index
  variables.explicatives <- liste$data.pars$variables.explicatives
  date.boundaries <- liste$data.pars$date.boundaries
  
  transformations <- swap.null(liste$data.pars$transformations,rep(default.transformations,length(index)+length(variables.explicatives)))
  scalingFactor <- swap.null(liste$data.pars$scalingFactor,rep(default.scalingFactor,length(index)+length(variables.explicatives)))
  freq.str <- swap.null(liste$data.pars$freq.str,default.freq.str)
  
  # Ajout des data.pars
  liste.modeles <- list(general.pars = list(
    load.model = load.model,
    test.Transformations = test.Transformations,
    type.a.rouler = type.a.rouler
  ),
  data.pars = list(
    index = index,
    variables.explicatives = variables.explicatives[!variables.explicatives %in% index],
    transformations = transformations,
    scalingFactor = scalingFactor,
    date.boundaries = date.boundaries, # Partagé par toutes les séries
    freq.str=freq.str)
  )
  if (type.a.rouler[1] == "ALL"){
    type.a.rouler <- names(liste)[!names(liste) %in% c("general.pars","data.pars")]
  }
  # print(liste.modeles)
  # 
  # cat(paste(names(liste)[!names(liste) %in% c("general.pars","data.pars")],sep=" "),"\n")
  
  for (item in Reduce(intersect, list(type.a.rouler,names(liste)[!names(liste) %in% c("general.pars","data.pars")]))){
    if (!item %in% admissible.types){
      type <- swap.null(liste[[item]]$type,
                        default.type)
    } else {
      type <- item
    }
    
    Transition.Type <- swap.null(liste[[item]]$Transition.Type,
                                 switch(item,
                                        "UNIVARIEE" = "Random",
                                        "MELANGE" = "Random",
                                        "HMM" = default.Transition.Type,
                                        NULL))
    
    index.name <- gsub(".","",index,fixed = TRUE)
    if (!is.null(Transition.Type)){
      item.name <- paste(paste(index.name,type,sep="_"),
                         Transition.Type, sep=".")
    } else {
      item.name <- paste(index.name,type, sep="_")
    }
    
    
    nbRegimes <- swap.null(liste[[item]]$K,1)
    
    for (k in 1:length(nbRegimes)){
      temp.liste.modeles <- list()
      K <- as.character(nbRegimes[k])
      
      mu0    <- liste$general.pars$mu0[[K]]
      sigma0 <- liste$general.pars$sigma0[[K]]
      gamma0 <- liste[[item]]$gamma0
      
      # Ajout du nombre de paramètres au nom de l'item
      item.name.K <- paste0(item.name,".K.",K)
      
      initial.Distribution <- liste[[item]]$initial.Distribution[[K]]
      
      tau <- swap.null(liste[[item]]$tau, 
                       default.tau)
      
      
      
      basic.list <- list(type=type,
                         Transition.Type = Transition.Type,
                         mu0 = mu0,
                         sigma0 = sigma0,
                         distribution=NULL)
      
      if (item == "DDMS"){
        nb.tau <- length(tau)
        for (i in 1:nb.tau){
          temp.liste.modeles[[paste0(item.name.K,".tau.",tau[i])]] <- 
            append(
              basic.list,
              list(
                tau = tau[i],
                gamma0 = gamma0[[as.character(tau[i])]]
              )
            )
        }
      } else if (item == "LOGIT"){
        variables.explicatives <- liste$data.pars$variables.explicatives
        nb.var.exp <- length(variables.explicatives)
        for (i in 1:nb.var.exp){
          
          variable.explicative <- 
            if (!variables.explicatives[i] %in% index){
              variables.explicatives[i]
            }
          temp.liste.modeles[[paste0(item.name.K,".",
                                     gsub(".","",
                                          variables.explicatives[i],
                                          fixed = TRUE))]] <- 
            append(
              basic.list,
              list(
                variable.explicative = variable.explicative,
                gamma0 = gamma0[[variables.explicatives[i]]]
              )
            )
        }
      } else {
        temp.liste.modeles[[item.name.K]] <- append(
          basic.list,
          list(gamma0 = gamma0[[K]])
        )
      }
      
      # cat("Liste prélim \n")
      # Hmisc::list.tree(temp.liste.modeles)
      # 
      nb.spec <- length(temp.liste.modeles)
      # cat("nb.spec =",nb.spec,"\n")
      
      for (i in 1:nb.spec){
        if (is.null(initial.Distribution)){
          liste.modeles[[names(temp.liste.modeles)[i]]] <- temp.liste.modeles[[i]]
        } else {
          for (init.dist in initial.Distribution){
            model.name <- paste0(names(temp.liste.modeles)[i],
                                 ".init.(",paste(init.dist, collapse="."),")")
            liste.modeles[[model.name]] <- append(
              temp.liste.modeles[[i]],
              list(initial.Distribution = init.dist)
            )
          } # End for initial.Distribution
        } # End if on initial.Distribution
      } #end for on temp.Modele
    }
  } # End for items dans liste initiale
  
  # cat("Liste finale \n")
  # Hmisc::list.tree(liste.modeles)
  return(liste.modeles)
}


