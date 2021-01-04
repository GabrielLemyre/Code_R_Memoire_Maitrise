# ————————————————————————————————————————————————————————————————————————————————————————————————————
# Hidden Markov Models
# MAIN - fichier d'entrainement et de traitement des résultats
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

if(tryCatch(!outside.run,error = function(err){return(TRUE)})){
  # Closing all connections
  closeAllConnections()
  # Vide la cache de R
  rm(list=ls())
  # Retrait de tous les graphiques ouvert
  graphics.off()
  
  outside.run <- FALSE
  ignored.files.OUTSIDE.RUN <- c()
} else {
  ignored.files.OUTSIDE.RUN <- "Parametrisations"
}

cat("outside.run =",outside.run,"\n")

# ————————————————————————————————————————————————————————
# ROULE TOUS LES MODÈLES SUR TOUS LES ÉCHANTILLONS
# ————————————————————————————————————————————————————————
# Fonction entraînant puis sauvegardant les modèles
# ————————————————————————————————————————————————————————


# ———————————————————————————————————————————————————————
# Changement du document de travail ----
# ———————————————————————————————————————————————————————
path <- '~/Documents/GitHub/Memoire_Maitrise/code/r'
setwd(path.expand(path)) # Setting Sourcing path

# profvis({
# ———————————————————————————————————————————————————————
# Wrapping functions ----
# ———————————————————————————————————————————————————————
# IMPORTATION DE TOUTES LES FONCTIONS EXTERNES
#   Contenues dans le document "path"
source("wrapper.R") 

# ———————————————————————————————————————————————————————
# Paramètres généraux ----
# ———————————————————————————————————————————————————————
# Heure et date de lancement de ce document
run.time.start <- cat.Sys.Time()

# Largeur de la console exortée
options("width"=85)

round.to <- 1
round.to.Pred.err <- 8
set.seed(5)

print.result.optim <- FALSE

send.email.when.done <- FALSE

lag.max.acf <- 200
n.sim <- 50

print.init = FALSE

train.test.prop <- 0.90

# Couleur à utiliser pour les graphiques
colors <- list(
  "2" = c("darkblue","red"),
  "3" = c("darkblue","gold1","red"),
  "4" = c("darkblue","green4","gold1","red")
)

# ———————————————————————————————————————————————————————
# Chemin pour exporter les résulats ----
# ———————————————————————————————————————————————————————
output.path <- '~/Documents/GitHub/Memoire_Maitrise'
dir.create(output.path, showWarnings = FALSE)

# ———————————————————————————————————————————————————————
# Chemin vers les donnnées ----
# ———————————————————————————————————————————————————————
data.path <- "/Users/gabriellemyre/Documents/GitHub/Memoire_Maitrise/DATA"

# ———————————————————————————————————————————————————————
# Chemin pour exporter les résultats ----
# ———————————————————————————————————————————————————————
# Chemin GÉNÉRAL vers les résultats
path.resultats <- file.path(paste(output.path,"__resultats",sep="/",collapse = ""))
dir.create(path.resultats, showWarnings = FALSE)

# ——————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# IMPORTATION DE TOUTES LES DONNÉES ----
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ——————————————————————————————————————————————————————————————————————————————
BD.complete <- Construction_BD_complete(data.path)

# ——————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# INFORMATIONS POUR ENTRAINEMENT ----
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ——————————————————————————————————————————————————————————————————————————————


if(tryCatch(outside.run,
            error = function(err){return(FALSE)})){
  liste.index <- list(
    run.model = run.model
  )
} else {
  liste.index <- list(
    SP.500 = SP.500
  )
  
}


# ——————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# ROUTINE COMPLÈTE SUR TOUS LES INDICES ET LEUR SPÉCIFICATIONS ----
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ——————————————————————————————————————————————————————————————————————————————

for (index.obj.prebuild in liste.index){ # Boucle sur tous les élément de la liste d'objet de modèle
  # Heure et date de lancement de ce document
  index.time.start <- cat.Sys.Time()
  
  # ——————————————————————————————————————————————————————————————————————————————
  # ////////////////////////////////////////////////////////////////////
  # PRÉPARATIONS GÉNÉRALES ———— assignation des variables générale et obtention des infos ----
  # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  # ——————————————————————————————————————————————————————————————————————————————
  # Obtention des variables d'information reçues avec l'"index.obj"
  info.names <- c("data.pars","general.pars")
  nb.info <- length(info.names)
  
  # Ajout des paramétrisations
  index.obj <- Build.liste.modeles(index.obj.prebuild)
  
  # Extraction du nom de l'indice à l'étude
  index.name <- as.character(index.obj$data.pars$index)
  
  # Liste des modèles à entraîner
  liste.modeles <- names(index.obj)[!names(index.obj) %in% info.names]
  
  index.obj$liste.modeles.converge <- c()
  
  # Give the index a name
  index.obj$name <- paste0(gsub(".","",index.obj$data.pars$index,fixed=TRUE),
                           " - (",index.obj$data.pars$freq.str,") - [",index.obj$data.pars$date.boundaries[1],", ",index.obj$data.pars$date.boundaries[2],"] - ",
                           run.time.start)
  kept.index.obj <- index.obj$name
  # ———————————————————————————————————————————————————————
  # Liste qui contiendra les résultats des tests statistiques ----
  # ———————————————————————————————————————————————————————
  
  var.generale <- c("mllk", "AIC.v", "AICc", "BIC.v","MAFE","MSFE","MASFE","MSSFE", "nb.param", "Temps", "Conv")
  in.sample.analysis.names <- c('Modèle - In Sample', var.generale)
  out.of.sample.analysis.names <- c('Modèle - Out of Sample', var.generale)
  retrait.colonnes <- c('Modèle - In Sample', 'Modèle - Out of Sample',"Erreur","Conv")
  
  index.obj$Analyse.InSample <- 
    data.frame(
      matrix(ncol=length(in.sample.analysis.names),
             nrow=0),
      stringsAsFactors = FALSE
    )
  
  index.obj$Analyse.OutofSample <- 
    data.frame(
      matrix(ncol=length(out.of.sample.analysis.names),
             nrow=0),
      stringsAsFactors = FALSE
    )
  
  # ———————————————————————————————————————————————————————
  # DÉFINITION DES CHEMINS POUR EXPORTATIONS ----
  # ———————————————————————————————————————————————————————
  # Chemin pour exporter les INDICES
  path.modele <- 
    file.path(paste(path.resultats,
                    gsub(".","",index.obj$data.pars$index,fixed=TRUE),sep="/",collapse = ""))
  
  dir.create(path.modele, showWarnings = FALSE)
  # Ajout path.modele à l'objet de modèle
  index.obj$path.modele <- path.modele
  
  # Chemin avec dates de l'échnatillon
  path.modele.date <- 
    file.path(paste(path.modele,
                    paste0("[ ",index.obj$data.pars$date.boundaries[1]," ~ ",index.obj$data.pars$date.boundaries[2]," ]"),sep="/",collapse = ""))
  
  dir.create(path.modele.date, showWarnings = FALSE)
  # Ajout path à l'objet de modèle
  index.obj$path.modele.date <- path.modele.date
  
  # Chemin pour exporter les MODÈLES (DATES)
  path.modele.ajd <- 
    file.path(paste(path.modele.date,
                    Sys.Date(),sep="/",collapse = ""))
  
  dir.create(path.modele.ajd, showWarnings = FALSE)
  # Ajout path à l'objet de modèle
  index.obj$path.modele.ajd <- path.modele.ajd
  
  # Chemin pour exporter les GRAPHIQUES (sub DATES)
  path.modele.date.graphs <- 
    file.path(paste(path.modele.ajd,
                    paste0(index.obj$name," - graphiques"),sep="/",collapse = ""))
  
  # Ajout path.modele.date.graphs à l'objet de modèle
  index.obj$path.modele.date.graphs <- path.modele.date.graphs
  
  # ———————————————————————————————————————————————————————
  # CONNECTING TO OUTPUT FILE FOR CONSOLE CONTENT----
  # ———————————————————————————————————————————————————————
  sink.file <- paste0(index.obj$path.modele.ajd,"/",index.obj$name,'.txt')
  print(paste0("Sink file : ",sink.file))
  sink(file=sink.file,append = TRUE)
  hequal()
  print(liste.modeles)
  hequal()
  
  # ———————————————————————————————————————————
  # Le modèle doit-il être importé 
  #   NULL, entraîne le modèle
  #   sinon, le modèle est importé (load.model doit contenir tout le path)
  # ———————————————————————————————————————————
  if (is.null(index.obj$general.pars$load.model)){
    model.loaded = FALSE
    # ——————————————————————————————————————————————————————————————————————————————
    # ////////////////////////////////////////////////////////////////////
    # PRÉPARATION DES DONNÉES ——— SPÉCIFIQUES À CET OBJET DE MODÈLE ---- 
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # ——————————————————————————————————————————————————————————————————————————————
    
    # ———————————————————————————————————————————
    # Préparation données pour modélisation et tests
    # ———————————————————————————————————————————
    
    index.obj$modele.data <- index.obj$data.pars %>% with({
      data.prep(BD.complete,
                data.names = c(index,variables.explicatives[!variables.explicatives %in% index]),
                start.date = date.boundaries[1],
                end.date = date.boundaries[2],
                freq.str = freq.str,
                transformations = transformations,
                scalingFactor = scalingFactor,
                date.col.name = "Date",
                train.test.ratio = train.test.prop
      )
    })
    
    
    
    
    # ——————————————————————————————————————————————————————————————————————————————
    # ////////////////////////////////////////////////////////////////////
    # FULL ROUTINE SUR TOUS LES MODÈLES DE L'"index.obj" EN COURS ----
    # ——————————————————————————————————————————————————————————————————————————————
    for (modele.info in liste.modeles){
      # MESSAGE DÉBUT DE TRAINING
      cat("\n\n")
      hline(10)
      cat("Début - time stamp :",as.character(cat.Sys.Time()),"\n")
      hbegin()
      
      
      # Construction du nom du modèle pour le tableau
      index.obj[[modele.info]]$nom.modele.row <- 
        tryCatch({
          index.obj[[modele.info]] %>% with({
            as.character(
              paste0(
                type,
                if (type %in% c("HMM")){
                  paste0("-",Transition.Type)
                },
                "(K=",length(mu0),
                if(type=="DDMS"){paste0(",tau=",tau)},")",
                if(!is.null(index.obj[[modele.info]]$initial.Distribution)){
                  paste0("-init=(",
                         paste0(as.character(index.obj[[modele.info]]$initial.Distribution),
                                collapse=", "),
                         ")")
                },
                if (!is.null(Transition.Type)){
                  paste0(
                    if (Transition.Type %in% c("UNIVARIEE","GAS")){
                      paste0("-init = (",
                             paste0(s(c((1-gamma0[2])/(2-gamma0[1]-gamma0[2]),
                                        (1-gamma0[1])/(2-gamma0[1]-gamma0[2])),2),collapse=", "),
                             ")")
                    },
                    if (Transition.Type=="LOGIT"){
                      paste0(" sur ",gsub(".","",index.obj$data.pars$index,fixed=TRUE),
                             " $\\sim$ ", 
                             as.character(swap.len0(gsub(".","",index.obj[[modele.info]]$variable.explicative,fixed=TRUE),paste0(index.obj$data.pars$index,"[-1]"))))
                    }
                  )
                } else {""}
              )
            )
          })  
        },
        error = function(error.message){
          err.message.builder(paste0("Erreur lors de la construction du nom \n       du modèle__",modele.info,"__ sur l'indice _",
                                     index.obj$data.pars$index," - (",index.obj$data.pars$freq.str,")_"," :\n"),
                              error.message)
        } )
      
      index.var.comb <- index.obj[[modele.info]]$nom.modele.row
      # index.var.comb <- paste0("__",modele.info,"__ sur l'indice _",
      #                          index.obj$data.pars$index," - (",index.obj$data.pars$freq.str,")_")
      
      
      # Test si le modèle existe
      tryCatch({
        model.spec <- index.obj[[modele.info]]
      },error = function (error.message) {
        err.message.builder(paste0("Erreur lors de l'importation \n       du modèle ",index.var.comb," :\n"),
                            error.message)
      })
      
      tryCatch({
        index.obj %>% with({
          # Test transformation des paramètres
          if (general.pars$test.Transformations){
            N2N.test <- index.obj[[modele.info]] %>% with({
              N2N.test(mu0,
                       sigma0,
                       gamma0,
                       type=type,
                       Transition.Type=Transition.Type)
            })
          }
        })
      },
      error = function (error.message){
        err.message.builder(paste0("Erreur lors du test des transformation des paramètres \n       du modèle ",index.var.comb," :\n"),
                            error.message)
      })
      
      # ————————————————————————————————————————————————————————————
      # ENTRAINEMENT DES MODÈLES ----
      # ————————————————————————————————————————————————————————————
      
      tryCatch({
        index.obj[[modele.info]] <- append(index.obj[[modele.info]],
                                           index.obj[[modele.info]] %>%  with({
                                             train.HMM(
                                               DATA = index.obj$modele.data$train.DATA, # BD préparée complète
                                               data.name = index.obj$data.pars$index,  # Index de l'objet de modèle
                                               exp.var.name = swap.len0(index.obj[[modele.info]]$variable.explicative,index.obj$data.pars$index),
                                               nom.modele.row = index.obj[[modele.info]]$nom.modele.row,
                                               mu0 = index.obj[[modele.info]]$mu0,
                                               sigma0 = index.obj[[modele.info]]$sigma0,
                                               gamma0 = index.obj[[modele.info]]$gamma0,
                                               type = index.obj[[modele.info]]$type,
                                               nbStepsBack = index.obj[[modele.info]]$tau,
                                               distribution = swap.null(index.obj[[modele.info]]$distribution, default="Normal"),
                                               Transition.Type = swap.null(index.obj[[modele.info]]$Transition.Type, default="Homogeneous"),
                                               initial.Distribution = swap.null(index.obj[[modele.info]]$initial.Distribution, default=NULL),
                                               print.result.optim=print.result.optim
                                             )
                                           })
        )
        # Ajout du modèle à la liste puisqu'il a convegé
        index.obj$liste.modeles.converge <- index.obj$liste.modeles.converge %>% append(modele.info)
      },
      error=function(error.message){
        err.message.builder(paste0("Erreur lors de l'optimisation \n       du modèle ",index.var.comb," :\n"),
                            error.message)
        index.obj[[modele.info]]$error <- as.character(error.message)
      },
      warning = function(warning.message){
        warn.message.builder(paste0("Erreur lors de l'optimisation \n       du modèle ",index.var.comb," :\n"),
                             warning.message)
        index.obj[[modele.info]]$warning <- as.character(warning.message)
      }) # FIN ENTRAINEMENT
      # ————————————————————————————————————————————————————————————
      # STATISTIQUES OUT.OF.SAMPLE ----
      # ————————————————————————————————————————————————————————————
      tryCatch({
        # Si le modèle a convergé
        if (!is.null(index.obj[[modele.info]]$HMM.Train)){
          tryCatch({
            
            # ————————————————————————————
            # Selection des données pour les tests
            data.var.test <- 
              as.numeric(unlist(index.obj$modele.data$full.DATA[index.obj$data.pars$index]))
            exp.var.test <- 
              tryCatch({
                if(!is.null(index.obj[[modele.info]]$variable.explicative)){
                  as.numeric(unlist(index.obj$modele.data$full.DATA[index.obj[[modele.info]]$variable.explicative]))
                }else{
                  NULL 
                }
              },
              error = function(err){
                return(NULL)
              })
            # ——————————————
            
            # ————————————————————————————
            # OOfS - Obtention des paramètres de travail
            model.spec <- index.obj[[modele.info]]
            index.obj[[modele.info]]$parvect <- 
              index.obj[[modele.info]] %>%  with({
                normal.HMM.N2W(
                  mu=HMM.Train$mu,
                  sigma=HMM.Train$sigma,
                  matGamma=HMM.Train$matGamma,
                  type=model.spec$type,
                  Transition.Type=model.spec$Transition.Type
                )})
            # ——————————————
            # ——————————————————————————————————————————————————————————————————————
            tryCatch({ # TRY FILTERING, bloquera tous ce qui est analyse OofS sinon
              # ————————————————————————————
              # OOfS - FULL FILTER
              out.of.sample.filtering.and.smoothing <- 
                index.obj[[modele.info]] %>% with({
                  # if (!type %in% c("UNIVARIEE")){
                  normal.HMM.KimFilter(
                    type=type,
                    mu=HMM.Train$mu,
                    sigma=HMM.Train$sigma,
                    matGamma=HMM.Train$matGamma,
                    DATA=data.var.test,
                    exp.var=exp.var.test,
                    nbStepsBack = index.obj[[modele.info]]$tau,
                    distribution = swap.null(index.obj[[modele.info]]$distribution, default="Normal"),
                    Transition.Type = swap.null(index.obj[[modele.info]]$Transition.Type, default="Homogeneous"),
                    initial.Distribution = swap.null(index.obj[[modele.info]]$initial.Distribution, default=NULL)
                  )
                })
              # ——————————————
              
              obs.index.OofS <- c((index.obj$modele.data$n.IS+1):index.obj$modele.data$n.total)
              index.obj[[modele.info]]$HMM.Test$first.OofS.date <- min(as.Date(index.obj$modele.data$test.DATES[1]),
                                                                       as.Date(index.obj$modele.data$real.end))
              
              index.obj[[modele.info]]$HMM.Test$diag.Gamma <- out.of.sample.filtering.and.smoothing$diag.Gamma
              
              # ——————————————
              # NOMBRE DE PARAMÈTRES LIBRES
              index.obj[[modele.info]]$nb.param <- 
                length(index.obj[[modele.info]]$parvect)
              # ——————————————
              
              # ——————————————
              # OOfS - LLK
              index.obj[[modele.info]]$HMM.Test$out.of.sample.LLK <- 
                out.of.sample.filtering.and.smoothing$llk - index.obj[[modele.info]]$HMM.Train$mllk
              # ——————————————
              
              
              # ——————————————
              # OOfS - MAFE
              index.obj[[modele.info]]$HMM.Test$out.of.sample.MAFE <- 
                mean(abs(out.of.sample.filtering.and.smoothing$pred.err.xt[obs.index.OofS]))
              # ——————————————
              
              
              # ——————————————
              # OOfS - MSFE
              index.obj[[modele.info]]$HMM.Test$out.of.sample.MSFE <- 
                mean(
                  (out.of.sample.filtering.and.smoothing$pred.err.xt[obs.index.OofS])^2
                )
              # ——————————————
              
              
              # ——————————————
              # OOfS - MASFE
              index.obj[[modele.info]]$HMM.Test$out.of.sample.MASFE <- 
                mean(abs(out.of.sample.filtering.and.smoothing$standard.pred.err.xt[obs.index.OofS]))
              # ——————————————
              
              
              # ——————————————
              # OOfS - MSSFE
              index.obj[[modele.info]]$HMM.Test$out.of.sample.MSSFE <- 
                mean((out.of.sample.filtering.and.smoothing$standard.pred.err.xt[obs.index.OofS])^2)
              # ——————————————
              
              # ——————————————
              # OOfS - FILTERED PROB
              index.obj[[modele.info]]$HMM.Test$out.of.sample.FILTERED.prob <- 
                out.of.sample.filtering.and.smoothing$p.ct.x1t
              # ——————————————
              
              # ——————————————
              # OOfS - PREDICTIVE PROB
              index.obj[[modele.info]]$HMM.Test$out.of.sample.PREDICTIVE.prob <- 
                out.of.sample.filtering.and.smoothing$p.ct.x1tm1
              # ——————————————
              
              # ——————————————
              # OOfS - PREDICTIVE PROB
              index.obj[[modele.info]]$HMM.Test$out.of.sample.SMOOTH.prob <- 
                out.of.sample.filtering.and.smoothing$p.ct.x1T
              # ——————————————
              
              # ——————————————
              # OOfS - AIC
              index.obj[[modele.info]]$HMM.Test$out.of.sample.AIC <- 
                index.obj[[modele.info]] %>% with({
                  2*(-HMM.Test$out.of.sample.LLK + nb.param)
                })
              # ——————————————
              
              # ——————————————
              # OOfS - AICc
              index.obj[[modele.info]]$HMM.Test$out.of.sample.AICc <- 
                index.obj[[modele.info]] %>% with({
                  HMM.Test$out.of.sample.AIC + (2*nb.param^2 + 2*nb.param)/(index.obj$modele.data$n.OofS - nb.param - 1)
                })
              # ——————————————
              
              # ——————————————
              # OOfS - BIC
              index.obj[[modele.info]]$HMM.Test$out.of.sample.BIC <- 
                index.obj[[modele.info]] %>% with({
                  nb.param*log(index.obj$modele.data$n.OofS) - 2 * HMM.Test$out.of.sample.LLK
                })
              # ——————————————
              
            },
            error=function(error.message){
              err.message.builder(paste0("Erreur lors du filtre pour analyse OofS \n       du modèle ",index.var.comb," :\n"),
                                  error.message)
            })
            
          },
          error=function(err.message){
            err.message.builder(paste0("Erreur intérieur lors du calcul des stats hors de l'échantillon \n       du modèle ",index.var.comb," :\n"),
                                err.message)
          })
        } else {
          hline()
          cat("Le modèle",modele.info,"n'a pas convergé\n")
          hline()
        }
      },
      error=function(err.message){
        err.message.builder(paste0("Erreur extérieur lors du calcul des stats hors de l'échantillon \n       du modèle ",index.var.comb," :\n"),
                            err.message)
        return(liste.index[[index.obj$data.pars$index]] <- index.obj)
        # return(index.obj=index.obj)
      })
      
      hend()
      cat("Fin  ",index.var.comb,"  time stamp :",as.character(cat.Sys.Time()),"\n")
      hline(10)
      cat("\n")
      
    } # Fin de l'entrainement et du calcul des statistiques hors de l'échantillon
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # ——————————————————————————————————————————————————————————————————————————————
  } else {
    # ——————————————————————————————————————————————————————————————————————————————
    # Si le modèle existe déjà, il est importé plutôt qu'entraîné
    # ——————————————————————————————————————————————————————————————————————————————
    index.obj <- readRDS(index.obj$general.pars$load.model)
    model.loaded = TRUE
    new.logit.sp500 <- readRDS("/Users/gabriellemyre/Documents/GitHub/Memoire_Maitrise/__resultats/SP500/[ 1999-12-31 ~ 2020-12-27 ]/2020-12-27/SP500 - (Q) - [1999-12-31, 2020-12-27] - 14h57m21s.rds")
    index.obj$`SP500_HMM.LOGIT.K.4.SP500.init.(0.1.0.0)` = new.logit.sp500$`SP500_HMM.LOGIT.K.4.SP500.init.(0.1.0.0)`
    
    # Ajout HHMM simplifié
    HHMM.simplifie.optim = readRDS("/Users/gabriellemyre/Documents/GitHub/Memoire_Maitrise/__resultats/__FINAUX__/SP500/[ 1999-12-31 ~ 2020-12-26 ]/SP500 - (Q) - [1999-12-31, 2020-12-29] - 17h37m16s.rds")
    index.obj$`SP500_HHMM.simplifie.K.4.init.(0.1.0.0)` = HHMM.simplifie.optim$`SP500_HHMM.simplifie.K.4.init.(0.1.0.0)`
    
    if (index.obj$data.pars$index == "SP.500"){
      index.obj$liste.modeles.converge = c(index.obj$liste.modeles.converge,"SP500_HHMM.simplifie.K.4.init.(0.1.0.0)")
      index.obj$liste.modeles.converge <- index.obj$liste.modeles.converge[c(1:8,15,14,9:13)]
    }
    
    # Mise à jour des paths
    tryCatch({
      index.obj$path.modele <- path.modele
      index.obj$path.modele.date <- path.modele.date
      index.obj$path.modele.ajd <- path.modele.ajd
      index.obj$path.modele.date.graphs <- path.modele.date.graphs
    },
    error = function(err){
      err.message.builder("Problème lors de la mise-à-jour des chemins pour exportation des résultats :\n",
                          err)
    })
  } # Fin boucle d'entraînement sur l'indice
  # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  # ——————————————————————————————————————————————————————————————————————————————
  
  # Si modèle télécharger, réinitialise tableau des résultats
  if (model.loaded ){
    index.obj$Analyse.InSample.prev <- index.obj$Analyse.InSample
    index.obj$Analyse.InSample <- 
      data.frame(
        matrix(ncol=length(in.sample.analysis.names),
               nrow=0),
        stringsAsFactors = FALSE
      )
    
    index.obj$Analyse.OutofSample.prev <- index.obj$Analyse.OutofSample
    index.obj$Analyse.OutofSample <- 
      data.frame(
        matrix(ncol=length(out.of.sample.analysis.names),
               nrow=0),
        stringsAsFactors = FALSE
      )
  }
  
  for (modele.info in index.obj$liste.modeles.converge){
    print(modele.info)
    Analyse.row.IS <- NULL
    Analyse.row.OofS <- NULL
    ex.var.name <- try_default(index.obj[[modele.info]]$variable.explicative, default=paste0(index.obj$data.pars$index,"[-1]"), quiet=TRUE) 
    
    # ROW in sample
    tryCatch({
      Analyse.row.IS <- index.obj[[modele.info]]$HMM.Train %>% with({
        data.frame(matrix(c(index.obj[[modele.info]]$nom.modele.row,
                            s(mllk,round.to), 
                            s(AIC.v,round.to),
                            s(AICc,round.to),
                            s(BIC.v,round.to), 
                            s(MAFE,round.to.Pred.err),
                            s(MSFE,round.to.Pred.err), 
                            s(MASFE,round.to.Pred.err),
                            s(MSSFE,round.to.Pred.err), 
                            index.obj[[modele.info]]$nb.param,
                            index.obj[[modele.info]]$timeSpent.String,
                            codeConv),
                          nrow=1),
                   stringsAsFactors = FALSE) # Time to run Full model
      })
      
      index.obj$Analyse.InSample <- index.obj$Analyse.InSample %>% rbind(Analyse.row.IS)  
    },
    error = function(err){
      err.message.builder(paste("Problème lors de la construction de la ligne Analyse.InSample pour",modele.info," :\n"),
                          err)
    },
    warning = function(warn.mess){
      warn.message.builder(paste("Warning lors de la construction de la ligne Analyse.OutofSample pour",modele.info," :\n"),
                           warn.mess)
    })
    
    # ROW out of sample
    tryCatch({
      Analyse.row.OofS <- index.obj[[modele.info]]$HMM.Test %>% with({
        data.frame(matrix(c(index.obj[[modele.info]]$nom.modele.row,
                            s(out.of.sample.LLK,round.to), 
                            s(out.of.sample.AIC,round.to),
                            s(out.of.sample.AICc,round.to),
                            s(out.of.sample.BIC,round.to), 
                            s(out.of.sample.MAFE,round.to.Pred.err),
                            s(out.of.sample.MSFE,round.to.Pred.err), 
                            s(out.of.sample.MASFE,round.to.Pred.err),
                            s(out.of.sample.MSSFE,round.to.Pred.err), 
                            index.obj[[modele.info]]$nb.param,
                            index.obj[[modele.info]]$timeSpent.String,
                            index.obj[[modele.info]]$HMM.Train$codeConv
        ),nrow=1),
        stringsAsFactors = FALSE) # Time to run Full model
      })
      
      # print(Analyse.row.OofS)
      index.obj$Analyse.OutofSample <- index.obj$Analyse.OutofSample %>% rbind(Analyse.row.OofS)
    },
    error = function(err){
      err.message.builder(paste("Problème lors de la construction de la ligne Analyse.OutofSample pour",modele.info," :\n"),
                          err)
    },
    warning = function(warn.mess){
      warn.message.builder(paste("Warning lors de la construction de la ligne Analyse.OutofSample pour",modele.info," :\n"),
                           warn.mess)
    })
  }
  # ————————————————————————————————————————————————————————————————————————————————————
  # Tente de consolider les résultats et retourne une table au format LaTeX et 
  #   une table de données au format R
  # ————————————————————————————————————————————————————————————————————————————————————
  tryCatch({
    # ———————————————————————————————————————————
    # Consolidation des résultats d'analyse et sauvegarde
    # ———————————————————————————————————————————
    tryCatch({
      colnames(index.obj$Analyse.InSample) <- in.sample.analysis.names
      tryCatch(index.obj$Analyse.InSample$Conv <- as.integer(as.character(index.obj$Analyse.InSample$Conv)),
               error = function(err.message) {
                 err.message.builder(paste0("Erreur lors de la conversion du code de convergence InS \n       du modèle ",paste0("__",modele.info,"__ sur l'indice _",index.obj$data.pars$index,"_")," :\n"),
                                     err.message)
               })
      
      colnames(index.obj$Analyse.OutofSample) <- out.of.sample.analysis.names
      tryCatch(index.obj$Analyse.OutofSample$Conv <- as.integer(index.obj$Analyse.OutofSample$Conv),
               error = function(err.message) {
                 err.message.builder(paste0("Erreur lors de la conversion du code de convergence OofS \n       du modèle ",paste0("__",modele.info,"__ sur l'indice _",index.obj$data.pars$index,"_")," :\n"),
                                     err.message)
               })
    },
    error = function(err){
      err.message.builder(paste("Problème lors de l'assignation des noms de colonnes pour résultats pour",modele.info," :\n"),
                          err)
    })
    # print(index.obj$Analyse.InSample)
    # cat("-----\n")
    # print(index.obj$Analyse.OutofSample)
    
    # ———————————————————————————————————————————
    # In sample analysis EXPORT LaTeX
    # ———————————————————————————————————————————
    # Sélection des variables à inclure dans le tableau
    variables.conservees.IS <- !colnames(index.obj$Analyse.InSample)%in%retrait.colonnes
    Analyse.matrix.in.sample.select <- matrix(unlist(index.obj$Analyse.InSample[,variables.conservees.IS]),
                                              ncol=dim(index.obj$Analyse.InSample[,variables.conservees.IS])[2])
    
    # Obtention de la ligne correspondant à un changement de paramétrisation
    double.hline.index.IS <- NULL
    split.model.names <- strsplit(as.character(index.obj$Analyse.InSample[,1]),
                                  "-",fixed=TRUE)
    for (i in 1:length(split.model.names)){
      current.name <- paste0(split.model.names[[i]][-length(split.model.names[[i]])],
                             collapse="")
      if (i>1){
        if (current.name!= previous.name){
          double.hline.index.IS <- double.hline.index.IS %>% c((i-1))
        }
      }
      previous.name <- current.name
    }
    
    # Obtention des indices correspondants aux max de LLK, AIC, AICc, BIC
    bestLLK.IS  <- which(as.numeric(as.character(index.obj$Analyse.InSample$mllk)) == max(as.numeric(as.character(index.obj$Analyse.InSample$mllk)), na.rm=TRUE))
    bestAIC.IS  <- which(as.numeric(as.character(index.obj$Analyse.InSample$AIC.v)) == min(as.numeric(as.character(index.obj$Analyse.InSample$AIC.v)), na.rm=TRUE))
    bestAICc.IS <- which(as.numeric(as.character(index.obj$Analyse.InSample$AICc)) == min(as.numeric(as.character(index.obj$Analyse.InSample$AICc)), na.rm=TRUE))
    bestBIC.IS  <- which(as.numeric(as.character(index.obj$Analyse.InSample$BIC.v)) == min(as.numeric(as.character(index.obj$Analyse.InSample$BIC.v)), na.rm=TRUE))
    bestMAFE.IS  <- which(as.numeric(as.character(index.obj$Analyse.InSample$MAFE)) == min(as.numeric(as.character(index.obj$Analyse.InSample$MAFE)), na.rm=TRUE))
    bestMSFE.IS  <- which(as.numeric(as.character(index.obj$Analyse.InSample$MSFE)) == min(as.numeric(as.character(index.obj$Analyse.InSample$MSFE)), na.rm=TRUE))
    bestMASFE.IS <- which(as.numeric(as.character(index.obj$Analyse.InSample$MASFE)) == min(as.numeric(as.character(index.obj$Analyse.InSample$MASFE)), na.rm=TRUE))
    bestMSSFE.IS  <- which(as.numeric(as.character(index.obj$Analyse.InSample$MSSFE)) == min(as.numeric(as.character(index.obj$Analyse.InSample$MSSFE)), na.rm=TRUE))
    
    color.mask.IS <- list(colors = c("yellow"))
    color.mask.IS$mask <- matrix(0, nrow=dim(Analyse.matrix.in.sample.select)[1], 
                                 ncol=dim(Analyse.matrix.in.sample.select)[2])
    color.mask.IS$mask[bestLLK.IS,1]  <- 1
    color.mask.IS$mask[bestAIC.IS,2]  <- 1
    color.mask.IS$mask[bestAICc.IS,3] <- 1
    color.mask.IS$mask[bestBIC.IS,4]  <- 1
    color.mask.IS$mask[bestMAFE.IS,5]  <- 1
    color.mask.IS$mask[bestMSFE.IS,6]  <- 1
    color.mask.IS$mask[bestMASFE.IS,7] <- 1
    color.mask.IS$mask[bestMSSFE.IS,8]  <- 1
    
    # ———————————————————————————————————————————
    # COMMANDE LATEX - In sample
    # ———————————————————————————————————————————
    Analyse.in.Sample.LaTeX <- Make.LaTeX.Table(
      R.Matrix.Object=Analyse.matrix.in.sample.select,
      minimal.table=TRUE,
      # n.dec=round.to,
      title=paste0("Résultats IN.SAMPLE sur ",index.obj$data.pars$index," pour ",
                   index.obj$modele.data$train.DATES[1],
                   " au ",
                   index.obj$modele.data$train.DATES[length(index.obj$modele.data$train.DATES)]),
      print.Cons=FALSE,
      # table.type="table",
      table.type="sidewaystable",
      Row.Pos="c",
      Col.Titles=colnames(index.obj$Analyse.InSample[,variables.conservees.IS]),
      Row.Titles=as.character(index.obj$Analyse.InSample[,1]),
      bold.col.title=TRUE,
      bold.row.title=TRUE,
      color.mask=color.mask.IS,
      left.line=TRUE,
      Top.Line=TRUE,
      hlines=TRUE,
      double.hline.index=double.hline.index.IS,
      contract.time=TRUE
    )
    
    # ———————————————————————————————————————————
    # Out of sample analysis EXPORT LaTeX
    # ———————————————————————————————————————————
    # Sélection des variables à inclure dans le tableau
    variables.conservees.OofS <- !colnames(index.obj$Analyse.OutofSample) %in% retrait.colonnes
    Analyse.OutofSample.select <- matrix(unlist(index.obj$Analyse.OutofSample[,variables.conservees.OofS]),
                                         ncol=dim(index.obj$Analyse.OutofSample[,variables.conservees.OofS])[2])
    
    # ———————————————————————————————————————————
    # Obtention de la ligne correspondant à un changement de paramétrisation
    #   Dans ces cas la ligne est doublé dans la table LaTeX
    #   résultante.
    # ———————————————————————————————————————————
    double.hline.index.OofS <- NULL
    split.model.names <- strsplit(as.character(index.obj$Analyse.OutofSample[,1]),
                                  "-",fixed=TRUE)
    for (i in 1:length(split.model.names)){
      current.name <- paste0(split.model.names[[i]][-length(split.model.names[[i]])],
                             collapse="")
      if (i>1){
        if (current.name!= previous.name){
          double.hline.index.OofS <- double.hline.index.OofS %>% c((i-1))
        }
      }
      previous.name <- current.name
    }
    
    # ———————————————————————————————————————————
    # Obtention des indices correspondants aux max de LLK, AIC, AICc, BIC.v
    # ———————————————————————————————————————————
    bestLLK.OofS  <- which(as.numeric(as.character(index.obj$Analyse.OutofSample$mllk)) == max(as.numeric(as.character(index.obj$Analyse.OutofSample$mllk)), na.rm=TRUE))
    bestAIC.OofS  <- which(as.numeric(as.character(index.obj$Analyse.OutofSample$AIC.v)) == min(as.numeric(as.character(index.obj$Analyse.OutofSample$AIC.v)), na.rm=TRUE))
    bestAICc.OofS <- which(as.numeric(as.character(index.obj$Analyse.OutofSample$AICc)) == min(as.numeric(as.character(index.obj$Analyse.OutofSample$AICc)), na.rm=TRUE))
    bestBIC.OofS  <- which(as.numeric(as.character(index.obj$Analyse.OutofSample$BIC.v)) == min(as.numeric(as.character(index.obj$Analyse.OutofSample$BIC.v)), na.rm=TRUE))
    bestMAFE.OofS  <- which(as.numeric(as.character(index.obj$Analyse.OutofSample$MAFE)) == min(as.numeric(as.character(index.obj$Analyse.OutofSample$MAFE)), na.rm=TRUE))
    bestMSFE.OofS  <- which(as.numeric(as.character(index.obj$Analyse.OutofSample$MSFE)) == min(as.numeric(as.character(index.obj$Analyse.OutofSample$MSFE)), na.rm=TRUE))
    bestMASFE.OofS <- which(as.numeric(as.character(index.obj$Analyse.OutofSample$MASFE)) == min(as.numeric(as.character(index.obj$Analyse.OutofSample$MASFE)), na.rm=TRUE))
    bestMSSFE.OofS  <- which(as.numeric(as.character(index.obj$Analyse.OutofSample$MSSFE)) == min(as.numeric(as.character(index.obj$Analyse.OutofSample$MSSFE)), na.rm=TRUE))
    
    color.mask.OofS <- list(colors = c("yellow"))
    color.mask.OofS$mask <- matrix(0, nrow=dim(Analyse.OutofSample.select)[1], 
                                   ncol=dim(Analyse.OutofSample.select)[2])
    color.mask.OofS$mask[bestLLK.OofS,1]  <- 1
    color.mask.OofS$mask[bestAIC.OofS,2]  <- 1
    color.mask.OofS$mask[bestAICc.OofS,3] <- 1
    color.mask.OofS$mask[bestBIC.OofS,4]  <- 1
    color.mask.OofS$mask[bestMAFE.OofS,5]  <- 1
    color.mask.OofS$mask[bestMSFE.OofS,6]  <- 1
    color.mask.OofS$mask[bestMASFE.OofS,7] <- 1
    color.mask.OofS$mask[bestMSSFE.OofS,8]  <- 1
    
    
    # ———————————————————————————————————————————
    # COMMANDE LATEX - out of sample
    # ———————————————————————————————————————————
    Analyse.out.of.Sample.LaTeX <- Make.LaTeX.Table(
      R.Matrix.Object=Analyse.OutofSample.select,
      minimal.table=TRUE,
      # n.dec=round.to,
      title=paste0("Résultats OUT.OF.SAMPLE pour ",
                   index.obj$modele.data$test.DATES[1],
                   " au ",
                   index.obj$modele.data$test.DATES[length(index.obj$modele.data$test.DATES)]),
      print.Cons=FALSE,
      # table.type="table",
      table.type="sidewaystable",
      Row.Pos="c",
      Col.Titles=colnames(index.obj$Analyse.OutofSample[,variables.conservees.OofS]),
      Row.Titles=as.character(index.obj$Analyse.OutofSample[,1]),
      bold.col.title=TRUE,
      bold.row.title=TRUE,
      color.mask=color.mask.OofS,
      left.line=TRUE,
      Top.Line=TRUE,
      hlines=TRUE,
      double.hline.index=double.hline.index.OofS,
      contract.time=TRUE
    )
    
    
    # Correction de la chaîne de caractères LaTeX
    string.finale.analyses <- gsub(
      as.character(
        "\\end{landscape}\n % ------------------------ \n% ------------------------ \n\\begin{landscape}\\noindent"
      ),
      "",
      as.character(paste0(Analyse.in.Sample.LaTeX,Analyse.out.of.Sample.LaTeX)),
      fixed = TRUE)
    
    tryCatch({ 
      # Sauvegarde de la chaîne de caractères LaTeX dans l'objet de modèle
      index.obj$string.finale.analyses <- string.finale.analyses
      # Sauvegarde de la chaîne de caractères LaTeX dans le Presse-papier
      Copie.Presse.Papier(string.finale.analyses)
      
      tex.table.file <- paste0(index.obj$path.modele.ajd,"/",kept.index.obj,".tex")
      file.test(tex.table.file) # Teste si document existe déjà et le crée sinon
      writeLines(string.finale.analyses, tex.table.file)
    },
    error = function (err.message){
      err.message.builder(paste0("Erreur lors de l'enregistrement de l'objet temporaire dans la liste finale pour modèle \n",
                                 paste0("__",modele.info,"__ sur l'indice _",index.obj$data.pars$index,"_")," :\n"),
                          err.message
      )
    },
    warning = function (warn.message){
      warn.message.builder(paste0("Warning lors de l'enregistrement de l'objet temporaire dans la liste finale pour modèle \n",
                                  paste0("__",modele.info,"__ sur l'indice _",index.obj$data.pars$index,"_")," :\n"),
                           warn.message
      )
    })
    
  },
  error=function(err.message){
    err.message.builder(paste0("Erreur lors de la consolidation des résultats -> LaTeX \n       du modèle ",paste0("__",modele.info,"__ sur l'indice _",index.obj$data.pars$index,"_")," :\n"),
                        err.message)
  }
  )# Fin consolidation des résultats
  # ————————————————————————————————————————————————————————————————————————————————————
  
  
  # ————————————————————————————————————————————————————————————————————————————————————
  # ////////////////////////////////////////////////////////////////////
  # CRÉATION DES GRAPHIQUES ----
  # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  # ————————————————————————————————————————————————————————————————————————————————————
  # Création du répertoire pour enregistrer les graphiques
  dir.create(index.obj$path.modele.date.graphs, showWarnings = FALSE)
  
  # ————————————————————————————————————————————————————————————————————————————————————
  # PROBABILITÉS DIVERSES
  # ————————————————————————————————————————————————————————————————————————————————————
  probs.filtrees <<- list()
  probs.lissees <<- list()
  for (modele.info in index.obj$liste.modeles.converge){
    if (!index.obj[[modele.info]]$type %in% c("UNIVARIEE","MELANGE")){

      # Building the name for the graphs
      if (print.init){
        model.name.graph = index.obj[[modele.info]]$nom.modele.row
      } else {
        split.sur = strsplit(index.obj[[modele.info]]$nom.modele.row," sur ",fixed=TRUE)
        split.init = strsplit(split.sur[[1]][1],"-init",fixed=TRUE)
        if (length(unlist(split.sur))>1){
          model.name.graph = paste0(split.init[[1]][1]," sur ",gsub("$\\sim$","~",split.sur[[1]][2], fixed=TRUE))
        } else {
          model.name.graph = split.init[[1]][1]
        }
      }

      index.var.comb <- paste0("__",modele.info,"__ sur l'indice _",index.obj$data.pars$index,"_")

      tryCatch({
        index.obj[[modele.info]] %>%
          graph_probabilites(titre = paste0("Diagramme des probabilitées pour le modèle \n",
                                            model.name.graph," - (",index.obj$data.pars$freq.str,")",
                                            "\nTrain :",
                                            index.obj$modele.data$train.DATES[1],
                                            " - ",
                                            index.obj$modele.data$train.DATES[length(index.obj$modele.data$train.DATES)],
                                            "\nTest  :",
                                            index.obj$modele.data$test.DATES[1],
                                            " - ",
                                            index.obj$modele.data$test.DATES[length(index.obj$modele.data$test.DATES)]),
                             full.DATA = index.obj$modele.data$full.DATA,
                             full.Dates = index.obj$modele.data$full.DATES,
                             index.name = index.obj$data.pars$index,
                             colors = colors,
                             graph.file.name = paste0(which(index.obj$liste.modeles.converge == modele.info),
                                                      " - ",
                                                      gsub(".","_",
                                                           model.name.graph,fixed=TRUE)),
                             where = index.obj$path.modele.date.graphs
          )
      },
      error = function (err.message){
        err.message.builder(
          paste0("Erreur lors de la création du graphique de probabilités filtrées \n       du modèle ",paste0("__",modele.info,"__ sur l'indice _",index.obj$data.pars$index,"_")," :\n"),
          err.message
        )
      },
      warning = function(warn.message){
        warn.message.builder(paste0("Warning lors de la création du graphique de probabilités filtrées \n       du modèle ",paste0("__",modele.info,"__ sur l'indice _",index.obj$data.pars$index,"_")," :\n"),
                             warn.message)
      })
    }
  }
  
  tryCatch({
    index.obj$probs.filtrees <<- probs.filtrees
    index.obj$probs.lissees <<- probs.lissees
  },
  error = function (err.message){
    err.message.builder(
      paste0("Erreur lors de la Sauvegarde des graphiques pour la série ",index.obj$data.pars$index," :\n"),
      err.message
    )
  },
  warning = function(warn.message){
    warn.message.builder(paste0("Warning lors de la Sauvegarde des graphiques pour la série ",index.obj$data.pars$index," :\n"),
                         warn.message)
  })
  
  
  # ————————————————————————————————————————————————————————————————————————————————————
  # AUTOCORRÉLATION EMPIRIQUE VS SIMULÉE
  # ————————————————————————————————————————————————————————————————————————————————————
  index.obj$simulated.acf =
    simulate_hmm(
      index.obj = index.obj,
      lag.max = lag.max.acf,
      graph.file.name = paste0("ACF - ",index.obj$data.pars$index," - ",cat.Sys.Time()),
      n = n.sim,
      where = index.obj$path.modele.date.graphs,
      print.init = print.init
    )

  # Sauvegarde du tableau des statistiques descriptives
  tex.desc_table.file <- paste0(index.obj$path.modele.date.graphs,"/Statistiques_descriptive_simulées_",kept.index.obj,".tex")
  file.test(tex.desc_table.file) # Teste si document existe déjà et le crée sinon
  writeLines(index.obj$simulated.acf$tab.desc, tex.desc_table.file)
  
  # ————————————————————————————————————————————————————————————————————————————————————
  # ////////////////////////////////////////////////////////////////////
  # SAUVEGARDE DES MODÈLES FINAUX ET DE LEURS CARACTÉRISTIQUES ----
  # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  # ————————————————————————————————————————————————————————————————————————————————————
  tryCatch({
    # Enregistre l'objet de modèle complet à la même position dans la liste initiale
    liste.index <- liste.index[!names(liste.index) %in% as.character(index.obj$data.pars$index)]
    liste.index[[as.character(index.obj$data.pars$index)]] <- index.obj
    
    # Sauvegarde le modèle sous le nom de l'indice pour usage interne
    assign(as.character(index.obj$data.pars$index),index.obj)
  },
  error=function(err.message){
    err.message.builder(paste0("Erreur lors de l'extraction du modèle \n",paste0("__",modele.info,"__ sur l'indice _",index.obj$data.pars$index,"_")," :\n"),
                        err.message)
  })
  
  # Sauvegarde de l'objet de modèle complet à path.modele
  tryCatch({save.Model(liste = liste.index[[as.character(index.obj$data.pars$index)]], 
                       model.name = kept.index.obj,
                       path = index.obj$path.modele.ajd,
                       add.timeDate = FALSE)},
           error=function(err.message){
             err.message.builder(paste0("Erreur lors de la sauvegarde du modèle \n",paste0("__",modele.info,"__ sur l'indice _",index.obj$data.pars$index,"_")," :\n"),
                                 err.message)
           },
           warning = function(warn.mess){
             warn.message.builder(paste0("Warning lors de la sauvegarde du modèle \n",paste0("__",modele.info,"__ sur l'indice _",index.obj$data.pars$index,"_")," :\n"),
                                  warn.mess)
           }
  )
  
  # ————————————————————————————————————————————————————————————————————————————————————
  # ////////////////////////////////////////////////////////////////////
  # ENVOI DES RÉSULTATS PAR COURRIEL ----
  # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  # ————————————————————————————————————————————————————————————————————————————————————
  tryCatch({
    if (send.email.when.done){
      cat("Sending final email")
      email.conv(adresse.email = adresse.email,
                 pwd.server = pwd.server,
                 host.server = host.server,
                 port.server = port.server,
                 message.email= list(Analyse.InSample=index.obj$Analyse.InSample %>% arrange(as.numeric(as.character(BIC.v))),
                                     Analyse.OutofSample=index.obj$Analyse.OutofSample %>% arrange(as.numeric(as.character(BIC.v)))),
                 sujet=paste(index.obj$data.pars$index,"- [",index.obj$data.pars$date.boundaries[1]," ~ ",index.obj$data.pars$date.boundaries[2],"] - Fin de l'entrainement --- début :",index.time.start,", fin :",cat.Sys.Time()," - ",index.obj$path.modele.ajd))
    }
  },
  error = function(err.message){
    err.message.builder(paste0("Erreur lors de l'envoi des résultats par courriel \npour l'indice _",index.obj$data.pars$index,"_ :\n"),
                        err.message)
  } )
  
  # Closing all connections
  closeAllConnections()
  # ————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
} # Fin : BOUCLE sur toutes les spécifications et tous les indices
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ——————————————————————————————————————————————————————————————————————————————
# })

# # Retrait de tous les graphiques ouvert
# graphics.off()
# Closing all connections
closeAllConnections()
