# ----------------------------------------------------------------------------------------------------
# WRAPPER
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# First version : 28 février 2020
# Last version  : 5 mars 2020
# ----------------------------------------------------------------------------------------------------

# --------------------------------------------------------
# WRAPPER FUNCTION
# --------------------------------------------------------
# Routine d'importation de toutes les fonctions et
#   packages necessaires à l'entrainement, aux analyses
#   et aux tests sur les modèles
# --------------------------------------------------------

# --------------------------------------------------------
# Modification à l'environnement global
# --------------------------------------------------------
options(max.print=1000000)

# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# LIBRAIRIES INTERNES
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# -———————————————————————————————————————————————————————————————————————————————————

# --------------------------------------------------------
# Sourcing private libraries
# R-Tools
# --------------------------------------------------------
path <- '~/Documents/GitHub/R-Tools'
setwd(path.expand(path)) # Setting Sourcing path

# Essai d'importer les fonctions de base, 
#   necessaires pour obtenir la fonction "sources_correctly"
file <- "BasicFunctions.R"
fn <- try(source(file), silent = TRUE)
if (inherits(fn, "try-error")){
  stop(paste("\n-------------------------\nFile :",file,"\nError :",fn))
}

# Autres fonctions externes standards
sources_correctly("RunTime.R")
sources_correctly("LaTeXTable.R")
# --------------------------------------------------------

# --------------------------------------------------------
# Changement du document de travail
# --------------------------------------------------------
path <- '~/Documents/GitHub/Memoire_Maitrise/code/r'
setwd(path.expand(path)) # Setting Sourcing path


# --------------------------------------------------------
# Sourcing private libraries
# Fonctions de bases pour l'entrainement du modèle
# --------------------------------------------------------
# Obtention de la liste de toutes les fonctions .R proches de wrapper.R
source.files.liste <- file.list(path,ignored.files.vector=c("__modelisation",
                                                            "wrapper.R","tests",
                                                            "resultats","PACKAGES"),print.inside.message = FALSE)
n.files <- length(source.files.liste)

# Sourcing all files in the given list
for (i in 1:n.files){
  hline()
  cat("sourcing",source.files.liste[i],"\n")
  sources_correctly(source.files.liste[i])
}


# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# IMPORTATION FONCTIONS FINALE DE MODÉLISATION ----
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ————————————————————————————————————————————————————————————————————————————————————
source.modelisation.liste <- file.list('~/Documents/GitHub/Memoire_Maitrise/code/r/__modelisation',
                                ignored.files.vector=c("main.R",ignored.files.OUTSIDE.RUN),
                                print.inside.message = FALSE)

n.files.mod <- length(source.modelisation.liste)

# Sourcing all files in the given list
for (i in 1:n.files.mod){
  hline()
  cat("sourcing",source.modelisation.liste[i],"\n")
  sources_correctly(source.modelisation.liste[i])
}