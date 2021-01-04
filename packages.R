# ----------------------------------------------------------------------------------------------------
# Markov Chains
# External packages
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# First version : 22 avril 2020
# Last version  : 22 avril 2020
# ----------------------------------------------------------------------------------------------------


# Lecture de fichier texte
# install.packages('plyr') # Installation du package
library(plyr)
# install.packages('dplyr') # Installation du package
library(dplyr)

# install.packages('scales') # Installation du package
library(scales)

# install.packages('latex2exp') # Installation du package
library(latex2exp)

# Utilisation potentielle du package MHSMM pour évaluer HMMs
# install.packages('mhsmm')
library(mhsmm)

# Package pour HMM with Generalized Autoregressive Score
# remove.packages("HMMGAS",
#                 lib="/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
# install.packages("~/Documents/GitHub/Memoire_Maitrise/code/r/PACKAGES/HMMGAS", repos = NULL, type="source")
# detach_package("HMMGAS")
# unloadNamespace("HMMGAS")
library("HMMGAS")

#load numDeriv package for computing numerical derivatives
#enables the use of functions grad and hessian
# install.packages('numDeriv')
library(numDeriv)

#load Rsolnp package for numerical optimization with solnp function
# install.packages('Rsolnp')
library(Rsolnp)

# install.packages('rootSolve')
library(rootSolve)

# install.packages('data.table')
library(data.table) # Permet d'utiliser rbindlist

# install.packages('matrixcalc') # To test weither matrix is positive definite
library(matrixcalc) # Permet d'utiliser rbindlist

# install.packages('magic') #combiner matrice par blocs
library(magic)

# install.packages('gridExtra') # add tables in plots
library(gridExtra)

# install.packages('grid') # add tables in plots
library(grid)

# install.packages('plotrix') # add tables in plots
library(plotrix)

# install.packages('timeDate') # add tables in plots
library(timeDate)

# install.packages('rootSolve') # Outil d'analyse de séries chronologiques
library("rootSolve", lib.loc="~/Library/R/3.1/library")

# install.packages('MASS') # permet d'entrainer une normale sans spécifier la llk explicitement
library(MASS) ## loading package MASS

# install.packages('igraph') # basics of network analysis and visualization
library(igraph)

# Instalation du package forecast permettant d'utiliser les fonctions ACF, PAcf et CCf
# install.packages('forecast', dependencies = TRUE) # Installation du package
library(forecast)

# Permet de créer des graphiques très complexes
# install.packages('ggplot2', dependencies = TRUE) # Installation du package
library(ggplot2)

# install.packages('ggfortify') # Installation du package
library("ggfortify")

# Instalation du package cowplot opermettant de combiner des graphiques
# install.packages('cowplot', dependencies = TRUE) # Installation du package
library(cowplot)

# Instalation du package base qui permet d'utiliser NormalizedPath
# install.packages('base', dependencies = TRUE) # Installation du package
library(base)

# Lecture de fichier texte
# install.packages('readtext') # Installation du package
library(readtext)

# -------------------------------------------
# LIBRAIRIES POUR LES ENVOIS AUTOMATIQUES DE COURRIELS
# -------------------------------------------

# 
# install.packages("emayili")
library(emayili)

# 
# install.packages("magrittr")
library(magrittr)

# 
# install.packages("htmltools")
library(htmltools)

# 
# install.packages("htmltools")
library(textutils)

# Tool for profiling code in R
library(profvis)
