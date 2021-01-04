# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# INFORMATION POUR PARAMÉTRISATION DU S&P TSX
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# Last version : 15 avril 2019
# Last version : 4 mars 2020
# ----------------------------------------------------------------------------------------------------

# --------------------------------------------------------
# CADUSD - paramétrisation
# --------------------------------------------------------
# 
# --------------------------------------------------------



# TEST SI DOCUMENT APPELÉ DE L'EXTÉRIEUR OU PAS,
#   Si appelé de l'extérieur, 
if(tryCatch(outside.run, error = function(err){return(TRUE)})){
  # Closing all connections
  closeAllConnections()
  # Vide la cache de R
  rm(list=ls())
  # Retrait de tous les graphiques ouvert
  graphics.off()
  outside.run <- TRUE
  ignored.files.OUTSIDE.RUN <- "Parametrisations"
} else {
  ignored.files.OUTSIDE.RUN <- c()
}

# Chemin vers le document qui roule le modèle
path.MAIN <- '~/Documents/GitHub/Memoire_Maitrise/code/r/__modelisation/main.R'

CADUSD=list(
  data.pars = list(index="CADUSD",
                   variables.explicatives = c("INDPRO.log",
                                              "CADUSD"),
                   transformations = c("diff","scale"),
                   scalingFactor = c(100,100),
                   date.boundaries = c(as.Date("2003-09-17"),as.Date("2020-12-26")), # Partagé par toutes les modèles
                   freq.str="Q"
  ),
  general.pars = list(
    type.a.rouler = c("GAS"),
    load.model = "/Users/gabriellemyre/Documents/GitHub/Memoire_Maitrise/__resultats/__FINAUX__/CADUSD/[ 2003-09-17 ~ 2020-12-26 ]/CADUSD - (Q) - [2003-09-17, 2020-12-26] - 23h36m00s.rds",
    # test.Transformations=TRUE,
    # 
    mu0 = list(
      "1" = c(0.03),
      "2" = c(0.035, -0.012),
      "3" = c(0.188,  0.024, -0.051),
      "4" = c(0.269, -0.047, -0.003,  0.073)
    ),
    # 
    sigma0 = list(
      "1" = c(0.584543),
      "2" = c(0.8997, 0.3974),
      "3" = c(1.1059, 0.3974, 0.8997),
      "4" = c(1.234, 0.417, 0.718, 0.236)
    )
  ),
  UNIVARIEE = list(),
  # 
  MELANGE = list(
    K = c(2,
          3,
          4),
    gamma0 = list("2" = c(0.25, 0.75),
                  "3" = c(0.109, 0.345, 0.545),
                  "4" = c(0.048, 0.460, 0.340, 0.152))
  ),
  # 
  HMM = list(
    K = c(2,
          3,
          4),
    Transition.Type = "Homogeneous",
    gamma0 = list(
      "2" = matrix(c(0.98, 0.02, 
                     0.02, 0.98),
                   ncol=2 , byrow=T),
      "3" = matrix(c(0.98, 0.01, 0.01,
                     0.01, 0.98, 0.01,
                     0.01, 0.01, 0.98),
                   byrow=T,nrow=3),
      "4" = matrix(c(0.97, 0.01, 0.01, 0.01, 
                     0.01, 0.97, 0.01, 0.01,
                     0.01, 0.01, 0.97, 0.01,
                     0.01, 0.01, 0.01, 0.97),
                   ncol=4 , byrow=T)
    ),
    initial.Distribution = list(
      "2" = list(c(0.25, 0.75)),
      "3" = list(c(0,1,0)),
      "4" = list(c(0, 1, 0, 0)))
  ),
  # 
  HHMM = list(
    K=4,
    gamma0 = list(
      "4" = matrix(
        c(0.97, 0.02, 0.01,
          0.97, 0.02, 0.01,
          0.97, 0.02, 0.01,
          0.97, 0.02, 0.01,
          0.9,       0.1,          0,
          0.9,       0.1,          0,
          0.6,       0.4,          0,
          0.6,       0.4,          0),
        ncol=3,byrow=T)
    ),
    initial.Distribution = list(
      "4" = list(
        c(0.051, 0.570, 0.373, 0.006)
      )
    )
  ),
  # 
  HHMM.simplifie = list(
    K=4,
    gamma0 = list(
      "4" = matrix(
        c(0.98, 0.02, # Matrice de tansition parents
          0.02, 0.98,
          0.98, 0.02,  # Matrice de tansition enfants 1
          0.02, 0.98,
          0.98, 0.02,  # Matrice de tansition enfants 2
          0.02, 0.98),
        ncol=2,byrow=T)
    ),
    initial.Distribution = list(
      "4" = list(
        c(0, 1, 0, 0)
      )
    )
  ),
  # 
  LOGIT = list(
    K=4,
    Transition.Type = "LOGIT",
    gamma0 = list(
      "SP.500" = c(11.9153937,  -3.2901744,  -3.2196457,
                   9.6929927,   4.6357821,   3.6808849,
                   4.8872555,   9.7381110,  -2.8478310,
                   -11.7162761,   5.8048835,   1.0552933,
                   3.3021530,   1.9993717,   1.3619955,
                   -0.4021961,   0.2381303,  -0.4848302,
                   -1.5521677,  -0.5938818,  -3.0357197,
                   -1.9905202,  -3.1908825,  -7.6027283),
      "INDPRO.log" = c(27.484,  23.343,  -4.746,  
                       25.823,  20.729,  20.574,
                       1.703,   6.783, -62.430, 
                       -13.096,  20.549,  19.142,
                       2.527,   3.125,  -0.060,   
                       1.380,  -1.645,   1.889,
                       -0.242,   2.200, -34.287,  
                       4.161,  -3.287,  -0.039),
      "CADUSD" = c(4.323,   0.457,  -8.830,   
                   8.834,   3.424,   2.446,
                   4.669,   9.886,  -5.177, 
                   -11.719,  -0.222,  -0.977,
                   0.479,   0.778,  -0.316,   
                   3.821,   1.758,   2.672,
                   -3.336,  -4.191,  -3.442,  
                   -1.991,  -7.316,  -4.868)
    ),
    initial.Distribution = list(
      "4" = list(
        # c(0.051, 0.570, 0.373, 0.006),
        c(0, 1, 0, 0)
        # c(0, 0, 1, 0)
      )
    )
  ),
  #
  GAS = list(
    K=2,
    Transition.Type = "GAS",
    gamma0 = list(
      "2" = c(0.25, 0.75, # Distribution initiale
              0.9, 0.9, # Diagonale de la matrice A
              0.9, 0.9)
    )
  ),
  # 
  DDMS = list(
    K = 4,
    tau = c(
      3,
      5,
      14
    ),
    gamma0 = list(
      "3" = c(
        5.55,  -7.09, -17.60,
        1.92,   9.55,  0.33,
        -6.72, -5.44, -34.49,
        -17.40,  35.18,   1.88,
        -28.12, -8.43,  5.69,
        19.15,  15.08,  17.69,
        10.07 ,  9.10, -15.38,
        -19.69, -32.79,   1.5
      ),
      "5" = c(
        8.96, -1.03, -2.92,
        4.09,  4.85,  2.05,
        4.96,  5.86,  3.87,
        -4.52,  4.13, 1.82,
        1.69, -0.80,  2.08,
        0.54,  0.26, -0.30,
        6.38,  5.789, -0.1,
        -1.01, -2.07,  1.5
      ),
      "14" = c(
        8.96, -1.03, -2.92,
        4.09,  4.85,  2.05,
        4.96,  5.86,  3.87,
        -4.52,  4.13, 1.82,
        1.69, -0.80,  2.08,
        0.54,  0.26, -0.30,
        6.38,  5.789, -0.1,
        -1.01, -2.07,  1.5
      )
    ),
    initial.Distribution = list(
      "4" = list(
        # c(0.051, 0.570, 0.373, 0.006),
        c(0, 1, 0, 0)
      )
    )
  )
)

if(outside.run){
  cat("Running", path.MAIN)
  run.model <- CADUSD
  source(path.MAIN)
}


