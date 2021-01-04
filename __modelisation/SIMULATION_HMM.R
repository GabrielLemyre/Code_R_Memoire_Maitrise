# ————————————————————————————————————————————————————————————————————————————————————————————————————
# Hidden Markov Models
# SIMULATION DE SÉQUENCES DE HMMs
# ————————————————————————————————————————————————————————————————————————————————————————————————————
# written
# Gabriel LEMYRE
# ————————————————————————————————————————————————————————————————————————————————————————————————————
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ————————————————————————————————————————————————————————————————————————————————————————————————————
# First version : 24 décembre, 2020
# Last version : 27 décembre, 2020
# ————————————————————————————————————————————————————————————————————————————————————————————————————

# ————————————————————————————————————————————————————
# CRÉATION D'UNE SÉQUENCE SIMULÉE DEPUIS UN HMM
# ————————————————————————————————————————————————————
# Fonction produisant une séquence de longueur T
# ————————————————————————————————————————————————————

simulate_hmm = function(index.obj,
                        n,
                        print.init = TRUE,
                        lag.max=40,
                        graph.file.name = "ACF",
                        where){
  # ————————————————————————————————————————————————————
  probs.filtrees.in <<- probs.filtrees
  print(length(probs.filtrees.in))
  probs.lissees.in <<- probs.lissees
  print(length(probs.lissees.in))
  
  row.names <- c()
  sim.sample.Ct.all = list()
  sim.sample.Xt.all = list()
  acf.empirique = list()
  acf.empirique.abs = list()
  acf.empirique.square = list()
  
  acf.SIM.uncentered.all = list()
  acf.SIM.abs.all = list()
  acf.SIM.square.all = list()
  des_stats.SIM.all = list()
  acf.SIM.levier.NEG.all = list()
  acf.SIM.levier.POS.all = list()
  acf.SIM.levier.DIFF.all = list()
  acf.SIM.Taylor.all = list()
  
  pdf.width <- 7.5 # En pouces
  pdf.height <- 9.75 # En pouces
  title.text.size <- 10 # Font size
  label.text.size <- 6
  element_line <- 0.1 # En cm
  plot.margins <- unit(c(0, 0.3, 0.2, 0.3), "cm")
  panel.margins <- unit(c(0, 0, 0, 0), "cm")
  step.x.axis <- 20
  # ————————————————————————————————————————————————————
  # Centered returns for the studied series
  index.data = as.vector(unlist(index.obj$modele.data$full.DATA[index.obj$data.pars$index]))
  centered.returns.Xt <- index.data-mean(index.data)
  
  
  # ————————————————————————————————————————————————————
  # AUTO-CORRELATION - EMPIRIQUE
  # ————————————————————————————————————————————————————
  # ACF empirique des log-rendements sans transformations
  acf.empirique.uncentered <- ggAcf(
    as.vector(unlist(index.obj$modele.data$full.DATA[index.obj$data.pars$index])),
    lag.max = as.numeric(lag.max),
    type = "correlation",
    plot = FALSE
  )[1:lag.max]
  
  # ACF empirique des log-rendements centrés au carré
  acf.empirique.square <- ggAcf(
    (centered.returns.Xt)^2,
    lag.max = as.numeric(lag.max),
    type = "correlation",
    plot = FALSE
  )[1:lag.max]
  
  # ACF empirique des log-rendements centrés en valeur absolue
  acf.empirique.abs <- ggAcf(
    abs(centered.returns.Xt),
    lag.max = as.numeric(lag.max),
    type = "correlation",
    plot = FALSE
  )[1:lag.max]
  
  
  # ————————————————————————————————————————————————————
  # STATISTIQUES DESCRIPTIVES - EMPIRIQUE
  # ————————————————————————————————————————————————————
  des_stats.empirique <- list("Moyenne" = mean(index.data),
                              "Écart-type" = sd(index.data),
                              "Asymétrie" = timeDate::skewness(index.data, method="moment"),
                              "Aplatissement" = timeDate::kurtosis(index.data, method="moment"), 
                              "Minimum" = min(index.data),
                              "Maximum" = max(index.data), 
                              "n" = length(index.data)
  )
  
  # ————————————————————————————————————————————————————
  # EFFET DE LEVIER - EMPIRIQUE
  # ————————————————————————————————————————————————————
  max.0 <- function(x){max(x,0)}
  negChoc.empirique <- sapply(-centered.returns.Xt,FUN=max.0)
  posChoc.empirique <- sapply(centered.returns.Xt,FUN=max.0)
  
  # Corrélation entre la valeur absolue des erreurs et un choc négatif
  acf.empirique.levier.NEG <- 
    ggCcf(
      abs(centered.returns.Xt), 
      negChoc.empirique, 
      lag.max = as.numeric(lag.max),
      type = "correlation",
      plot = FALSE
    )[1:lag.max]
  
  # Corrélation entre la valeur absolue des erreurs et un choc positif
  acf.empirique.levier.POS <- 
    ggCcf(
      abs(centered.returns.Xt), 
      posChoc.empirique, 
      lag.max = as.numeric(lag.max),
      type = "correlation",
      plot = FALSE
    )[1:lag.max]
  
  
  acf.empirique.levier.DIFF = acf.empirique.levier.NEG
  acf.empirique.levier.DIFF$acf = acf.empirique.levier.NEG$acf - acf.empirique.levier.POS$acf
  
  # ————————————————————————————————————————————————————
  # EFFET DE TAYLOR - EMPIRIQUE
  # ————————————————————————————————————————————————————
  acf.empirique.Taylor = acf.empirique.abs
  acf.empirique.Taylor$acf = acf.empirique.abs$acf - acf.empirique.square$acf
  
  
  # ————————————————————————————————————————————————————
  # Dataframe des autocorrélations empiriques
  # ————————————————————————————————————————————————————
  acf.empirique.df <- data.frame(
    lag = acf.empirique.abs$lag,
    acf.empirique.uncentered = acf.empirique.uncentered$acf,
    acf.empirique.square = acf.empirique.abs$acf,
    acf.empirique.abs = acf.empirique.abs$acf,
    acf.empirique.levier.DIFF = acf.empirique.levier.DIFF$acf,
    acf.empirique.Taylor = acf.empirique.Taylor$acf
  )
  
  max.acf.empirique.uncentered <- max(acf.empirique.uncentered$acf)
  min.acf.empirique.uncentered <- min(acf.empirique.uncentered$acf)
  max.acf.empirique.square <- max(acf.empirique.square$acf)
  min.acf.empirique.square <- min(acf.empirique.square$acf)
  max.acf.empirique.abs <- max(acf.empirique.abs$acf)
  min.acf.empirique.abs <- min(acf.empirique.abs$acf)
  max.acf.empirique.levier.DIFF <- max(acf.empirique.levier.DIFF$acf)
  min.acf.empirique.levier.DIFF <- min(acf.empirique.levier.DIFF$acf)
  max.acf.empirique.Taylor <- max(acf.empirique.Taylor$acf)
  min.acf.empirique.Taylor <- min(acf.empirique.Taylor$acf)
  
  # ————————————————————————————————————————————————————
  # FUNCTION TO CALCULATE RANDOM INTEGRAL VALUE
  # ————————————————————————————————————————————————————
  random_int_hmm = function(nbRegimes,prob){
    sample(c(1:nbRegimes), size = 1, replace = TRUE, prob = prob)
  }
  # ————————————————————————————————————————————————————
  # BOUCLE SUR TOUS LES MODÈLES POUR CALCULER LES STATS
  # ————————————————————————————————————————————————————
  for (modele.info in index.obj$liste.modeles.converge){
    
    if (!is.null(index.obj[[modele.info]]$HMM.Train) & !swap.null(index.obj[[modele.info]]$Transition.Type,"Homogeneous") %in% c("GAS") & !index.obj[[modele.info]]$type %in% c("UNIVARIEE", "MELANGE")){
      # ——————————————————————————
      # Building the name for the graphs
      # ——————————————————————————
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
      
      # —————————————————————————————————————
      cat("Simulating :",model.name.graph,"\n")
      # —————————————————————————————————————
      
      # Sauvegarde le nom pour mettre sur le graphique final
      row.names <- c(row.names,model.name.graph)
      
      # Initialisation de la matrice pour recevoir la matrice de transition
      Gamma.Built = matrix()
      
      # Extraction des paramètres du modèle
      mu = index.obj[[modele.info]]$HMM.Train$mu
      sigma = index.obj[[modele.info]]$HMM.Train$sigma
      Gamma.v = index.obj[[modele.info]]$HMM.Train$matGamma
      
      type = index.obj[[modele.info]]$type
      Transition.Type = swap.null(index.obj[[modele.info]]$Transition.Type,"Homogeneous")
      nbStepsBack = index.obj[[modele.info]]$tau
      
      nbRegimes = length(mu)
      
      # Initialisation états équiprobables
      dist.initiale.diag.Gamma = rep(1,nbRegimes)/nbRegimes
      
      # ————————————————————————————————————————————————————
      # CONSTRUCTION DE GAMMA
      # ————————————————————————————————————————————————————
      if (type %in% c("HHMM","HHMM.simplifie","DDMS")){
        Gamma.Built.obj = Gamma.Build(prob.i=Gamma.v,
                                      type=type,
                                      Transition.Type=Transition.Type,
                                      nbStepsBack=nbStepsBack)
        
        Gamma.Built = Gamma.Built.obj$matGamma
        
        if (type %in% c("DDMS")){
          nbRegimes = nbRegimes*nbStepsBack
          mu = mu %x% rep(1,nbStepsBack)
          sigma = sigma %x% rep(1,nbStepsBack)
          dist.initiale.diag.Gamma = dist.initiale.diag.Gamma %x% c(1,rep(0,nbStepsBack-1)) # Kronecker product
        }
      } else if (!Transition.Type %in% c("LOGIT")){
        Gamma.Built = Gamma.v
      }
      
      
      # ————————————————————————————————————————————————————
      # SIMULATION
      # ————————————————————————————————————————————————————
      sim.sample.Ct <- rep(NA,n)
      sim.sample.Xt = rep(NA,n)
      centered.returns.SIM = rep(NA,n)
      
      # ——————————
      # t = 1
      # ——————————
      sim.sample.Ct[1] = random_int_hmm(nbRegimes,dist.initiale.diag.Gamma)
      sim.sample.Xt[1] = rnorm(1,mean = mu[sim.sample.Ct[1]], sd = sigma[sim.sample.Ct[1]])
      centered.returns.SIM[1] = sim.sample.Xt[1]-mu[sim.sample.Ct[1]]
      if (Transition.Type %in% c("LOGIT")){
        Gamma.Built.obj = Gamma.Build(prob.i = Gamma.v,
                                      Transition.Type=Transition.Type,
                                      exp.var=sim.sample.Xt[1],
                                      type=type)
        Gamma.Built = Gamma.Built.obj$matGamma
      }
      
      # ——————————
      # t = 2, ..., n
      # ——————————
      for (t in 2:n){
        sim.sample.Ct[t] = random_int_hmm(nbRegimes,Gamma.Built[sim.sample.Ct[t-1],])
        sim.sample.Xt[t] = rnorm(1,mean = mu[sim.sample.Ct[t]], sd = sigma[sim.sample.Ct[t]])
        centered.returns.SIM[t] = sim.sample.Xt[t]-mu[sim.sample.Ct[t]]
        
        if (Transition.Type %in% c("LOGIT")){
          Gamma.Built.obj = Gamma.Build(prob.i = Gamma.v,
                                        Transition.Type=Transition.Type,
                                        exp.var=sim.sample.Xt[t],
                                        type=type)
          Gamma.Built = Gamma.Built.obj$matGamma
        }
      }
      
      # Si DDMS, on rassemble les états simulés par paquets de longueurs
      #   tau "nbStepsBack" puisque correspondent au même état de la 
      #   chaîne d'origine
      if (type %in% c("DDMS")){
        sim.sample.Ct = ceiling(sim.sample.Ct/nbStepsBack)
      }
      
      # ENREGISTREMENT DES SÉQUENCES SIMULÉES
      sim.sample.Ct.all[[as.character(modele.info)]] = sim.sample.Ct
      sim.sample.Xt.all[[as.character(modele.info)]] = sim.sample.Xt
      
      # ————————————————————————————————————————————————————
      # AUTO-CORRELATION - SIMULÉES
      # ————————————————————————————————————————————————————
      # centered.returns.SIM <- as.vector(sim.sample.Xt)-mean(sim.sample.Xt)
      
      # ——————————————————————————
      # SANS TRANSFORMATION - SIMULÉES
      # ——————————————————————————
      acf.SIM <- 
        ggAcf(
          as.vector(sim.sample.Xt),
          lag.max = as.numeric(lag.max),
          type = "correlation",
          plot = FALSE
        )[1:lag.max]
      
      acf.SIM.uncentered.all[[as.character(modele.info)]] = acf.SIM
      
      # ——————————————————————————
      # ABSOLUTE VALUE - SIMULÉES
      # ——————————————————————————
      acf.SIM.abs <- 
        ggAcf(
          abs(centered.returns.SIM),
          lag.max = as.numeric(lag.max),
          type = "correlation",
          plot = FALSE
        )[1:lag.max]
      
      acf.SIM.abs.all[[as.character(modele.info)]] = acf.SIM.abs
      
      # ——————————————————————————
      # SQUARED - SIMULÉES
      # ——————————————————————————
      acf.SIM.square <- 
        ggAcf(
          (centered.returns.SIM)^2,
          lag.max = as.numeric(lag.max),
          type = "correlation",
          plot = FALSE
        )[1:lag.max]
      
      acf.SIM.square.all[[as.character(modele.info)]] = acf.SIM.square
      
      
      # ————————————————————————————————————————————————————
      # STATISTIQUES DESCRIPTIVES - SIMULÉES
      # ————————————————————————————————————————————————————
      des_stats.SIM <- list("Moyenne" = mean(sim.sample.Xt),
                            "Écart-type" = sd(sim.sample.Xt),
                            "Asymétrie" = timeDate::skewness(sim.sample.Xt, method="moment"),
                            "Aplatissement" = timeDate::kurtosis(sim.sample.Xt, method="moment"), 
                            "Minimum" = min(sim.sample.Xt),
                            "Maximum" = max(sim.sample.Xt), 
                            "n" = length(sim.sample.Xt)
      )
      
      des_stats.SIM.all[[as.character(modele.info)]] = des_stats.SIM
      
      # ————————————————————————————————————————————————————
      # EFFET DE LEVIER - SIMULÉES
      # ————————————————————————————————————————————————————
      max.0 <- function(x){max(x,0)}
      negChoc <- sapply(-centered.returns.SIM,FUN=max.0)
      posChoc <- sapply(centered.returns.SIM,FUN=max.0)
      
      # Corrélation entre la valeur absolue des erreurs et un choc négatif
      acf.SIM.levier.NEG <- 
        ggCcf(
          abs(centered.returns.SIM), 
          negChoc, 
          lag.max = as.numeric(lag.max),
          type = "correlation",
          # calc.ci = FALSE,
          plot = FALSE
        )[1:lag.max]
      
      acf.SIM.levier.NEG.all[[as.character(modele.info)]] = acf.SIM.levier.NEG
      
      # Corrélation entre la valeur absolue des erreurs et un choc positif
      acf.SIM.levier.POS <- 
        ggCcf(
          abs(centered.returns.SIM), 
          posChoc, 
          lag.max = as.numeric(lag.max),
          type = "correlation",
          # calc.ci = FALSE,
          plot = FALSE
        )[1:lag.max]
      
      acf.SIM.levier.POS.all[[as.character(modele.info)]] = acf.SIM.levier.POS
      
      acf.SIM.levier.DIFF = acf.SIM.levier.NEG
      acf.SIM.levier.DIFF$acf = acf.SIM.levier.NEG$acf - acf.SIM.levier.POS$acf
      acf.SIM.levier.DIFF.all[[as.character(modele.info)]] = acf.SIM.levier.DIFF
      
      # ————————————————————————————————————————————————————
      # EFFET DE TAYLOR - SIMULÉES
      # ————————————————————————————————————————————————————
      acf.SIM.Taylor = acf.SIM.abs
      acf.SIM.Taylor$acf = acf.SIM.abs$acf - acf.SIM.square$acf
      acf.SIM.Taylor.all[[as.character(modele.info)]] = acf.SIM.Taylor
      
      
      probs.filtrees
      # ————————————————————————————————————————————————————
      # PROBABILITÉS FILTRÉES
      # ————————————————————————————————————————————————————
      
      
      
      
      # ————————————————————————————————————————————————————
      # FIN DE L'ÉTUDE DU MODÈLE
      # ————————————————————————————————————————————————————
    } else {
      cat("Le modèle",modele.info,"n'est pas simulé.\n")
    }
  } # FIN BOUCLE SUR TOUS LES MODÈLES
  # ————————————————————————————————————————————————————
  
  # ————————————————————————————————————————————————————
  # SANS TRANSFORMATION
  # ————————————————————————————————————————————————————
  # Combinaison des plots de ACF
  final.plot.uncentered <- autoplot(acf.SIM.uncentered.all, ncol = 1)
  # final.plot <- new('ggmultiplot', plots = acf.sim.sample.Ct.all, ncol = 1)

  for (i in 1:length(acf.SIM.uncentered.all)) {
    final.plot.uncentered[i] <-
      final.plot.uncentered[i] +
      ggtitle(row.names[i]) +
      # geom_line(acf.empirique.df, mapping = aes(x = lag, y = acf)) +
      stat_smooth(acf.empirique.df, method = "loess", formula = y ~ x, size = 0.1, mapping = aes(x = lag, y = acf.empirique.uncentered)) +
      scale_y_continuous(limits = c(1.5*min(min.acf.empirique.uncentered,-0.01),1.01*max.acf.empirique.uncentered)) +
      scale_x_discrete(limits = c(1,seq.int(from = 0, to = lag.max, by = step.x.axis)[-1]),
                       breaks = c(1,seq.int(from = 0, to = lag.max, by = step.x.axis)[-1]),
                       labels = as.character(c(1,seq.int(from = 0, to = lag.max, by = step.x.axis)[-1]))
      ) +
      theme(plot.title = element_text(size = title.text.size, face = "bold")) +
      theme(axis.text=element_text(size=label.text.size),
            axis.title=element_text(size=label.text.size,face="bold"))

    if (i != length(acf.SIM.uncentered.all)){final.plot.uncentered[i] = final.plot.uncentered[i] + xlab("")}
  }

  final.plot.uncentered.fin =
    gridExtra::grid.arrange(
      grobs = final.plot.uncentered@plots,
      top = textGrob(paste0("Autocorrélation des log-rendements\nSimulations depuis ",index.name),
                     gp=gpar(fontsize=1.2*title.text.size,fontface = 'bold')),
      ncol=1
    )
  final.plot.uncentered.fin.bg <- cowplot::ggdraw(final.plot.uncentered.fin) +
    theme(plot.background = element_rect(fill="white", color = NA))

  pdf(paste(where,paste0(graph.file.name,"_1_uncentered.pdf"),sep="/"),
      width = pdf.width,
      height = pdf.height)
  plot(final.plot.uncentered.fin.bg)
  dev.off()


  # ————————————————————————————————————————————————————
  # ABSOLUTE VALUE
  # ————————————————————————————————————————————————————
  # Combinaison des plots de ACF
  final.plot.absolute <- autoplot(acf.SIM.abs.all, ncol = 1)
  # final.plot <- new('ggmultiplot', plots = acf.sim.sample.Ct.all, ncol = 1)

  for (i in 1:length(acf.SIM.abs.all)) {
    final.plot.absolute[i] <-
      final.plot.absolute[i] +
      ggtitle(row.names[i]) +
      # geom_line(acf.empirique.df, mapping = aes(x = lag, y = acf)) +
      stat_smooth(acf.empirique.df, method = "loess", formula = y ~ x, size = 0.1, mapping = aes(x = lag, y = acf.empirique.abs)) +
      scale_y_continuous(limits = c(1.5*min(min.acf.empirique.abs,-0.01),1.01*max.acf.empirique.abs)) +
      scale_x_discrete(limits = c(1,seq.int(from = 0, to = lag.max, by = step.x.axis)[-1]),
                       breaks = c(1,seq.int(from = 0, to = lag.max, by = step.x.axis)[-1]),
                       labels = as.character(c(1,seq.int(from = 0, to = lag.max, by = step.x.axis)[-1]))
      ) +
      theme(plot.title = element_text(size = title.text.size, face = "bold")) +
      theme(axis.text=element_text(size=label.text.size),
            axis.title=element_text(size=label.text.size,face="bold"))

    if (i != length(acf.SIM.abs.all)){final.plot.absolute[i] = final.plot.absolute[i] + xlab("")}
  }

  final.plot.absolute.fin =
    gridExtra::grid.arrange(
      grobs = final.plot.absolute@plots,
      top = textGrob(paste0("Autocorrélation des log-rendements centrés en valeur absolue\nSimulations depuis ",index.name),
                     gp=gpar(fontsize=1.2*title.text.size,fontface = 'bold')),
      ncol=1
    )
  final.plot.absolute.fin.bg <- cowplot::ggdraw(final.plot.absolute.fin) +
    theme(plot.background = element_rect(fill="white", color = NA))

  pdf(paste(where,paste0(graph.file.name,"_2_absolute_residues.pdf"),sep="/"),
      width = pdf.width,
      height = pdf.height)
  plot(final.plot.absolute.fin.bg)
  dev.off()

  # ————————————————————————————————————————————————————
  # SQUARED
  # ————————————————————————————————————————————————————
  # Combinaison des plots de ACF
  final.plot.squared <- autoplot(acf.SIM.square.all, ncol = 1)
  # final.plot <- new('ggmultiplot', plots = acf.sim.sample.Ct.all, ncol = 1)

  for (i in 1:length(acf.SIM.square.all)) {
    final.plot.squared[i] <-
      final.plot.squared[i] +
      ggtitle(row.names[i]) +
      # geom_line(acf.empirique.df, mapping = aes(x = lag, y = acf)) +
      stat_smooth(acf.empirique.df, method = "loess", formula = y ~ x, size = 0.1, mapping = aes(x = lag, y = acf.empirique.square)) +
      scale_y_continuous(limits = c(1.5*min(min.acf.empirique.square,-0.01),1.01*max.acf.empirique.square)) +
      scale_x_discrete(limits = c(1,seq.int(from = 0, to = lag.max, by = step.x.axis)[-1]),
                       breaks = c(1,seq.int(from = 0, to = lag.max, by = step.x.axis)[-1]),
                       labels = as.character(c(1,seq.int(from = 0, to = lag.max, by = step.x.axis)[-1]))
      ) +
      theme(plot.title = element_text(size = title.text.size, face = "bold")) +
      theme(axis.text=element_text(size=label.text.size),
            axis.title=element_text(size=label.text.size,face="bold"))

    if (i != length(acf.SIM.square.all)){final.plot.squared[i] = final.plot.squared[i] + xlab("")}
  }

  final.plot.squared.fin =
    gridExtra::grid.arrange(
      grobs = final.plot.squared@plots,
      top = textGrob(paste0("Autocorrélation des log-rendements centrés au carré\nSimulations depuis ",index.name),
                     gp=gpar(fontsize=1.2*title.text.size,fontface = 'bold')),
      ncol=1
    )
  final.plot.squared.fin.bg <- cowplot::ggdraw(final.plot.squared.fin) +
    theme(plot.background = element_rect(fill="white", color = NA))

  pdf(paste(where,paste0(graph.file.name,"_3_squared_residues.pdf"),sep="/"),
      width = pdf.width,
      height = pdf.height)
  plot(final.plot.squared.fin.bg)
  dev.off()


  # ————————————————————————————————————————————————————
  # STATISTIQUES DESCRIPTIVES
  # ————————————————————————————————————————————————————
  n.modele <- length(des_stats.SIM.all)
  des_stats_final = matrix(c(index.name,as.vector(round(as.numeric(as.character(unlist(des_stats.empirique))),3))),
                           ncol=8,
                           nrow=1)
  for (i in 1:n.modele) {
    des_stats_final = rbind(des_stats_final,
                            matrix(c(row.names[i],as.vector(round(as.numeric(as.character(unlist(des_stats.SIM.all[[i]]))),3))),nrow=1))
  }
  des_stats_final <- as.data.frame(des_stats_final)
  colnames(des_stats_final) = c("","mu","sigma","asymétrie","aplatissment","Min","Max","n")

  cat("\n")
  hequal()
  cat("Analyse par simulation des statistiques descriptives\n")
  hline()
  print(des_stats_final)
  hequal()

  des_stats_final.for.LaTeX = matrix(unlist(des_stats_final),
                                     ncol=dim(des_stats_final)[2])

  tab.desc <- Make.LaTeX.Table(des_stats_final.for.LaTeX[,-1],
                               minimal.table=TRUE,
                               Col.Titles = colnames(des_stats_final),
                               Row.Titles=as.character(des_stats_final[,1]),
                               Cross.Lines = T,
                               n.dec = 3,
                               Row.Pos = 'c',
                               title=paste("Statistiques descriptives sur ",index.name,sep=""),
                               print.Cons=FALSE,
                               copy.CB=TRUE)
  # ————————————————————————————————————————————————————
  # EFFET DE LEVIER
  # ————————————————————————————————————————————————————
  # Combinaison des plots de ACF
  final.plot.LEVIER <- autoplot(acf.SIM.levier.DIFF.all, ncol = 1)
  # final.plot <- new('ggmultiplot', plots = acf.sim.sample.Ct.all, ncol = 1)

  for (i in 1:length(acf.SIM.levier.DIFF.all)) {
    final.plot.LEVIER[i] <-
      final.plot.LEVIER[i] +
      ggtitle(row.names[i]) +
      # geom_line(acf.empirique.df, mapping = aes(x = lag, y = acf)) +
      stat_smooth(acf.empirique.df, method = "loess", formula = y ~ x, size = 0.1, mapping = aes(x = lag, y = acf.empirique.levier.DIFF)) +
      scale_y_continuous(limits = c(1.5*min(min.acf.empirique.levier.DIFF,-0.01),1.01*max.acf.empirique.levier.DIFF)) +
      scale_x_discrete(limits = c(1,seq.int(from = 0, to = lag.max, by = step.x.axis)[-1]),
                       breaks = c(1,seq.int(from = 0, to = lag.max, by = step.x.axis)[-1]),
                       labels = as.character(c(1,seq.int(from = 0, to = lag.max, by = step.x.axis)[-1]))
      ) +
      theme(plot.title = element_text(size = title.text.size, face = "bold")) +
      theme(axis.text=element_text(size=label.text.size),
            axis.title=element_text(size=label.text.size,face="bold"))

    if (i != length(acf.SIM.levier.DIFF.all)){final.plot.LEVIER[i] = final.plot.LEVIER[i] + xlab("")}
  }

  final.plot.LEVIER.fin =
    gridExtra::grid.arrange(
      grobs = final.plot.LEVIER@plots,
      top = textGrob(paste0("Effet de Levier \nSimulations depuis ",index.name),
                     gp=gpar(fontsize=1.2*title.text.size,fontface = 'bold')),
      ncol=1
    )
  final.plot.LEVIER.fin.bg <- cowplot::ggdraw(final.plot.LEVIER.fin) +
    theme(plot.background = element_rect(fill="white", color = NA))

  pdf(paste(where,paste0(graph.file.name,"_4_effet_levier.pdf"),sep="/"),
      width = pdf.width,
      height = pdf.height)
  plot(final.plot.LEVIER.fin.bg)
  dev.off()


  # ————————————————————————————————————————————————————
  # EFFET DE TAYLOR
  # ————————————————————————————————————————————————————
  # Combinaison des plots de ACF
  final.plot.TAYLOR <- autoplot(acf.SIM.Taylor.all, ncol = 1)

  for (i in 1:length(acf.SIM.Taylor.all)) {
    final.plot.TAYLOR[i] <-
      final.plot.TAYLOR[i] +
      ggtitle(row.names[i]) +
      # geom_line(acf.empirique.df, mapping = aes(x = lag, y = acf)) +
      stat_smooth(acf.empirique.df, method = "loess", formula = y ~ x, size = 0.1, mapping = aes(x = lag, y = acf.empirique.Taylor)) +
      scale_y_continuous(limits = c(1.5*min(min.acf.empirique.Taylor,-0.01),1.01*max.acf.empirique.Taylor)) +
      scale_x_discrete(limits = c(1,seq.int(from = 0, to = lag.max, by = step.x.axis)[-1]),
                       breaks = c(1,seq.int(from = 0, to = lag.max, by = step.x.axis)[-1]),
                       labels = as.character(c(1,seq.int(from = 0, to = lag.max, by = step.x.axis)[-1]))
      ) +
      theme(plot.title = element_text(size = title.text.size, face = "bold")) +
      theme(axis.text=element_text(size=label.text.size),
            axis.title=element_text(size=label.text.size,face="bold"))

    if (i != length(acf.SIM.Taylor.all)){final.plot.TAYLOR[i] = final.plot.TAYLOR[i] + xlab("")}
  }

  final.plot.TAYLOR.fin =
    gridExtra::grid.arrange(
      grobs = final.plot.TAYLOR@plots,
      top = textGrob(paste0("Effet de Taylor \nSimulations depuis ",index.name),
                     gp=gpar(fontsize=1.2*title.text.size,fontface = 'bold')),
      ncol=1
    )
  final.plot.TAYLOR.fin.bg <- cowplot::ggdraw(final.plot.TAYLOR.fin) +
    theme(plot.background = element_rect(fill="white", color = NA))

  pdf(paste(where,paste0(graph.file.name,"_5_effet_Taylor.pdf"),sep="/"),
      width = pdf.width,
      height = pdf.height)
  plot(final.plot.TAYLOR.fin.bg)
  dev.off()
  
  #  Nom modèle GAS
  # split.sur = strsplit(index.obj$SP500_HMM.GAS.K.2$nom.modele.row," sur ",fixed=TRUE)
  # split.init = strsplit(split.sur[[1]][1],"-init",fixed=TRUE)
  # GAS.name.graph = split.init[[1]][1]
  # row.names <- c(row.names[1:7],GAS.name.graph,row.names[8:length(row.names)])
  # ————————————————————————————————————————————————————
  # PROBABILITÉS FILTRÉES
  # ————————————————————————————————————————————————————
  # Combinaison des plots de ACF
  # final.plot.filtrees <- autoplot(probs.filtrees, ncol = 1)
  # 
  # for (i in 1:length(probs.filtrees)) {
  #   final.plot.filtrees[i] <-
  #     final.plot.filtrees[i] + 
  #     ggtitle(row.names[i]) +
  #     theme(plot.title = element_text(size = title.text.size, face = "bold")) +
  #     theme(axis.text=element_text(size=label.text.size),
  #           axis.title=element_text(size=label.text.size,face="bold"))
  #   
  #   if (i != length(probs.filtrees)){
  #     final.plot.filtrees[i] = final.plot.filtrees[i] + xlab("")+ 
  #       theme(axis.text.x = element_blank(),
  #             axis.title.x = element_blank(),
  #             axis.ticks.x = element_blank()
  #       )
  #   }
  # }
  # 
  # final.plot.filtrees.fin = 
  #   gridExtra::grid.arrange(
  #     grobs = final.plot.filtrees@plots, 
  #     top = textGrob(paste0("Probabilités filtrées sur ",index.name), 
  #                    gp=gpar(fontsize=1.2*title.text.size,fontface = 'bold')),
  #     ncol=1,
  #     heights = c(1,1,1,1,1,1,1,1,1,1,1.25)
  #   )
  # final.plot.filtrees.fin.bg <- cowplot::ggdraw(final.plot.filtrees.fin) + 
  #   theme(plot.background = element_rect(fill="white", color = NA))
  # 
  # pdf(paste(where,paste0(graph.file.name,"probabilites_filtrees.pdf"),sep="/"),
  #     width = pdf.width,
  #     height = pdf.height)
  # plot(final.plot.filtrees.fin.bg)
  # dev.off()
  
  
  # ————————————————————————————————————————————————————
  # PROBABILITÉS LISSÉES
  # ————————————————————————————————————————————————————
  # Combinaison des plots de ACF
  # final.plot.lissees <- autoplot(probs.lissees, ncol = 1)
  # 
  # for (i in 1:length(probs.lissees)) {
  #   final.plot.lissees[i] <-
  #     final.plot.lissees[i] +
  #     ggtitle(row.names[i]) +
  #     theme(plot.title = element_text(size = title.text.size, face = "bold")) +
  #     theme(axis.text=element_text(size=label.text.size),
  #           axis.title=element_text(size=label.text.size,face="bold")) +
  #     theme(legend.position="none")
  #   
  #   if (i != length(probs.lissees)){
  #     final.plot.lissees[i] = final.plot.lissees[i] + xlab("")+ 
  #       theme(axis.text.x = element_blank(),
  #             axis.title.x = element_blank(),
  #             axis.ticks.x = element_blank()
  #       )
  #   }
  # }
  # 
  # final.plot.lissees.fin =
  #   gridExtra::grid.arrange(
  #     grobs = final.plot.lissees@plots,
  #     top = textGrob(paste0("Probabilités lissées sur ",index.name),
  #                    gp=gpar(fontsize=1.2*title.text.size,fontface = 'bold')),
  #     ncol=1,
  #     heights = c(1,1,1,1,1,1,1,1,1,1,1.25)
  #   )
  # final.plot.lissees.fin.bg <- cowplot::ggdraw(final.plot.lissees.fin) +
  #   theme(plot.background = element_rect(fill="white", color = NA))
  # 
  # pdf(paste(where,paste0(graph.file.name,"probabilites_lissees.pdf"),sep="/"),
  #     width = pdf.width,
  #     height = pdf.height)
  # plot(final.plot.lissees.fin.bg)
  # dev.off()
  
  
  # ————————————————————————————————————————————————————
  return(
    list(
      tab.desc=tab.desc,
      sim.sample.Ct.all = sim.sample.Ct.all,
      sim.sample.Xt.all = sim.sample.Xt.all,
      acf.SIM.uncentered.all = acf.SIM.uncentered.all,
      acf.SIM.abs.all = acf.SIM.abs.all,
      acf.SIM.square.all = acf.SIM.square.all,
      des_stats.SIM.all = des_stats.SIM.all,
      acf.SIM.levier.NEG.all = acf.SIM.levier.NEG.all,
      acf.SIM.levier.POS.all = acf.SIM.levier.POS.all,
      acf.SIM.levier.DIFF.all = acf.SIM.levier.DIFF.all,
      acf.SIM.Taylor.all = acf.SIM.Taylor.all,
      final.plot.uncentered = final.plot.uncentered,
      final.plot.absolute = final.plot.absolute,
      final.plot.squared = final.plot.squared,
      final.plot.LEVIER = final.plot.LEVIER,
      final.plot.TAYLOR = final.plot.TAYLOR
    )
  )
}
# ————————————————————————————————————————————————————
# ——————————————————————
# 
# n = 100
# simulated.acf = simulate_hmm(index.obj,
#                              lag.max=200,
#                              print.init = FALSE,
#                              graph.file.name = paste0("ACF - ",cat.Sys.Time()),
#                              n= n,
#                              where="/Users/gabriellemyre/Documents/GitHub/Memoire_Maitrise/__resultats/__FINAUX__")
# simulated.acf$acf.sim.sample.Ct.all
# cat(simulated.acf$tab.desc)
# ————————————————————————————————————————————————————
# 