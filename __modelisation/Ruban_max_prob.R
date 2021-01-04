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

ruban_max_prob_state = function(index.obj,
                        print.init = TRUE,
                        graph.file.name = "ACF",
                        where){
      # ————————————————————————————————————————————————————
      row.names <- c()
      ruban.max.prob = list()
      
      pdf.width <- 7.5 # En pouces
      pdf.height <- 9.75 # En pouces
      title.text.size <- 10 # Font size
      label.text.size <- 6
      element_line <- 0.1 # En cm
      plot.margins <- unit(c(0, 0.3, 0.2, 0.3), "cm")
      panel.margins <- unit(c(0, 0, 0, 0), "cm")
      # ————————————————————————————————————————————————————
      
      # ————————————————————————————————————————————————————
      # FUNCTION TO CALCULATE RANDOM INTEGRAL VALUE
      # ————————————————————————————————————————————————————
      state.prob.max = function(nbRegimes,prob){
            sample(c(1:nbRegimes), size = 1, replace = TRUE, prob = prob)
      }
      # ————————————————————————————————————————————————————
      # BOUCLE SUR TOUS LES MODÈLES POUR CALCULER LES STATS
      # ————————————————————————————————————————————————————
      for (modele.info in index.obj$liste.modeles.converge){
            if (!is.null(index.obj[[modele.info]]$HMM.Train) & !index.obj[[modele.info]]$type %in% c("UNIVARIEE", "MELANGE")){
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
                  
                  # ————————————————————————————————————————————————————
                  # EFFET DE TAYLOR - SIMULÉES
                  # ————————————————————————————————————————————————————
                  acf.SIM.Taylor = acf.SIM.abs
                  acf.SIM.Taylor$acf = acf.SIM.abs$acf - acf.SIM.square$acf
                  acf.SIM.Taylor.all[[as.character(modele.info)]] = acf.SIM.Taylor
                  
                  # ————————————————————————————————————————————————————
                  # FIN DE L'ÉTUDE DU MODÈLE
                  # ————————————————————————————————————————————————————
            } else {
                  cat("Le modèle",modele.info,"n'est pas simulé.\n")
            }
      } # FIN BOUCLE SUR TOUS LES MODÈLES
      # ————————————————————————————————————————————————————
      
      
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
                  scale_x_discrete(limits = c(1,lag.max)) +
                  theme(plot.title = element_text(size = title.text.size, face = "bold")) +
                  theme(axis.text=element_text(size=label.text.size),
                        axis.title=element_text(size=label.text.size,face="bold"))
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
      
      # ————————————————————————————————————————————————————
      return(
            list(
                  ruban.max.prob = ruban.max.prob
            )
      )
}
# ————————————————————————————————————————————————————
# ——————————————————————
# 
# n = 1000
# simulated.acf = simulate_hmm(index.obj,
#                              lag.max=200,
#                              print.init = FALSE,
#                              graph.file.name = paste0("ACF - ",cat.Sys.Time()),
#                              n= n,
#                              where="/Users/gabriellemyre/Documents/GitHub/Memoire_Maitrise/__resultats/__Résumé résultats pour Maciej")
# simulated.acf$acf.sim.sample.Ct.all
# cat(simulated.acf$tab.desc)
# ————————————————————————————————————————————————————
# 