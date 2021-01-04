# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Graphique des probabilités filtrées et lissées
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# First version : january 8th, 2020
# Last version : february 13th, 2020
# ----------------------------------------------------------------------------------------------------

# --------------------------------------------------------
# CRÉATION ET EXPORTATION D'UN GRAPHIQUE AFFICHANT LES PROBABILITÉS FILTRÉES 
#   IN AND AOUT OF SAMPLE 
# --------------------------------------------------------
# Fonction produisant des graphiques
# --------------------------------------------------------

graph_probabilites <- function(objet.modele,
                               titre,
                               full.Dates,
                               full.DATA,
                               index.name,
                               test.Dates = NULL,
                               colors,
                               graph.file.name = "Graph",
                               where){
  
  # IS.probs = index.obj[[modele.info]]$filtered.Prob,
  # IS.Dates = index.obj$modele.data$full.Dates,
  # OofS.probs = index.obj[[modele.info]]$HMM.Test$out.of.sample.FILTERED.prob,
  # OofS.Dates = index.obj$modele.data$test.DATES,
  # ———————————————————————————————————————————————————————
  # Paramètres graphiques ----
  # ———————————————————————————————————————————————————————
  # Nombre de points à afficher sur l'axe horizontal
  nbTicks=100
  
  hline()
  cat("Création du graphique",titre,"\n")
  
  # PDF DIMENSIONS FOR OUTPUTED FILE
  pdf.width <- 14 # En cm
  pdf.height <- 16 # En cm
  text.size <- 6 # Font size
  legend.text.size <- 8
  element_line <- 0.1 # En cm
  plot.margins <- unit(c(0, 0.3, 0.2, 0.3), "cm")
  panel.margins <- unit(c(0, 0, 0, 0), "cm")
  
  # Options ggplot
  alpha=1
  size=0.1
  vline.size=0.5
  vline.color <- "white"
  vlinetype <- "longdash"
  hline.color <- "black"
  
  # ———————————————————————————————————————————————————————
  # VARIABLES pour GRAPHS ----
  # ———————————————————————————————————————————————————————
  # n.IS <- length(train.Dates)
  # cat('n.IS =',n.IS,"\n")
  n.OofS <- length(full.Dates)
  first.OofS.date <- objet.modele$HMM.Test$first.OofS.date
  
  # ticks.step <- floor(n.OofS/nbTicks)
  # x.ticks <- seq(ticks.step,n.OofS,ticks.step)
  
  # DATA SETS
  index.DATA <- unlist(full.DATA[index.name])
  
  variable.explicative = objet.modele$variable.explicative
  
  # nbRegimes depuis la taille de mu0
  nbRegimes <- length(objet.modele$mu0)
  
  # DIAGONALE DE GAMMA
  OofS.diag.Gamma <- tryCatch({
    as.matrix(objet.modele$HMM.Test$diag.Gamma[1:nbRegimes,])
  },
  error = function(err){
    return(NULL)
  },
  warning = function(warn){
    return(NULL)
  })
  
  # OUT OF SAMPLE
  OofS.probs.filtered <- as.matrix(objet.modele$HMM.Test$out.of.sample.FILTERED.prob)
  OofS.probs.smooth <- as.matrix(objet.modele$HMM.Test$out.of.sample.SMOOTH.prob)
  
  # OofS.probs.filtered[is.na(OofS.probs.filtered)] <- 0
  # OofS.probs.smooth[is.na(OofS.probs.smooth)] <- 0
  
  # COLLAPSE PAR RÉGIME SI DDMS
  if (objet.modele$type=="DDMS"){
    matrice.de.somme <- diag(nbRegimes) %x% t(rep(1,swap.null(objet.modele$tau,1)))
    # IS.probs.filtered <- matrice.de.somme %*% IS.probs.filtered
    # IS.probs.smooth <- matrice.de.somme %*% IS.probs.smooth
    OofS.probs.filtered <- matrice.de.somme %*% OofS.probs.filtered
    OofS.probs.smooth <- matrice.de.somme %*% OofS.probs.smooth
  }
  
  # Noms à donner au groupe dans le graph
  regime.names <- paste0("mu = ",round(objet.modele$HMM.Train$mu,3),",\n",
                         "sigma = ",round(objet.modele$HMM.Train$sigma,3))
  
  # cat("sigma =",objet.modele$HMM.Train$sigma,"\n")
  
  # ORDERING THE ROWS
  temp.for.order <- data.frame(idx = c(1:nbRegimes),
                               sigma = objet.modele$HMM.Train$sigma,
                               stringsAsFactors = FALSE) %>% arrange(sigma)
  
  idx <- as.numeric(temp.for.order$idx)
  # cat("idx =",idx,"\n")
  
  regime.names.order <- regime.names[idx]
  
  # IS.probs.filtered <- matrice.de.somme %*% IS.probs.filtered
  # IS.probs.smooth <- matrice.de.somme %*% IS.probs.smooth
  # OofS.probs.filtered <- matrice.de.somme %*% OofS.probs.filtered
  # OofS.probs.smooth <- matrice.de.somme %*% OofS.probs.smooth
  
  name.IS <- "In Sample"
  name.OofS <- "Out of Sample"
  #  Essai de faire le graph
  tryCatch({
    # Matrice de somme cumulée sur les lignes
    matrice.somme.cumule <- lower.tri(diag(nbRegimes)) + diag(nbRegimes)
    
    # # Assignation des groupes aux observations
    # groups.IS <- matrice.somme.cumule %*% matrix(1,
    #                                              ncol=dim(IS.probs.filtered)[2],
    #                                              nrow=dim(IS.probs.filtered)[1])
    
    # Assignation des groupes aux observations
    groups.OofS <- matrice.somme.cumule %*% matrix(1,
                                                   ncol=dim(OofS.probs.filtered)[2],
                                                   nrow=dim(OofS.probs.filtered)[1])
    # print(groups.OofS)
    # cat("Dates :",typeof(.0001),"\n")
    # cat("Filtered :",length(as.numeric(as.character(as.vector(OofS.probs.filtered)))),"\n")
    # cat("Smooth :",length(as.numeric(as.character(as.vector(OofS.probs.filtered)))),"\n")
    
    plot.data <- data.frame(
      index.DATA = as.numeric(rep(index.DATA, each=nbRegimes)),
      dates = full.Dates[as.numeric(rep(c(1:n.OofS), each = nbRegimes))],
      filtered = as.numeric(as.character(OofS.probs.filtered)),
      smooth = as.numeric(as.character(OofS.probs.smooth)),
      group = as.vector(regime.names[groups.OofS])
    )
    # print(head(plot.data))
    
    
    # plot(x=full.Dates[c(1:n.OofS)],y=OofS.probs.filtered[1,],type="l")
    # plot(x=full.Dates[c(1:n.OofS)],y=OofS.probs.filtered[2,],type="l")
    
    
    if (!is.null(OofS.diag.Gamma)){
      plot.data$diag.Gamma = as.numeric(as.character(OofS.diag.Gamma))
      # pdf.height <- 5/4*pdf.height
    }
    
    if (!is.null(variable.explicative)){
      if (variable.explicative!=index.name){
        var.exp.DATA = rep(as.numeric(as.character(unlist(full.DATA[variable.explicative]))), each=nbRegimes)
        
        plot.data$var.exp.DATA = var.exp.DATA
        # pdf.height = 4/3*pdf.height
        
        add.var.exp = TRUE
      } else {
        add.var.exp = FALSE
      }
    } else {
      add.var.exp = FALSE
    }
    
    plot.data$dates <- as.Date(as.character(plot.data$dates))
    plot.data$group <- factor(plot.data$group, levels = regime.names.order )
    # plot.data$group <- factor(plot.data$group, levels = regime.names.order )
    
    
    # SÉRIE À L'ÉTUDE
    index.plot <<-  
      ggplot(plot.data, aes(x=dates, y=index.DATA))+
      geom_line(size=size) + 
      geom_hline(yintercept = 0, 
                 size = size, 
                 color = hline.color) +
      geom_vline(xintercept = first.OofS.date, 
                 size = vline.size, 
                 color = vline.color, 
                 linetype=vlinetype) +
      xlab("") + 
      # xlim(1,n.OofS) +
      ylab(index.name) + 
      theme(legend.position="none") +
      theme(plot.margin = plot.margins,
            panel.spacing = panel.margins,
            text = element_text(size = text.size), 
            element_line(size = element_line),
            axis.text.x = element_text(angle=90, hjust=0.5, vjust=0.5,size = text.size),
            axis.title.x = element_blank()
            # axis.title.x = element_blank(),
            # axis.ticks.x = element_blank(),
            # panel.grid = element_blank(),
            # panel.border = element_blank()
      ) + 
      scale_x_date(breaks = "1 year", 
                   limits = as.Date(full.Dates[c(1,n.OofS)]), 
                   labels=date_format("%Y")) +
      scale_y_continuous(expand = c(0, 0))
    
    # PROBABILITÉS FILTRÉES
    Graph.filtered <-  
      ggplot(plot.data, aes(x=dates, y=filtered, 
                            fill=group)) + 
      geom_area() + 
      geom_hline(yintercept = 0.5, 
                 size = size, 
                 color = hline.color) +
      geom_vline(xintercept = first.OofS.date, 
                 size = vline.size, 
                 color = hline.color, 
                 linetype=vlinetype) +
      ylab(unname(TeX("$P\\[\\,C_t\\, |\\, X_{1:t}\\,\\]"))) + 
      scale_fill_manual(values = alpha(colors[[as.character(nbRegimes)]][1:nbRegimes], alpha)) +
      theme(legend.position="none") + 
      theme(plot.margin = plot.margins,
            panel.spacing = panel.margins,
            text = element_text(size = text.size),
            element_line(size = element_line),
            axis.text.x = element_text(angle=90, hjust=0.5, vjust=0.5,size = text.size),
            axis.title.x = element_blank()
            #       # axis.ticks.x = element_blank(),
            #       # panel.grid = element_blank(),
            #       # panel.border = element_blank()
      ) + 
      scale_x_date(breaks = "1 year", 
                   limits = as.Date(full.Dates[c(1,n.OofS)]), 
                   labels=date_format("%Y")) +
      scale_y_continuous(expand = c(0, 0), limits = c(-0.001,1.001))
    
    probs.filtrees[[as.character(objet.modele$nom.modele.row)]] <<- Graph.filtered
    
    tryCatch({
      # PROBABILITÉS LISSÉES
      Graph.smooth <-  
        ggplot(plot.data, aes(x=dates, y=smooth, 
                              fill=group)) + 
        geom_area() + 
        geom_hline(yintercept = 0.5, 
                   size = size, 
                   color = hline.color) +
        geom_vline(xintercept = first.OofS.date, 
                   size = vline.size, 
                   color = hline.color, 
                   linetype=vlinetype) +
        ylab(unname(TeX("$P\\[\\,C_t\\, |\\, X_{1:T}\\,\\]"))) + 
        scale_fill_manual(values = alpha(colors[[as.character(nbRegimes)]][1:nbRegimes], alpha)) +
        theme(legend.position="bottom",
              plot.margin = plot.margins,
              panel.spacing = panel.margins,
              text = element_text(size = text.size),
              legend.text=element_text(size=legend.text.size),
              element_line(size = element_line),
              axis.text.x = element_text(angle=90, hjust=0.5, vjust=0.5,size = text.size)
              # axis.title.x = element_blank(),
              # panel.grid = element_blank(),
              # panel.border = element_blank()
        ) + 
        labs(fill = "Distributions :") + 
        scale_x_date(breaks = "1 year", 
                     limits = as.Date(full.Dates[c(1,n.OofS)]), 
                     labels=date_format("%Y")) +
        scale_y_continuous(expand = c(0, 0), limits = c(-0.001,1.001))
      
      
      probs.lissees[[as.character(objet.modele$nom.modele.row)]] <<- Graph.smooth
      
    },
    warning = function(warn){
      warn.message.builder("Warning :\n",warn)
    },
    error = function(err){
      err.message.builder("Warning :\n",err)
    })
    
    if (!is.null(OofS.diag.Gamma)){
      # DIAGONALE DE LA MATRICE GAMMA
      diag.Gamma.plot <-  
        ggplot(plot.data, aes(x=dates, y=diag.Gamma, colour=group)) + 
        geom_line() + 
        xlab("") + 
        ylab("Diagonale de Gamma") + 
        geom_hline(yintercept = 0.5, 
                   size = size, 
                   color = hline.color) +
        geom_vline(xintercept = first.OofS.date, 
                   size = vline.size, 
                   color = hline.color, 
                   linetype = vlinetype) +
        scale_colour_manual(values = alpha(colors[[as.character(nbRegimes)]][1:nbRegimes], alpha)) +
        theme(legend.position="none",
              plot.margin = plot.margins,
              panel.spacing = panel.margins,
              text = element_text(size = text.size),
              legend.text=element_text(size=legend.text.size),
              element_line(size = element_line),
              axis.text.x = element_text(angle=90, hjust=0.5, vjust=0.5,size = text.size),
              axis.title.x = element_blank()
              # panel.grid = element_blank(),
              # panel.border = element_blank()
        ) + 
        scale_x_date(breaks = "1 year", 
                     limits = as.Date(full.Dates[c(1,n.OofS)]), 
                     labels=date_format("%Y")) +
        scale_y_continuous(expand = c(0, 0), limits = c(-0.001,1.001))
      
      if (add.var.exp){
        # SÉRIE EXPLICATIVE
        var.exp.plot <-  
          ggplot(plot.data, aes(x=dates, y=var.exp.DATA))+
          geom_line(size=size) + 
          geom_hline(yintercept = 0, 
                     size = size, 
                     color = hline.color) +
          geom_vline(xintercept = first.OofS.date, 
                     size = vline.size, 
                     color = vline.color, 
                     linetype=vlinetype) +
          xlab("") + 
          ylab(variable.explicative) + 
          theme(legend.position="none") +
          theme(plot.margin = plot.margins,
                panel.spacing = panel.margins,
                text = element_text(size = text.size), 
                element_line(size = element_line),
                axis.text.x = element_text(angle=90, hjust=0.5, vjust=0.5,size = text.size),
                axis.title.x = element_blank()
                # axis.title.x = element_blank(),
                # axis.ticks.x = element_blank(),
                # panel.grid = element_blank(),
                # panel.border = element_blank()
          ) + 
          scale_x_date(breaks = "1 year", 
                       limits = as.Date(full.Dates[c(1,n.OofS)]), 
                       labels=date_format("%Y")) +
          scale_y_continuous(expand = c(0, 0))
        
        # FULL PLOT
        full.plot <- plot_grid(index.plot,
                               var.exp.plot,
                               diag.Gamma.plot,
                               Graph.filtered,
                               Graph.smooth, 
                               # labels = c(index.name,variable.explicative,'FILTER', 'SMOOTH'),
                               # label_size = 9,
                               align="v",
                               rel_heights = c(1,1,1,1,1.75),
                               ncol = 1)
      } else {
        # FULL PLOT
        full.plot <- plot_grid(index.plot,
                               diag.Gamma.plot,
                               Graph.filtered,
                               Graph.smooth, 
                               # labels = c(index.name,variables.explicatives,'FILTER', 'SMOOTH'),
                               # label_size = 9,
                               align="v",
                               rel_heights = c(1,1,1,1.5),
                               ncol = 1)
      }
      
    } else {
      # FULL PLOT
      full.plot <- plot_grid(index.plot,
                             Graph.filtered,
                             Graph.smooth, 
                             # labels = c(index.name,'FILTER', 'SMOOTH'),
                             # label_size = 9,
                             align="v",
                             rel_heights = c(1,1,1.3),
                             ncol = 1)
    }
    
    titre.GGPLOT <- ggdraw() + draw_label(titre, size = text.size,fontface='bold')
    # Ajout du titre
    final.plot <- plot_grid(titre.GGPLOT, 
                            full.plot, 
                            ncol=1, 
                            rel_heights=c(0.15, 1)) # rel_heights values control title margins
    
    
    # plot(final.plot)
    
    tryCatch({
      ggsave(
        filename = paste0(graph.file.name,".pdf"),
        plot = final.plot,
        device = "pdf",
        path = where,
        scale = 1,
        width = pdf.width,
        height = pdf.height,
        units = c("cm"),
        dpi = 300,
        limitsize = TRUE
      )
    },
    error = function(err){
      err.message.builder(paste("Erreur dans graphique_probabilites.R lors de la sauvegarde du graphique\n      ",
                                graph.file.name," :"),
                          err)
    },
    warning = function(warn){
      warn.message.builder(paste("Warning dans graphique_probabilites.R lors de la sauvegarde du graphique pour :",graph.file.name,"\n      "),
                           warn)
    })
  },
  error = function (err){
    # Retourne l'erreur sinon
    err.message.builder(paste("Erreur lors de la création du graphique du filtre d'Hamilton pour :",graph.file.name,"\n"),
                        err)
  },
  warning = function(warn){
    warn.message.builder("Warning lors de la création du graphique du filtre d'Hamilton\n",
                         warn)
  })
  
  
}


