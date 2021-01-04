# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Unconstraining parameters
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# First version : january 8th, 2020
# Last version : 5 mars 2020
# ----------------------------------------------------------------------------------------------------

# --------------------------------------------------------
# TEST DE NORMALITÉ
# --------------------------------------------------------
# Fonctions permettant de tester la normalité de l'échantillon
#     En premier lieu on trace la densité des observations,
#     puis on trace le graphique des quantiles de la loi normale
#     en comparaison avec ceux de notre échantillon.
# --------------------------------------------------------

tests.Normalite <- function(DATA.df, index,GraphPath, descriptive.stat.console=FALSE){
  # --------------------------------------------------------
  # Paramètres graphiques
  # --------------------------------------------------------
  # Nombre de points à afficher sur l'axe horizontal
  nbTicks=10
  
  # JPEG DIMENSIONS FOR OUTPUTED FILE
  image.width <- 1250
  image.heigth <- 666
  
  # Quantité pour taille et options affichage graphique
  cex.axis   <- 1.5
  cex.mtext  <- 2
  cex.legend <- 2.5
  cex.pch    <- 1

  # --------------------
  hline()
  for (i in 1:length(index)) {
    index.name <- index[i]
    # Unpacking du data.frame
    dates <- rownames(DATA.df)
    DATA <- DATA.df[,index.name]
    
    # Moyenne et Écart-type échantillonale pour densité normale univariée
    mu <- mean(DATA)
    sigma <- sd(DATA)
    # Données standardisées
    z.norm<-(DATA-mu)/sigma
    
    #descriptive statistics
    des_stats <- list("Moyenne" = mu,
                      "Écart-type" = sigma,
                      "Asymétrie" = timeDate::skewness(DATA, method="moment"),
                      "Aplatissement" = timeDate::kurtosis(DATA, method="moment"), 
                      "Minimum" = min(DATA),
                      "Maximum" = max(DATA), 
                      "n" = length(DATA)
    )
    
    # Sauvegarde sous forme de liste et de tableaux LaTeX
    # Résultat copié au presse-papier
    tab.desc <- Make.LaTeX.Table(matrix(des_stats,nrow=1),
                                 Col.Titles = names(des_stats),
                                 Cross.Lines = T, 
                                 Row.Pos = 'c',
                                 title=paste("Statistiques descriptives sur ",index.name,sep=""),
                                 print.Cons=descriptive.stat.console,
                                 copy.CB=TRUE)
    
    
    # Plot the normality tests and export resulting plot
    jpeg(paste(GraphPath,"/normality/Normality_Histogramme_",index.name,".png",sep=""), width = image.width, height = image.heigth)
    layout(matrix(c(1,1,1,1,2,3,2,3),4,2,byrow=TRUE))
    par(mar = rep(6, 4)) # Set the margin on all sides to 2
    
    # --------------------------------------------------------
    # GRAPHIQUE : Densité paramétrique, non-paramétrique et histogramme
    # --------------------------------------------------------
    borne.xlim.hist <- max(abs(min(DATA)),abs(max(DATA)))
    # Estimation de la densité des données et histogramme plus courbe normale
    hist_plot <- hist(DATA, freq=FALSE, breaks=300,
                      plot=TRUE, col="light gray", border="white",
                      xlim=c(-borne.xlim.hist,borne.xlim.hist), xaxs="i", xlab="", ylab="",main="",axes=FALSE)
    # Ajout des axes, orientation du texte et ajustement de l'écart en
    #   l'axe et les étiquettes
    axis(side=1, cex.axis=cex.axis, padj=-0.2)
    axis(side=2, cex.axis=cex.axis, las=1)
    # Ajout de la densité non-paramétrique
    lines(density(DATA), col="black", lwd=1)
    # Ajout de la courbe d'une loi normale avec paramètre empirique
    curve(dnorm(x, mean=mu, sd=sigma),
          col="black", lwd=1, lty=2, add=TRUE)
    # Ajout d'un titre sur le graphique, modification de son emplacement,
    #   de sa taille et de sa police
    mtext(paste(index.name,"\nDensité des rendements",sep=""), side=3, line=0.5, font=4, las=0, cex=cex.mtext)
    # Ajout d'une légende sur le graphique
    legend("topleft", inset=0.05, legend=c("Non-paramétrique", "Normale"),
           lty=c(1,2), lwd=1, pch=NA, pt.cex=cex.pch, col=c("black","black"),
           bty="n", cex=cex.legend, ncol=1, x.intersp=1, y.intersp=1.5)
    
    
    # --------------------------------------------------------
    # GRAPHIQUE : Comparaison des quantilles de la loi normale et de l'échantillon
    # --------------------------------------------------------
    # Comparaison des quantiles de l'échantillon par rapport à ceux d'une normale
    qqnorm(z.norm, main="") ## drawing the QQplot 
    mtext("Q-Q Plot", side=3, line=2, font=4, las=0, cex=cex.mtext)
    abline(0,1) ## drawing a 45-degree reference line
    
    
    # --------------------------------------------------------
    # GRAPHIQUE : Box Plot de l'échantillon
    # --------------------------------------------------------
    # box plot
    boxplot(DATA,
            border="black",
            lwd=1,
            varwidth = FALSE, xaxt="n", yaxt="n",xlab="",ylab="")
    title(main = "Box Plot",
          cex.main=cex.mtext*1.5)
    title(xlab="log-rendements",
          line=1)
    axis(side=2, las=1)
    title(ylab="Observations",
          line=2)
    
    # Fermeture de l'objet de graphique
    dev.off()
    cat("Graphique sauvegardé :",paste(GraphPath,"/normality/Normality_Histogramme_",index.name,".png",sep=""),"\n")
    hline()
  }

  return(list(des_stats=des_stats,
              tab.desc=tab.desc,
              hist_plot=hist_plot))
}

