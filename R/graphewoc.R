# ==========================================================================
# Program: graphewoc.R
# Description: contains some graph function utils to modele_complet_saison
# Author: Anne Bisson
# Creation: 13/01/2017
# Update: 18/01/2017 (add of graphical function for surface changes)
# Update: 30/01/2017 (bug correction)
# ==========================================================================

#' Stock of one simulation
#'
#' @param outputSimu a dataframe from one simulation of model_complet_saison
#' @param nbnT number of year simulated in model_complet season (length of colomn time/365)
#' @param np number of parcelle of bush in simulation of modele_complet season
#' @param title character title of graph
#' @param cbPalette colour choosed (particular format, see Rcolorchart)
#'
#' @return a ggplot graph format pdf
#' @export
#'
#' @examples
#' # not run
plotStock <- function(outputSimu,nbnT,np, title, cbPalette){
  labels <- c( P1 = "Plants \n (bush)", R1 = "Dead-Roots \n (Bush)",
               D1 = "Organic \n (Bush)", M1 = "Nutrients \n (Bush)",
               PC  = "Plants \n (Compound)", DC = "Organic \n (Compound)",
               MC = "Nutrients \n (Compound)")
  outputSimu <- as.data.frame(outputSimu) %>%
    dplyr::select(time,P1,R1,D1,M1,PC,DC,MC) %>%
    dplyr::filter(time > (nbnT-np)*365+1)  %>%
    reshape2::melt(id=c("time")) %>%
    dplyr::mutate(time = time/365)%>%
    dplyr::mutate(Type = substring(variable,1,1), Place = substring(variable,2,3))
  PlotGG <- ggplot2::ggplot(data=outputSimu,ggplot2::aes(x = time, y = value, colour= Place)) + ggplot2::geom_line(size=1) +
    ggplot2::xlab("time (years)") + ggplot2::ylab("N (in kgN.ha-1)") +
    ggplot2::facet_grid(variable ~  . , scales = "free_y", labeller=ggplot2::labeller(variable = labels))  +
    ggplot2::theme(legend.position='none',
          strip.text.y = ggplot2::element_text(size=13),                   # taille legend facet_grid
          axis.text=ggplot2::element_text(size=12),                        # taille label axe x
          axis.title=ggplot2::element_text(size=14),                       # taille legend axe x
          legend.title = ggplot2::element_text(size = 14),
          legend.key.size = ggplot2::unit(0.8, "cm"),
          legend.text = ggplot2::element_text(size = 14)) +
    ggplot2::scale_colour_manual(values=cbPalette) # To use for fills, add
  ggplot2::ggsave(paste0(title,".pdf"), width = 20, height = 25, units = "cm")
  return(PlotGG)
}

#' Rendement et production in function of cycle
#'
#' @param datain a dataframe
#' @param cbPalette color choosed (3)
#' @param titleGraph character title of graph
#' @param rend_chosed character rendement variable chosed
#' @param prod_chosed character production variable chosed
#'
#' @return a ggplot2 graph
#' @export
#'
#' @examples
#' # not run
plotRdtProd.cycle <-function(datain,cbPalette,titleGraph, rend_chosed, prod_chosed){
  datain2 <- dplyr::select_(datain, "np","nf", "frac", "livestock", "id", "airB",
                     rend = rend_chosed,  prod = prod_chosed) %>%
    dplyr::mutate(livestock = as.factor(livestock), id = as.factor(id), np = as.factor(np)) %>%
    reshape2::melt(id=c("np","nf","frac", "livestock","id","airB"))
  labelsFacet <- c(rend =  "Yield (in kgN/ha)", prod = "Production (in KgN)")
  # rajouter tous les paramètres qui changeront à l'appel de la fonction
  ggip <- ggplot2::ggplot(datain2) + ggplot2::aes(x = 1-frac, y = value, colour = np, ymin =0) +
    ggplot2::geom_line(size = 1, aes(group = id, linetype = livestock)) +
    ggplot2::xlim(0,1) +
    #ggtitle("Yields (in kgN/ha cultivated) and production (in kgN) in whole agro-ecosystem for different lengths of cycle and duration of fallow, with and without livestock") +
    ggplot2::xlab("Fractions of cycle with fallows") + ggplot2::ylab("") +
    ggplot2::facet_grid(variable ~ airB,  scales = "free", labeller=ggplot2::labeller(variable = labelsFacet)) +
    ggplot2::labs(colour = "Length of cycle") +
    ggplot2::theme(legend.position='bottom', legend.box = "horizontal",              # position de la legende
          strip.text.y = ggplot2::element_text(size=15),            # taille legend facet_grid
          axis.text=ggplot2::element_text(size=12),                 # taille label axe x
          axis.title=ggplot2::element_text(size=15),                # taille legend axe x
          legend.title = ggplot2::element_text(size = 14),
          legend.key.size = ggplot2::unit(0.8, "cm"),
          legend.text = ggplot2::element_text(size = 14)) +
    ggplot2::scale_colour_manual(values=cbPalette) # To use for fills
  ggplot2::ggsave(titleGraph, width = 15, height = 15, units = "cm")
  return(ggip)
}

#' Production in function of cycle and cattle
#'
#' @param datain a dataframe
#' @param cbPalette chosen color (3)
#' @param titleGraph character title of graph
#' @param prod_chosed character production variable chosed
#' @param airB_chosed factor
#'
#' @return a ggplot2 graph
#' @export
#'
#' @examples
#' # not run
plotLiveProd.cycle <-function(datain,cbPalette,titleGraph, prod_chosed, airB_chosed){
  datain2 <- dplyr::filter(datain, airB == airB_chosed) %>%
    dplyr::mutate(livestock = as.factor(paste0("L",livestock)),  np = as.factor(np)) %>%
    dplyr::select_("np","nf", "frac","livestock" , "id", prod = prod_chosed) %>%
    reshape2::melt(id=c("np","nf","frac", "livestock","id"))
  labels_bl <- c(L0 = "Production (in KgN) \n without livestock", L410 = "Production (in KgN) \n with livestock")
  ggip <- ggplot2::ggplot(datain2) + ggplot2::aes(x = 1-frac, y = value, colour = np, ymin =0) +
    ggplot2::geom_line(size = 1, aes(group = id, linetype = livestock)) +
    ggplot2::xlim(0,1) +
    #ggtitle("Yields (in kgN/ha cultivated) and production (in kgN) in whole agro-ecosystem for different lengths of cycle and duration of fallow, with and without livestock") +
    ggplot2::xlab("Fractions of cycle with fallows") + ggplot2::ylab("") +
    ggplot2::facet_grid(livestock ~ . ,  scales = "fixed" , labeller=labeller(livestock = labels_bl)) +
    ggplot2::labs(colour = "Length of cycle") +
    ggplot2::theme(legend.position='bottom', legend.box = "horizontal",              # position of legend
          strip.text.y = ggplot2::element_text(size=15),            # size legend facet_grid
          axis.text=ggplot2::element_text(size=12),                 # size label axis x
          axis.title=ggplot2::element_text(size=15),                # size legend axis x
          legend.title = ggplot2::element_text(size = 14),
          legend.key.size = ggplot2::unit(0.8, "cm"),
          legend.text = ggplot2::element_text(size = 14)) +
    ggplot2::scale_colour_manual(values=cbPalette) # To use for fills
  ggplot2::ggsave(titleGraph, width = 15, height = 15, units = "cm")
  return(ggip)
}

#' Rendement et production in function of surface of compound and bush
#'
#' @param datain a dataframe
#' @param cbPalette color choosed (3)
#' @param titleGraph character title of graph
#' @param livestock_chosed number
#'
#' @return a ggplot2 graph
#' @export
#'
#' @examples
#' # not run
plotRdtProd.surface <-function(datain,cbPalette,titleGraph, livestock_chosed){
  datain <- dplyr::filter(datain, livestock == livestock_chosed) %>%
    dplyr::select(rend_comp,  prod_comp,rend_bush, prod_bush,rend_tot, prod_tot, aireCa) %>%
    reshape2::melt(id=c("aireCa")) %>%
    dplyr::mutate(Type = substr(variable,1,4), Ring = substr(variable,6,9))
  labelsFacet <- c(rend =  "Yield (in kgN/ha)", prod = "Production (in KgN)")
  # rajouter tous les paramètres qui changeront à l'appel de la fonction
  ggip <- ggplot2::ggplot(datain) + ggplot2::aes(x = aireCa, y = value, ymin =0) +
    ggplot2::geom_line(size = 2, ggplot2::aes(colour = Ring)) +   ggplot2::xlim(0,100) +
    #ggtitle("Yield in kgN/ha and Production in kgN in Compound Ring")+
    ggplot2::xlab("% of Compound area over Bush + Compound area ") + ggplot2::ylab("") +
    ggplot2::facet_grid(Type ~ . ,  scales = "free", labeller=ggplot2::labeller(Type = labelsFacet)) +
    ggplot2::theme(legend.position='bottom',                        # position de la legende
          strip.text.y = ggplot2::element_text(size=15),            # taille legend facet_grid
          axis.text=ggplot2::element_text(size=12),                 # taille label axe x
          axis.title=ggplot2::element_text(size=15),                # taille legend axe x
          legend.title = ggplot2::element_text(size = 14),
          legend.key.size = ggplot2::unit(0.8, "cm"),
          legend.text = ggplot2::element_text(size = 14)) +
    ggplot2::scale_colour_manual(values=cbPalette)  # To use for fills
  ggplot2::ggsave(titleGraph, width = 15, height = 15, units = "cm")
  return(ggip)
}


#------------------------- End of file --------------------------------------------
