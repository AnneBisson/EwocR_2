# ==========================================================================
# Program: modelewoc.R
# Description: contains model for Bisson et al 2017
# Author:Anne Bisson
# Creation: 20/01/2017
# Update: 30/01/2017 bug correction and add function for figure of the paper
# Update: 14/01/2019 correction of the model and add function for figure in Bisson and al 2019

# ==========================================================================
# SAISON DE VEGETATION: second membre du systeme d'equations
# ==========================================================================

#' systems of differential equations modeling the behavior
#' of the agro-ecosystem during the growing season
#'
#' @param t numeric vector, duration of simulation
#' @param x numeric vector, initial  condition of state variables
#' @param parametre list, list of initialization of parameters models
#' @details xxxx
#'
#' @return list, contains vector and dataframe
#' @export
#'
#' @examples
#' # not run
GrowingSeason<-function(t,x,parametre) {
  with(as.list(c(x,parametre)), {

    # Initialisation
    # **************
    dx                  = rep(0,(length(x)));

    ODBf <- ODBc <- OMBc <- OMBf <- ODC <- OMC  <- IPc <- IPf <- food_in_fa <-
      food_in_crop <- food_in_case <- food_ou_fa <- food_ou_crop <- food_ou_case <- MCB <- MCf <- 0;

    #  Bush subsystem
    # ******************
    # sumPfb = sum of stocks of nitrogen in plants of follow subunits of the bush subsystem
    # seq(1,4*(np-1)+1,4) is the vector of index of Pi variables, i.e. plants of each plots of bush ring
    if (nf==0){
      sumPfb  = 0
    } else {
      sumPfb              = (1/nf)*sum(x[seq(1,4*(np-1)+1,4)]*Fall)
    }

    for(i in 1:np) { # loop on plots de l'aureole de Brousse

      # We get the values of the state from which the derivatives will be calculated
      # second member of the equation system
      P                   = x[4*(i-1)+1];  # Plants
      R                   = x[4*(i-1)+2];  # Woody Roots (variable D)
      D                   = x[4*(i-1)+3];  # Organic Matter (variable O)
      M                   = x[4*(i-1)+4];  # Inorganic Matter(variable I)

      if(nf==0){
        phi_j_wet =0
      } else {
        phi_j_wet = zeta_wet*gmax*((co_fabu*(1/nf)*P)/(KK + co_fabu*sumPfb))*kapa*psi
      }


      if(Fall[i]==1){# fallow plot
        # Correspondance with notations of the paper:
        # bJ=u3; KJ=K3; P=P3; M=M3; ipJ=ip3; aireB/np=\alpha_3; co_fabu=b3; KK=KB;
        # R=D3; r0=d; D=O3; m=m3; odJ=od3; fu=\nu; tau=h; lambdK=\lambda_k; omJ=om3;
        # derivation of P
        dx[4*(i-1)+1]       = bJ*(1- P/KJ)*P*M + ipJ*(1- P/KJ)*P - (np/aireB)*phi_j_wet;
        # derivation of R
        dx[4*(i-1)+2]       = - r0*R;
        # derivation of D
        dx[4*(i-1)+3]       = - m*D + id- odJ*D + r0*R + (1-fu)*lambdK*(tau*(np/aireB)*phi_j_wet + (1-tau)*kapa*rho*(np/(aireB*nf)));
        # derivation of M
        dx[4*(i-1)+4]       = im + m*D - bJ*(1 - P/KJ)*M*P - omJ*M + fu*lambdK*(tau*(np/aireB)*phi_j_wet + (1-tau)*kapa*rho*(np/(aireB*nf)));

        # others outpus
        food_in_fa          = food_in_fa + phi_j_wet;
        food_ou_fa          = food_ou_fa + lambdK*(tau*phi_j_wet) + lambdK*(1-tau)*kapa*rho*(1/nf)
        ODBf                = ODBf +  odJ*D;
        OMBf                = OMBf +  omJ*M;
        IPf                 = IPf + ipJ*(1- P/KJ)*P;
        MCf                 = MCf + m*D;

      } else{# crop plot
        # Correspondance with notations of the paper:
        # bB=u2; K[i]=K2; P=P2; M=M2; ipB=ip2; r0=d; R=D2; m=m2; D=O2; odB=od2; omB=om2;
        # derivation of P
        dx[4*(i-1)+1]       = bB*(1 - P/K[i])*M*P + ipB*P*(1 - P/K[i]);
        # derivation of R
        dx[4*(i-1)+2]       = -r0*R;
        # derivation of D
        dx[4*(i-1)+3]       = r0*R + id - m*D - odB*D;
        # derivation of M
        dx[4*(i-1)+4]       = im + m*D - bB*(1 - P/K[i])*M*P - omB*M;

        # others outputs
        ODBc                = ODBc +  odB*D;
        OMBc                = OMBc +  omB*M;
        IPc                 = IPc + ipB*(1 - P/K[i])*P;
        MCB                 = MCB + m*D;
      }
    }

    # Compound ring
    # ******************
    # We get the values of the state on the halo of the box, values from which will be calculated the
    # Derivatives (the second member of the equation system)
    PC                  = x[np*4+1];
    DC                  = x[np*4+2];
    MC                  = x[np*4+3];
    # Correspondance with notations of the paper:
    # bC=u1; KC=K1; PC=P1; omC=om1; MC=M1;
    # m=m1; DC=O1; odC=od1; lambd=\lambda_V; aireC=\alpha_1; V0=V((n-1)T+\tau)
    # derivation of P
    dx[np*4+1]          = bC*(1 - PC/KC)*MC*PC;
    # derivation of D
    dx[np*4+2]          = - m*DC - odC*DC + id + (lambd*V0/365)/aireC;
    # derivation of M
    dx[np*4+3]          = m*DC - bC*(1 - PC/KC)*MC*PC - omC*MC + im;

    # other outputs
    # **************
    dx[np*4+4]          = omC*MC; #OMC
    dx[np*4+5]          = OMBc/(np-nf); #OMBc
    dx[np*4+6]          = OMBf/nf; #OMBf
    dx[np*4+7]          = odC*DC; #ODC
    dx[np*4+8]          = ODBc/(np-nf); #ODBc
    dx[np*4+9]          = ODBf/nf; #OMBf
    dx[np*4+10]         = IPf/nf; #fixation Fllw
    dx[np*4+11]         = IPc/(np-nf); #fixation crop
    dx[np*4+12]         = food_in_fa #fallow intake cattle
    dx[np*4+13]         = 0 # compound intake cattle
    dx[np*4+14]         = 0 # bush filed intake cattle
    dx[np*4+15]         = kapa*rho - food_in_fa #savanna
    dx[np*4+16]         = m*DC  #miCA
    dx[np*4+17]         = MCf/max(nf,1) #mif
    dx[np*4+18]         = MCB/(np-nf) #mic
    dx[np*4+19]         = food_ou_fa; #fallow dung
    dx[np*4+20]         = 0 #compound dung
    dx[np*4+21]         = 0 #bush field dung
    return(list(dx))
  })
}


# ==========================================================================
#           SAISON SECHE: second membre du systeme d'equations
# ==========================================================================
#' systems of differential equations modeling the behavior of the agro-ecosystem during the dry season
#'
#' @param t numeric vector, duration of simulation
#' @param x numeric vector, initial  condition of state variables
#' @param parametre list, list of initialization of parameters 'models
#' @details xxxx
#'
#' @return list, contains dx
#' @export
#'
#' @examples
#' # not run
#' # test<-DrySeason(t=,x=,parametre=)
DrySeason<-function(t,x,parametre) {
  with(as.list(c(x,parametre)), {

    # Initialisation
    # **************
    dx                  = rep(0,(length(x)));

    OMBc <- OMBf <- food_in_fa <- food_in_crop <- sumPfb  <- food_in_case <- food_ou_fa <- food_ou_crop <- food_ou_case <- MCB <- MCf <- 0;

    # Aureole de Brousse
    # ******************
    # sumPfb = Sum of Nitrogen in Fallow Plants of Bush ring
    # sumPfc = Sum of Nitrogen in Crop Plants of Bush ring
    # seq(1,4*(np-1)+1,4) est le vecteur des indices des variables Pi, c'est à dire les plantes,
    # de chaque parcelle de l'aureole de Brousse
    if (nf==0){
      sumPfb            = 0
    } else {
      sumPfb            = (1/nf)*sum(x[seq(1,4*(np-1)+1,4)]*Fall)
    }

    if (np-nf==0){
      sumPfc            = 0
    } else {
      sumPfc            = (1/(np-nf))*sum(x[seq(1,4*(np-1)+1,4)]*as.numeric(Fall== 0))
    }


    Ptot = ((nf/np)*aireB*co_fabu*sumPfb + (1-nf/np)*aireB*co_chbu*sumPfc + aireC*co_case*x[np*4+1])/(aireB+aireC)

    for(i in 1:np) {# loop on each subunits of the bush subsystem

      # On recupere les valeurs de l'etat à partir desquelles vont etre calculees les derivees (le second membre du système
      # d'equations)
      P                   = x[4*(i-1)+1];  # Plantes
      R                   = x[4*(i-1)+2];  # Residus (variable R du document pdf)
      D                   = x[4*(i-1)+3];  # Detritus (variable O du document pdf)
      M                   = x[4*(i-1)+4];  # azote Mineral (variable I du document pdf)


      if(Fall[i]==1){# parcelle en fallow

        phi_bf_dry = zeta_dry*gmax*(((aireB*co_fabu*P)/(aireB+aireC))/(KK + Ptot))*kapa*psi

        # Correspondance with notations of the paper:
        # ccJ=c3; aireB/np=\alpha_3; co_chbu=b2; co_case=b1; co_fabu=b3; P=P3; R=D3; D=O3; r0=d; m=m3; fu=\nu; tau=h; KK=KB;
        # derivation of P
        dx[4*(i-1)+1]       = - ccJ*P  - (np/aireB)*phi_bf_dry;
        # derivation of R
        dx[4*(i-1)+2]       = -r0*R;
        # derivation of D
        dx[4*(i-1)+3]       =  ccJ*P + id  - m*D + r0*R + lambdK*(1-fu)*tau*(np/aireB)*phi_bf_dry;
        # derivation of M
        dx[4*(i-1)+4]       =  im + m*D - omJ*M + lambdK*fu*tau*(np/aireB)*phi_bf_dry;

        # !!! pas verif
        food_in_fa          = food_in_fa + phi_bf_dry;
        food_ou_fa          = food_ou_fa + tau*lambdK*phi_bf_dry;
        MCf                 = MCf + m*D;
        OMBf                = OMBf +  omJ*M;

      }else{#parcelle en culture

        phi_bc_dry = zeta_dry*gmax*(((aireB*co_chbu*P)/(aireB+aireC))/(KK + Ptot))*kapa*psi

        # Correspondance with notations of the paper:
        # ccB=c2; P=P2; R=D2; D=O2; aireB/np=\alpha_2; co_chbu=b2; co_case=b1; co_fabu=b3; r0=d; m=m2; fu=nu; tau=h; KK=KB;
        # derivation of P
        dx[4*(i-1)+1]       = - ccB*P - (np/aireB)*phi_bc_dry;
        # derivation of R
        dx[4*(i-1)+2]       = - r0*R;
        # derivation of D
        dx[4*(i-1)+3]       = ccB*P + id - m*D  + r0*R + lambdK*(1-fu)*tau*(np/aireB)*phi_bc_dry;
        # derivation of M
        dx[4*(i-1)+4]       = im + m*D - omB*M + lambdK*fu*tau*(np/aireB)*phi_bc_dry;

        # !!! pas verif
        food_in_crop        = food_in_crop + phi_bc_dry;
        food_ou_crop        = food_ou_crop + tau*lambdK*phi_bc_dry;
        MCB                 = MCB + m*D;
        OMBc                = OMBc +  omB*M;
      }
    }
    # Compound ring
    # ******************
    # We recover the state values on the Compound ring, values from which the derivatives will be calculated
    # (the second member of the system of equations)
    PC                  = x[np*4+1];
    DC                  = x[np*4+2];
    MC                  = x[np*4+3];

    phi_c_dry = zeta_dry*gmax*(((aireC*co_case*PC)/(aireB+aireC))/(KK + Ptot))*kapa*psi
    # Correspondance with notations of the paper:
    # ccC=c1; aireC=\alpha_1; co_case=b1; KK=KB; co_fabu=b3 ; co_chbu=b2 ; m=m1; PC=P1; DC=O1; lambd=\lambda_V
    # fu=\nu; tau=h
    # derivation of P
    dx[np*4+1]          = - ccC*PC - (1/aireC)*phi_c_dry;
    # derivation of D
    dx[np*4+2]          = ccC*PC  - m*DC + id  + (lambd*V0)/(aireC*365) + lambdK*(1-fu)*( tau*(1/aireC)*phi_c_dry + (1-tau)*(kapa*rho/aireC));
    # derivation of M
    dx[np*4+3]          = m*DC + im - omC*MC + lambdK*fu*(tau*(1/aireC)*phi_c_dry + (1-tau)*(kapa*rho/aireC));

    # Others Outputs
    # **************
    dx[np*4+4]          = omC*MC; #OMC
    dx[np*4+5]          = OMBc/(np-nf); #OMBc
    dx[np*4+6]          = OMBf/nf; #OMBf
    dx[np*4+7]          = 0; #ODC
    dx[np*4+8]          = 0; #ODBc
    dx[np*4+9]          = 0; #OMBf
    dx[np*4+10]         = 0; #fixation crop
    dx[np*4+11]         = 0; #fixation Fllw
    dx[np*4+12]         = food_in_fa; #
    dx[np*4+13]         = food_in_crop ; #mean P in bush crop
    dx[np*4+14]         = phi_c_dry
    dx[np*4+15]         = kapa*rho - food_in_fa - food_in_crop - phi_c_dry;  # sava cattle
    dx[np*4+16]         = m*DC  #miCA
    dx[np*4+17]         = MCf/nf #miF
    dx[np*4+18]         = MCB/(np-nf)  # miC
    dx[np*4+19]         = food_ou_fa; #fallow dung
    dx[np*4+20]         = food_ou_crop; #bush fields dung
    dx[np*4+21]         = lambdK*(tau*phi_c_dry + (1-tau)*(kapa*rho)); #compound dung
    if (is.nan(x[1])==T){break}
    return(list(dx));
  })
}


#' Function modeling the behavior of the agro-ecosystem during nbnT years
#'
#' @param nbnT duration of simulation (in years)
#' @param np length of a cycle (in years)
#' @param nf duration of fallow (in years) with nf < np
#' @param parms num
#' @details It take as argument the numbers of years to simulate (nbnT),
#'    the number of parcels in bush ring (np),
#'    the number of parcels in Fllws in bush ring (nf)
#'    the list of parms (see the end of document)
#'
#' @return list
#' @export
#'
#' @examples
#' # not run
CompleteModelSeason <- function(nbnT, np,nf, parms) {

  # Creation d'une liste de parametres pour la simulation
  # *****************************************************
  parms2              <- list()
  # nnumper of year per cycle = number of plot in bush ring
  parms2[["np"]]      <- np
  # Number of year of fallow in a cycle = number of fallow plot in the same time.
  # Constant number from a year to an other
  parms2[["nf"]]      <- nf
  # vecteur des capacite d'accueil de chaque parcelle de l'aureole de Brousse: on l'initiale à la valeur de capacite d'accueil
  # maximale qui est parms[["KBmax"]]
  parms2[["K"]]       <- rep(parms[["KBmax"]],np)
  # vecteur d'indice caracterisant l'etat de chaque parcelle de brousse: 0 pour "en culture", et 1 pour "en fallow"
  # c'est un vecteur de meme taille que le nombre de parcelle dans l'aureole de brousse, c'est à dire de longueur np
  parms2[["Fall"]]    <-  c(rep(0,np-nf),rep(1,nf))
  # Harvest that will become household waste
  parms2[["V0"]]      <- 0
  # aire de l'Compound ring: aire totale - aire de l'aureole de Brousse
  parms2[["aireC"]]   <- parms[["aireT"]] - parms[["aireB"]]

  # value of the initial condition of P at the beginning of the growing season
  sigma               <- 1

  # initial condition for the simulation
  # **************************************
  # initialisation: there are np*4+21 variables d'etat dans le modele:
  # - np parcelles de brousse avec 4 variables d'etat pour chaque parcelle:
  #       P --> plantes
  #       R --> residus, le D du document pdf
  #       D --> detritus, le O du document pdf
  #       M --> azote minerale
  # - 1 parcelle pour l'Compound ring avec 3 variables d'etat
  #
  # On initialise au valeurs suivantes:
  # P=sigma, R = 0, D = 1000, M = 0
  n0                  = rep(0,np*4+21);
  # plot of bush ring
  for(i in 1:np) {
    n0[4*(i-1)+1]       = sigma;
    n0[4*(i-1)+2]       = 0;
    n0[4*(i-1)+3]       = 1000;
    n0[4*(i-1)+4]       = 0;}
  # plot of Compound ring
  n0[np*4+1]          = sigma;
  n0[np*4+2]          = 1000;
  n0[np*4+3]          = 0;

  # initialisation des sorties de la fonction
  # dans output on stockera les variables d'etat
  # Dans V_harvest, on stockera la valeur annuelle de la recolte stockee dans le village
  # (recolte issue de la brousse=fallow+culture, et de la case)
  # Dans v_harvest_compound, on stockera la partie de la recolte stockee dans le village qui est issue de l'Compound ring uniquement
  # Dans V_exp, on stockera la partie de la recolte issue de la Brousse cultivee qui est exportee en dehors du village

  V_harvest <- V_harvest_compound <- V_exp <- output <- C_clearing <-  NULL

  # Simulation of the nbnT years
  # ******************************
  for(j in 0:(nbnT-1)) { # loop sur les years de simulation. Comme on commence à 0, on endi a nbnT-1.

    # an = numero du dernier jour de l'annee j-1:
    # l'annee 0 commence au jour 1 et endit au jour 365
    # l'annee 1 commence au jour 366 et endit au jour 2*365
    # ... etc
    # l'annee j commence au jour 365*j+1 et endit au jour 365*(j+1)
    an                  = 365*j

    ## growing season
    # *************************************
    # tsaisV = vecteur de temps correspondant à la saison de vegetation de l'annee j.
    # La saison de vegetation de l'annee j commence au jour 365*j+1 et endi au jour 365*j+lgthss, lgthss etant
    # la longueur (en jours) de la saison de vegetation au cours d'une annee
    tsaisV              = c((an+1):(an+parms[["lgthss"]]));  # le pas de temps du vecteur est le jour
    # tsaisV=seq((an+1),(an+parms[["lgthss"]]), by = 0.1);  # le pas de temps du vecteur est le dixième de jour

    # simulation du système pendant la saison de vegetation de l'annee j
    # n0 est la condition initiale
    # tsaisV est le vecteur temps
    # GrowingSeason est la fonction de second membre du système d'equations pendant la saison de vegetation
    # c(parms,parms2) sont les paramètres du modèle
    # on utilise une schema numerique de Runge Kutta d'ordre 4
    solV                = deSolve::ode(n0,tsaisV,GrowingSeason,c(parms,parms2), method = "rk4");

    ## Recoltes
    # *************************************
    # on recupère les valeurs des variables d'etat en end de saison de vegetation (à an+parms[["lgthss"])
    # attention, la premiere colonne de solV est le temps.
    lgsv                <- as.vector(solV[dim(solV)[1],])

    # VV = recolte stockee dans le village (recolte issue de la brousse=fallow+culture, et de la case)
    # VVC = partie de la recolte stockee dans le village qui est issue de l'Compound ring
    # VVE = partie de la recolte issue de la Brousse cultivee qui est exportee en dehors du village
    # CC  = clearing, partie perdue lors du passage jachere ->  crop
    VV <- VVC <- VVE <- CC <- 0

    # initialisation du vecteur de conditions initiales CI2 pour la dry season
    CI2                 <- rep(0,np*4+21);
    # Sur chacune des aureoles, c'est la même chose:
    # - sur le pool de plantes, on enlève ce qui est recolte.  parms2[["gmB"]] et parms2[["gmC"]] etant les
    # pourcentages de recolte on plots de brousse cultivee et de case
    # - sur les autres pools, la recolte ne change rien: on reprend donc les valeurs de la end de la saison de vegetation

    # aureole de Brousse
    for(i in 1:np) { # loop on plots de l'aureole de Brousse

      if(parms2[["Fall"]][i]==1){ # Parcelles de Brousse en jachere
        # sur le pool de plantes, on enlève ce qui est recolte.
        # remarque: there is un decalage d'un indice 1 entre CI2 et lgsv car le premier element de lgsv est le temps
        CI2[4*(i-1)+1]      = lgsv[4*(i-1)+2] ;
        # sur les autres pools, la recolte ne change rien: on reprend donc les valeurs de la end de la saison de vegetation
        CI2[4*(i-1)+2]      = lgsv[4*(i-1)+3];
        CI2[4*(i-1)+3]      = lgsv[4*(i-1)+4];
        CI2[4*(i-1)+4]      = lgsv[4*(i-1)+5];


      }else{ # Parcelles de Brousse en culture
        # sur le pool de plantes, on enlève ce qui est recolte. parms2[["gmB"]] etant le pourcentage de recolte sur les
        # parcelles de brousse en culture
        # remarque: there is un decalage d'un indice 1 entre CI2 et lgsv car le premier element de lgsv est le temps
        CI2[4*(i-1)+1]      = (1 - parms[["gmB"]]) * lgsv[4*(i-1)+2];
        # sur les autres pools, la recolte ne change rien: on reprend donc les valeurs de la end de la saison de vegetation
        CI2[4*(i-1)+2]      = lgsv[4*(i-1)+3];
        CI2[4*(i-1)+3]      = lgsv[4*(i-1)+4];
        CI2[4*(i-1)+4]      = lgsv[4*(i-1)+5];
        # partie de la recolte (fallow+culture+case) qui est stockee dans le village
        # parms[["exp"]] est le pourcentage de la recolte qui est exportee en dehors du village
        VV                  = VV +  (parms[["aireB"]]/np)*parms[["gmB"]]*(1-parms[["exp"]])*lgsv[4*(i-1)+2];
        # partie de la recolte issue de la brousse cultivee qui est exportee en dehors du village
        VVE                 = VVE + (parms[["aireB"]]/np)*parms[["gmB"]]*(parms[["exp"]])*lgsv[4*(i-1)+2];
      }
    }
    # Compound ring
    # sur le pool de plantes, on enlève ce qui est recolté. parms2[["gmC"]] etant le pourcentage de recolte sur les
    # parcelles de case
    # remarque: there is un decalage d'un indice 1 entre CI2 et lgsv car le premier element de lgsv est le temps
    CI2[np*4+1]         = (1 - parms[["gmC"]])*lgsv[np*4+2];

    # Pour toutes les autres variables d'etat, la recolte ne change rien: on reprend les valeurs de la end de la saison
    # de vegetation
    # remarque: there is un decalage d'un indice 1 entre CI2 et lgsv car le premier element de lgsv est le temps
    for (indice in 2:21){
      CI2[np*4+indice]    = lgsv[np*4+indice+1];
    }

    # partie de la recolte (fallow+culture+case) qui est stockee dans le village
    VV                  = VV + parms[["gmC"]]*lgsv[np*4+2]*parms2[["aireC"]] ;
    parms2[["V0"]]      = VV

    # artie de la recolte stockee dans le village qui est issue de l'Compound ring
    VVC                 = parms[["gmC"]]*lgsv[np*4+2]*parms2[["aireC"]];

    ##  dry season
    # *************************************
    # tsaisS = vecteur de temps correspondant à la dry season de l'annee j.
    # La dry season de l'annee j commence au jour 365*j+lgthss+1 et endi au jour 365*(j+1), lgthss etant
    # la longueur (en jours) de la saison de vegetation au cours d'une annee
    tsaisS              = c((an+parms[["lgthss"]]+1):(an+365));   # le pas de temps du vecteur est le jour
    # tsaisS=seq((an+parms[["lgthss"]]+1),(an+365), by = 0.1);    # le pas de temps du vecteur est le dixième de jour

    # simulation du système pendant la saison de vegetation de l'annee j
    # CI2 est la condition initiale
    # tsaisS est le vecteur temps
    # DrySeason est la fonction de second membre du système d'equations pendant la dry season
    # c(parms,parms2) sont les paramètres du modèle
    # on utilise une schema numerique de Runge Kutta d'ordre 4
    solS                = deSolve::ode(CI2,tsaisS,DrySeason,c(parms,parms2), method = "rk4");

    ## seedling/rotate
    # *************************************
    # on recupère les valeurs des variables d'etat en end de saison de sèche (à an+365)
    lgss                <- as.vector(solS[dim(solS)[1],]) # values at (n+1)T
    # attention, la premiere colonne de solS est le temps.

    # update du vecteur d'indice caracterisant l'etat de chaque parcelle de brousse:
    # 0 pour "en culture", et 1 pour "en fallow"
    # on decale tous les indices de 1 de sorte qu'une et une seule parcelle va passer de la culture à la fallow
    # et une et une seule parcelle va passer de la fallow à la culture
    Fall                <- parms2[["Fall"]]
    jach_temp           <- c(Fall[length(Fall)],Fall[1:length(Fall)-1])
    parms2[["Fall"]]    <- jach_temp

    # initialisation du vecteur de conditions initiales CII pour la saison de vegetation
    CII                 <- rep(0,np*4+21);
    # remarque: there is un decalage d'un indice 1 entre CII et lgss car le premier element de lgss est le temps

    # aureole de brousse
    for(i in 1:np) { # loop on the subunits of the bush subsystem

      # update of the carrying capacity of each cropland subunit of the bush subsystem depending on the initial stock
      # of I which is given by  lgss[4*(i-1)+5]
      parms2[["K"]][i]    <- min(parms[["parKB"]]*lgss[4*(i-1)+5],parms[["KBmax"]]) ## new capacity for each bush field

      if(parms2[["Fall"]][i]==1){# parcelle en fallow
        # Concerning the stock of plants in fallows, we get back the value at the end of the dry season
        # except in the case where the stock is too weak, in this case we make a sowing
        CII[4*(i-1)+1]      = max(lgss[4*(i-1)+2],sigma/100);
        # the other pools do not change: the values of the end of the dry season are used
        CII[4*(i-1)+2]      = lgss[4*(i-1)+3];
        CII[4*(i-1)+3]      = lgss[4*(i-1)+4];
        CII[4*(i-1)+4]      = lgss[4*(i-1)+5];

      } else {# cropland subunits
        if (jach_temp[i] - Fall[i] < 0 ){ # rotation: the subunit was in state of fallow the year before
          # semis pour les plantes
          CII[4*(i-1)+1]      = sigma;
          # on rajoute le stock de racines de la jachere au stock de residus: pendant la jachere ces racines sont comprises
          # dans le pool de plantes P
          CII[4*(i-1)+2]      = lgss[4*(i-1)+3]+ parms[["root"]]*lgss[4*(i-1)+2];
          # les autres pools ne changent pas: on reprend les valeurs de la end de la dry season
          CII[4*(i-1)+3]      = lgss[4*(i-1)+4];
          CII[4*(i-1)+4]      = lgss[4*(i-1)+5];
          CC                  = (1 - parms[["root"]])*lgss[4*(i-1)+2]*(parms[["aireB"]]/np)

        } else { # parcelle qui etait dejà en culture l'annee d'avant
          # semis pour les plantes
          CII[4*(i-1)+1]      = sigma;
          # les autres pools ne changent pas: on reprend les valeurs de la end de la dry season
          CII[4*(i-1)+2]      = lgss[4*(i-1)+3];
          CII[4*(i-1)+3]      = lgss[4*(i-1)+4];
          CII[4*(i-1)+4]      = lgss[4*(i-1)+5];
        }
      }
    }
    # Compound ring
    # semis pour les plantes
    CII[np*4+1]         = sigma;
    # les autres pools ne changent pas: on reprend les valeurs de la end de la dry season
    CII[np*4+2]         = lgss[np*4+3];
    CII[np*4+3]         = lgss[np*4+4];

    # update de la condition initiale
    n0                  <- CII

    # stockage des sorties
    output              <- rbind(output,rbind(solV,solS))

    # stockage des recoltes
    # VV = recolte stockee dans le village (recolte issue de la brousse=fallow+culture, et de la case)
    # VVC = partie de la recolte stockee dans le village qui est issue de l'Compound ring
    # VVE = partie de la recolte issue de la Brousse cultivee qui est exportee en dehors du village
    V_harvest           <- c(V_harvest,VV)
    V_harvest_compound  <- c(V_harvest_compound,VVC)
    V_exp               <- c(V_exp, VVE)
    C_clearing          <- c(C_clearing,CC)
  }
  colnames(output) <-  c("time", c(paste0(rep(c("P","R","D","M"),np),rep(1:np, each = 4))),
                         "PC","DC","MC","OM_Case","OM_Gdnt","OM_Fllw","OD_Case","OD_Gdnt","OD_Fllw","IP_Fllw","IP_Gdnt",
                         "food_in_Fllw","food_in_Gdnt","food_in_case","food_in_sava",
                         "Mi_Case", "Mi_Fllw", "Mi_Gdnt","food_ou_Fllw","food_ou_Gdnt","food_ou_case")
  return(list(output = output, harvest = V_harvest, compound_harvest = V_harvest_compound, exp_harvest = V_exp, Clearing = C_clearing))
}


# ==========================================================================
#     Fonction function_AS
# ==========================================================================
#' function_AS: cette fonction renvoie des sorties ponctuelles du modele definies a la fin du code
#'
#' @param param_AS parametres dont on fait l'analyse de sensibilité (matrix nbcol = nb de parametres, nbrow = nb de jeux de parametres)
#' @param npp longueur de rotation
#' @param nff durée en jachère (nff<npp)
#' @param nbnTt nombre d'années simulées
#'
#' @return sorties d'interet pour chaque jeu de parametres
#' @export
#'
#' @examples
#' # not run
function_AS <- function(param_AS,npp,nff,nbnTt){
  # CONDITIONS DE SIMULATION
  # nb d'years simulees
  nbnT = nbnTt ## nb years simulees
  np <- npp
  nf <- nff
  # initialisation pour la sauvegarde des sorties ponctuelles pour chaque jeu de parametres
  sorties <- matrix(0, nrow=1, ncol=31)
  # loop des scenarios du design AS
  for (i in 1:nrow(param_AS)) {
    tim <- paste0(i, "/", 1)
    print(tim)
    # STRUCTURE & PARAMETRES DU MODELE
    # redefinitions des parametres qui changent
    for (j in 1:length(parms)){
      parms[[j]] = param_AS[i,j]
    }

    # SIMULATIONS
    stock_fluxes <- model_complet_saison(nbnT,np,nf,parms)
    #attach(stock_fluxes)
    check <- as.data.frame(stock_fluxes$output) %>%  filter(time > (nbnT-np)*365) %>%
      mutate(time = time - max(time) + 365) %>%   filter(time == 365) %>%
      select(OM_Case,OM_Gdnt,OM_Fllw, OD_Case,OD_Gdnt,OD_Fllw,
             IP_Gdnt,IP_Fllw, Mi_Case,Mi_Fllw,Mi_Gdnt,
             P1,D1,M1,P10,D10,M10,PC,DC,MC,
             food_in_Fllw,food_in_Gdnt,food_in_case,food_in_sava,
             food_ou_Fllw,food_ou_Gdnt,food_ou_case)
    Vec_interet <- c(check,stock_fluxes$harvest[nbnT],stock_fluxes$compound_harvest[nbnT],
                     stock_fluxes$exp_harvest[nbnT],
                     stock_fluxes$harvest[nbnT]-stock_fluxes$compound_harvest[nbnT]+stock_fluxes$exp_harvest[nbnT])
    #detach(stock_fluxes)
    # Sorties d'interet
    sorties[i,]  <- as.numeric(Vec_interet) ;
    head(sorties)
  }# end loop scenarios AS
  return(sorties)
} # end fonction du modele


#' modele_complet_season in function of different configuration of the agro-ecosystem
#' (length of cycle and duration of fallow)
#'
#' @param nbc vector, length of cycle tested : C'est aussi le nb de parcelles dans l'aureole de brousse,
# ce qui permet d'echelonner les mises en fallow et en culture. Each year, one plot is mise
# en jachere (passage culture --> fallow) et une autre qui est mise en culture (passage fallow --> culture)
#' @param livestock vector, livestock charge tested in TLU
#' @param aireBush vector, aire of bush, compound area is then calculated with bush area.
#' @param parms, list of parms used for model_complet_season
#' @param nbnT, number, nb of years simulated
#'
#' @return dataframe usable by graphic function plotRdtProd.cycle
#' @export
#'
#' @examples
#' # not run
FallowCycle <- function(nbc,livestock,aireBush,parms,nbnT){
  print("Warning ! This function is likely to last several hours or even days")
  outCycle <- NULL
  for (ba in 1:length(aireBush)) {
    parms[["aireB"]] <- aireBush[ba]
    for (bet in 1:length(livestock)) {
      parms[["kapa"]] <- livestock[bet]
      for (lg in 1:length(nbc)){
        np <- nbc[lg]
        for(nf in 0:np){ ## fraction de l'aureole de brousse en jachere
          print(paste0("nf/np = ",nf,"/",np))
          stock_Fluxes <- CompleteModelSeason(nbnT,np,nf,parms)
          check <- as.data.frame(stock_Fluxes$output) %>%  dplyr::filter(time > (nbnT-np)*365) %>%
            dplyr::mutate(time = time - max(time) + 365) %>%   dplyr::filter(time == 365) %>%
            dplyr::select(OM_Case,OM_Gdnt,OM_Fllw, OD_Case,OD_Gdnt,OD_Fllw,
                          IP_Gdnt,IP_Fllw, Mi_Case,Mi_Fllw,Mi_Gdnt,
                          P1,D1,M1,P10,D10,M10,PC,DC,MC,
                          food_in_Fllw,food_in_Gdnt,food_in_case,food_in_sava,
                          food_ou_Fllw,food_ou_Gdnt,food_ou_case)
          Vec_interet <-
            c(check,stock_Fluxes$harvest[nbnT],stock_Fluxes$compound_harvest[nbnT],stock_Fluxes$exp_harvest[nbnT],
              stock_Fluxes$harvest[nbnT]-stock_Fluxes$compound_harvest[nbnT]+stock_Fluxes$exp_harvest[nbnT])
          # Sorties d'interet
          outCycle <- rbind(outCycle,c(parms[["aireB"]],parms[["kapa"]] ,np,nf, as.numeric(Vec_interet)))
        }
      }
    }
  }
  colnames(outCycle) <-
    c("aireB","livestock","np","nf","OM_Case", "OM_Gdnt", "OM_Fllw","OD_Case", "OD_Gdnt", "OD_Fllw",
      "IP_Gdnt", "IP_Fllw","Mi_Case", "Mi_Fllw", "Mi_Gdnt", "P1" , "D1", "M1","P10" , "D10", "M10",
      "PC", "DC", "MC", "food_in_Fllw", "food_in_Gdnt", "food_in_case", "food_in_sava",
      "food_ou_Fllw", "food_ou_Gdnt", "food_ou_case",
      "harvest", "compound_harvest", "export_harvest", "bush_harvest");
  outCycle[is.na(outCycle)] <- 0
  outCycle <- dplyr::mutate(as.data.frame(outCycle),
                            frac       = (np-nf)/np,
                            prod_bush  = bush_harvest/parms[["gmB"]],
                            prod_comp  = compound_harvest/parms[["gmC"]],
                            prod_tot   = prod_bush+prod_comp,
                            total_harv = bush_harvest + compound_harvest,
                            rend_comp  = prod_comp /(parms[["aireT"]] - aireB),
                            rend_bush  = prod_bush /(aireB*frac),
                            rend_tot   = (prod_comp +prod_bush)/(parms[["aireT"]] - (1-frac)*aireB),
                            airB       = as.factor(paste0("%",(aireB/parms[["aireT"]])*100)),
                            aireB      = as.factor(aireB),
                            livestock  = as.factor(livestock),
                            np         = as.factor(np),
                            nf         = as.factor(nf),
                            id         = paste0(livestock,np))
  return(outCycle)
}



#'  CompleteModelSeason in function of different configuration of the agro-ecosystem
#'  (surface of Bush:Compound Ratio)
#'
#' @param np length of a cycle (in years)
#' @param nf duration of fallow (in years) with nf < np
#' @param livestock herbivory pressure (in Tropical livestock Units)
#' @param aireBush surface of bush ring
#' @param parms parms
#' @param nbnT duration of simulation (in years)
#'
#' @return dataframe usable by graphic function plotRdtProd.surface
#' @export
#'
#' @examples
#' # not run
SurfaceVariation <- function(np,nf,livestock,aireBush,parms,nbnT){
  print("Warning ! This function is likely to last several hours or even days")
  OutSurface <- NULL
  for (ba in 1:length(aireBush)) {
    parms[["aireB"]] = aireBush[ba]
    print("aireB",ba)
    for (bet in 1:length(livestock)) {
      parms[["kapa"]] = livestock[bet]
      stock_fluxes <- CompleteModelSeason(nbnT,np,nf,parms)
      check <- as.data.frame(stock_fluxes$output) %>%  dplyr::filter(time > (nbnT-np)*365) %>%
        dplyr::mutate(time = time - max(time) + 365) %>%   dplyr::filter(time == 365) %>%
        dplyr::select(OM_Case,OM_Gdnt,OM_Fllw, OD_Case,OD_Gdnt,OD_Fllw,
                      IP_Gdnt,IP_Fllw, Mi_Case,Mi_Fllw,Mi_Gdnt,
                      P1,D1,M1,P10,D10,M10,PC,DC,MC,
                      food_in_Fllw,food_in_Gdnt,food_in_case,food_in_sava,
                      food_ou_Fllw,food_ou_Gdnt,food_ou_case)
      Vec_interet <- c(check,stock_fluxes$harvest[nbnT],stock_fluxes$compound_harvest[nbnT],stock_fluxes$exp_harvest[nbnT],
                       stock_fluxes$harvest[nbnT]-stock_fluxes$compound_harvest[nbnT]+stock_fluxes$exp_harvest[nbnT])
      # Sorties d'interet
      OutSurface <- rbind(OutSurface,c(parms[["aireB"]],parms[["kapa"]] ,np,nf, as.numeric(Vec_interet)))

    }
  }
  colnames(OutSurface) <- c("aireB","livestock","np","nf",
                            "OM_Case", "OM_Gdnt", "OM_Fllw", "OD_Case", "OD_Gdnt", "OD_Fllw",
                            "IP_Gdnt", "IP_Fllw", "Mi_Case", "Mi_Fllw", "Mi_Gdnt",
                            "P1" , "D1", "M1","P10" , "D10", "M10", "PC", "DC", "MC",
                            "food_in_Fllw", "food_in_Gdnt", "food_in_case", "food_in_sava",
                            "food_ou_Fllw", "food_ou_Gdnt", "food_ou_case",
                            "harvest", "compound_harvest", "export_harvest", "bush_harvest");
  OutSurface[is.na(OutSurface)] <- 0
  OutSurface <- dplyr::mutate(as.data.frame(OutSurface),
                              frac       = (np-nf)/np,
                              total_harv = bush_harvest + compound_harvest,
                              prod_bush  = bush_harvest/parms[["gmB"]],
                              prod_comp  = compound_harvest/parms[["gmC"]],
                              prod_tot   = prod_bush+prod_comp,
                              rend_comp  = prod_comp /(parms[["aireT"]] - aireB),
                              rend_bush  = prod_bush /(aireB*frac),
                              rend_tot   = (prod_comp +prod_bush)/(parms[["aireT"]] - (1-frac)*aireB),
                              livestock  = as.factor(livestock),
                              airB       = paste0("%",(aireB/parms[["aireT"]])*100),
                              testbush   = prod_bush*25/((aireB/parms[["aireT"]])*100),
                              airB       = as.factor(airB),
                              aireBu     = (aireB/parms[["aireT"]])*100,
                              aireCa     = 100 - aireBu)
  return(OutSurface)
}



#------------------------- End of file --------------------------------------------

