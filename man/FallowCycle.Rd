% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelewoc.R
\name{FallowCycle}
\alias{FallowCycle}
\title{modele_complet_season in function of different configuration of the agro-ecosystem
(length of cycle and duration of fallow)}
\usage{
FallowCycle(nbc, livestock, aireBush, parms, nbnT)
}
\arguments{
\item{nbc}{vector, length of cycle tested : C'est aussi le nb de parcelles dans l'aureole de brousse,}

\item{livestock}{vector, livestock charge tested in TLU}

\item{aireBush}{vector, aire of bush, compound area is then calculated with bush area.}

\item{parms, }{list of parms used for model_complet_season}

\item{nbnT, }{number, nb of years simulated}
}
\value{
dataframe usable by graphic function plotRdtProd.cycle
}
\description{
modele_complet_season in function of different configuration of the agro-ecosystem
(length of cycle and duration of fallow)
}
\examples{
# not run
}
