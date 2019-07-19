#' Outputs of several simulations
#'
#' A dataset containing Outputs of several simulations obtained with the parameters parms in which only change
#' the herbivory pressure (livestock is 0 or 410)
#' the length of cycle (10, 20, 50 and 100 years long cycle)
#' and the proportion of fallow during cycle.
#'
#' @format A data frame with 368 rows and 47 variables:
#' \describe{
#' \item{\code{aireB}}{a numeric vector}
#' \item{\code{kapa}}{a numeric vector}
#' \item{\code{np}}{a numeric vector}
#' \item{\code{nf}}{a numeric vector}
#' \item{\code{OM_Case}}{a numeric vector}
#' \item{\code{OM_Gdnt}}{a numeric vector}
#' \item{\code{OM_Fllw}}{a numeric vector}
#' \item{\code{OD_Case}}{a numeric vector}
#' \item{\code{OD_Gdnt}}{a numeric vector}
#' \item{\code{OD_Fllw}}{a numeric vector}
#' \item{\code{IP_Gdnt}}{a numeric vector}
#' \item{\code{IP_Fllw}}{a numeric vector}
#' \item{\code{Mi_Case}}{a numeric vector}
#' \item{\code{Mi_Fllw}}{a numeric vector}
#' \item{\code{Mi_Gdnt}}{a numeric vector}
#' \item{\code{P1}}{a numeric vector}
#' \item{\code{D1}}{a numeric vector}
#' \item{\code{M1}}{a numeric vector}
#' \item{\code{P10}}{a numeric vector}
#' \item{\code{D10}}{a numeric vector}
#' \item{\code{M10}}{a numeric vector}
#' \item{\code{PC}}{a numeric vector}
#' \item{\code{DC}}{a numeric vector}
#' \item{\code{MC}}{a numeric vector}
#' \item{\code{food_in_Fllw}}{a numeric vector}
#' \item{\code{food_in_Gdnt}}{a numeric vector}
#' \item{\code{food_in_case}}{a numeric vector}
#' \item{\code{food_in_sava}}{a numeric vector}
#' \item{\code{food_ou_Fllw}}{a numeric vector}
#' \item{\code{food_ou_Gdnt}}{a numeric vector}
#' \item{\code{food_ou_case}}{a numeric vector}
#' \item{\code{harvest}}{a numeric vector}
#' \item{\code{compound_harvest}}{a numeric vector}
#' \item{\code{export_harvest}}{a numeric vector}
#' \item{\code{bush_harvest}}{a numeric vector}
#' \item{\code{frac}}{a numeric vector}
#' \item{\code{total_harv}}{a numeric vector}
#' \item{\code{prod_bush}}{a numeric vector}
#' \item{\code{prod_comp}}{a numeric vector}
#' \item{\code{prod_tot}}{a numeric vector}
#' \item{\code{rend_comp}}{a numeric vector}
#' \item{\code{rend_bush}}{a numeric vector}
#' \item{\code{rend_tot}}{a numeric vector}
#' \item{\code{Livestock}}{a factor with levels \code{0} \code{410}}
#' \item{\code{airB}}{a factor with levels \code{75 \%}}
#' \item{\code{testbush}}{a numeric vector}
#' \item{\code{id}}{a character vector}
#'}
"cycleData"

#' parameters
#'
#' A list used in modelling
#'
#' @format A list of 39 parameters for the ewoc model
#' \describe{
#' \item{\code{KC}}{a numeric}
#' \item{\code{KJ}}{a numeric}
#' \item{\code{KBmax}}{a numeric}
#' \item{\code{ipB}}{a numeric}
#' \item{\code{ipJ}}{a numeric}
#' \item{\code{bC}}{a numeric}
#' \item{\code{bB}}{a numeric}
#' \item{\code{bJ}}{a numeric}
#' \item{\code{ccC}}{a numeric}
#' \item{\code{ccB}}{a numeric}
#' \item{\code{ccJ}}{a numeric}
#' \item{\code{m}}{a numeric}
#' \item{\code{odC}}{a numeric}
#' \item{\code{odB}}{a numeric}
#' \item{\code{odJ}}{a numeric}
#' \item{\code{omC}}{a numeric}
#' \item{\code{omB}}{a numeric}
#' \item{\code{omJ}}{a numeric}
#' \item{\code{rho}}{a numeric}
#' \item{\code{kapa}}{a numeric}
#' \item{\code{fu}}{a numeric}
#' \item{\code{KK}}{a numeric}
#' \item{\code{lambdK}}{a numeric}
#' \item{\code{tau}}{a numeric}
#' \item{\code{co_case}}{a numeric}
#' \item{\code{co_chbu}}{a numeric}
#' \item{\code{co_fabu}}{a numeric}
#' \item{\code{aireB}}{a numeric}
#' \item{\code{aireT}}{a numeric}
#' \item{\code{im}}{a numeric}
#' \item{\code{id}}{a numeric}
#' \item{\code{lambd}}{a numeric}
#' \item{\code{gmC}}{a numeric}
#' \item{\code{gmB}}{a numeric}
#' \item{\code{exp}}{a numeric}
#' \item{\code{lgthss}}{a numeric}
#' \item{\code{r0}}{a numeric}
#' \item{\code{parKB}}{a numeric}
#' \item{\code{root}}{a numeric}
#' \item{\code{psi}}{a numeric}
#' \item{\code{zeta_dry}}{a numeric}
#' \item{\code{zeta_wet}}{a numeric}
#' \item{\code{gmax}}{a numeric}
#'}
"parms"

