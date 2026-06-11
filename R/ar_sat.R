#' Calculate argon saturation concentration
#'
#' Computes the equilibrium dissolved argon (Ar) concentration in water at a
#' given temperature, salinity, and barometric pressure using the equation of
#' Hamme and Emerson (2004, Eqns. 1-2, Table 4). Vapor pressure of water is
#' calculated via the Antoine equation (Stull 1947) and used to apply a
#' pressure correction to the saturation concentration.
#'
#' @param t water temperature in degrees Celsius
#' @param sal salinity in per mille (ppt)
#' @param baromet.press.atm barometric pressure in atmospheres
#'
#' @return equilibrium dissolved Ar concentration in micromoles per
#'   kilogram (umol/kg)
#'
#' @references
#'
#' Hamme, R. C. & Emerson, S. R. (2004). \emph{The solubility of neon,
#' nitrogen and argon in distilled water and seawater}, Deep-Sea Research I,
#' 51(11), 1517-1528.
#'
#' Stull, D. R. (1947). \emph{Vapor Pressure of Pure Substances. Organic and
#' Inorganic Compounds.} Industrial & Engineering Chemistry, 39(4), 517-540.
#' doi: 10.1021/ie50448a022
#'
#' @examples
#' # Ar saturation at 20 degrees C, fresh water, 1 atm
#' ar_sat(t = 20, sal = 0, baromet.press.atm = 1)
#'
#' @export
ar_sat <- function(t, sal, baromet.press.atm) {
  # Antoine equation: vapor pressure of water [bar], valid -18 to 100 C
  # (Stull 1947)
  vapor.press <- exp(4.6543 - (1435.264 / ((t + 273.15) + -64.848)))
  vapor.press <- vapor.press * 0.98692 # [bar] to [atm]

  # pressure correction: (P - Pw) / (1 atm - Pw)
  press.corr <- (baromet.press.atm - vapor.press) / (1 - vapor.press)

  # Coefficients [umol/kg] (Hamme and Emerson 2004, Table 4)
  A0 <- 2.79150
  A1 <- 3.17609
  A2 <- 4.13116
  A3 <- 4.90379
  B0 <- -6.96233e-3
  B1 <- -7.66670e-3
  B2 <- -1.16888e-2

  # Scaled temperature (Hamme and Emerson 2004, Eqn. 2)
  TS <- log((298.15 - t) / (273.15 + t))
  S <- sal

  # Saturation concentration (Hamme and Emerson 2004, Eqn. 1)
  lnAr.sat <- A0 +
    A1 * TS +
    A2 * TS^2 +
    A3 * TS^3 +
    S * (B0 + B1 * TS + B2 * TS^2)
  Ar.sat <- exp(lnAr.sat)

  Ar.sat * press.corr
}
