#' Calculate oxygen saturation concentration
#'
#' Computes the equilibrium dissolved oxygen concentration in water at a given
#' temperature, salinity, and barometric pressure using the combined fit
#' equation of Garcia and Gordon (1992, Eqn. 8). Vapor pressure of water is
#' calculated via the Antoine equation (Stull 1947) and used to apply a
#' pressure correction to the saturation concentration.
#'
#' @param t water temperature in degrees Celsius
#' @param sal salinity in per mille (ppt)
#' @param baromet.press.atm barometric pressure in atmospheres
#'
#' @return equilibrium dissolved oxygen concentration in micromoles per
#'   kilogram (umol/kg)
#'
#' @references
#'
#' Garcia, H., and L. Gordon (1992), \emph{Oxygen solubility in seawater:
#' Better fitting equations}, Limnology and Oceanography, 37(6).
#'
#' Stull, D. R. (1947). \emph{Vapor Pressure of Pure Substances. Organic and
#' Inorganic Compounds.} Industrial & Engineering Chemistry, 39(4), 517-540.
#' doi: 10.1021/ie50448a022
#'
#' @examples
#' # O2 saturation at 20 degrees C, fresh water, 1 atm
#' o2_sat(t = 20, sal = 0, baromet.press.atm = 1)
#'
#' @export
o2_sat <- function(t, sal, baromet.press.atm) {
  # Antoine equation: vapor pressure of water [bar], valid -18 to 100 C
  # (Stull 1947)
  vapor.press <- exp(4.6543 - (1435.264 / ((t + 273.15) + -64.848)))
  vapor.press <- vapor.press * 0.98692 # [bar] to [atm]

  # pressure correction: (P - Pw) / (1 atm - Pw)
  press.corr <- (baromet.press.atm - vapor.press) / (1 - vapor.press)

  # Combined fit coefficients [umol/kg] (Garcia and Gordon 1992, Table 1)
  A0 <- 5.80818
  A1 <- 3.20684
  A2 <- 4.1189
  A3 <- 4.93845
  A4 <- 1.01567
  A5 <- 1.41575
  B0 <- -7.01211 * 10^-3
  B1 <- -7.25958 * 10^-3
  B2 <- -7.93334 * 10^-3
  B3 <- -5.54491 * 10^-3
  C0 <- -1.32412 * 10^-7

  # Scaled temperature (Garcia and Gordon 1992, Eqn. 8)
  TS <- log((298.15 - t) / (273.15 + t))
  S <- sal

  lnO2.sat <- A0 +
    A1 * TS +
    A2 * TS^2 +
    A3 * TS^2 +
    A3 * TS^3 +
    A4 * TS^4 +
    A5 * TS^5 +
    S * (B0 + B1 * TS + B2 * TS^2 + B3 * TS^3) +
    C0 * S^2
  O2.sat <- exp(lnO2.sat)

  O2.sat * press.corr
}
