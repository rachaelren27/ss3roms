#' Region Ocean Modeling System (ROMS) data
#'
#' A data set provided by Kristin N. Marshall based on
#' ROMS output for the California Current.
#'
#' @format A data frame with one row per year and
#'   a column for every variable plus an initial column noting the `year`.
#'   Values of `-999` are used instead of NA because Stock Synthesis will
#'   predict a value for a given year if provided data and has a negative fleet.
#'   These initial years with missing data for upwelling and eddy kinetic energy
#'   are because they represent the preconditioning period.
#' \describe{
#'   \item{year}{annual time step}
#'   \item{UWpre.a}{upwelling}
#'   \item{LSTyolk}{long-shore transport}
#'   \item{EKEpre.ms.c}{eddy kinetic energy}
#' }
"ROMS"
