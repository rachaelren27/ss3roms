#' Add a catch-per-unit-effort fleet
#'
#' Add a fleet that includes catch-per-unit-effort (CPUE) data
#' to input files for a Stock Synthesis model.
#'
#' @details # HARD-CODED
#' Search for the term `HARD-CODED` to find variables that are
#' defined in the code rather than available as input arguments.
#' For example, the timing of a survey is hard coded to 1.
#' Future development could specify these variables using
#' input arguments or
#' calculations.
#'
#' @param datlist A list created by [r4ss::SS_readdat()].
#' @param ctllist A list created by [r4ss::SS_readctl()].
#' @param data A matrix or data frame with the following data:
#'   * years,
#'   * season,
#'   * catch or catch-per-unit-effort data, and
#'   * uncertainty about the data.
#'   The matrix will be referenced only by its column number,
#'   not by column name, so order matters but column names do not.
#'   All values in the matrix are numeric.
#'   Currently, you can only provide catch **or** catch-per-unit-effort data.
#'   If you want both, add a fleet for one type and then augment
#'   `CPUE` or `catch` in the returned list. Be careful with this though because
#'   you might also have to change the type of fleet in `fleetinfo` of
#'   the returned list. See fleettype.
#' @param fleetname A character string supplying the name for the new fleet.
#' @param fleettype A character value that corresponds to the name in `datlist`
#'   for the object of interest.
#'   The default is to include the information as a survey, i.e., `"CPUE"`,
#'   rather than as a fishery, i.e., `"catch"`.
#' @param units An integer supplying the units type.
#' @param minageselected The minimum age that is fully selected by the fleet.
#'   The default is age-0 fish and the entry should always be an integer.
#'   This parameter will be used for age-based selectivity type 11.
#' @param maxageselected The maximum age that is fully selected by the fleet.
#'   The default is age-0 fish and the entry should always be an integer.
#'   This parameter will be used for age-based selectivity type 11.
#' @importFrom rlang :=
#' @author Kelli F. Johnson
#' @export
#' @return An invisible list of lists is returned with the following two items:
#' * datlist and
#' * ctllist.
#' @seealso
#' * [r4ss::SS_writectl()]
#' * [r4ss::SS_writedat()]
#'
#' @examples
#' \dontrun{
#' # Read in files
#' dat <- r4ss::SS_readdat(
#'   file = system.file(
#'     "extdata", "models", "PacificHake", "hake_data.ss",
#'     package = "ss3roms"
#'   ),
#'   verbose = FALSE
#' )
#' ctl <- r4ss::SS_readctl(
#'   file = system.file(
#'     "extdata", "models", "PacificHake", "hake_control.ss",
#'     package = "ss3roms"
#'   ),
#'   use_datlist = TRUE,
#'   datlist = dat,
#'   verbose = FALSE
#' )
#' data(ROMS)
#' newlists <- add_fleet(
#'    datlist = dat,
#'    ctllist = ctl,
#'    data = data.frame(
#'      year = ROMS[["year"]],
#'      seas = 7,
#'      obs = ROMS[["UWpre.a"]],
#'      se_log = 0.01
#'    ),
#'    fleetname = "Age0_upwelling",
#'    fleettype = "CPUE",
#'  )
#'  fs::dir_copy(
#'    path = system.file("extdata", "models", "PacificHake",
#'      package = "ss3roms"
#'    ),
#'    new_path = file.path("test"),
#'    overwrite = TRUE
#'  )
#'  r4ss::SS_writectl(
#'    ctllist = newlists[["ctllist"]],
#'    outfile = file.path(
#'      "test",
#'      basename(newlists[["ctllist"]][["sourcefile"]])
#'    ),
#'    overwrite = TRUE,
#'    verbose = FALSE
#'  )
#'  r4ss::SS_writedat(
#'    datlist = newlists[["datlist"]],
#'    outfile = file.path(
#'      "test",
#'      basename(newlists[["datlist"]][["sourcefile"]])
#'    ),
#'    overwrite = TRUE,
#'    verbose = FALSE
#'  )
#' }
add_fleet <- function(datlist,
                      ctllist,
                      data,
                      fleetname,
                      fleettype = c("CPUE", "catch"),
                      units = 0,
                      minageselected = 0,
                      maxageselected = 0) {
  # Input checks
  stopifnot(NCOL(data) == 4)
  # Determine some info from the arguments
  fleettype <- match.arg(fleettype, several.ok = FALSE)
  datlistfleetcolname <- switch(fleettype,
    "catch" = "fleet",
    "CPUE" = "index"
  )
  
  # Change data file
  datlist[["Nfleets"]] <- datlist[["Nfleets"]] + 1
  if (fleettype == "CPUE") {
    datlist[["Nsurveys"]] <- datlist[["Nsurveys"]] + 1
  }
  datlist[["fleetnames"]][datlist[["Nfleets"]]] <- fleetname
  datlist[["fleetinfo"]] <- dplyr::bind_rows(
    datlist[["fleetinfo"]],
    data.frame(
      # 1 is fishery, 2 is discard, 3 is survey
      "type" = ifelse(fleettype == "CPUE", 3, 1),
      # HARD-CODED surveytiming, area, and need_catch_mult
      "surveytiming" = 1,
      "area" = 1,
      "units" = units,
      "need_catch_mult" = 0,
      "fleetname" = fleetname
    )
  )
  datlist[["surveytiming"]] <- datlist[["fleetinfo"]][, "surveytiming"]
  datlist[["CPUEinfo"]][datlist[["Nfleets"]], ] <- c(
    "Fleet" = datlist[["Nfleets"]],
    "Units" = units,
    # HARD-CODED Errtype and SD_Report
    "Errtype" = 0,
    "SD_Report" = 0
  )
  row.names(datlist[["CPUEinfo"]]) <- datlist[["fleetnames"]]
  datlist[[fleettype]] <- datlist[[fleettype]] %>%
  dplyr::bind_rows(
    (data %>%
      `colnames<-`(colnames(datlist[[fleettype]])[-3]) %>%
      dplyr::mutate(!!datlistfleetcolname := 
        ifelse(obs == -999, -1, 1) * datlist[["Nfleets"]]) %>%
      dplyr::relocate(!!datlistfleetcolname, .after = "seas")
    )
  )
  # Do not turn on Dirichlet-multinomial parameter for composition data
  datlist[["len_info"]][datlist[["Nfleets"]], ] <- c(-1, 0.001, 0, 0, 0, 0, 0.001)
  row.names(datlist[["len_info"]]) <- datlist[["fleetnames"]]
  datlist[["age_info"]][datlist[["Nfleets"]], ] <- c(-1, 0.001, 0, 0, 0, 0, 0.001)
  row.names(datlist[["age_info"]]) <- datlist[["fleetnames"]]
  
  # Alter the control file to add selectivity parameters for new survey
  ctllist[["Nfleets"]] <- datlist[["Nfleets"]]
  ctllist[["fleetnames"]] <- datlist[["fleetnames"]]
  
  ctllist[["Q_options"]][datlist[["Nsurveys"]], ] <-
    ctllist[["Q_options"]][datlist[["Nsurveys"]] - 1, ]
  ctllist[["Q_options"]][datlist[["Nsurveys"]], "fleet"] <- datlist[["Nfleets"]]
  row.names(ctllist[["Q_options"]]) <-
    datlist[["fleetnames"]][ctllist[["Q_options"]][["fleet"]]]
  ctllist[["Q_parms"]][(datlist[["Nsurveys"]] * 2) - 1:0, ] <-
    ctllist[["Q_parms"]][NROW(ctllist[["Q_parms"]]) - 1:0, ]
  row.names(ctllist[["Q_parms"]])[NROW(ctllist[["Q_parms"]]) - 1:0] <- gsub(
    "\\([0-9]+\\)", 
    glue::glue('({datlist[["Nfleets"]]})'),
    row.names(ctllist[["Q_parms"]])[1:2]
  )
  # HARD-CODED Turn off length-based selectivity 
  ctllist[["size_selex_types"]][datlist[["Nfleets"]], ] <- rep(0, 4)
  row.names(ctllist[["size_selex_types"]]) <- datlist[["fleetnames"]]
  # HARD-CODED Set age-based selectivity to use
  # pattern = 11
  # discard = 0, i.e., no discarding
  # male = 0, i.e., no male-specific selectivity parameters
  # special = 0, i.e., no special selectivity information
  ctllist[["age_selex_types"]][datlist[["Nfleets"]], ] <- rep(0, 4)
  row.names(ctllist[["age_selex_types"]]) <- datlist[["fleetnames"]]
  if(units < 30) {
    # HARD-CODED: for units < 30 age-selectivity type is 11
    ctllist[["age_selex_types"]][datlist[["Nfleets"]], 1] <- 11
    ctllist[["age_selex_parms"]][NROW(ctllist[["age_selex_parms"]]) + 1:2, ] <- data.frame(
      LO = 0, HI = maxageselected + 1, INIT = c(minageselected, maxageselected),
      PRIOR = -1, PR_SD = 0.01, PR_type = 0,
      PHASE = -1,
      "env_var&link" = 0, dev_link = 0, dev_minyr = 0, dev_maxyr = 0, dev_PH = 0,
      Block = 0, Block_Fxn = 0
    ) %>%
      magrittr::set_rownames(
        glue::glue('AgeSel_P_{1:2}_',
                   '{datlist[["fleetnames"]][datlist[["Nfleets"]]]}',
                   '({datlist[["Nfleets"]]})'
        )
      )
  }  
  # return the lists that will need to be written to the disk
  invisible(list(datlist = datlist, ctllist = ctllist))
}
