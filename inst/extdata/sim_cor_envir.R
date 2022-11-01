library(here)
library(ss3roms)
library(ggplot2)
library(r4ss)

data(ROMS)
ROMS <- dplyr::slice(ROMS, -1)

dat <- SS_readdat(
  file = system.file(
    "extdata", "models", "PacificHake", "hake_data.ss",
    package = "ss3roms"
  ),
  verbose = FALSE
)

ctl <- SS_readctl(
  file = system.file(
    "extdata", "models", "PacificHake", "hake_control.ss",
    package = "ss3roms"
  ),
  use_datlist = TRUE,
  datlist = dat,
  verbose = FALSE,
  version = 3.30
)

# --- get recruitment deviations from base model -------------------------------
peel <- 14
term.year <- 2020 - peel

temp <- SSgetoutput(dir = c(here('inst/extdata/models/PacificHake'),
                                  here('inst/extdata/models', 
                                       'retrospectives', 
                                       paste0('retro-',peel))),
                          forecast = FALSE) %>% 
          SSsummarize() 

base.rec.devs <- temp$recdevs
base.rec.devs.sub <- base.rec.devs %>% dplyr::filter(Yr >= 1981 & Yr <= 2010)

# base terminal year rec dev
base.rec.dev <- base.rec.devs %>% dplyr::filter(Yr == term.year) %>%
  dplyr::pull(replist1)

# base retro terminal year rec dev
base.err <- abs(base.rec.devs %>% dplyr::filter(Yr == term.year) %>%
                  dplyr::pull(replist2) - base.rec.dev)


# --- create folders -----------------------------------------------------------
# naming convention: paste0('test_rand', <correlation index>, '-', <seed index>)
# ex: 'test_rand1-4' correlation level #1, random seed #4
create_folder <- function(dir.name) {
    dir.create(here(dir.name))
    
    # copy ss.exe file
    list.of.files <- list.files(here('inst/extdata/models',
                                     'retrospectives', 
                                     paste0('retro-', peel)),
                                pattern = "ss.exe")
    # file.create(here(dir.name, 'ss.exe'))
    file.copy(from = here('inst/extdata/models', 
                          'retrospectives', 
                          paste0('retro-', peel), 
                          list.of.files),
              to = here(dir.name, list.of.files), overwrite = TRUE)
}


# --- generate simulated random data ---------------------------------------------
simcor <- function (x, ymean=0, ysd=0.5, correlation=0) {
  n <- length(x)
  y <- rnorm(n)
  z <- correlation * scale(x)[,1] + sqrt(1 - correlation^2) *
    scale(resid(lm(y ~ x)))[,1]
  yresult <- ymean + ysd * z
  yresult
}

simautocorr <- function (rec.devs, ymean=0, ysd=0.5, corr=0) {
    rho <- arima(base.rec.devs$replist1, order = c(1,0,0))$coef[[1]]
}

# --- fit retrospectives to correlated environmental data ----------------------
sim_fit_retro <- function(seed.ind, corr, corr.ind){
  s <- seed.ind*46
  set.seed(s)
  
  ROMS <- ROMS %>% dplyr::mutate(rand = simcor(base.rec.devs.sub$replist1,
                                                   correlation = corr,
                                               ymean = mean(base.rec.devs$replist1),
                                               ysd = sd(base.rec.devs$replist1)))

  newlists <- add_fleet(
    datlist = dat,
    ctllist = ctl,
    data = data.frame(
      year = ROMS[["year"]],
      seas = 7,
      obs = exp(ROMS$rand),
      se_log = 0.05
    ),
    fleetname = "env",
    fleettype = "CPUE",
    units = 31
  )

  dirname <- paste0('test_rand', corr.ind, '-', seed.ind)
  create_folder(dirname)

  copy_SS_inputs(
    dir.old = system.file("extdata", "models", "PacificHake",
                          package = "ss3roms"
    ),
    dir.new = file.path(here(dirname)),
    overwrite = TRUE
  )

  SS_writectl(
    ctllist = newlists[["ctllist"]],
    outfile = file.path(
      here(dirname),
      basename(newlists[["ctllist"]][["sourcefile"]])
    ),
    overwrite = TRUE,
    verbose = FALSE
  )
  SS_writedat(
    datlist = newlists[["datlist"]],
    outfile = file.path(
      here(dirname),
      basename(newlists[["datlist"]][["sourcefile"]])
    ),
    overwrite = TRUE,
    verbose = FALSE
  )

  # fit retrospective
  SS_doRetro(masterdir = here(dirname),
                   oldsubdir = '',
                   years = -peel,
                   extras = "-nohess"
  )
}

# --- calculate environmentally-linked errors ----------------------------------
future::plan("multisession", workers = 11)

  base.errs <- c()
  env.errs <- c()
  avg.se <- c()
  dirs <- c()
  num.seed <- 50
  ind <- 1
  

  for (corr in c(0, 0.25, 0.5, 0.75, 0.9)) {
    1:num.seed %>%
      furrr::future_map(sim_fit_retro,
                        corr = corr,
                        corr.ind = ind,
                        .options = furrr::furrr_options(seed = T))

    for(s in 1:num.seed) {
      dirs[s] <- here(paste0('test_rand', ind, '-', s),
                      'retrospectives', paste0('retro-', peel))
    }

    temp <- SSgetoutput(dirvec = dirs,
                              forecast = FALSE) %>%
      SSsummarize()

    rec.devs <- temp$recdevs
    # added.se <- temp$indices %>% dplyr::filter(Fleet == 4 & Yr == term.year)

    # get environmentally-linked recruitment deviation errors and avg se
    for(s in 1:num.seed) {
      env.errs[(ind - 1)*num.seed + s] <- abs(rec.devs %>% 
                                                dplyr::filter(Yr == term.year) %>%
                           dplyr::pull(paste0('replist', s)) - base.rec.dev)
      
      # avg.se[ind] <- mean(added.se$SE - added.se$SE_input)
    }
    
    ind <- ind + 1
  }

  
# --- plot errors --------------------------------------------------------------
corrs <- rep(c(0, 0.25,0.5,0.75,0.9), each = num.seed)
errs <- as.data.frame(cbind(env.errs, corrs))

errs.plot <- ggplot(data = errs, aes(x = corrs, y = env.errs, group = corrs)) + 
  geom_boxplot() +
  xlab("correlation") + 
  ylab("absolute error") + 
  ylim(0,2) +
  geom_hline(yintercept = base.err, color = "red")
