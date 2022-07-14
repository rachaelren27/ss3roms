library(here)
library(ss3roms)
library(ggplot2)

data(ROMS)
ROMS <- dplyr::slice(ROMS, -1)

# --- get recruitment deviations from base model -------------------------------
temp <- r4ss::SSgetoutput(dir = c(here('inst/extdata/models/PacificHake'),
                                  here('inst/extdata/models', 'retrospectives', 'retro-15')),
                          forecast = FALSE) %>% 
          r4ss::SSsummarize() 

rec.devs <- temp$recdevs
rec.devs.sub <- rec.devs %>% dplyr::filter(Yr >= 1981 & Yr <= 2010)

peel <- 15
term.year <- 2021 - peel

# base terminal year rec dev
base.rec.dev <- rec.devs %>% dplyr::filter(Yr == term.year) %>%
  dplyr::pull(replist1)

# base retro terminal year rec dev
base.err <- abs(rec.devs %>% dplyr::filter(Yr == term.year) %>%
                  dplyr::pull(replist2) - base.rec.dev)


# --- create simulated random data ---------------------------------------------
simcor <- function (x, ymean=0, ysd=1, correlation=0) {
  n <- length(x)
  y <- rnorm(n)
  z <- correlation * scale(x)[,1] + sqrt(1 - correlation^2) *
    scale(resid(lm(y ~ x)))[,1]
  yresult <- ymean + ysd * z
  yresult
}


# --- fit retrospectives to correlated environmental data ----------------------
sim_fit_retro <- function(i, corr, corr_ind){
  s <- i*46
  set.seed(s)
  
  ROMS <- ROMS %>% dplyr::mutate(rand = simcor(rec.devs.sub$replist3,
                                                   correlation = corr))
  
  newlists <- add_fleet(
    datlist = dat,
    ctllist = ctl,
    data = data.frame(
      year = ROMS[["year"]],
      seas = 7,
      obs = exp(ROMS$rand),
      se_log = 0.01
    ),
    fleetname = "env",
    fleettype = "CPUE", 
    units = 31
  )
  
  dirname <- paste0('test_rand', corr_ind, '-', i)
  
  r4ss::copy_SS_inputs(
    dir.old = system.file("extdata", "models", "PacificHake",
                          package = "ss3roms"
    ),
    dir.new = file.path(here(dirname)),
    overwrite = TRUE
  )
  
  r4ss::SS_writectl(
    ctllist = newlists[["ctllist"]],
    outfile = file.path(
      here(dirname),
      basename(newlists[["ctllist"]][["sourcefile"]])
    ),
    overwrite = TRUE,
    verbose = FALSE
  )
  r4ss::SS_writedat(
    datlist = newlists[["datlist"]],
    outfile = file.path(
      here(dirname),
      basename(newlists[["datlist"]][["sourcefile"]])
    ),
    overwrite = TRUE,
    verbose = FALSE
  )
  
  # fit retrospective
  r4ss::SS_doRetro(masterdir = here(dirname),
                   oldsubdir = '', 
                   years = -peel
  )
}

# --- calculate environmentally-linked errors ----------------------------------
future::plan("multisession", workers = 12)

  base.errs <- c()
  env.errs <- c()
  dirs <- c()
  ind <- 1
  for (corr in c(0.25, 0.5, 0.75, 0.9)) {
    # 1:10 %>%
    #   furrr::future_map(sim_fit_retro,
    #                     corr = corr,
    #                     corr_ind = ind,
    #                     .options = furrr::furrr_options(seed = T))

    for(s in 1:10) {
      dirs[s] <- here(paste0('test_rand', ind, '-', s),
                      'retrospectives', 'retro-15')
    }
  
    temp <- r4ss::SSgetoutput(dirvec = dirs,
                              forecast = FALSE) %>%
      r4ss::SSsummarize()
    
    rec.devs <- temp$recdevs
    
    # get environmentally-linked recruitment deviation errors
    for(i in 1:10) {
      env.errs[(ind - 1)*10 + i] <- abs(rec.devs %>% dplyr::filter(Yr == term.year) %>%
                           dplyr::pull(paste0('replist', i)) - base.rec.dev)
    }
    
    ind <- ind + 1
  }
  
# --- plot errors --------------------------------------------------------------
corrs <- rep(c(0.25,0.5,0.75,0.9), each = 10)
errs <- as.data.frame(cbind(env.errs, corrs))

errs.plot <- ggplot(data = errs, aes(x = corrs, y = env.errs, group = corrs)) + 
  geom_boxplot() + 
  xlab("correlation") + 
  ylab("error") + 
  geom_hline(yintercept = base.err, color = "red")