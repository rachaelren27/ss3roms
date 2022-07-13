library(here)
library(ss3roms)
library(ggplot2)

### Plot of environmental drivers vs. rec devs
data(ROMS)
ROMS <- dplyr::slice(ROMS, -1)

# add fake environmental drivers
set.seed(1500)
# ROMS <- ROMS %>% dplyr::mutate(test = rnorm(nrow(ROMS)))

simcor <- function (x, ymean=0, ysd=1, correlation=0) {
  n <- length(x)
  y <- rnorm(n)
  z <- correlation * scale(x)[,1] + sqrt(1 - correlation^2) *
    scale(resid(lm(y ~ x)))[,1]
  yresult <- ymean + ysd * z
  yresult
}

rec.devs.sub <- rec.devs %>% dplyr::filter(Yr >= 1981 & Yr <= 2010)
ROMS <- ROMS %>% dplyr::mutate(rand0.25 = simcor(rec.devs.sub$replist3,
                                                 correlation = 0.25))

temp <- r4ss::SSgetoutput(dirvec = here('inst/extdata/models/PacificHake'))
bifurcation <- readr::read_csv(here('data-raw', 'bifurcation_index.csv'), comment = '#') %>%
  dplyr::mutate(year = year + 1) # bifurcation index impacts preconditioning phase, recruitment year - 1
ROMS.big <- dplyr::slice(ROMS, -1) %>%
  dplyr::left_join(temp$replist1$recruit[,c('Yr', 'dev')], by = c('year' = 'Yr')) %>%
  dplyr::right_join(bifurcation) %>%
  dplyr::arrange(year) %>%
  dplyr::rename(`Bifurcation index precond.` = bifurcation_index,
                `Eddy kinetic energy\nprecond.` = EKEpre.ms.c,
                `Longshore transport yolk` = LSTyolk,
                `Upwelling precond.` = UWpre.a)

# Adjust control and data files
dat <- r4ss::SS_readdat(
  file = system.file(
    "extdata", "models", "PacificHake", "hake_data.ss",
    package = "ss3roms"
  ),
  verbose = FALSE
)
ctl <- r4ss::SS_readctl(
  file = system.file(
    "extdata", "models", "PacificHake", "hake_control.ss",
    package = "ss3roms"
  ),
  use_datlist = TRUE,
  datlist = dat,
  verbose = FALSE,
  version = 3.30
)

newlists <- add_fleet(
  datlist = dat,
  ctllist = ctl,
  data = data.frame(
    year = ROMS[["year"]],
    seas = 7,
    obs = exp(ROMS$rand0.25),
    se_log = 0.01
  ),
  fleetname = "env",
  fleettype = "CPUE",
  units = 31
)
dirname <- 'test_rand2'

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


base.err <- c()
env.err <- c()
j <- 1
for(peel in 5:15){
    # get retrospectives
    term.year <- 2021 - peel
    
    temp <- r4ss::SSgetoutput(dirvec = c(here('inst/extdata/models', 'retrospectives', paste0('retro-', peel)),
                                         here('test_rand1', 'retrospectives', paste0('retro-', peel)),
                                         here('inst/extdata/models/PacificHake')),
                              forecast = FALSE) %>%
      r4ss::SSsummarize()
    
    rec.devs <- temp$recdevs
    
    # get base terminal year recruitment deviation
    base.err[j] <- abs(rec.devs %>% dplyr::filter(Yr == term.year) %>% 
      dplyr::pull(replist1) - rec.devs %>% dplyr::filter(Yr == term.year) %>% 
      dplyr::pull(replist3))
    
    # get environmentally-linked terminal year recruitment deviation
    env.err[j] <- abs(rec.devs %>% dplyr::filter(Yr == term.year) %>% 
      dplyr::pull(replist2) - rec.devs %>% dplyr::filter(Yr == term.year) %>% 
        dplyr::pull(replist3))
    
    j <- j + 1
}

# # graph base vs environmentally-linked terminal data by number of years peeled
# term.years <- 2021 - 5:15
# term.err <- as.data.frame(cbind(term.years, base.err, env.err))
# 
# comp.term.err <- ggplot(data = term.err, aes(x = term.years)) +
#   geom_line(aes(y = base.err, col = "Base Error")) +
#   geom_line(aes(y = env.err, col = "Envir Error")) + 
#   xlab("Terminal Year") +
#   ylab("Recruitment Deviation Error")

# # graph recruitment deviation by year
# names(rec.devs)[1:3] <- c(paste('Age 1', term.year, 'retro'),
#                           paste('Env', term.year, 'retro'),
#                           '2021 Age 1')
# rec.devs[rec.devs$Yr > 2021-peel, 1:2] <- NA
# plot.dat <- rec.devs %>%
#   tidyr::pivot_longer(cols = 1:3, names_to = 'model', values_to = 'rec.dev') %>%
#   dplyr::filter(!grepl('Early|Late|Fore', Label),
#                 Yr > 1990) %>%
#   dplyr::mutate(model = factor(model))
# 
# big.plot <- plot.dat %>%
#   ggplot() +
#   geom_line(aes(x = Yr, y = rec.dev, col = model)) +
#   # geom_line(aes(x = Yr, y = age1.2021)) +
#   # ggsidekick::theme_sleek(18) +
#   labs(x = 'Year', y = 'Recruitment deviation') +
#   geom_hline(yintercept = 0, col = 'red') +
#   # scale_color_manual(values = inauguration::inauguration('inauguration_2021', 3), drop = FALSE) +
#   NULL
# 
# png('assess_skill.png', width = 11, height = 4.5, units = 'in', res = 2000)
# big.plot +
#   geom_point(aes(x = Yr, y = rec.dev, col = model), cex = 2,
#              data = dplyr::filter(plot.dat, Yr==term.year))
# dev.off()
# 
# rand <- simcor(rec.devs.sub$replist3, correlation = 0.25)
# ROMS <- as.data.frame(cbind(ROMS, rand))
# colnames(ROMS)[ncol(ROMS)] <- 'rand0.25'