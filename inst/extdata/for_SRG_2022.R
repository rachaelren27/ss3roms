library(here)
library(ss3roms)
library(ggplot2)

### Plot of environmental drivers vs. rec devs
data(ROMS)
ROMS <- dplyr::slice(ROMS, -1)

# add fake environmental driver
set.seed(1400)

# ROMS <- ROMS %>% dplyr::mutate(test = rnorm(nrow(ROMS)))

simcor <- function (x, ymean=0, ysd=1, correlation=0) {
  n <- length(x)
  y <- rnorm(n)
  z <- correlation * scale(x)[,1] + sqrt(1 - correlation^2) * 
    scale(resid(lm(y ~ x)))[,1]
  yresult <- ymean + ysd * z
  yresult
}

temp <- r4ss::SSgetoutput(dir = here('inst/extdata/models/Sablefish'),
                          forecast = FALSE) %>% 
  r4ss::SSsummarize() 
rec.devs <- temp$recdevs

rand <- simcor(rec.devs$replist1,
               correlation = 0.25,
               ymean = mean(rec.devs$replist1),
               ysd = sd(rec.devs$replist1))

rand <- faux::rnorm_pre(rec.devs$replist1, mu = mean(rec.devs$replist1), sd = sd(rec.devs$replist1), r = 0.25)

# rand <- rec.devs.sub$replist1

temp <- r4ss::SSgetoutput(dirvec = here('inst/extdata/models/PetraleSole'))
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

ROMS.plot <- ROMS.big %>%
  tidyr::pivot_longer(cols = -c('year', 'dev'), names_to = 'indicator', values_to = 'value') %>%
  dplyr::filter(!is.na(value)) 
ROMS.cor <- ROMS.plot %>%
  dplyr::group_by(indicator) %>%
  dplyr::summarize(pearson.cor = cor(value, dev, use = 'pairwise.complete.obs'))

png('correlations.png', width = 12, height = 4, units = 'in', res = 2000)
ROMS.plot %>%
  ggplot() +
  geom_point(aes(x = value, y = dev, col = year)) +
  geom_text(x = Inf, y = Inf, hjust = 1, vjust = 1, size = 5,
            aes(label = paste('r =', round(pearson.cor, 2))), 
            data = ROMS.cor) +
  facet_wrap(~indicator, scales = 'free_x', nrow = 1) +
  labs(x = 'Environmental index', y = 'Recruitment deviation') +
  # ggsidekick::theme_sleek() +
  scale_x_continuous(n.breaks = 4) +
  fishualize::scale_color_fish(option = 'Ostracion_whitleyi')
  #scale_color_gradientn(colors = calecopal('chaparral3',5))
dev.off()

### Retrospectives looking at including env driver

# Adjust control and data files
dat <- r4ss::SS_readdat(
  file = here("inst", "extdata", "models", "Sablefish", "data.ss"),
  verbose = FALSE
)

ctl <- r4ss::SS_readctl(
  file = here("inst", "extdata", "models", "Sablefish", "control.ss"),
  use_datlist = TRUE,
  datlist = dat,
  verbose = FALSE,
)

newlists <- add_fleet(
  datlist = dat,
  ctllist = ctl,
  data = data.frame(
    year = rec.devs$Yr,
    seas = 7,
    obs = exp(rand),
    se_log = 0.01
  ),
  fleetname = "env",
  fleettype = "CPUE", 
  units = 31
)
dirname <- 'test_rand'
                  
r4ss::copy_SS_inputs(
  dir.old = here("inst", "extdata", "models", "Sablefish"),
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

# Run SS in MLE mode
r4ss::run_SS_models(dirvec = here('inst/extdata/models/Sablefish'), 
                    skipfinished = FALSE, 
                    model = here('inst/extdata/models/Sablefish/ss'))

# Do retrospectives
peel <- 15

r4ss::SS_doRetro(
  masterdir = here('inst/extdata/models'),
  oldsubdir = 'Sablefish',
  newsubdir = "s_retrospectives",
  years = -peel,
  extras = '-nohess'
)

r4ss::SS_doRetro(masterdir = here(dirname),
                 oldsubdir = '', 
                 years = -peel,
                 newsubdir = "s_retrospectives",
                 extras = '-nohess'
)

# # --- comparing retrospectives -------------------------------------------------
# # read in output
# retroModels <- r4ss::SSgetoutput(dirvec = here(dirname, "retrospectives",
#                                          paste0("retro-", 10:15)),
#                            forecast = FALSE) %>% 
#   r4ss::SSsummarize()
# 
# # set the ending year of each model in the set
# endyrvec <- 2020 - (10:15)
# # make comparison plot
# compare_retro <- r4ss::SSplotComparisons(retroModels, endyrvec = endyrvec, new = FALSE)
# 
# # make Squid Plot of recdev retrospectives
# par(mfrow = c(2, 1))
# # first scaled relative to most recent estimate
# scaled_squid <- r4ss::SSplotRetroRecruits(retroModels,
#                     endyrvec = endyrvec, cohorts = 2005:2010,
#                     relative = TRUE, legend = FALSE
# )

 
temp <- r4ss::SSgetoutput(dirvec = c(here('inst/extdata/models', 's_retrospectives', paste0('retro-', peel)),
                                     here(dirname, 's_retrospectives', paste0('retro-', peel)),
                                   here('inst/extdata/models/Sablefish')),
                          forecast = FALSE) %>%
  r4ss::SSsummarize()

rec.devs <- temp$recdevs

# --- calculate CI -------------------------------------------------------------
# rec.devs.lower <- temp$recdevsLower
# rec.devs.lower[rec.devs.lower$Yr > 2015-peel, 2] <- NA
# lowerCI <- rec.devs.lower$replist2
#   
# rec.devs.upper <- temp$recdevsUpper
# rec.devs.upper[rec.devs.upper$Yr > 2015-peel, 2] <- NA
# upperCI <- rec.devs.upper$replist2
# 
# rec.devs <- cbind(rec.devs, lowerCI, upperCI)


# names(rec.devs)[1:6] <- c(paste0('retro', 10:15))
# names(rec.devs)[7] <- 'age1.2021'
# for(ii in 1:nrow(rec.devs)) {
#   for(jj in 1:6) {
#     if(!is.na(rec.devs$retro14[ii])){
#       if(rec.devs$retro14[ii] > 2020-(10:15)[jj]) {
#         rec.devs[ii,jj] <- NA
#       }
#     }
#   }
# }

# --- plot rec devs ------------------------------------------------------------
names(rec.devs)[1:3] <- c(paste('Age 1', 2020 - peel, 'retro'),
                          paste('Env', 2020 - peel, 'retro'),
                          '2020 Age 1')
# names(rec.devs)[1] <- paste('Age 1', 2020 - peel, 'retro')
rec.devs[rec.devs$Yr > 2020-peel, 1:2] <- NA
plot.dat <- rec.devs %>%
  tidyr::pivot_longer(cols = 1:3, names_to = 'model', values_to = 'rec.dev') %>%
  dplyr::filter(
                Yr > 1990 & Yr <= 2020) %>%
  dplyr::mutate(model = factor(model))

# plot.dat <- plot.dat %>% dplyr::filter(model != '2020 Age 1')

big.plot <- plot.dat %>%
  ggplot() +
  ylim(-4, 4) + 
  geom_line(aes(x = Yr, y = rec.dev, col = model)) +
  # scale_color_manual(values=c("#00BA38","#619CFF")) + 
  # geom_line(aes(x = Yr, y = age1.2021)) +
  # ggsidekick::theme_sleek(18) +
  # geom_ribbon(aes(x = Yr,
  #                 y = rec.dev,
  #                 ymin = lowerCI,
  #                 ymax = upperCI),
  #             alpha = 0.1) + 
  labs(x = 'Year', y = 'Recruitment deviation') +
  geom_hline(yintercept = 0, col = 'red') +
  # scale_color_manual(values = inauguration::inauguration('inauguration_2021', 3), drop = FALSE) +
  NULL



# png('assess_skill.png', width = 11, height = 4.5, units = 'in', res = 2000)
assess.skill <- big.plot +
  geom_point(aes(x = Yr, y = rec.dev, col = model), cex = 2,
             data = dplyr::filter(plot.dat, Yr==2020 - peel))
# dev.off()

png('assess_skill0.png', width = 11, height = 4.5, units = 'in', res = 2000)
big.plot %+%
  dplyr::filter(plot.dat, model == '2021 Age 1') +
  geom_point(aes(x = Yr, y = rec.dev, col = model), cex = 2,
             data = dplyr::filter(plot.dat, Yr == 2021 - peel, model == '2021 Age 1'))
  
dev.off()

temp %>%
  r4ss::SSsummarize() %>%
  r4ss::SSplotRetroRecruits(cohorts = 2000:2005, endyrvec = 2021 - (5:15))



qplot(x = year, y = bifurcation_index, geom = 'line', data = bifurcation)


