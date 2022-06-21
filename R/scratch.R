library(here)
library(ss3roms)
library(ggplot2)
library(future)

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
data(ROMS)
ROMS <- dplyr::slice(ROMS, -1)
temp <- r4ss::SSgetoutput(dirvec = 'inst/extdata/models/PacificHake') 
bifurcation <- readr::read_csv(here('data-raw', 'bifurcation_index.csv'), comment = '#') %>%
   dplyr::mutate(year = year + 1) # bifurcation index impacts preconditioning phase, recruitment year - 1
ROMS.big <- dplyr::slice(ROMS, -1) %>%
   dplyr::left_join(temp$replist$recruit[,c('Yr', 'dev')], by = c('year' = 'Yr')) %>%
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
   ggsidekick::theme_sleek(18) +
   scale_x_continuous(n.breaks = 4)
dev.off()

newlists <- add_fleet(
   datlist = dat,
   ctllist = ctl,
   data = data.frame(
      year = ROMS[["year"]],
      seas = 7,
      obs = exp(-ROMS$EKEpre.ms.c),
#      obs = exp(ROMS[["EKEpre.ms.c"]]),
#      obs = exp(1.1 * ROMS$dev),
      se_log = 0.01
   ),
   fleetname = "env",
   fleettype = "CPUE", 
   units = 31
)
dirname <- 'test_negEKE'
#dirname <- 'test_posEKE'

fs::dir_copy(
   path = system.file("extdata", "models", "PacificHake",
                      package = "ss3roms"
   ),
   new_path = file.path(dirname),
   overwrite = TRUE
)

r4ss::copy_SS_inputs(
   dir.old = system.file("extdata", "models", "PacificHake",
                                           package = "ss3roms"
   ),
   dir.new = file.path(dirname),
   overwrite = TRUE
)

r4ss::SS_writectl(
   ctllist = newlists[["ctllist"]],
   outfile = file.path(
      dirname,
      basename(newlists[["ctllist"]][["sourcefile"]])
   ),
   overwrite = TRUE,
   verbose = FALSE
)
r4ss::SS_writedat(
   datlist = newlists[["datlist"]],
   outfile = file.path(
      dirname,
      basename(newlists[["datlist"]][["sourcefile"]])
   ),
   overwrite = TRUE,
   verbose = FALSE
)

plan(multisession, workers = 3)
#future_map(
purrr::map(
   .x = 'test_negEKE', #c('test_negEKE', 'test_posEKE'), 
   .f = r4ss::run_SS_models, 
   model = 'inst/extdata/bin/Windows64/ss', 
   skipfinished = FALSE
)

r4ss::run_SS_models(dirvec = 'test_negEKE', skipfinished = FALSE, model = 'inst/extdata/bin/Windows64/ss')

r4ss::run_SS_models(dirvec = c('test_negEKE', 'test_posEKE'),
                    model = 'inst/extdata/bin/Windows64/ss',
                    skipfinished = FALSE)

temp <- r4ss::SSgetoutput(dirvec = c('test_negEKE', 'inst/extdata/models/PacificHake')) %>%
   r4ss::SSsummarize() %>% 
   r4ss::SSplotComparisons(subplots = 11, legendlabels = c('exp(-EKE)', '2021 Age1'))
abline(v = c(1980.5, 2010.5))


r4ss::SSgetoutput(dirvec = c('test_posEKE', 'test_negEKE')) %>%
   r4ss::SSsummarize() %>%
   r4ss::SSplotComparisons(subplots = 13, legendlabels = c('+EKE', '-EKE', '2021 Age1'))

dplyr::right_join(temp$recdevs, ROMS, by = c('Yr' = 'year')) %>%
   tidyr::pivot_longer(cols = c('replist3', 'EKEpre.ms.c'), names_to = 'name') %>%
   ggplot(aes(x = Yr, y = value)) +
   geom_point() +
   geom_line()+
   xlab('Yr') +
   facet_wrap(~name, scales = 'free_y')

r4ss::SS_doRetro(
   masterdir = here('inst/extdata/models'),
   oldsubdir = 'PacificHake',
   years = -(10:15)
)

dplyr::left_join(temp$indices, temp$recdevs) %>% 
   dplyr::filter(Fleet_name == "Acoustic_Survey", Obs != 1) %>% View

r4ss::SS_output('test_negEKE') %>%
   r4ss::SS_plots(html = TRUE)

r4ss::SS_doRetro(masterdir = here('test_negEKE'),
                 oldsubdir = '', 
                 years = -(10:15)
)

temp <- r4ss::SSgetoutput(dirvec = paste0(here('test_negEKE', 'retrospectives', 'retro'), -(10:15))) %>%
   r4ss::SSsummarize()

temp %>%
   r4ss::SSsummarize() %>%
   r4ss::SSplotRetroRecruits(cohorts = 2000:2010, endyrvec = 2020 - 10:15)
