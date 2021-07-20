# KILIFI
# First load multi_res_fit_sigma_revision.RData from supplementary data of:
# 1. Lipsitch M, et al. (2012) Estimating rates of carriage acquisition and clearance
# and competitive ability for pneumococcal serotypes in Kenya with a Markov transition
# model. Epidemiology 23(4):510–519.

library(data.table)
library(ggplot)

data = as.data.table(t(fit$res))
data[, serotype := c(colnames(fit$res))]
data[, variable := substr(serotype, 1, 1)]
data[, variable := factor(variable, levels = unique(variable))]
data[, serotype := substr(serotype, 3, nchar(serotype))]
data[, serotype := factor(serotype, levels = unique(serotype))]

data[, inv.estimate := 1/estimate]
data[, inv.lower := 1/upper]
data[, inv.upper := 1/lower]

# try to recreate Fig. 2, as a check...
ggplot(data[variable == "b"]) + geom_pointrange(aes(x = serotype, y = inv.estimate, ymin = inv.lower, ymax = inv.upper)) # Good
ggplot(data[variable == "a"]) + geom_pointrange(aes(x = serotype, y = estimate, ymin = lower, ymax = upper)) # Good
ggplot(data[variable == "r"]) + geom_pointrange(aes(x = serotype, y = estimate, ymin = lower, ymax = upper)) + scale_y_log10() # Good

# calculate incidence-weighted time to clearance
weighted.mean(data[variable == "b"]$inv.estimate, data[variable == "a"]$estimate) # 71.4 days; clearance rate of 0.426 per month

fwrite(data, "~/Dropbox/Elegy/Calculations/kilifi.txt", sep = "\t")

# check carriage and resistance prevalence
q = fread("~/Dropbox/Elegy/Runs/M6-HighCarriage/K_spp_twi.txt", skip = "trial\t")
q[trial == "treatment;1;0;-1", mean(S+SR+R)]
q[trial == "treatment;1;0;-1", mean(Rfrac)]

# calculate interv_strength for a given clearance factor and susceptibility (transmission) factor
dv = 0.656
sv = 3.49
cat(paste(round(c(0,  1 - sv * (1 - seq(0.1, 0.9, by = 0.1)),  1 - (1 - seq(0.9, 0.1, by = -0.1))/dv), 5), collapse = ", "))

# Targets: 90.0% carriage prevalence and 81.4% resistance prevalence (Kobayashi et al 2017) -- Kibera (Nairobi) and Lwak.
# Duration of carriage: 71.4 days estimated from Lipsitch et al 2012 (above)