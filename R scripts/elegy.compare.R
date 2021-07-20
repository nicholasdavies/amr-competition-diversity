library(ggplot2)
library(data.table)

# compare ODE versus root-finding runs

compare = function(stem, vars)
{
    r0 = fread(paste0(stem, "_root_false.txt"), skip = "trial\t");
    r1 = fread(paste0(stem, "_root_true.txt"),  skip = "trial\t");
    
    rr = data.table(v0 = unlist(r0[, ..vars], use.names = F),
                    v1 = unlist(r1[, ..vars], use.names = F),
                    variable = rep(vars, each = nrow(r0)),
                    trial = 1:nrow(r0));
    
    #ggplot(rr) + geom_point(aes(x = v0, y = v1)) + facet_wrap(~variable, scales = "free") + labs(x = "ODE", y = "Root")
    ggplot(rr) + geom_point(aes(x = trial, y = v0 - v1)) + facet_wrap(~variable, scales = "free") + labs(x = "trial", y = "ODE - Root")
}


compare("~/Dropbox/Elegy/Runs/M4-Check/spp_whc_tau", c("X1", "S1", "Sr1", "Rs1", "R1"))


r = fread(paste0("~/Dropbox/Elegy/Runs/M4-Check/spp_whc_tau", "_root_true.txt"), skip = "trial\t");
ggplot(r) + geom_histogram(aes(k), bins = 30) + xlim(0.25,1)


# compare Elegy to Ode
ele = fread("~/Dropbox/Elegy/Runs/M1-MCMC/spp_neu_oldschool.txt", skip = "trial\t");
ele = fread("~/Dropbox/Elegy/Runs/M1-MCMC/spp_neu_oldschool_1thread.txt", skip = "trial\t");
ele = fread("~/Dropbox/Elegy/Runs/M1-MCMC/spp_neu_oldschool_normal.txt", skip = "trial\t");

cutoff = 0
ggplot(ele[trial >= cutoff]) + geom_histogram(aes(beta), bins = 30)
ggplot(ele[trial >= cutoff]) + geom_histogram(aes(c), bins = 30)
ggplot(ele[trial >= cutoff]) + geom_histogram(aes(k), bins = 30)
ggplot(ele[trial >= cutoff]) + geom_histogram(aes(l_noise), bins = 30)
ggplot(ele[trial >= cutoff]) + geom_histogram(aes(SR), bins = 30)
ggplot(ele[trial >= cutoff]) + geom_histogram(aes(ll), bins = 30)

ggplot(ele) + geom_point(aes(trial, Rfrac26))

ggplot(ele[trial >= -1950]) + geom_line(aes(x = trial, y = ll, colour = chain, group = chain)) + ylim(-130,-100)

ggplot() + geom_freqpoly(data = ele[trial < -1000 & trial > -1800], aes(x = ll), colour = "red") + geom_freqpoly(data = ele[trial >= -1000], aes(x = ll), colour = "blue")

ggplot() + geom_freqpoly(data = ele[(trial %% 100) < 50], aes(x = ll), colour = "red")  + geom_freqpoly(data = ele[trial %% 100 >= 50], aes(x = ll), colour = "blue")

q = rnorm(10000)
p = dbeta(rbeta(10000, 0.9, 0.9), 0.9, 0.9, log = T)
hist(p, breaks = 100)


# compare to just random sampling
ran = fread("~/Dropbox/Elegy/Runs/M1-MCMC/spp_neu_oldschool_random.txt", skip = "trial\t");
ran$weight = exp(ran$ll + 100)
max(ran$ll)
ran$type = ran$R26 > 0.01
ggplot(ran) + geom_histogram(aes(ll, weight = weight)) + xlim(-115,-100) + facet_grid(type~.)
ggplot(ran) + geom_histogram(aes(beta, weight = weight)) + facet_grid(type~.)
ggplot(ran) + geom_histogram(aes(k, weight = weight)) + facet_grid(type~.)
ggplot(ran) + geom_histogram(aes(c, weight = weight)) + facet_grid(type~.)
ggplot(ran) + geom_histogram(aes(l_noise, weight = weight)) + facet_grid(type~.)
ggplot(ran) + geom_histogram(aes(R26, weight = weight), bins = 100) + facet_grid(type~.)

ggplot(ran) + geom_histogram(aes(l_noise, weight = weight))

ggplot(ran) + geom_bin2d(aes(c, l_noise, weight = weight))

summary(ran$R26)



# investigate DEMCMC
dn = fread("~/Dropbox/Modules/MCMC/out.txt")
dn = dn[trial >= -1000000]

norm = function(m, x) exp(-(x-m) * (x-m))
fake = data.table(x = seq(-5, 10, by = 0.05))
fake$y = norm(0, fake$x) + 0.5 * norm(5, fake$x);

ggplot(dn) + geom_density2d(aes(x = trial, y = x))
ggplot() + 
    geom_freqpoly(data = dn[trial < -500000], aes(x, colour = "pre"), bins = 100) + 
    geom_freqpoly(data = dn[trial >= -500000], aes(x, colour = "post"), bins = 100) +
    geom_line(data = fake, aes(x, 285000 * y), linetype = "22")

ggplot(dn) + geom_bin2d(aes(trial, x))

ggplot(dn[chain == 0 & trial > -10000]) + geom_line(aes(x = trial, y = x))
