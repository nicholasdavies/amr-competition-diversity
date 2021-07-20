# Find potential for coexistence under various models
require(data.table)
require(ggplot2)
require(smoothr)
require(gsl)
library(ks)

regionrate = function(shape, r, n)
{
    X0 = gamma_inc_Q(shape + 1, qgamma(r/n, shape, 1));
    X1 = 0;
    if (r < n - 1) {
        X1 = gamma_inc_Q(shape + 1, qgamma((r+1)/n, shape, 1));
    }
    return (n * (X0 - X1));
}

inv_regionrate = function(q, r, n)
{
    return (sapply(q, function(q) uniroot(function(x) regionrate(x, r, n) - q, c(0, 500))$root, USE.NAMES = F))
}

do = function(model, scenario.c, scenario.tau, cn, xn, xv, yn = NA, yv = NA)
{
    # Assemble x/y specification for command line
    xyspec = paste0("-", xn, " ", xv);
    if (!is.na(yn)) {
        xyspec = paste0(xyspec, " -", yn, " ", yv);
    }
    
    # Find appropriate cost value
    cost = as.numeric(system(paste("~/Dropbox/Elegy/elegy", "~/Dropbox/Elegy/Runs/M3-Potential/config.cfg", scenario.c, xyspec), intern = T));
    costspec = paste0("-", cn, " ", cost);
    
    # Given above cost value, find where Rfrac lifts off x-axis
    tau = as.numeric(system(paste("~/Dropbox/Elegy/elegy", "~/Dropbox/Elegy/Runs/M3-Potential/config.cfg", scenario.tau, xyspec, costspec), intern = T));
    
    return (data.table(model = model, xvar = xn, x = xv, yvar = yn, y = yv, costvar = cn, cost = cost, tau0 = tau, pot = (0.1 - tau) / 0.1));
}

smooth = function(x, rows, cols, centre = 0.2)
{
    nx = x;
    neighbour = (1 - centre) / 4;
    for (rr in 1:rows) {
        for (cc in 1:cols) {
            weight = centre;
            i = (rr - 1) * cols + cc;
            v = x[i] * centre; 
            
            if (rr > 1)    { weight = weight + neighbour; v = v + x[i - cols] * neighbour; }
            if (rr < rows) { weight = weight + neighbour; v = v + x[i + cols] * neighbour; }
            if (cc > 1)    { weight = weight + neighbour; v = v + x[i - 1] * neighbour; }
            if (cc < cols) { weight = weight + neighbour; v = v + x[i + 1] * neighbour; }
            
            nx[i] = v / weight;
        }
    }
    return (nx);
}

smooth2 = function(x, rows, cols, centre = 0.2)
{
    nx = x;
    neighbour = (1 - centre) / 4;
    for (rr in 1:rows) {
        for (cc in 1:cols) {
            i = (rr - 1) * cols + cc;
            
            v = x[i] * centre; 
            v = v + (if (rr > 1)    x[i - cols] else x[i]) * neighbour;
            v = v + (if (rr < rows) x[i + cols] else x[i]) * neighbour;
            v = v + (if (cc > 1)    x[i - 1]    else x[i]) * neighbour;
            v = v + (if (cc < cols) x[i + 1]    else x[i]) * neighbour;
            
            nx[i] = v;
        }
    }
    return (nx);
}

# Within-host competition A
potential1 = NULL;
for (k in seq(0.5, 10, by = 0.1))
{
    cat(".");
    potential1 = rbind(potential1, do("neu", 2, 3, "c", "k", k));
}

# Within-host competition B
potential2 = NULL;
for (k in seq(0.5, 10, by = 0.1))
{
    cat(".");
    potential2 = rbind(potential2, do("whc", 4, 5, "b", "k", k));
}

# Pathogen variability
potential3 = NULL;
for (days in seq(48, 200, by = 4))
{
    delta = 1 - 365.25 / (0.65 * days * 12);
    for (mult in seq(1.2, 10, by = 0.2))
    {
        a = log(mult)/log(26/25);
        cat(".");
        potential3 = rbind(potential3, cbind(do("dty", 6, 7, "c", "delta", delta, "a", a), data.frame(days = days, mult = mult)));
    }
}

# Population heterogeneity
potential4 = NULL;
for (g in seq(0.5, 0.99, by = 0.01))
{
    for (top10 in seq(1.1, 5.0, by = 0.1))
    {
        cat(".");
        shape = inv_regionrate(top10, 9, 10);
        potential4 = rbind(potential4, cbind(do("twi", 8, 9, "c", "g", g, "shape", shape), data.frame(top10 = top10)));
    }
}

# PLOTTING
theme_set(theme_cowplot(font_size = 7, font_family = "Helvetica", line_size = 0.25))
pl1 = ggplot(potential1) + geom_line(aes(x = x, y = pot)) + labs(x = "Efficiency of co-colonisation, k", y = "Potential for coexistence", title = "Treatment competition")

pl2 = ggplot(potential2) + geom_line(aes(x = x, y = pot)) + labs(x = "Efficiency of co-colonisation, k", y = "Potential for coexistence", title = "Growth competition")

potential3$ipot = floor(10 * potential3$pot)
pl3 = ggplot(potential3) + 
    geom_raster(aes(x = days, y = mult, fill = ipot)) + 
    geom_contour(aes(x = days, y = mult, z = ipot), colour = "#ffffff") +
    labs(x = "Duration of carriage of longest-lasting serotype (days)", y = "Beta multiplier for vanishingly rare strain", title = "Pathogen diversity")
pl3
ggsave("~/Dropbox/Elegy/Submission/STM Response/Figures/S2-pd-bg.pdf", width = 8, height = 6, units = "cm")

potential4$ipot = floor(10 * potential4$pot)
pl4 = ggplot(potential4) + 
    geom_raster(aes(x = x, y = top10, fill = pot)) + 
    labs(x = "Assortativity of subpopulations, g", y = "Relative consumption by top 10%", title = "Treatment diversity") + 
    ggplot2:::limits(c(0, 1), "fill")

plot_grid(pl1, pl2, pl3, pl4, ncol = 1)
ggsave("~/Dropbox/Elegy/Figures/5-coexistence-new.pdf", width = 8, height = 24, units = "cm", device = cairo_pdf);


fwrite(potential1, file = "~/Dropbox/Elegy/potential1.txt", sep = "\t")
fwrite(potential2, file = "~/Dropbox/Elegy/potential2.txt", sep = "\t")
fwrite(potential3, file = "~/Dropbox/Elegy/potential3_new.txt", sep = "\t")
fwrite(potential4, file = "~/Dropbox/Elegy/potential4_new.txt", sep = "\t")

potential1 = fread("~/Dropbox/Elegy/Runs/M3-Potential/potential1.txt")
potential2 = fread("~/Dropbox/Elegy/Runs/M3-Potential/potential2.txt")
potential3 = fread("~/Dropbox/Elegy/Runs/M3-Potential/potential3_new.txt")
potential3[, days := 365.25 / ((1 - x) * 0.65 * 12)];
potential4 = fread("~/Dropbox/Elegy/Runs/M3-Potential/potential4_new.txt")
potential4[, top10 := round(regionrate(y, 9, 10), 1)]

regionrate(0.5, 9, 10)
regionrate(1, 9, 10)
regionrate(2, 9, 10)
regionrate(5, 9, 10)
regionrate(10, 9, 10)
regionrate(30, 9, 10)
regionrate(300, 9, 10)

# fake d-types calc?
u = 0.7
tau_star = 0.1
rho_star = 0.3

alpha = 2
delta = 0.5

potential3approx = NULL;
for (alpha in seq(1.1, 10, by = 0.1))
{
    for (delta in seq(0, 1, by = 0.01))
    {
        b = log(alpha);
        q = -log(1 - rho_star * (1 - exp(-b))) / b;
        uq = u * (1 - delta + 2 * q * delta);
        c_1mc = tau_star / uq;
        tau_hat = c_1mc * u * (1 - delta);

        pot = (tau_star - tau_hat) / tau_star;
        
        potential3approx = rbind(potential3approx, data.table(
            model = "dty_approx",
            xvar = "delta", x = delta,
            yvar = "alpha", y = alpha,
            costvar = "c", cost = c_1mc / (1 + c_1mc),
            tau0 = tau_hat, pot = pot))
    }
}

ggplot(potential3approx) + geom_raster(aes(x = x, y = y, fill = pot)) + geom_contour(aes(x = x, y = y, z = pot), breaks = seq(0, 1, by = 0.1)) + labs(x = "Clearance rate range, δ", y = "Abundance range, α", title = "Pathogen variability")
# + annotate("text", x = labs3x, y = labs3y, label = labs3, size = 2, colour = "white")


# overlay of posteriors
d_td = fread("~/Dropbox/Elegy/Runs/M1-MCMC/spp_cm_twi.txt", skip = "trial\t")[trial >= 0]
d_td$rel10 = regionrate(d_td$shape, 9, 10)

ggplot(d_td) +
    geom_density2d(aes(x = g, y = rel10)) +
    xlim(0.5, 1) +
    ylim(1, 5)
ggsave("~/Dropbox/Elegy/Submission/STM Response/Figures/S2-td.pdf", width = 8, height = 6, units = "cm")

d_pd = fread("~/Dropbox/Elegy/Runs/M1-MCMC/spp_cm_dty.txt", skip = "trial\t")[trial >= 0]

ggplot(d_pd) +
    geom_density2d(aes(x = delta, y = a)) +
    xlim(0, 1) +
    ylim(10, 50)
ggsave("~/Dropbox/Elegy/Submission/STM Response/Figures/S2-pd.pdf", width = 8, height = 6, units = "cm")

d_tc = fread("~/Dropbox/Elegy/Runs/M1-MCMC/spp_cm_neu.txt", skip = "trial\t")[trial >= 0]
d_tc[, as.list(hdi(k), 0.9)]

d_gc = fread("~/Dropbox/Elegy/Runs/M1-MCMC/spp_cm_whc.txt", skip = "trial\t")[trial >= 0]
d_gc[, as.list(hdi(k), 0.9)]

d_td[, sum(g>0.995)/.N]

