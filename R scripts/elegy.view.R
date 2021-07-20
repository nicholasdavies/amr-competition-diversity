library(ggplot2)
library(ggstance)
library(ggridges)
library(cowplot)
library(reshape2)
library(coda)
library(binom)
library(plyr)
library(data.table)
library(HDInterval)
library(msm)
library(loo)
library(cowplot)
library(Cairo)
library(gsl)
library(grid)
library(stringr)

Path = "~/Dropbox/Elegy/Runs/";
Colours = rev(hcl(h = seq(15, 375, length = 5), l = 65, c = 100)[1:4])

#
# FIT PLOTTING
#

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
    return (sapply(q, function(q) uniroot(function(x) regionrate(x, r, n) - q, c(0, 150))$root, USE.NAMES = F))
}

loaddata = function(filename, incomplete = FALSE, prop = NA, dty.short.frac = 0.5)
{
    # load data
    if (incomplete) {
        data = read.table(paste(Path, filename, ".txt", sep = ""), header = TRUE, fill = T)
    } else {
        data = fread(paste(Path, filename, ".txt", sep = ""), skip = "trial\t", data.table = F, sep = "\t", fill = T)
    }

    # remove burn-in
    if (is.numeric(data$trial)) {
        if (max(data$trial) >= 0) {
            data = data[data$trial >= 0, ];
            if (!is.na(prop)) {
                data = data[data$trial >= quantile(data$trial, prop, names = F), ];
            }
        } else {
            if (is.na(prop)) {
                prop = 0.5;
            }
            data = data[data$trial >= quantile(data$trial, prop, names = F), ];
        }
    } else {
        data = data[round(nrow(data) * prop):nrow(data),];
    }
    
    # add convenience columns
    data$carriage = data$S + data$R + data$SR
    data$sigma = data$l_noise
    data$dual = data$SR / data$carriage
    if ("shape" %in% names(data)) {
        data$top10 = regionrate(data$shape, 9, 10);
    }
    if ("delta" %in% names(data)) {
        longest = (365.25 / 12) / (u * (1 - data$delta));
        shortest = (365.25 / 12) / (u * (1 + data$delta));
        data$extreme = cbind(longest, shortest)[cbind(1:nrow(data), rbinom(nrow(data), 1, dty.short.frac) + 1)];
    }
    if ("a" %in% names(data)) {
        data$strength = (1 + 1/25)^data$a;
    }

    return (data)
}

rbind0 = function(a, b) {
    if (is.null(a)) {
        return (b);
    } else {
        return (rbind(a, b));
    }
}

lighten = function(cols, a = .5)
{
    sapply(cols, function(h) rgb(t(col2rgb(h) * a + 255 * (1-a)), maxColorValue = 255), USE.NAMES = F)
}

# deviance information criterion
DIC <- function(data)
{
    pD <- var(-2 * data$ll)
    Db <- mean(-2 * data$ll)
    return (list(DIC = pD + Db, pD = pD, D = Db))
}

# estimated Bayesian information criterion
eBIC <- function(data)
{
    first <- which(colnames(data) == "ll") + 1
    last <- which(colnames(data) == "S") - 1
    n <- length(r)
    k <- last - first + 1
    ml <- max(data$ll)
    
    return(list(BIC = (log(n) - log(2*pi)) * k - 2 * ml, n = n, k = k))
}

# Akaike information criterion
AIC <- function(data)
{
    first <- which(colnames(data) == "ll") + 1
    last <- which(colnames(data) == "S") - 1
    k <- last - first + 1
    
    return (2*k - 2 * max(data$ll))
}

#
# WAIC
#

# truncated normal-binomial point density
norm.binom = function(x, mu, sigma, n, r)
{
    dnorm(x, mu, sigma) * dbinom(r, n, x)
}

# overall likelihood for normal-binomial density
point.likelihood.normal.trunc = function(model.rho, model.sigma, data.n, data.r)
{
    a = integrate(norm.binom, 0, 1, model.rho, model.sigma, data.n, data.r)$value
    b = pnorm(1, model.rho, model.sigma) - pnorm(0, model.rho, model.sigma)
    return (a / b)
}

# convert data set into set of log pointwise probabilities
# if supplied, per-country weights must have a mean of 1
get.lpp = function(d, model = "binomial", weights = NA)
{
    lpp = matrix(0, nrow(d), length(r))
    col0 = which(names(d) == "Rfrac0")
    colj = which(names(d) == "Rfrac1") - col0
    for (row in 1:nrow(d))
    {
        for (col in 1:length(r))
        {
            rho = d[row, col0 + colj * (col - 1)]
            if (model == "binomial") {
                lpp[row, col] = dbinom(r[col], n[col], rho, log = T);
            } else if (model == "normal.trunc") {
                lpp[row, col] = log(point.likelihood.normal.trunc(rho, d$l_noise[row], n[col], r[col]));
            } else {
                stop("invalid model");
            }
            
            if (all(!is.na(weights))) {
                lpp[row, col] = lpp[row, col] * weights[col];
            }
        }
    }
    return (lpp)
}

lppd = function(lpp)
{
    sum(log(colMeans(exp(lpp))))
}

pWAIC2 = function(lpp)
{
    sum(apply(lpp, 2, function(col) var(col)))
}

WAIC = function(lpp)
{
    -2 * (lppd(lpp) - pWAIC2(lpp))
}

# summarize likelihoods for a set of MCMC chains
summarise.likelihoods <- function(data.list, names)
{
    AICs = lapply(data.list, function(a) { return (AIC(a)) })
    DICs = lapply(data.list, function(a) { return (DIC(a)) })
    eBICs = lapply(data.list, function(a) { return (eBIC(a)) })
    LLD = lapply(data.list, function(a) { return (hdi(-2 * a$ll, credMass = 0.95)) })
    cat(paste(c("Model   ", names, "\n"), collapse="\t"))
    cat(paste(c("AIC     ", AICs, "\n"), collapse="\t"))
    cat(paste(c("DIC     ", DICs, "\n"), collapse="\t"))
    cat(paste(c("eBIC    ", eBICs, "\n"), collapse="\t"))
    cat(paste(c("HDI95(D)", LLD, "\n"), collapse="\t"))
}

# takes a long time...
summarise.WAICs = function(data.list, names, model, weights = NA)
{
    cat("\nWAIC and LOO-CV\n")
    waics = vector("list", length(data.list))
    loos = vector("list", length(data.list))
    for (i in 1:length(data.list))
    {
        lpp = get.lpp(data.list[[i]], model, weights)
        name = names[[i]]
        
        myWAIC = WAIC(lpp)
        theirWAIC = waic(lpp)
        loo_relative_eff = relative_eff(exp(lpp), data.list[[i]]$chain + 1)
        lool = loo(lpp, r_eff = loo_relative_eff)
        cat(paste0("\n", name, "\n\nWAIC:\n"))
        cat(paste0(myWAIC, "\n"))
        cat("\nloo-calculated WAIC:\n")
        print(theirWAIC)
        cat("\nLOOIC:\n")
        print(lool)
        cat("\n")
        
        loos[[i]] = lool
        waics[[i]] = theirWAIC
    }
    cat("\nComparison of WAICs:\n")
    print(compare(x = waics))
    cat("\nComparison of LOOICs:\n")
    print(compare(x = loos))
}

figname = function(stem)
{
    return (paste0("~/Dropbox/Elegy/Submission/STM Response/", stem, "-", filename, postfix, ".pdf"));
}

diagnose = function(data.list, name.list)
{
    for (i in 1:length(data.list))
    {
        x = data.list[[i]]
        chains = x$chain
        first = which(names(x) == "ll")
        last = which(names(x) == "S") - 1
        x = x[, first:last]
        M = mcmc.list(lapply(split(x, chains), mcmc))
        es = effectiveSize(M)
        gd = gelman.diag(M)
    
        cat(paste(c("Model", name.list[[i]], "\n"), collapse="\t"))
        cat(paste(c("Parameter", names(es), "\n"), collapse="\t"))
        cat(paste(c("Eff. size", es, "\n"), collapse="\t"))
        cat(paste(c("R est.", gd$psrf[,1], "\n"), collapse="\t"))
        cat(paste(c("R 95% UCI", gd$psrf[,2], "\n"), collapse="\t"))
        cat(paste(c("MPSRF", gd$mpsrf, "\n"), collapse="\t"))
    }
}


# combine compartments to country level
combine_compartments = function(d, nc, ng, h)
{
    h = rep_len(h, nc * ng);
    l0 = which(names(d) == "epochs") + 1;
    l1 = which(names(d) == "Rfrac0");
    cols = names(d)[l0:l1];
    cols = substr(cols, 1, nchar(cols) - 1);
    if (any(cols != c("X", "S", "R", "Rfrac"))) {
        stop("combine_compartments is only set up to work with the 'simple' base model.");
    }
    
    dd = d[, 1:(l0 - 1)];
    
    for (i in 1:nc) {
        idx = ((i - 1) * ng):(i * ng - 1);
        xc = paste0("X", idx);
        sc = paste0("S", idx);
        rc = paste0("R", idx);
        newcols = data.frame(x = rowSums(t(t(d[, xc]) * h[idx + 1])),
                             s = rowSums(t(t(d[, sc]) * h[idx + 1])),
                             r = rowSums(t(t(d[, rc]) * h[idx + 1])));
        newcols$rfrac = newcols$r / (1 - newcols$x);
        newcols$rrftop = d[, paste0("Rfrac", i * ng - 1)] / newcols$rfrac;
        names(newcols) = paste0(c("X", "S", "R", "Rfrac", "rrftop"), i - 1);
        dd = cbind(dd, newcols);
    }
    
    r0 = which(names(dd) == "rrftop0");
    r1 = which(names(dd) == "rrftop1") - r0;
    rrftoprand = sample(0:(nc - 1), nrow(dd), replace = T);
    rrftoprand = dd[cbind(1:nrow(dd), r0 + r1 * rrftoprand)];
    names(rrftoprand) = "rrftoprand";
    
    return (cbind(dd, d[, (l0 + nc * ng * (l1 - l0 + 1)):ncol(d)], rrftoprand));
}

# get distribution of treatment rates and resistance prevalences
tau_factors = function(dlist, namevec, z0vec, n_countries, all_varying = T, override_z = NA)
{
    D = NULL;
    for (i in 1:length(dlist))
    {
        if (all_varying) {
            z0 = which(names(dlist[[i]]) == z0vec[i]);
            zs = melt(t(dlist[[i]][, z0:(z0 + n_countries - 1)]))$value;
        } else if (!is.na(override_z)) {
            zs = rep(override_z, each = nrow(dlist[[i]]) * n_countries);
        } else {
            zs = rep(dlist[[i]][[z0vec[i]]], each = n_countries);
        }
        d = data.table(model = namevec[i], country = rep(1:n_countries), z = zs);
        d[, tau := z * tau];
        rfrac0 = which(names(dlist[[i]]) == "Rfrac0");
        rfracD = which(names(dlist[[i]]) == "Rfrac1") - rfrac0;
        d[, rfrac := melt(t(dlist[[i]][, seq(rfrac0, rfrac0 + (n_countries - 1) * rfracD, by = rfracD)]))$value];
        D = rbind0(D, d);
    }
    return (D);
}

# Posterior plot as ridgelines
ridges = function(dlist, dnames, dlongnames, vars, descs, priors, multipliers = 1, nrow = 1, trim = 0, legpos = "bottom", legdir = "horizontal", legjust = c(0.5, 0.5), legks = 0.2)
{
    # prepare data table for posteriors
    cat("Preparing posteriors...\n");
    D = NULL;
    for (i in 1:length(dlist))
    {
        d = dlist[[i]][, names(dlist[[i]]) %in% vars];
        d$model = dnames[i];
        D = rbind.fill(D, d);
    }
    D = as.data.table(melt(D, id.vars = "model"));
    D$model = factor(D$model, levels = c("prior", dnames));
    D = D[, .SD[((value >= quantile(value, probs = trim/2, na.rm = T)) & (value <= quantile(value, probs = 1 - trim/2, na.rm = T)))], by = .(model, variable)]
    D$value[D$value == 0] = NA;
    D$variable = factor(descs[match(D$variable, vars)], levels = descs);
    
    for (i in 1:length(descs))
    {
        D[variable == descs[i], value := value * multipliers[(i - 1) %% length(multipliers) + 1]];
    }
    
    # prepare data table for priors
    cat("Preparing priors...\n");
    P = NULL;
    fvars = levels(D$variable);
    for (j in 1:length(fvars))
    {
        xmin = min(D[variable == fvars[j], value], na.rm = T);
        xmax = max(D[variable == fvars[j], value], na.rm = T);
        if (!is.finite(xmin) | !is.finite(xmax)) next;
        extra = (xmax - xmin) * 0.25;
        xmin = xmin - extra;
        xmax = xmax + extra;
        xval = seq(xmin, xmax, length.out = 101);
        base = length(dlist) + 1;
        yval = priors[[j]](xval);
        nonzero = yval > 0;
        yval = base + 2 * yval / max(yval);
        p = data.table(x = xval, ymin = rep(base, length(xval)), ymax = yval, yval = yval, variable = fvars[j], model = factor("prior", levels = c("prior", dnames)));
        p$ymin[1:10] = NA;
        p$ymax[1:10] = NA;
        p$ymin[91:101] = NA;
        p$ymax[91:101] = NA;
        p = p[nonzero, ];
        P = rbind.fill(P, p);
    }
    P$variable = factor(P$variable, levels = descs);
    
    # colour scale
    colour_palette = c("#dddddd", Colours);
    names(colour_palette) = c("prior", dnames);

    cat("Plotting...\n");
    ggplot(D) + 
        geom_line(data = P, aes(x = x, y = yval), colour = "#000000", size = 0.3, linetype = "22") + 
        geom_ribbon(data = P, aes(x = x, ymin = ymin, ymax = ymax, fill = model)) + 
        geom_line(data = P, aes(x = x, y = ymax), colour = "#000000", size = 0.3) + 
        geom_density_ridges(aes(x = value, y = model, fill = model), 
                            scale = max(1, length(dlist) - 1.2), alpha = 0.5, size = 0.3, stat = "binline", bins = 25) + 
        facet_wrap(~variable, scales = "free_x", strip.position = "top", nrow = nrow) + 
        theme(axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
              strip.placement = "outside", strip.background = element_blank(), 
              panel.spacing.x = unit(0.2, "cm"), panel.spacing.y = unit(0.15, "cm"), panel.background = element_blank(), plot.background = element_blank(),
              legend.position = legpos, legend.direction = legdir, 
              legend.justification = legjust, legend.title = element_blank(), legend.key.size = unit(legks, "cm")) + 
        labs(x = NULL, y = NULL) + 
        scale_y_discrete(limits = unique(rev(D$model))) + 
        coord_cartesian(ylim = c(0.6, max(P$ymax + 0.25, na.rm = T)), expand = F) + 
        scale_fill_manual(breaks = c("prior", dnames), labels = paste0("  ", c("Prior", dlongnames), "  "), values = colour_palette) + 
        scale_colour_manual(breaks = c("prior", dnames), values = c("#666666", rep("#000000", length(dlist))))
}

# Prior-generating functions
prnull    = function()            return (function(x) rep(NA, length(x)));
prunif    = function(a, b)        return (function(x) dunif(x, a, b));
prgamma   = function(k, th)       return (function(x) dgamma(x, shape = k, scale = th));
prbeta    = function(a, b, s = 1) return (function(x) dbeta(x / s, a, b));
prnorm    = function(mu, sd)      return (function(x) dnorm(x, mu, sd));
prlnorm   = function(mu, sd)      return (function(x) dlnorm(x, mu, sd));
prhist    = function(hist)        return (function(x) sapply(x, function(y) hist$density[which.min(abs(hist$mids - y))], USE.NAMES = F))

# Visual check of chain settling on a stable distribution
inspect_convergence = function(d, nrow)
{
    first = which(colnames(d) == "ll") + 1;
    last = which(colnames(d) == "S") - 1;
    d$trialbin = floor(d$trial / ((max(d$trial) + 1) / 30));
    d = melt(d[, c(first:last, which(colnames(d) == "trialbin"))], id.vars = "trialbin");
    ggplot(d) + geom_density_ridges(aes(x = value, y = as.factor(trialbin), fill = as.factor(trialbin)), colour = "black", size = 0.1) + facet_wrap(~variable, scales = "free_x", nrow = nrow) + theme(legend.position = "none", axis.text.y = element_text(size = 3))
}


#
# INTERVENTION PLOTTING
#

# load intervention set
loadint = function(prefix, models, postfix, dlongnames, cut = NA)
{
    d = NULL;
    for (m in 1:length(models))
    {
        e = fread(paste0(Path, prefix, models[m], postfix), skip = "\t");
        if (!is.na(cut[m])) {
            e = e[floor(cut[m] * nrow(e) + 1):nrow(e),];
        }
        e[, c("iname", "icoverage", "istrength", "itime") := tstrsplit(trial, split = ';', fixed = T)];
        e[, icoverage := as.numeric(icoverage)];
        e[, istrength := as.numeric(istrength)];
        e[, itime := as.numeric(itime)];
        e[, model := models[m]];
        e[, Y := S + R + SR];
        d = rbind0(d, e[, .(Rfrac, Y, iname, icoverage, istrength, itime, model)]);
    }
    d[, model := factor(model, levels = models, labels = dlongnames)];
    d[, iname := factor(iname, levels = c("none", "vaccine", "clearance", "treatment", "hightreatment_between", "hightreatment_within"))];
    d[, R := Y * Rfrac];
    return (d)
}

# load per-country intervention set
loadint_percountry = function(prefix, models, postfix, dlongnames, n_countries, n_regions, cut = NA)
{
    d = NULL;
    for (m in 1:length(models))
    {
        e = fread(paste0(Path, prefix, models[m], postfix), skip = "\t");
        if (!is.na(cut[m])) {
            e = e[floor(cut[m] * nrow(e) + 1):nrow(e),];
        }
        e[, c("iname", "icoverage", "istrength", "itime") := tstrsplit(trial, split = ';', fixed = T)];
        e[, icoverage := as.numeric(icoverage)];
        e[, istrength := as.numeric(istrength)];
        e[, itime := as.numeric(itime)];
        e[, model := models[m]];
        for (co in 0:(n_countries - 1))
        {
            e[, country := co];
            if (models[m] %like% "twi") {
                Rfc0 = which(names(e) == paste0("Rfrac", co * n_regions));
                Rfc1 = which(names(e) == paste0("Rfrac", (co + 1) * n_regions - 1));
                Rf = e[, seq(Rfc0, Rfc1, length.out = n_regions), with = F];
                Xc0 = which(names(e) == paste0("X", co * n_regions));
                Xc1 = which(names(e) == paste0("X", (co + 1) * n_regions - 1));
                Y = 1 - e[, seq(Xc0, Xc1, length.out = n_regions), with = F];
                e$Y = rowMeans(Y);
                e$Rfrac = rowSums(Y * Rf) / rowSums(Y);
            } else {
                e$Rfrac = e[[paste0("Rfrac", co)]];
                e$Y = 1 - e[[paste0("X", co)]];
            }
            d = rbind0(d, e[, .(country, Rfrac, Y, iname, icoverage, istrength, itime, model)]);
        }
    }
    d[, model := factor(model, levels = models, labels = dlongnames)];
    d[, iname := factor(iname, levels = c("none", "vaccine", "clearance", "treatment", "hightreatment_between", "hightreatment_within"))];
    d[, R := Y * Rfrac];
    return (d)
}

# Establish a synonym in intervention set
synonym = function(d, iname_from, istrength_from, iname_to, istrength_to)
{
    s = d[iname == iname_from & istrength == istrength_from];
    s$iname = iname_to;
    s$istrength = istrength_to;
    return (rbind(d, s));
}

# prettify intervention names
pretty_names = function(d)
{
    d$iname = factor(d$iname, levels = c("none", "vaccine", "clearance", "treatment", "hightreatment_between", "hightreatment_within"),
                     labels = c("None", "Acquisition-blocking vaccine", "Clearance-accelerating vaccine", "Reduced treatment", "Reduced treatment (targeted)", "Reduced treatment (targeted within-country)"));
    return (d);
}

# process as differential
diffint = function(interv, baseline)
{
    baseline = baseline[, .(Rfrac, Y, R)];
    setnames(baseline, c("Rfrac", "Y", "R"), c("Rfrac0", "Y0", "R0"));
    return (cbind(interv, baseline[rep(1:nrow(baseline), each = nrow(interv) / nrow(baseline))]))
}

# how much of intervention B is equivalent to intervention C at strength S, in terms of outcome x
baseint = function(d, metric, i_compare, i_strengths, x, interp.method = "linear")
{
    # get points for interpolation
    pts = split(metric[[x]], ceiling(seq_along(metric[[x]]) / uniqueN(metric$istrength)))
    
    # for each intervention strength
    D = NULL;
    dd = d[iname == i_compare];
    dd = rbind(d[iname == "treatment" & istrength == 0.0], dd);
    dd$iname = i_compare;
    
    for (i_str in i_strengths)
    {
        # extract comparator dataset
        co = dd[istrength == i_str];
        
        # get interpolated points
        xi = mapply(function(pts, xi) approxfun(pts, seq_along(pts), interp.method, rule = 2)(xi), pts = pts, xi = co[[x]]);

        # transform into levels
        eq = approx(seq_along(unique(metric$istrength)), unique(metric$istrength), xi, rule = 2);
        
        D = rbind0(D, data.table(x = eq$x, y = eq$y, iname = i_compare, istrength = i_str, model = co$model));
    }
    return (D);
}

htitle = function(text, size = 6, x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5)
{
    ggdraw() + draw_label(text, x = x, y = y, size = size, hjust = hjust, vjust = vjust)
}




# # # # # # # # # # # #
#                     #
#  FIGURE GENERATION  #
#                     #
# # # # # # # # # # # #

source("~/Dropbox/Elegy/R scripts/elegy.select.R")
select("spp");

# set base theme
theme_set(theme_cowplot(font_size = 7, font_family = "Helvetica", line_size = 0.25))

# ALL FIGURES
cat("Fit figures...\n");
cat("Loading 1/4...\n");
d_neu = loaddata(sprintf("M1-MCMC/%s_cm_neu%s", filename, postfix));
cat("Loading 2/4...\n"); 
d_whc = loaddata(sprintf("M1-MCMC/%s_cm_whc%s", filename, postfix)); 
cat("Loading 3/4...\n"); 
d_dty = loaddata(sprintf("M1-MCMC/%s_cm_dty%s", filename, postfix), dty.short.frac = 0.145); 
cat("Loading 4/4...\n"); 
d_twi = loaddata(sprintf("M1-MCMC/%s_cm_twi%s", filename, postfix)); 
d_twi = combine_compartments(d_twi, length(r), 10, rep(0.1, 10)); 

dnames = c("twi", "dty", "neu", "whc");
dlongnames = c("Treatment diversity", "Pathogen diversity", "Treatment competition", "Growth competition");
dabbrevnames = c("Treatment div.", "Pathogen div.", "Treatment comp.", "Growth comp.");

d_neu$model = factor("neu", levels = dnames);
d_whc$model = factor("whc", levels = dnames);
d_dty$model = factor("dty", levels = dnames);
d_twi$model = factor("twi", levels = dnames);

d_neu$l_noise = d_neu$sigma = 0.06
d_whc$l_noise = d_whc$sigma = 0.06
d_dty$l_noise = d_dty$sigma = 0.06
d_twi$l_noise = d_twi$sigma = 0.06

dlist = list(d_twi, d_dty, d_neu, d_whc);

# # inspect convergence of chains
# cat("Convergence...\n");
# pl1 = inspect_convergence(d_neu, 2) + labs(title = dlongnames[1], y = "time (binned)")
# pl2 = inspect_convergence(d_whc, 2) + labs(title = dlongnames[2], y = "time (binned)")
# pl3 = inspect_convergence(d_dty, 2) + labs(title = dlongnames[3], y = "time (binned)")
# pl4 = inspect_convergence(d_twi, 2) + labs(title = dlongnames[4], y = "time (binned)")
# FS2 = plot_grid(pl1, pl2, pl3, pl4, nrow = 4)
# ggsave(figname("S2-converge"), FS2, width = 26, height = 20, units = "cm", device = cairo_pdf);

# generate posterior
cat("Posterior...\n");
F1p = ridges(dlist, dnames, dabbrevnames,
 c(     "c", "b",     "beta",          "g", "shape",        "a",   "delta",    "k", 
   "X0", "Rfrac", "carriage", "rrftoprand", "top10", "strength", "extreme", "dual"),
 c("Transmission\ncost of resistance, c", "Growth rate of\nsensitive strain, b", "Transmission rate, β\n(per month)", "Assortativity of\nsubpopulations, g", "Shape of treatment-\nrate distribution, κ", "Power of diversifying\nselection, a", "Range of\nclearance rates, δ", "Relative rate of\nco-colonisation, k",
   "dummy", "\nResistance\nfrequency (%)", "\nCarriage\nprevalence (%)", "Resistant carriage\nin top 10%,\nrelative to mean", "Treatment rate\nin top 10%,\nrelative to mean", "β multiplier for\nvanishingly\nrare strains", "Natural carriage\nduration of\nsubtypes (days)", "Fraction of dual\ncarriers among all\ncarriers (%)"),
 list(prbeta(1.5, 8.5), prgamma(2, 0.5), prgamma(beta.shape, beta.scale), prbeta(10, 1.5), prgamma(4, 2), prgamma(2, 5), prbeta(20, 25), prnorm(1, 0.5), prnull(), prnull(), prnull(), prnull(), prnull(), prnull(), prnull(), prnull()),
 multipliers = c(1, 1, 1, 1, 1, 1, 1, 1,
     1, 100, 100, 100, 1, 1, 1, 1),
 2, 0.004, c(-0.012, 0.005), "vertical", c(0,0), 0.345)
G1p = plot_to_gtable(F1p)
G1p = gtable_remove_grobs(G1p, c("panel-2-1", "strip-t-1-2", "axis-b-1-2")) # unclear what the logic in naming is here . . .

ggsave(figname("2-posterior"),
       ggdraw(G1p) + 
            draw_text("b", 1/178, 48/50, 10, 0, 0, fontface = "bold") +
            draw_text("shortest", 0.767, 0.36, 6, 0.5, 0.5) +
            draw_text("longest", 0.821, 0.36, 6, 0.5, 0.5) +
            draw_grob(rectGrob(gp = gpar(fill = "#ffffff", col = NA)), x = 0.778, y = 0.13, width = 0.008, height = 0.1),
        width = 18.4, height = 7, units = "cm");

# get distribution of treatment rates and rfracs
cat("Model fits...\n");
if (postfix == "") {
    tx = tau_factors(dlist, dnames, "z", length(r), F, 1);
} else {
    tx = tau_factors(dlist, dnames, c("z0", "z0", "z0", "z0"), length(r))
}

# get empirical
empirical = binom.confint(r, n, methods = "bayes", conf.level = 0.95, type = "highest", prior.shape1 = 1, prior.shape2 = 1);
# plot as ribbon showing quantiles, with estimated consumption on x axis
qu = tx[, as.list(c(mean(tau), hdi(rfrac, credMass = 0.95), hdi(rfrac, credMass = 0.5))), by = .(model, country)];
names(qu) = c("model", "country", "avg_tau", "q1", "q4", "q2", "q3")
qu[, erf_lower := rep(empirical$lower, 4)];
qu[, erf_upper := rep(empirical$upper, 4)];
qu[, erf_mean := rep(empirical$mean, 4)];
qu$model = factor(qu$model, levels = dnames, labels = dlongnames);

if (filename == "spp") {
    waics = data.frame(label = c("WAIC: 246.2 ± 67.3", "244.2 ± 66.4", "248.3 ± 67.4", "247.6 ± 67.2"), model = as.factor(dlongnames), x = 0.32, y = 58);
} else if (filename == "ecp") {
    waics = data.frame(label = c("WAIC: 449.4", "559.8", "775.0", "527.9"), model = as.factor(dlongnames), x = 0.164, y = 88);
}

F1f = ggplot(qu) + 
    geom_ribbon(aes(x = 12 * avg_tau, ymin = 100 * q1, ymax = 100 * q4, fill = model), alpha = 0.5) + 
    geom_ribbon(aes(x = 12 * avg_tau, ymin = 100 * q2, ymax = 100 * q3, fill = model), alpha = 1) + 
    geom_linerange(aes(x = 12 * avg_tau, ymin = 100 * erf_lower, ymax = 100 * erf_upper), size = 0.3) + 
    geom_point(aes(x = 12 * avg_tau, y = 100 * erf_mean), size = 0.1) + 
    labs(x = "Antibiotic consumption rate (courses per person per year)", y = "Resistance\nfrequency (%)") + 
    facet_wrap(~model, nrow = 1) + 
    expand_limits(y = 0) +
    geom_text(data = waics, aes(x = x, y = y, label = label), size = 6 / ggplot2:::.pt, hjust = 0) +
    scale_fill_manual(values = Colours) +
    theme(legend.position = "none", strip.background = element_blank(), 
          strip.text = element_text(margin = margin(t = 0, r = 0, b = 2, l = 0)),
          axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)))
ggsave(figname("2-fit"), ggdraw(F1f) + 
            draw_text("a", 1/178, 33/35, 10, 0, 0, fontface = "bold"), width = 18.4, height = 3.5, units = "cm", device = cairo_pdf);

#
# INTERVENTIONS FIGURES
#

# load and split data set
cat("Loading interventions...\n");
master = loadint(paste0("M2-Interventions/", filename, "_"), c("cm_twi", "cm_dty", "cm_neu", "cm_whc"), paste0(postfix, ".txt"), dlongnames);
master = synonym(master, "treatment", 1, "hightreatment_between", 1);
master = synonym(master, "vaccine", 1, "clearance", 1);
interv = master[istrength >= 0 & istrength <= 1.0];
metric = master[iname == "treatment"];
baseline = master[iname == "treatment" & istrength == 0];
baseline$iname = "none";

nones = interv[istrength == 0];
nones$iname = "vaccine"; interv = rbind(interv, nones);
nones$iname = "clearance"; interv = rbind(interv, nones);
nones$iname = "treatment"; interv = rbind(interv, nones);
nones$iname = "hightreatment_between"; interv = rbind(interv, nones);

# Intervention plots
cat("Intervention plots...\n");
test = diffint(interv, baseline)
test$plotx = test$istrength - 0.0375 + 0.015 * as.numeric(test$model)
test[Y < 0.001, Rfrac := NA]
test[is.nan(Y), Y := 0]

test$value = test$Y * 100;
test$variable = "Carriage prevalence (%)";
triplicate = test;
test$value = test$Rfrac * 100;
test$variable = "Resistance frequency (%)";
triplicate = rbind(triplicate, test);
test$value = test$R * 100;
test$variable = "Resistant carriage (%)";
triplicate = rbind(triplicate, test);
triplicate = triplicate[iname != "hightreatment_between"];
triplicate = pretty_names(triplicate);

hdit = triplicate[, as.list(c(hdi(value, credMass = 0.95), list(m = mean(value, na.rm = T)))), by = .(variable, plotx, iname, icoverage, istrength, model)]
names(hdit)[7:9] = c("lo95", "hi95", "mean")
if (filename == "spp") {
    hdit = hdit[!(iname %like% "vaccine") | (istrength <= 0.6)];
}

F2 = ggplot(hdit) +
    geom_point(aes(x = 100 * plotx, y = mean, colour = model), size = 0.4) + 
    geom_linerange(aes(x = 100 * plotx, ymin = lo95, ymax = hi95, colour = model), size = 0.3) + 
    facet_grid(variable ~ iname, switch = "y", scales = "free") + 
    labs(x = " ", y = NULL) + 
    scale_x_continuous(breaks = seq(0, 100, by = 10), expand = c(.02, .05)) + 
    scale_colour_manual(values = Colours) +
    theme(panel.background = element_rect(fill = "#f4f4f4"), 
          panel.grid.major.y = element_line(colour = "#ffffff"), 
          strip.placement = "outside", 
          strip.background = element_rect(fill = NA), 
          legend.position = c(0.67, 0.7), 
          legend.direction = "vertical", 
          legend.justification = c(0, 0), 
          legend.title = element_blank(),
          legend.key.height = unit(7, "pt"),
          legend.key.width = unit(4, "pt"),
          axis.title.x = element_text(margin = margin(t = -2, r = 0, b = 0, l = 0, unit = "pt")),
          strip.text.y = element_text(margin = margin(t = 0, r = 1, b = 0, l = 2, unit = "pt")),
          plot.margin = margin(0, 2, 2, -1, "pt"),
          panel.spacing.x = unit(2, "pt"))

ggsave(figname("3-interv"), 
       ggdraw(F2) + 
          draw_text(c("a", "b", "c"), 0/114, c(78/80, 52/80, 27/80), 10, 0, 0, fontface = "bold") +
          draw_text(c("Vaccine efficacy (%)", "Treatment-rate reduction (%)"), c(22/58, 49/58), 1/200, 6, 0.5, 0),
       width = 12.5, height = 10, units = "cm");

hdit2 = triplicate[, as.list(c(hdi(value, credMass = 0.95), hdi(value, credMass = 0.5), list(m = mean(value, na.rm = T)))), by = .(variable, plotx, iname, icoverage, istrength, model)]
names(hdit2)[7:11] = c("lo95", "hi95", "lo50", "hi50", "mean")
if (filename == "spp") {
    hdit2 = hdit2[!(iname %like% "vaccine") | (istrength <= 0.6)];
}

ggplot(hdit2) +
    geom_ribbon(aes(x = 100 * istrength, ymin = lo50, ymax = hi50, fill = model), colour = NA, alpha = 1.0, size = 0.3) + 
    geom_ribbon(aes(x = 100 * istrength, ymin = lo95, ymax = hi95, fill = model), colour = NA, alpha = 0.3, size = 0.3) + 
    geom_line(aes(x = 100 * istrength, y = mean, colour = model), size = 0.3) + 
    facet_grid(variable ~ iname, switch = "y", scales = "free") + 
    labs(x = " ", y = NULL) + 
    scale_x_continuous(breaks = seq(0, 100, by = 10), expand = c(.02, .05)) + 
    scale_colour_manual(values = Colours) +
    scale_fill_manual(values = Colours) +
    theme(panel.background = element_rect(fill = "#f4f4f4"), 
          panel.grid.major.y = element_line(colour = "#ffffff"), 
          strip.placement = "outside", 
          strip.background = element_rect(fill = NA), 
          legend.position = c(0.67, 0.68), 
          legend.direction = "vertical", 
          legend.justification = c(0, 0), 
          legend.title = element_blank(),
          legend.key.height = unit(7, "pt"),
          legend.key.width = unit(4, "pt"),
          axis.title.x = element_text(margin = margin(t = -2, r = 0, b = 0, l = 0, unit = "pt")),
          strip.text.y = element_text(margin = margin(t = 0, r = 1, b = 0, l = 2, unit = "pt")),
          plot.margin = margin(0, 2, 2, -1, "pt"),
          panel.spacing.x = unit(2, "pt"))

# Equivalence plots
cat("Loading interventions...\n");
master = loadint(paste0("M7-Variation/vacc_", filename, "_"), c("cm_twi", "cm_dty", "cm_neu", "cm_whc"), paste0(postfix, ".txt"), dlongnames);
interv = master[istrength >= 0 & istrength <= 1.0];
metric = master[iname == "treatment"];
baseline = master[iname == "treatment" & istrength == 0];
baseline$iname = "none";

cat("Equivalence plots...\n");
master[R < 0.001, R := 0]
metric[R < 0.001, R := 0]
eqV = baseint(master, metric, "vaccine", seq(0.0, 0.6, by = 0.02), "R", "linear");
eqC = baseint(master, metric, "clearance", seq(0.0, 0.6, by = 0.02), "R", "linear");
eq = rbind(eqV, eqC)
eq$rx = 100 * eq$y
eq2 = eq[, .(rx), by = .(iname, istrength, model)]
names(eq2) = c("iname", "istrength", "model", "m")
eq2 = pretty_names(eq2)

eq2[duplicated(m) & m > 0, m := 100, by = .(iname, model)]
eq2[, m2 := m - shift(m), by = .(iname, model)]

F3a = ggplot() +
    annotate("ribbon", x = c(0, 60), ymin = c(8.8, 8.8), ymax = c(23.1, 23.1), fill = "#dddddd") + 
    annotate("line", x = c(0, 60), y = c(15, 15), linetype = "22", size = 0.2) +
    geom_line(data = eq2[m2 != 0], aes(x = 100 * istrength, y = m, colour = model, group = model)) + 
    facet_wrap(~iname, nrow = 1, scales = "free") + 
    coord_cartesian(expand = F, xlim = c(0, 60), ylim = c(-42, 100)) +
    scale_x_continuous(breaks = 10*(0:6)) +
    scale_y_continuous(breaks = 20*(-2:5)) +
    scale_colour_manual(values = Colours, guide = guide_legend(nrow = 2)) +
    labs(x = "Vaccine efficacy (%)", y = "Equivalent reduction\nin antibiotic use (%)") + 
    theme(strip.placement = "outside", strip.background = element_rect(fill = NA), 
          legend.position = "bottom", legend.direction = "horizontal",
          legend.title = element_blank(), 
          legend.text = element_text(margin = margin(t = -15, b = -15, unit = "pt")),
          legend.key.height = unit(6, "pt"), legend.key.width = unit(8, "pt"),
          axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
          panel.spacing.x = unit(6, "pt"))

F3a

# Specific results
qq = eq2[iname == "Acquisition-blocking vaccine" & model == "Treatment competition" & istrength %in% c(0.10, 0.12)]$m
10 + 2 * (15 - qq[1]) / (qq[2] - qq[1])

qq = eq2[iname == "Acquisition-blocking vaccine" & model == "Growth competition" & istrength %in% c(0.5, 0.52)]$m
50 + 2 * (15 - qq[1]) / (qq[2] - qq[1])

qq = eq2[iname == "Clearance-accelerating vaccine" & model == "Treatment competition" & istrength %in% c(0.06, 0.08)]$m
6 + 2 * (15 - qq[1]) / (qq[2] - qq[1])

qq = eq2[iname == "Clearance-accelerating vaccine" & model == "Growth competition" & istrength %in% c(0.50, 0.52)]$m
50 + 2 * (15 - qq[1]) / (qq[2] - qq[1])

# Time for vaccine to be effective
vt = loadint(paste0("M4-VaccineTime/", filename, "_cm_"), c("twi", "dty", "neu", "whc"), paste0(postfix, ".txt"), dlongnames);
vt$trial = rep(c(seq(1, 1000), seq(1, 1000), seq(1, 1000), seq(1, 1000)), each = 55)
vts = vt[itime >= 0, as.list(c(hdi(R, credMass = 0.95), mean(R))) , by = .(model, itime, istrength, iname)];
names(vts)[5:7] = c("R0", "R1", "R");
vth = vt[itime < 0 & iname != "none", .(R = mean(R)), by = .(model, itime, istrength, iname)];
vt0 = vt[itime == 0.0001 & iname != "none", .(R = mean(R)), by = .(model, itime, istrength, iname)];

F3b = ggplot(pretty_names(vts)) + 
    geom_ribbon(aes(x = itime / 12, ymin = 100 * R0, ymax = 100 * R1, fill = model), alpha = 0.3) + 
    geom_line(aes(x = itime / 12, y = 100 * R, colour = model)) + 
    geom_segment(data = pretty_names(vt0), aes(y = 100 * R, yend = 100 * R, colour = model), x = -0.4, xend = 0) +
    facet_wrap(~iname, scales = "free", nrow = 1) + 
    coord_cartesian(expand = F, ylim = c(0, 16), xlim = c(-0.5, 20)) + 
    labs(x = "Time after intervention (years)", y = "Resistant carriage (%)") + 
    scale_x_continuous(breaks = 2*(0:10)) + 
    scale_y_continuous(breaks = 2*(0:8)) +
    scale_colour_manual(values = Colours) +
    scale_fill_manual(values = Colours) +
    theme(strip.placement = "outside", strip.background = element_rect(fill = NA), 
          legend.position = "none",
          axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
          panel.spacing.x = unit(6, "pt"))
F3b

merge(vts[itime == 240], vth, by = c("model", "istrength", "iname"))[, abs(R.y - R.x), by = c("model", "istrength", "iname")]

# Per-country
cat("Loading per-country...\n");
master_pc = loadint_percountry(paste0("M2-Interventions/", filename, "_cm_"), c("twi", "dty", "neu", "whc"), paste0(postfix, ".txt"), dlongnames, length(r), 10);

baseline_pc = master_pc[iname == "treatment" & istrength == 0];
baseline_pc$iname = "none";
interv_pc = master_pc[istrength > 0 & istrength <= 1.0];
test_pc = diffint(interv_pc, baseline_pc)

# Countries already ordered from lowest to highest consumption
test_pc$cc = countrycodes[test_pc$country + 1]
test_pc$cc = factor(test_pc$cc, levels = countrycodes)

# s. pneumo pneumonia cases change in both R and Y
# see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3700032/ which has per-country estimates of s. pneumo caused pneumonia.
cat("Per-country plot...\n");
test_pc$model = factor(test_pc$model, levels = dlongnames)

tpc = pretty_names(test_pc[iname == "vaccine" & istrength == 0.3])
tpcY = tpc[, mean(disconv * (Y - Y0)), by = .(model, cc)]
F3c1 = ggplot(tpc) + 
    geom_hline(aes(yintercept = 0), colour = "#666666", size = 0.2) + 
    geom_point(data = tpcY, aes(y = V1, x = cc), colour = "black", size = 0.9, shape = "diamond open", stroke = 0.2) + 
    geom_violin(aes(y = disconv * (R - R0), x = cc, fill = model), scale = "width", colour = NA) + 
    facet_wrap(~model, nrow = 1) + 
    labs(x = NULL, y = NULL) + 
    coord_flip(ylim = c(-30, 17)) +
    scale_fill_manual(values = Colours) +
    theme(legend.position = "none", axis.text.y.left = element_text(size = 5),
        panel.background = element_rect(fill = "#f4f4f4"), 
        panel.grid.major.x = element_line(colour = "#ffffff"),
        panel.grid.major.y = element_line(colour = "#ffffff"),
        panel.spacing.x = unit(4, "pt"),
        plot.margin = margin(1, 2, 1, 2, "pt"),
        strip.text = element_text(size = 6),
        strip.background = element_blank(), strip.placement = "outside")

tpc = pretty_names(test_pc[iname == "clearance" & istrength == 0.3])
tpcY = tpc[, mean(disconv * (Y - Y0)), by = .(model, cc)]
F3c2 = ggplot(tpc) + 
    geom_hline(aes(yintercept = 0), colour = "#666666", size = 0.2) + 
    geom_point(data = tpcY, aes(y = V1, x = cc), colour = "black", size = 0.9, shape = "diamond open", stroke = 0.2) + 
    geom_violin(aes(y = disconv * (R - R0), x = cc, fill = model), scale = "width", colour = NA) + 
    facet_wrap(~model, nrow = 1) + 
    labs(x = NULL, y = NULL) + 
    coord_flip(ylim = c(-30, 17)) +
    scale_fill_manual(values = Colours) +
    theme(legend.position = "none", axis.text.y.left = element_text(size = 5),
        panel.background = element_rect(fill = "#f4f4f4"), 
        panel.grid.major.x = element_line(colour = "#ffffff"),
        panel.grid.major.y = element_line(colour = "#ffffff"),
        panel.spacing.x = unit(4, "pt"),
        plot.margin = margin(1, 2, 1, 2, "pt"),
        strip.text = element_text(size = 6),
        strip.background = element_blank(), strip.placement = "outside")

F3 = plot_grid(plot_grid(F3a, F3b, nrow = 2, align = "v", axis = "b", labels = c("a", "b"), label_size = 10, rel_heights = c(1.2, 1)),
          plot_grid(htitle("Acquisition-blocking vaccine"), F3c1, htitle("Clearance-accelerating vaccine"), F3c2,
              htitle(bquote("Change in" ~ .(disname) ~ "cases (per" ~ .(discnt) ~ .(disgroup) ~ "per year)")),
              nrow = 5, rel_heights = c(0.28, 4.58, 0.28, 4.58, 0.28),
              labels = c("c", "", "", "", "", "", ""), label_size = 10),
          ncol = 2, rel_widths = c(1.05, 1.2))

ggsave(figname("4-policy"), F3, width = 18.4, height = 10, units = "cm", device = cairo_pdf);

#
# ALTERNATIVE SETTING PLOT
#

# load and split data set
select("spp")
dlongnames = c("Resistant", "Res. B", "Res. D", "Res. T");
master = loadint(paste0("M6-HighCarriage/K_", filename, "_"), c("twi", "dty", "neu", "whc"), paste0(postfix, ".txt"), dlongnames);
master$istrength = rep(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1), 4000);

interv = master[istrength >= 0 & istrength <= 1.0];
metric = master[iname == "treatment"];
baseline = master[iname == "treatment" & istrength == 0];
baseline$iname = "none";

nones = interv[istrength == 0];
nones$iname = "vaccine"; interv = rbind(interv, nones);
nones$iname = "clearance"; interv = rbind(interv, nones);

# Intervention plots
test = diffint(interv, baseline)
test = test[iname != "treatment"]

test$plotx = test$istrength - 0.0375 + 0.015 * as.numeric(test$model)
test$value = test$Y * 100;
test$variable = "Carriage prevalence (%)";
carr = pretty_names(test);

test$value = test$R * 100;
test$variable = "Resistant carriage (%)";
rcarr = pretty_names(test);

hdicarr = carr[, mean(value, na.rm = T), by = .(variable, plotx, iname, icoverage, istrength, model)]
names(hdicarr)[7] = "mean"

hdircarr = rcarr[, as.list(c(hdi(value, credMass = 0.95), list(m = mean(value, na.rm = T)))), by = .(variable, plotx, iname, icoverage, istrength, model)]
names(hdircarr)[7:9] = c("lo95", "hi95", "mean")

F4 = ggplot() +
    geom_point(data = hdicarr, aes(x = 100 * plotx, y = mean, colour = model, shape = "Total"), size = 0.5, stroke = 0.2) +
    geom_point(data = hdircarr, aes(x = 100 * plotx, y = mean, colour = model), size = 0.05) + 
    geom_linerange(data = hdircarr, aes(x = 100 * plotx, ymin = lo95, ymax = hi95, colour = model), size = 0.3) + 
    facet_wrap(~iname, nrow = 1) + 
    labs(x = "Vaccine efficacy (%)", y = "Carriage prevalence (%)") + 
    scale_x_continuous(breaks = seq(0, 100, by = 10)) + 
    scale_shape_manual(values = "diamond open") +
    scale_colour_manual(breaks = c("Total", "Resistant"), values = Colours) +
    guides(shape = guide_legend(order = 1), colour = guide_legend(order = 2, override.aes = list(colour = "#888888"))) +
    theme(panel.background = element_rect(fill = "#f4f4f4"), 
          panel.grid.major.y = element_line(colour = "#ffffff"), 
          strip.placement = "outside", 
          strip.background = element_rect(fill = NA), 
          legend.position = c(0, 0.32), 
          legend.direction = "vertical", 
          legend.justification = c(0, 0), 
          legend.title = element_blank(),
          legend.key.height = unit(7, "pt"),
          legend.key.width = unit(4, "pt"),
          legend.spacing.y = unit(-3, "pt"),
          axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)))

F4D = ggdraw(F4) +
    draw_grob(rectGrob(gp = gpar(fill = Colours[1], col = NA)), x = 0.112, y = 0.37, width = 0.01, height = 0.02) +
    draw_grob(rectGrob(gp = gpar(fill = Colours[2], col = NA)), x = 0.112, y = 0.31, width = 0.01, height = 0.02) +
    draw_grob(rectGrob(gp = gpar(fill = Colours[3], col = NA)), x = 0.112, y = 0.25, width = 0.01, height = 0.02) +
    draw_grob(rectGrob(gp = gpar(fill = Colours[4], col = NA)), x = 0.112, y = 0.19, width = 0.01, height = 0.02) +
    draw_text("Treatment diversity",   x = 0.138, y = 0.37, hjust = 0, vjust = 0.1, size = 6) +
    draw_text("Pathogen diversity",    x = 0.138, y = 0.31, hjust = 0, vjust = 0.1, size = 6) +
    draw_text("Treatment competition", x = 0.138, y = 0.25, hjust = 0, vjust = 0.1, size = 6) +
    draw_text("Growth competition",    x = 0.138, y = 0.19, hjust = 0, vjust = 0.1, size = 6) 

ggsave(figname("5-highburden"), F4D, width = 9, height = 5, units = "cm", device = cairo_pdf);


#
# DIAGNOSTICS
#

stop();

cat("Diagnostics...\n");
pweights = populations / mean(populations);
sink(paste0("~/Dropbox/Elegy/Submission/STM Response/diagnostics-", filename, postfix, "-cm.txt"));
if (postfix %like% "tau") {
    summarise.WAICs(dlist, dnames, "binomial", pweights)
} else {
    summarise.WAICs(dlist, dnames, "normal.trunc", pweights);
}
summarise.likelihoods(dlist, dnames);
diagnose(dlist, dnames);
sink(NULL);


#
# VERIFYING D-TYPES MODEL
#
library(stringr)
d = fread("~/Dropbox/Elegy/Runs/M5-Analysis/dtypes.txt", skip = "trial\t");
d = melt(d, id.vars = 1:10)
d = d[!(variable %like% "X") & !(variable %like% "Rfrac")]
the_d = str_locate(d$variable, "d")[, 1]

d[, compartment := str_sub(variable, 1, 1)]
d[, dtype := as.numeric(str_sub(variable, 2, the_d - 1))]
d[, country := as.numeric(str_sub(variable, the_d + 1))]
drows = 2 * 25 * 27
d[, trial := rep(1:(nrow(d) / drows), drows)]

abundance = d[, .(value = sum(value)), by = .(country, dtype, trial)]
dp1 = ggplot(abundance) + geom_violin(aes(x = as.factor(dtype), y = 100 * value), fill = "grey", scale = "width") + theme(legend.position = "none") + labs(x = "D-type", y = "Carriage prevalence (%)")

rfrac = d[compartment == "R", .(value = sum(value)), by = .(country, dtype, trial)]
rfrac$value = rfrac$value / abundance$value
rfrac = rfrac[, .(value = mean(value)), by = .(country, dtype)]
dp2 = ggplot(rfrac) + geom_line(aes(x = dtype, y = 100 * value, colour = country, group = country)) + labs(x = "D-type", y = "Resistance frequency (%)", colour = "Country ID") + scale_x_continuous(limits = c(0, 24), breaks = 0:24, expand = expand_scale(0, 0.6))

plot_grid(dp1, dp2, nrow = 2, align = "v", axis = "tblr")

ggsave("~/Dropbox/Elegy/Figures/Supplementary/S-dtypes.pdf", width = 10, height = 12, units = "cm")



#
# IMPACT OF PRIORS
#

compare = function(posterior, modelname, sweepname)
{
    d = fread(paste0("~/Dropbox/Elegy/Runs/M7-Variation/sweep_spp_cm_", modelname, "_", sweepname, ".txt"), skip = "trial\t");
    d = d[, .SD, .SDcols = c("lp", sweepname)];
    d[, nlp := exp(lp - max(lp))];
    firstbin = d[[sweepname]][1];
    bw = d[[sweepname]][2] - firstbin;
    
    ggplot() + 
        geom_histogram(data = posterior, aes_string(sweepname, "stat(ndensity)"), 
                       fill = "#ccddff", bins = 25, center = firstbin) +
        geom_col(data = d, aes_string(sweepname, "nlp"), width = bw / 5, fill = "#229955") +
        labs(y = "") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
}

p1a = compare(d_neu, "neu", "beta")
p1b = compare(d_neu, "neu", "c")
p1c = compare(d_neu, "neu", "k")

p2a = compare(d_whc, "whc", "beta")
p2b = compare(d_whc, "whc", "b")
p2c = compare(d_whc, "whc", "k")

p3a = compare(d_dty, "dty", "beta")#
p3b = compare(d_dty, "dty", "c")#?
p3c = compare(d_dty, "dty", "a")
p3d = compare(d_dty, "dty", "delta")##

p4a = compare(d_twi, "twi", "beta")#?
p4b = compare(d_twi, "twi", "c")
p4c = compare(d_twi, "twi", "g")#~
p4d = compare(d_twi, "twi", "shape")

pg = plot_grid(p4a, p4b, p4c, p4d,
               p3a, p3b, p3c, p3d,
               p1a, p1b, p1c, p1a,
               p2a, p2b, p2c, p1a,
          ncol = 4, nrow = 4)

ggsave("~/Dropbox/Elegy/Submission/STM Response/Figures/pprior.pdf", pg, width = 10, height = 10, units = "cm")

