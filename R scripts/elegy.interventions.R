library(data.table)
library(ggplot2)
library(ggstance)
library(cowplot)

IP = "~/Documents/Mechanisms/Elegy/Runs/M2-Interventions/";

rbind0 = function(a, b)
{
    if (is.null(a)) {
        return (b)
    } else {
        return (rbind(a, b))
    }
}

alpha.to.colour = function(plot, attributes = c("colour", "fill"), mode = "lighten")
{
    q = ggplot_build(plot)
    for (i in 1:length(q$data))
    {
        for (attribute in attributes)
        {
            if (all(c(attribute, "alpha") %in% names(q$data[[i]])) & all(is.numeric(q$data[[i]]$alpha)))
            {
                col = col2rgb(unlist(q$data[[i]][[attribute]]))
                alpha = unlist(q$data[[i]][["alpha"]])
                if (mode == "lighten") {
                    col = col * alpha + 255 * (1 - alpha)
                } else if (mode == "darken") {
                    col = col * alpha
                } else if (mode == "luminance") {
                    col = col * abs(alpha - 1) + ifelse(alpha > 0.5, 510 * (1 - alpha), 0)
                }
                q$data[[i]][[attribute]] = rgb(t(col), maxColorValue = 255)
                q$data[[i]][["alpha"]] = 1
            }
        }
    }
    plot(ggplot_gtable(q))
}

# load intervention set
loadint = function(prefix, models, postfix, dlongnames, cut = NA)
{
    d = NULL;
    for (m in 1:length(models))
    {
        e = fread(paste0(IP, prefix, models[m], postfix), skip = "\t");
        if (!is.na(cut[m])) {
            e = e[floor(cut[m] * nrow(e) + 1):nrow(e),];
        }
        e[, c("iname", "icoverage", "istrength", "itime") := tstrsplit(trial, split = ';', fixed = T)];
        e[, icoverage := as.numeric(icoverage)];
        e[, istrength := as.numeric(istrength)];
        e[, itime := as.numeric(itime)];
        e[, model := models[m]];
        e[, Y := S + R + SR];
        d = rbind0(d, e[, .(Rfrac, Y, iname, icoverage, istrength, model)]);
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
        e = fread(paste0(IP, prefix, models[m], postfix), skip = "\t");
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
            d = rbind0(d, e[, .(country, Rfrac, Y, iname, icoverage, istrength, model)]);
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
        labels = c("None", "Transmission-blocking vaccine", "Clearance-hastening vaccine", "Reduced treatment", "Reduced treatment (targeted)", "Reduced treatment (targeted within-country)"));
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
baseint = function(d, metric, i_compare, i_strengths, x, interp.method = "hyman", tol = 2e-6)
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
        pts = lapply(pts, function(x) { x[x < tol] = 0; return (x) });
        xi = mapply(function(pts, xi) splinefun(pmax(pts, tol), seq_along(pts), interp.method)(xi), pts = pts, xi = co[[x]]);
        
        # transform into levels
        eq = approx(seq_along(unique(metric$istrength)), unique(metric$istrength), xi, rule = 2);
        
        D = rbind0(D, data.table(x = eq$x, y = eq$y, iname = i_compare, istrength = i_str, model = co$model));
    }
    return (D);
}

check.monotone = function(L)
{
    for (pts in L)
    {
        pts = pmax(pts, 2e-6);
        if (is.unsorted(pts) & is.unsorted(rev(pts))) {
            print(pts);
        }
    }
}

htitle = function(text, size = 6, x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5)
{
    ggdraw() + draw_label(text, x = x, y = y, size = size, hjust = hjust, vjust = vjust)
}

# make selection
source("~/Dropbox/Elegy/elegy.select.R")
select("ecp_tau")

# load and split data set
dlongnames = c("Within-host competition A", "Within-host competition B", "Pathogen variability", "Population heterogeneity");

# master = loadint("spp_", c("neu_tau_long", "whc_tau_long", "dty_tau", "twi_tau"), ".txt", dlongnames, cut = c(0.5, 0.5, 0, 0));
master = loadint(paste0(filename, "_"), c("neu", "whc", "dty", "twi"), paste0(postfix, ".txt"), dlongnames, cut = cuts);
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

test = diffint(interv, baseline)
test$plotx = test$istrength - 0.045 + 0.015 * as.numeric(test$model)

theme_set(theme_cowplot(font_size = 6, font_family = "Helvetica", line_size = 0.25))

test$value = test$Y * 100; test$variable = "Carriage prevalence (%)";
triplicate = test;
test$value = test$Rfrac * 100; test$variable = "Resistance frequency (%)";
triplicate = rbind(triplicate, test);
test$value = test$R * 100; test$variable = "Resistant carriage (%)";
triplicate = pretty_names(rbind(triplicate, test));

summt = triplicate[, .(quantile(value, 0.1, na.rm = T), quantile(value, 0.5, na.rm = T), quantile(value, 0.9, na.rm = T)), by = .(variable, plotx, iname, icoverage, istrength, model)]
hdit = triplicate[, as.list(c(hdi(value, credMass = 0.9), hdi(value, credMass = 0.5), list(m = mean(value, na.rm = T)))), by = .(variable, plotx, iname, icoverage, istrength, model)]
names(hdit)[7:11] = c("lo90", "hi90", "lo50", "hi50", "mean")

# ggplot(triplicate) + geom_violin(aes(x = 100 * plotx, y = value, fill = model, group = plotx), colour = NA, scale = "width") + facet_grid(variable ~ iname, switch = "y", scales = "free_y") + labs(x = "Intervention efficacy", y = NULL) + scale_x_continuous(breaks = seq(0, 100, by = 10)) + theme(panel.background = element_rect(fill = "#f4f4f4"), panel.grid.major.y = element_line(colour = "#ffffff"), strip.placement = "outside", strip.background = element_rect(fill = NA), legend.position = "bottom", legend.direction = "horizontal", legend.justification = c(0.5, 0.5), legend.title = element_blank())
# 
# ggsave("~/Dropbox/Elegy/Figures/2-interventions-violins-B.pdf", width = 17.8, height = 10, units = "cm");

ggplot(hdit) + geom_point(aes(x = 100 * plotx, y = mean, colour = model), size = 0.05) + geom_linerange(aes(x = 100 * plotx, ymin = lo90, ymax = hi90, colour = model), size = 0.3) + facet_grid(variable ~ iname, switch = "y", scales = "free_y") + labs(x = "Intervention efficacy", y = NULL) + scale_x_continuous(breaks = seq(0, 100, by = 10)) + theme(panel.background = element_rect(fill = "#f4f4f4"), panel.grid.major.y = element_line(colour = "#ffffff"), strip.placement = "outside", strip.background = element_rect(fill = NA), legend.position = "bottom", legend.direction = "horizontal", legend.justification = c(0.5, 0.5), legend.title = element_blank())

ggsave(figname("2-interv"), width = 17.8, height = 10, units = "cm");

# ggplot(hdit) + geom_linerange(aes(x = 100 * plotx, ymin = lo50, ymax = hi50, colour = model), size = 0.5) + geom_linerange(aes(x = 100 * plotx, ymin = lo90, ymax = hi90, colour = model), size = 0.2) + facet_grid(variable ~ iname, switch = "y", scales = "free_y") + labs(x = "Intervention efficacy", y = NULL) + scale_x_continuous(breaks = seq(0, 100, by = 10)) + theme(panel.background = element_rect(fill = "#f4f4f4"), panel.grid.major.y = element_line(colour = "#ffffff"), strip.placement = "outside", strip.background = element_rect(fill = NA), legend.position = "bottom", legend.direction = "horizontal", legend.justification = c(0.5, 0.5), legend.title = element_blank())
# 
# ggsave("~/Dropbox/Elegy/Figures/2-interventions-lineribbon-B.pdf", width = 17.8, height = 10, units = "cm");



# Equivalence plots
eqV = baseint(master, metric, "vaccine", c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), "R", "hyman", 2e-6);
eqC = baseint(master, metric, "clearance", c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), "R", "hyman", 2e-6);
eqH = baseint(master, metric, "hightreatment_between", c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), "R", "hyman", 2e-6);
eq = rbind(eqV, eqC, eqH)
eq$rx = 100 * eq$y
eq2 = eq[, as.list(c(quantile(rx, c(0.05, 0.25, 0.75, 0.95), na.rm = T), list(m = mean(rx)))), by = .(iname, istrength, model)]

names(eq2) = c("iname", "istrength", "model", "q5", "q25", "q75", "q95", "m")

#ggplot(pretty_names(eq2), aes(x = istrength)) + annotate("ribbon", x = c(0,1), ymin = c(8.8, 8.8), ymax = c(23.1, 23.1), fill = "#dddddd") + geom_ribbon(aes(ymin = q5, ymax = q95, fill = model), alpha = 0.5) + geom_ribbon(aes(ymin = q25, ymax = q75, fill = model)) + facet_grid(iname ~ model) + labs(x = "Intervention efficacy", y = "Prescriptions equivalent (%)") + theme(legend.position = "none", strip.placement = "outside", strip.background = element_rect(fill = NA))

ggplot(pretty_names(eq2), aes(x = istrength)) +
    annotate("ribbon", x = c(0,1), ymin = c(8.8, 8.8), ymax = c(23.1, 23.1), fill = "#dddddd") + 
    annotate("line", x = c(0, 1), y = c(15, 15), linetype = "22", size = 0.2) +
    geom_line(aes(y = m, colour = model, group = model)) + 
    facet_wrap(~iname, ncol = 1, scales = "free_x") + 
    labs(x = "Intervention efficacy", y = "Prescriptions equivalent (%)") + 
    theme(strip.placement = "outside", strip.background = element_rect(fill = NA), 
        legend.position = c(0.21, 0.94), legend.direction = "vertical", legend.justification = c(0.5, 0.5), legend.title = element_blank(), 
        legend.text = element_text(margin = margin(t = -15, b = -15, unit = "pt")), legend.key.height = unit(8, "pt"))

ggsave(figname("3-equiv"), width = 8.7, height = 12, units = "cm");

ggplot(pretty_names(eq2), aes(x = istrength)) + 
    annotate("ribbon", x = c(0,1), ymin = c(8.8, 8.8), ymax = c(23.1, 23.1), fill = "#dddddd") + 
    annotate("line", x = c(0, 1), y = c(15, 15), linetype = "22", size = 0.2) +
    geom_line(aes(y = m, colour = model, group = model)) + 
    facet_wrap(~iname, nrow = 1, scales = "free_x") + 
    labs(x = "Intervention efficacy", y = "Prescriptions equivalent (%)") + 
    theme(strip.placement = "outside", strip.background = element_rect(fill = NA), legend.position = "none")

ggsave(figname("3-equivh"), width = 17.8, height = 5, units = "cm");



# per-country
master_pc = loadint_percountry(paste0(filename, "_"), c("neu", "whc", "dty", "twi"), paste0(postfix, ".txt"), dlongnames, length(r), 10, cut = cuts);

baseline_pc = master_pc[iname == "treatment" & istrength == 0];
baseline_pc$iname = "none";
interv_pc = master_pc[istrength > 0 & istrength <= 1.0];
test_pc = diffint(interv_pc, baseline_pc)

# per-country for within-country targeted reductions
wmaster_pc = loadint_percountry(paste0("within_", filename, "_"), "twi", ".txt", "Population heterogeneity", length(r), 10, cut = 0);
wbaseline_pc = wmaster_pc[iname == "hightreatment_within" & istrength == 0];
wbaseline_pc$iname = "none";
winterv_pc = wmaster_pc[istrength > 0 & istrength <= 1.0];
wtest_pc = diffint(winterv_pc, wbaseline_pc)

# order by model-estimated resistance prevalence -- after rerunning with correct Rfrac, check this is right
model_rfrac = test_pc[, mean(Rfrac0), by = country]$V1
test_pc$cc = countrycodes[test_pc$country + 1]
test_pc$cc = factor(test_pc$cc, levels = countrycodes[order(model_rfrac)])
wtest_pc$cc = countrycodes[wtest_pc$country + 1]
wtest_pc$cc = factor(wtest_pc$cc, levels = countrycodes[order(model_rfrac)])

# s. pneumo pneumonia cases change in both R and Y
# see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3700032/ which has per-country estimates of s. pneumo caused pneumonia. UK: 0.6 per 1000 child-years
plot_grid(plot_grid(ggdraw(), htitle(bquote("Per-country impact of 10% efficacy" ~ italic(.(species)) ~ "vaccine")), htitle("Per-country impact of 10% reduction in penicillin use"),
                nrow = 1, rel_widths = c(1, 8.8, 8.8)),
    ggplot(pretty_names(test_pc[istrength == 0.1])) + 
        geom_vline(aes(xintercept = 0), colour = "#666666", size = 0.2) + 
        geom_violinh(aes(x = disconv * (Y - Y0), y = cc), fill = "grey", scale = "width", colour = NA) + 
        geom_violinh(aes(x = disconv * (R - R0), y = cc, fill = model), scale = "width", colour = NA, alpha = 0.75) + 
        facet_grid(model ~ iname, scales = "free") + 
        labs(x = bquote("Change in" ~ .(disname) ~ "cases (per" ~ 10^6 ~ .(disgroup) ~ "per year)"), y = "Country") + 
        theme(legend.position = "none", axis.text.y.left = element_text(size = 5), panel.background = element_rect(fill = "#f4f4f4"), 
            panel.grid.major.x = element_line(colour = "#ffffff"), strip.background = element_blank()),
    nrow = 2, rel_heights = c(0.5, 19.5))

ggsave(figname("4-country"), width = 17.8, height = 20, units = "cm");

