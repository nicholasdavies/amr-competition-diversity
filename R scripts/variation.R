library(data.table)
library(ggplot2)
library(cowplot)
library(binom)
library(Cairo)
library(stringr)
library(jcolors)

Mod1 = "Treatment diversity";
Mod2 = "Pathogen diversity";
Mod3 = "Treatment competition";
Mod4 = "Growth competition";
Colours = rev(hcl(h = seq(15, 375, length = 5), l = 65, c = 100)[1:4])

loadvar = function(model, variable, varcol)
{
    if (variable == "S") {
        v = fread(paste0("~/Dropbox/Elegy/Runs/M7-Variation/spp_cm_", model, ".txt"), skip = "trial\t");
    } else {
        v = fread(paste0("~/Dropbox/Elegy/Runs/M7-Variation/spp_cm_", model, "_", variable, ".txt"), skip = "trial\t");
    }

    if (model == "twi") {
        v2 = v[, 1:which(names(v) == "epochs")];
        for (ci in 0:26)
        {
            groups = (ci * 10):(ci * 10 + 9);
            set(v2, NULL, paste0("X", ci), rowSums(v[, which(names(v) %in% paste0("X", groups)), with = F]) / 10);
            set(v2, NULL, paste0("S", ci), rowSums(v[, which(names(v) %in% paste0("S", groups)), with = F]) / 10);
            set(v2, NULL, paste0("R", ci), rowSums(v[, which(names(v) %in% paste0("R", groups)), with = F]) / 10);
            set(v2, NULL, paste0("Rfrac", ci), v2[[paste0("R", ci)]] / (v2[[paste0("S", ci)]] + v2[[paste0("R", ci)]]));
        }
        v = v2;
    }

    v = melt(v, id.vars = "lp", measure.vars = patterns(paste0("^Rfrac[0-9]|^X[0-9]|^", varcol, "[0-9]")));
    v = data.table("variable" = v[variable %like% paste0("^", varcol, "[0-9]"), value], 
              "Rfrac" = v[variable %like% "^Rfrac[0-9]", value],
              "X" = v[variable %like% "^X[0-9]", value]);
    v$tau = c(4.3,  4.5,  4.9,  6.0,  6.1,  6.2,  6.6,  6.9,  7.0,  7.2,  7.7,  8.7,  9.6,  9.8, 10.0, 10.2, 10.5, 11.3, 11.3, 11.5, 12.1, 12.6, 13.6, 13.9, 14.6, 15.0, 15.8);
    v$r   = c( 19,    2,    0,   69,    0,   16,   10,    8,   33,   31,   70,    0,    5,    3,   33,   34,    6,   25,   33,   76,  192,    2,  143,    3,   43,  226,    5);
    v$n   = c(939,   72,   64,  522,   31,  313,  616,  205,  146, 1028, 1744,   13,   32,   42,  195, 1030,   21,  137,  202,  435,  860,   35, 1511,   67,  291,  663,   15);
    v$cc  = c('NL', 'DE', 'EE', 'FI', 'LV', 'AT', 'NO', 'CZ', 'HU', 'SE', 'GB', 'MT', 'BG', 'IS', 'SI', 'DK', 'PL', 'HR', 'PT', 'IE', 'ES', 'LU', 'BE', 'LT', 'IT', 'FR', 'CY');
    names(v)[1] = variable;
    v$cc = factor(v$cc, levels = v$cc);
    v$model = model;
    
    v = cbind(v, binom.confint(v$r, v$n, methods = "bayes", prior.shape1 = 1, prior.shape2 = 1)[, 4:6]);
    v[, mode := r / n];
    return (v)
}

plotvar = function(model, variable, col, varcol = NULL, title = NULL, adj = 0)
{
    if (is.null(varcol)) varcol = variable;
    v = loadvar(model, variable, varcol);
    
    p1 = ggplot(v) + 
        geom_pointrange(aes(x = cc, y = 100 * mean, ymin = 100 * lower, ymax = 100 * upper), fatten = 0.2) + 
        geom_line(aes(x = cc, y = 100 * Rfrac, group = 1), colour = col) +
        labs(x = NULL, y = "Resistance\nfrequency (%)") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5));
    
    if (adj > 0) {
        v[[variable]] = adj / v[[variable]];
    }
    
    p2 = ggplot(v) + 
        geom_point(aes_string(x = "cc", y = variable), size = 1) +
        labs(x = NULL, y = ifelse(is.null(title), variable, title)) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
              panel.grid.major.y = element_line(colour = "white"), panel.background = element_rect(fill = "grey95"));
    
    if (is.null(title)) {
        plot_grid(p1, p2, ncol = 1, align = "v")
    } else {
        plot_grid(ggdraw() + draw_text(str_replace(title, "\n", " "), size = 6), p1, p2, 
                  ncol = 1, align = "v", rel_heights = c(1, 5, 5))
    }
}

loadvacc = function(model, variable = NULL)
{
    v = fread(paste0(paste(c("~/Dropbox/Elegy/Runs/M7-Variation/vacc_spp_cm", model, variable), collapse = "_"), ".txt"), skip = "trial\t");

    if (model == "twi") {
        v2 = v[, 1:which(names(v) == "epochs")];
        for (ci in 0:26)
        {
            groups = (ci * 10):(ci * 10 + 9);
            set(v2, NULL, paste0("X", ci), rowSums(v[, which(names(v) %in% paste0("X", groups)), with = F]) / 10);
            set(v2, NULL, paste0("S", ci), rowSums(v[, which(names(v) %in% paste0("S", groups)), with = F]) / 10);
            set(v2, NULL, paste0("R", ci), rowSums(v[, which(names(v) %in% paste0("R", groups)), with = F]) / 10);
            set(v2, NULL, paste0("Rfrac", ci), v2[[paste0("R", ci)]] / (v2[[paste0("S", ci)]] + v2[[paste0("R", ci)]]));
        }
        v = v2;
    }

    v[, c("iname", "icoverage", "istrength", "itime") := tstrsplit(trial, ";", fixed = T, type.convert = T)]
    return (v[, S:itime])
}

loadvacc2 = function(model, variables)
{
    v = cbind(loadvacc(model), variable = "base");
    for (variable in variables) {
        v = rbind(v, cbind(loadvacc(model, variable), variable = ifelse(variable == "bigG", "f", ifelse(variable == "g" & model != "twi", "f", ifelse(variable == "shape", "kappa", variable)))));
    }
    v[, X := 1 - S - R - SR];
    setcolorder(v, c("X", "Rfrac"));
    
    v = melt(v, id.vars = c("iname", "icoverage", "istrength", "itime", "variable"),
             measure.vars = patterns("^Rfrac[0-9]*$", "^X[0-9]*$"),
             variable.name = "country", value.name = c("Rfrac", "X"));
    
    ccs = c('Europe', 'NL', 'DE', 'EE', 'FI', 'LV', 'AT', 'NO', 'CZ', 'HU', 'SE', 'GB', 'MT', 'BG', 'IS', 'SI', 'DK', 'PL', 'HR', 'PT', 'IE', 'ES', 'LU', 'BE', 'LT', 'IT', 'FR', 'CY');
    v$cc = factor(ccs[v$country], levels = ccs);
    v[iname == "vaccine", iname := "Acquisition-blocking"];
    v[iname == "clearance", iname := "Clearance-accelerating"];
    
    v$model = model;
    
    return (v)
}

plotvacc = function(model, variables, title)
{
    v = loadvacc2(model, variables);

    colours = c("base" = "#aaaaaa", "beta" = "red", "c" = "orange", "b" = "orange",
                "u" = "#aaaa00", "k" = "green", "g" = "#000088", "z" = "#008888",
                "a" = "#0088ff", "delta" = "#880088", "f" = "blue", "kappa" = "purple");
    
    v$base = ifelse(v$variable == "base", "yes", "no");
    
    v$variable = factor(v$variable, levels = c("base", "c", "b", "beta", "u", "f", "z", "k", "a", "delta", "g", "kappa"), ordered = T);
    
    p1 = ggplot(v[iname == "Acquisition-blocking" & X < 0.999]) + 
        geom_line(aes(x = istrength * 100, y = Rfrac * (1 - X) * 100,
                      group = interaction(variable, iname), colour = variable, linetype = base, size = base)) + 
        facet_wrap(~cc, nrow = 4, scales = "free_y") + 
        scale_colour_manual(values = colours, guide = guide_legend(nrow = 1)) +
        scale_linetype_manual(values = c("yes" = "solid", "no" = "dashed"), guide = "none") +
        scale_size_manual(values = c("yes" = 1, "no" = 0.3), guide = "none") +
        labs(y = "Resistant carriage (%)", x = "Vaccine efficacy (%)", 
             colour = "Parameter", linetype = "Vaccine type", title = "Acquisition-blocking vaccine") + 
        theme(legend.position = "none")

    p2 = ggplot(v[iname == "Clearance-accelerating" & X < 0.999]) + 
        geom_line(aes(x = istrength * 100, y = Rfrac * (1 - X) * 100,
                      group = interaction(variable, iname), colour = variable, linetype = base, size = base)) + 
        facet_wrap(~cc, nrow = 4, scales = "free_y") + 
        scale_colour_manual(values = colours, guide = guide_legend(nrow = 1)) +
        scale_linetype_manual(values = c("yes" = "solid", "no" = "dashed"), guide = "none") +
        scale_size_manual(values = c("yes" = 1, "no" = 0.3), guide = "none") +
        labs(y = "Resistant carriage (%)", x = "Vaccine efficacy (%)", 
             colour = "Parameter", linetype = "Vaccine type", title = "Clearance-accelerating vaccine") + 
        theme(legend.position = "bottom", legend.justification = "center")
    
    plot_grid(ggdraw() + draw_text(title, size = 10), p1, p2, nrow = 3, rel_heights = c(1, 10, 11))
}

plotvacc_example = function(vacctype, showylabel)
{
    v1 = loadvacc2("twi", c("c", "z"));
    v2 = loadvacc2("dty", c("c", "delta", "z"));
    v3 = loadvacc2("neu", c("beta", "c", "u", "k", "z"));
    v4 = loadvacc2("whc", c("b", "z"));
    
    v = rbind(v1, v2, v3, v4);
    
    v[model == "twi", model := Mod1];
    v[model == "dty", model := Mod2];
    v[model == "neu", model := Mod3];
    v[model == "whc", model := Mod4];
    v$model = factor(v$model, levels = unique(v$model));
    
    v[variable == "delta", variable := "δ"];
    v[variable == "beta", variable := "β"];

    v$variable = factor(v$variable, levels = c("base", "c", "b", "β", "u", "f", "z", "k", "a", "δ", "g", "κ"), ordered = T);
    v$isbase = ifelse(v$variable == "base", "yes", "no");
    
    v[is.nan(Rfrac), Rfrac := 0];
    v[is.nan(X), X := 0];
    
    ggplot(v[cc == "Europe" & iname != "treatment" & istrength <= 0.6 & iname == vacctype]) + 
        geom_line(aes(x = istrength * 100, y = Rfrac * (1 - X) * 100,
                      group = interaction(variable, iname), colour = variable, linetype = isbase, size = isbase)) + 
        facet_wrap(~model, scales = "free", ncol = 1) + 
        labs(y = ifelse(showylabel, "Resistant carriage (%)", ""), x = "Vaccine efficacy (%)", 
             colour = NULL, linetype = NULL, title = vacctype) +
        scale_linetype_manual(values = c("yes" = "solid", "no" = "21"), guide = "none") +
        scale_size_manual(values = c("yes" = 0.8, "no" = 0.3), guide = "none") +
        scale_colour_manual(values = c("base" = "#b0b0b0", "β" = "#943CB4", "δ" = "#3CA437", "c" = "#5B6DC8", "b" = "#5B6DC8", "k" = "#3CA437", "u" = "#6ACDC5", "z" = "#D33B44"), guide = "none") +
        theme(legend.position = c(0, 0.98), strip.background = element_blank(), strip.text = element_blank(),
            legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.41, "cm"), plot.title = element_text(face = "plain", size = 6))
}

ploterror = function(mv)
{
    v = NULL;
    for (model in names(mv)) {
        for (variable in mv[[model]]) {
            if (model == "twi") {
                if (variable %in% c("beta", "c", "u")) {
                    w = loadvar(model, variable, "x");
                } else if (variable == "bigG") {
                    w = loadvar(model, "bigG", "G");
                } else {
                    w = loadvar(model, variable, variable);
                }
            } else {
                w = loadvar(model, variable, variable);
            }
            names(w)[1] = "value";
            if (variable == "g") {
                w$variable = "f";
            } else if (variable == "bigG") {
                w$variable = "g";
            } else if (variable == "beta") {
                w$variable = "β";
            } else if (variable == "delta") {
                w$variable = "δ";
            } else if (variable == "shape") {
                w$variable = "κ";
            } else {
                w$variable = variable;
            }
            v = rbind(v, w);
        }
    }
    
    v[, error := abs(mean - Rfrac)];
    v = v[, .(ME = max(100*error)), by = .(model, variable)];
    variables = c("c", "b", "β", "u", "f", "z", "k", "a", "δ", "g", "κ");
    v[, variable := factor(variable, levels = rev(variables))];
    v[, model := factor(model, levels = unique(model))];
    v[, ypos := as.integer(variable)];# + as.integer(model) / 5 - 0.5];
    
    invalidcl = "#999999";
    thecolours = c(
        "neu.c" = "#5B6DC8", "neu.β" = "#943CB4", "neu.u" = "#6ACDC5", "neu.z" = "#D33B44", "neu.k" = "#3CA437", "neu.f" = invalidcl,
        "whc.b" = "#5B6DC8", "whc.β" = invalidcl, "whc.u" = invalidcl, "whc.z" = "#D33B44", "whc.k" = invalidcl, "whc.f" = invalidcl, 
        "dty.c" = "#5B6DC8", "dty.β" = invalidcl, "dty.u" = invalidcl, "dty.z" = "#D33B44", "dty.a" = invalidcl, "dty.δ" = "#3CA437", "dty.f" = invalidcl,
        "twi.c" = "#5B6DC8", "twi.β" = invalidcl, "twi.u" = invalidcl, "twi.z" = "#D33B44", "twi.g" = invalidcl, "twi.κ" = invalidcl, "twi.f" = invalidcl);
    
    ggplot(v) + 
        geom_point(aes(x = ME, y = variable, colour = interaction(model, variable)), shape = "diamond") +
        theme(panel.grid.major.y = element_line(colour = "white"), panel.background = element_rect(fill = "grey95"), 
            legend.position = "none", strip.text = element_blank(), strip.background = element_blank()) +
        xlim(0, max(v$ME)) +
        facet_wrap(~model, ncol = 1, drop = T, scales = "free") +
        scale_colour_manual(values = thecolours) +
        labs(y = NULL, x = "Predictive error (%)", colour = NULL)
}


# Plot setup
theme_set(theme_cowplot(font_size = 7, font_family = "Helvetica", line_size = 0.25))


# Fits
p = plot_grid(ggdraw() + draw_text(Mod3, size = 10), plot_grid(
    plotvar("neu", "beta", "red", title = "Transmission rate,\nbeta (per month)"),
    plotvar("neu", "c", "orange", title = "Transmission cost\nof resistance, c"),
    plotvar("neu", "u", "#888800", title = "Clearance rate, u\n(per month)"),
    plotvar("neu", "k", "green", title = "Relative rate of\nco-colonisation, k"),
    plotvar("neu", "g", "blue", title = "Assortativity of\ncountries, f"),
    plotvar("neu", "z", "#008888", title = "Defined daily doses\nper course, z", adj = 5), nrow = 2), ncol = 1, rel_heights = c(0.5, 11.5))
ggsave("~/Dropbox/Elegy/Figures/Supplementary/variation-neu.pdf", p, width = 20, height = 12, units = "cm")

p = plot_grid(ggdraw() + draw_text(Mod4, size = 10), plot_grid(
    plotvar("whc", "beta", "red", title = "Transmission rate,\nbeta (per month)"),
    plotvar("whc", "b", "orange", title = "Growth rate of\nsensitive strain, b"),
    plotvar("whc", "u", "#888800", title = "Clearance rate, u\n(per month)"),
    plotvar("whc", "k", "green", title = "Relative rate of\nco-colonisation, k"),
    plotvar("whc", "g", "blue", title = "Assortativity of\ncountries, f"),
    plotvar("whc", "z", "#008888", title = "Defined daily doses\nper course, z", adj = 5), nrow = 2), ncol = 1, rel_heights = c(0.5, 11.5))
ggsave("~/Dropbox/Elegy/Figures/Supplementary/variation-whc.pdf", p, width = 20, height = 12, units = "cm")

p = plot_grid(ggdraw() + draw_text(Mod2, size = 10), plot_grid(
    plotvar("dty", "beta", "red", title = "Transmission rate,\nbeta (per month)"),
    plotvar("dty", "c", "orange", title = "Transmission cost\nof resistance, c"),
    plotvar("dty", "u", "#888800", title = "Clearance rate, u\n(per month)"),
    plotvar("dty", "a", "#0088ff", title = "Power of diversifying\nselection, a"),
    plotvar("dty", "delta", "#888888", title = "Range of clearance\nrates, delta"),
    plotvar("dty", "g", "blue", title = "Assortativity of\ncountries, f"),
    plotvar("dty", "z", "#008888", title = "Defined daily doses\nper course, z", adj = 5), nrow = 3), ncol = 1, rel_heights = c(0.5, 16.75))
ggsave("~/Dropbox/Elegy/Figures/Supplementary/variation-dty.pdf", p, width = 20, height = 17.25, units = "cm")

p = plot_grid(ggdraw() + draw_text(Mod1, size = 10), plot_grid(
    plotvar("twi", "beta", "red", "x", title = "Transmission rate,\nbeta (per month)"),
    plotvar("twi", "c", "orange", "x", title = "Transmission cost\nof resistance, c"),
    plotvar("twi", "u", "#888800", "x", title = "Clearance rate, u\n(per month)"),
    plotvar("twi", "g", "#8800ff", title = "Assortativity of\nsubpopulations, g"),
    plotvar("twi", "shape", "purple", title = "Shape of treatment\nrate distribution, kappa"),
    plotvar("twi", "bigG", "blue", "G", title = "Assortativity of\ncountries, f"),
    plotvar("twi", "z", "#008888", title = "Defined daily doses\nper course, z", adj = 5), nrow = 3), ncol = 1, rel_heights = c(0.5, 16.75))
ggsave("~/Dropbox/Elegy/Figures/Supplementary/variation-twi.pdf", p, width = 20, height = 17.25, units = "cm")



# Impact of vaccination
p = plotvacc("neu", c("g"), Mod3)
ggsave("~/Dropbox/Elegy/Figures/Supplementary/S-variation-neu-other.pdf", p, width = 20, height = 20, units = "cm")

p = plotvacc("whc", c("beta", "u", "k", "g"), Mod4)
ggsave("~/Dropbox/Elegy/Figures/Supplementary/S-variation-whc-other.pdf", p, width = 20, height = 20, units = "cm")

p = plotvacc("dty", c("beta", "u", "a", "g"), Mod2)
ggsave("~/Dropbox/Elegy/Figures/Supplementary/S-variation-dty-other.pdf", p, width = 20, height = 20, units = "cm")

p = plotvacc("twi", c("beta", "u", "g", "bigG", "shape"), Mod1)
ggsave("~/Dropbox/Elegy/Figures/Supplementary/S-variation-twi-other.pdf", p, width = 20, height = 20, units = "cm")


p = plotvacc("neu", c("beta", "c", "u", "k", "z"), Mod3)
ggsave("~/Dropbox/Elegy/Figures/Supplementary/S-variation-neu-fit.pdf", p, width = 20, height = 20, units = "cm")

p = plotvacc("whc", c("b", "z"), Mod4)
ggsave("~/Dropbox/Elegy/Figures/Supplementary/S-variation-whc-fit.pdf", p, width = 20, height = 20, units = "cm")

p = plotvacc("dty", c("c", "delta", "z"), Mod2)
ggsave("~/Dropbox/Elegy/Figures/Supplementary/S-variation-dty-fit.pdf", p, width = 20, height = 20, units = "cm")

p = plotvacc("twi", c("c", "z"), Mod1)
ggsave("~/Dropbox/Elegy/Figures/Supplementary/S-variation-twi-fit.pdf", p, width = 20, height = 20, units = "cm")


# Example figure
v = rbind(loadvar("twi", "beta", "x"),
    loadvar("dty", "beta", "beta"),
    loadvar("neu", "beta", "beta"),
    loadvar("whc", "beta", "beta"));

v$Y = pmin(0.899, pmax(0.101, 1 - v$X));
v$Ybinned = floor(v$Y * 10) / 10 + 0.05;
for (m in unique(v$model)) {
    for (y in v[model == m, unique(Ybinned)]) {
        v[Ybinned == y & model == m, stack := 1:.N]
    }
}

v[model == "twi", model := Mod1]
v[model == "dty", model := Mod2]
v[model == "neu", model := Mod3]
v[model == "whc", model := Mod4]
v$model = factor(v$model, levels = c(Mod1, Mod2, Mod3, Mod4));
    
p1 = ggplot(v) + 
    geom_pointrange(aes(x = cc, y = 100 * mean, ymin = 100 * lower, ymax = 100 * upper), fatten = 0.2) + 
    geom_line(aes(x = cc, y = 100 * Rfrac, group = model, colour = model)) +
    geom_text(data = data.table(model = factor(c(Mod1, Mod2, Mod3, Mod4), levels = c(Mod1, Mod2, Mod3, Mod4))), aes(label = model), 
        x = 1, y = 55, hjust = 0, size = 6 / ggplot2:::.pt) +
    labs(x = NULL, y = "Resistance frequency (%)", colour = NULL) +
    facet_wrap(~model, ncol = 1, scales = "free") +
    scale_colour_manual(values = Colours) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none", strip.background = element_blank(), strip.text = element_blank());
    
p2 = ploterror(list(
               twi = c("beta", "c", "u", "g", "bigG", "shape", "z"),
               dty = c("beta", "c", "u", "a", "delta", "g", "z"),
               neu = c("beta", "c", "u", "k", "g", "z"),
               whc = c("beta", "b", "u", "k", "g", "z"))) + geom_vline(xintercept = 2.5, linetype = "22", size = 0.3, colour = "#666666")

p3 = plotvacc_example("Acquisition-blocking", T)
p4 = plotvacc_example("Clearance-accelerating", F)

plot_grid(p1, p2, p3, p4, ncol = 4, align = "h", rel_widths = c(2, 1, 1, 1), 
    labels = c("a", "b", "c", ""), label_size = 10, label_x = c(0, 0, 0, 0)) # 6.79 by 3.94 inches

ggsave("~/Dropbox/Elegy/Submission/STM Response/Figures/6-variation.pdf", width = 18.4, height = 10, units = "cm", device = cairo_pdf())

library(eurostat)

q = fread("~/Dropbox/ResistanceDrivers/Data/tps00001.txt", header = T, check.names = T)[, c(1,11)]
names(q) = c("RegionName", "pop")
q = rbind(q, data.table(RegionName = "Norway", pop = 5258000))
q = rbind(q, data.table(RegionName = "Iceland", pop = 338349))


res = fread("~/Dropbox/ResistanceDrivers/Data/ecdc-resistance.csv", na.strings = "-")[Population == "Streptococcus pneumoniae|Penicillins" & Indicator == "Non-susceptible (I and R) isolates proportion  "]
complete = res[, sum(is.na(NumValue)) == 0, by = RegionCode][V1 == T, RegionCode]
res = res[RegionCode %in% complete & RegionName != "Belgium"]
res = merge(res, q, by = "RegionName")

rsmoothed = res
rsmoothed = rbind(rsmoothed[, .(RegionCode, NumValue, Time, pop)], rsmoothed[, .(RegionCode, NumValue, Time = Time + 1, pop)], rsmoothed[, .(RegionCode, NumValue, Time = Time - 1, pop)])

europe = rsmoothed[, weighted.mean(NumValue, weights = pop, na.rm = T), by = Time]
europe

plpl2 = ggplot(res[RegionCode %in% complete]) + geom_line(aes(Time, NumValue, colour = RegionCode, group = RegionCode)) + geom_line(data = europe, aes(Time, V1), size = 2, linetype = "32") + scale_x_continuous(breaks = seq(2005, 2017, by = 2), limits = c(2005, 2017)) + labs(x = "Year", y = "Frequency of\npenicillin non-susceptibility (%)", colour = "Country")


rN = fread("~/Dropbox/ResistanceDrivers/Data/ecdc-resistance.csv", na.strings = "-")[Population == "Streptococcus pneumoniae|Penicillins" & Indicator == "Total tested isolates", .(Time, RegionName, N = NumValue)]
rR = fread("~/Dropbox/ResistanceDrivers/Data/ecdc-resistance.csv", na.strings = "-")[Population == "Streptococcus pneumoniae|Penicillins" & Indicator == "Non-susceptible (I and R) isolates", .(Time, RegionName, R = NumValue)]
res2 = merge(rN, rR)
res2 = cbind(res2, res2[, binom.confint(R, N, methods = "bayes", prior.shape1 = 1, prior.shape2 = 1)])

cons = fread("~/Dropbox/ResistanceDrivers/Data/ecdc-drugs.txt")[sector == "AC" & atc == "J01C"]
rc = merge(res2, cons, by.x = c("RegionName", "Time"), by.y = c("country", "year"))

rc[RegionName == "Belgium" & Time > 2007, N := 0]
rc = rc[N > 0]

# quantify linear trend
summary(lm(mean ~ ddd.per.thousand, data = rc))

lm(mean ~ ddd.per.thousand, data = rc)$coefficients[[1]]
plpl1 = ggplot(rc, aes(x = ddd.per.thousand, y = 100 * mean)) + geom_pointrange(aes(ymin = 100 * lower, ymax = 100 * upper, colour = as.factor(Time)), fatten = 0.3, size = 0.3) + geom_smooth(method = "lm") + labs(x = "Defined daily doses per 1000 people per day, J01C, primary care", y = "Frequency of\npenicillin non-susceptibility (%)", colour = "Year") + scale_colour_manual(values = hsv(h = seq(0, 1, length.out = 13), s = 1, v = seq(0.5, 1, length.out = 13)), guide = guide_legend(ncol = 2))

# trend in trends?
rc0 = rc[Time <= 2009]
rc0 = merge(rc0, rc0[, .(Completeness = .N), by = RegionName], by = "RegionName")
rc0 = rc0[Completeness == 5]

trends = NULL
for (y in 2005:2009) {
    model = lm(mean ~ ddd.per.thousand, data = rc0[Time == y]);
    trends = rbind(trends, data.table(intercept = model$coefficients[[1]], 
        slope = model$coefficients[[2]], mu = rc0[Time == y, mean(mean)], year = y))
}

summary(lm(mu ~ year, data = trends))
summary(lm(slope ~ year, data = trends))

plpl3 = ggplot(rc0, aes(x = ddd.per.thousand, y = 100 * mean)) + 
    geom_pointrange(aes(ymin = 100 * lower, ymax = 100 * upper, colour = as.factor(Time)), fatten = 0.3, size = 0.3) + 
    geom_smooth(method = "lm") + labs(x = "Defined daily doses per 1000 people per day, J01C, primary care", y = "Frequency of\npenicillin non-susceptibility (%)", colour = "Year") + 
    scale_colour_manual(values = hsv(h = seq(0, 1, length.out = 13), s = 1, v = seq(0.5, 1, length.out = 13)), guide = guide_legend(ncol = 2)) +
    theme(legend.position = "none") +
    facet_wrap(~Time, nrow = 1)

plot_grid(plpl1, plpl2, plpl3, nrow = 3, align = "v", labels = c("a", "b", "c"), label_size = 10, rel_heights = c(2,2,1.5))
ggsave("~/Dropbox/Elegy/Submission/STM Response/Figures/salient.patterns.pdf", width = 17.25, height = 16, units = "cm")

ggplot(rc[RegionName != "Belgium" & N > 10], aes(x = ddd.per.thousand, y = mean)) + geom_pointrange(aes(ymin = lower, ymax = upper)) + geom_smooth(method = "lm") + scale_y_log10() + scale_x_log10()

rc2 = rc[, .(ddd = mean(ddd.per.thousand, na.rm = T), N = sum(N), R = sum(R)), by = .(RegionName)]
rc2 = cbind(rc2, rc2[, binom.confint(R, N, methods = "bayes", prior.shape1 = 1, prior.shape2 = 1)])

ggplot(rc2, aes(x = ddd, y = mean)) + geom_pointrange(aes(ymin = lower, ymax = upper)) + geom_smooth(method = "lm")
ggplot(rc2, aes(x = ddd, y = mean)) + geom_pointrange(aes(ymin = lower, ymax = upper)) + geom_smooth(method = "lm") + scale_y_log10() + scale_x_log10()
