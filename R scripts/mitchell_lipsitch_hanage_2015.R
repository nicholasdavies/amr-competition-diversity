# Mitchell, Lipsitch, Hanage model

library(deSolve)
library(ggplot2)
library(ggthemes)

mlh = function(time, state, parameters) { with(as.list(c(state, parameters)),
{
    U = (1 - f) - IVs - INs - IRs - IVNs - IVRs - INRs
    V = f       - IVv - INv - IRv - IVNv - IVRv - INRv
    
    LV = beta * (IVs + IVv + 0.5 * (IVNs + IVNv + IVRs + IVRv))
    LN = beta * (INs + INv + 0.5 * (IVNs + IVNv + INRs + INRv))
    LR = a*beta*(IRs + IRv + 0.5 * (IVRs + IVRv + INRs + INRv))
    
    dIVs  = LV*U + k*LV/2*IVRs + k*LV/2*IVNs - k*LN*IVs - k*LR*IVs - (d+u)*IVs
    dINs  = LN*U + k*LN/2*IVNs + k*LN/2*INRs - k*LV*INs - k*LR*INs - (d+u)*INs
    dIRs  = LR*U + k*LR/2*IVRs + k*LR/2*INRs - k*LN*IRs - k*LV*IRs - (d+u)*IRs
    
    dIVNs = k*LN*IVs + k*LV*INs + k*LN/2*IVRs + k*LV/2*INRs - (k*LV/2 + k*LN/2 + k*LR)*IVNs - (d+u)*IVNs
    dIVRs = k*LR*IVs + k*LV*IRs + k*LR/2*IVNs + k*LV/2*INRs - (k*LV/2 + k*LR/2 + k*LN)*IVRs - (d+u)*IVRs
    dINRs = k*LR*INs + k*LN*IRs + k*LR/2*IVNs + k*LN/2*IVRs - (k*LN/2 + k*LR/2 + k*LV)*INRs - (d+u)*INRs
    
    LV = v * LV

    dIVv  = LV*V + k*LV/2*IVRv + k*LV/2*IVNv - k*LN*IVv - k*LR*IVv - (d+u)*IVv
    dINv  = LN*V + k*LN/2*IVNv + k*LN/2*INRv - k*LV*INv - k*LR*INv - (d+u)*INv
    dIRv  = LR*V + k*LR/2*IVRv + k*LR/2*INRv - k*LN*IRv - k*LV*IRv - (d+u)*IRv
    
    dIVNv = k*LN*IVv + k*LV*INv + k*LN/2*IVRv + k*LV/2*INRv - (k*LV/2 + k*LN/2 + k*LR)*IVNv - (d+u)*IVNv
    dIVRv = k*LR*IVv + k*LV*IRv + k*LR/2*IVNv + k*LV/2*INRv - (k*LV/2 + k*LR/2 + k*LN)*IVRv - (d+u)*IVRv
    dINRv = k*LR*INv + k*LN*IRv + k*LR/2*IVNv + k*LN/2*IVRv - (k*LN/2 + k*LR/2 + k*LV)*INRv - (d+u)*INRv
    
    return (list(c(dIVs, dINs, dIRs, dIVNs, dIVRs, dINRs, dIVv, dINv, dIRv, dIVNv, dIVRv, dINRv)))
}) }

f = 0.5
parameters = c(beta = 1.5, d = 1, u = 1/60, a = 1.01, k = 0.5, v = 0.95, f = f)
init = c(IVs = (1-f)/3, INs = (1-f)/3, IRs = (1-f)/3, IVNs = 0, IVRs = 0, INRs = 0,
    IVv = f/3, INv = f/3, IRv = f/3, IVNv = 0, IVRv = 0, INRv = 0)

sol = as.data.table(ode(y = init, times = seq(0, 2000, by = 10), func = mlh, parms = parameters, method = "lsode", maxsteps = 50000));
sol = melt(sol, id.vars = "time", variable.factor = F);
sol[, vstatus := ifelse(variable %like% "s$", "unvacc", "vacc")];
sol[, compartment := substr(variable, 1, nchar(variable) - 1)];

ggplot(sol) + geom_line(aes(x = time, y = value, colour = compartment)) + facet_grid(compartment ~ vstatus, scales = "free") + geom_hline(aes(yintercept = 0))
