library(data.table)
library(wpp2019)
library(countrycode)

data(pop)
popul = data.table(name = pop$name, pop = pop$`2005`)
popul[, cc := countrycode(name, origin = "country.name", destination = "iso3c")]

travel = fread("~/Dropbox/Elegy/R scripts/tour_dem_tnw_1_Data.csv")
travel[, cc := countrycode(GEO, origin = "country.name", destination = "iso3c")]
travel[, Value := as.numeric(gsub(" |:", "", Value))]

d = merge(popul,
    travel[TIME == 2007 & 
        PARTNER == "Outbound" &
        PURPOSE == "Total" &
        DURATION == "1 night or over",],
    by = "cc")

d[, proportion.abroad := as.numeric(Value) / (365 * pop * 1000)]
d[, weighted.mean(proportion.abroad, pop, na.rm = T)]
