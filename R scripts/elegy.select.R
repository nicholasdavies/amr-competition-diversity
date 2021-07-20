# set appropriate variables for different data sets
select = function(set)
{
    if (set == "spp" | set == "spp_tau")
    {
        filename <<- "spp";
        if (set == "spp") {
            postfix <<- "";
        } else {
            postfix <<- "_tau";
        }
        species <<- "S. pneumoniae";
        countrycodes <<- c("NL", "DE", "EE", "FI", "LV", "AT", "NO", "CZ", "HU", "SE", "GB", "MT", "BG", "IS", "SI", "DK", "PL", "HR", "PT", "IE", "ES", "LU", "BE", "LT", "IT", "FR", "CY");
        courses <<- 5;
        tau <<- 365.25 / (12000 * courses) * c(4.3,  4.5,  4.9,  6.0,  6.1,  6.2,  6.6,  6.9,  7.0,  7.2,  7.7,  8.7,  9.6,  9.8, 10.0, 10.2, 10.5, 11.3, 11.3, 11.5, 12.1, 12.6, 13.6, 13.9, 14.6, 15.0, 15.8);
        r <<- c( 19,    2,    0,   69,    0,   16,   10,    8,   33,   31,   70,    0,    5,    3,   33,   34,    6,   25,   33,   76,  192,    2,  143,    3,   43,  226,    5);
        n <<- c(939,   72,   64,  522,   31,  313,  616,  205,  146, 1028, 1744,   13,   32,   42,  195, 1030,   21,  137,  202,  435,  860,   35, 1511,   67,  291,  663,   15);
        populations <<- c(16381696, 82266372, 1340680, 5288720, 2200325, 8295487, 4709153, 10298828, 10055780, 9148092, 61322463, 406724, 7545338, 311566, 2018122, 5461438, 38120560, 4312749, 10542964, 4398942, 45226803, 479993, 10625700, 3231294, 58438310, 63826129, 767125);
        u <<- 0.65
        disconv <<- 122
        discnt <<- "100,000"
        disname <<- "severe pneumococcal pneumonia"
        disgroup <<- "children under 5"
        beta.shape <<- 5
        beta.scale <<- 0.35
        carriage.mean <<- 0.5
        carriage.sd <<- 0.002
    }
    else if (set == "ecp" | set == "ecp_tau")
    {
        filename <<- "ecp";
        if (set == "ecp") {
            postfix <<- "";
        } else {
            postfix <<- "_tau";
        }
        species <<- "E. coli";
        countrycodes <<- c("NL", "DE", "EE", "FI", "SE", "NO", "LV", "AT", "HU", "BG", "CZ", "MT", "SK", "GB", "LT", "SI", "PL", "DK", "IS", "HR", "PT", "CY", "LU", "ES", "GR", "IE", "IT", "BE", "FR", "RO");
        courses <<- 10;
        tau <<- 365.25 / (12000 * courses) * c(4.4,  4.5,  4.8,  6.3,  6.3,  6.4,  6.5,  6.6,  7.1,  8.2,  8.3,  8.8,  8.9,  8.9,  9.6,  9.6, 10.4, 10.7, 10.7, 11.9, 12.2, 13.3, 13.4, 14.5, 14.6, 15.5, 15.5, 16.5, 18.8, 18.8);
        r <<- c(2540, 3974,   93,  889,  135, 1511,  103, 2433, 1193,   95, 1721,  132,  551, 3366,  347,  727,  224, 2080,   77,  576, 2990,   84,  209, 4104,  605, 1752, 2280, 1552, 6239,  189);
        n <<- c(5376, 8053,  196, 2472,  396, 3299,  192, 4880, 1970,  143, 3172,  238,  878, 5117,  582, 1326,  346, 4594,  173, 1042, 5177,  123,  347, 6427, 1079, 2646, 3385, 2674,10946,  259);
        u <<- 0.4;
        populations <<- c(16939923, 81686611, 1315407, 5479531, 9799186, 5190239, 1977527, 8642699, 9843028, 7177991, 10546059, 445053, 5423801, 65128861, 2904910, 2063531, 37986412, 5683483, 330815, 4207993, 10358076, 847664, 569604, 46444832, 10820883, 4701957, 60730582, 11274196, 66593366, 19815481);
        disconv <<- 1000
        discnt <<- "1,000"
        disname <<- "ExPEC carriage"
        disgroup <<- "people"
        beta.shape <<- 5
        beta.scale <<- 0.32
        carriage.mean <<- 0.75
        carriage.sd <<- 0.002
    }
}

# zero out parameters - needed to declare them
filename = "null"
postfix = "null"
species = "null"
countrycodes = c("null")
r = c(0)
n = c(0)
u = 0
populations = c(0)
tau = c(0)
courses = 0
disconv = 0
disexp = 0
disname = "null"
disgroup = "null"
beta.shape = 0
beta.scale = 0
carriage.mean = 0
carriage.sd = 0
