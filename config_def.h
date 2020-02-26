// config_def.h
// Contains definitions of all config parameters to be loaded from the
// configuration file or the command line.

PARAMETER ( int,            model,              1 );        // 1 = single carriage; 2 = mixed carriage; 3 = D-type
PARAMETER ( int,            groups,             1 );        // number of subgroups
PARAMETER ( int,            group_span,         1 );        // Average every X groups for evaluating with mcmc
PARAMETER ( int,            group_connectivity, 1 );        // Every adjacent X groups are connected
PARAMETER ( bool,           separate_groups,    true );     // Where groups can be separated and solved individually, do so.
PARAMETER ( bool,           skip_contact,       true );     // when group_connectivity is 1, skip multiplying by contact matrix
PARAMETER ( bool,           vaccinate,          false );    // whether to track separate vaccinated groups
PARAMETER ( int,            whc_steps,          2 );        // for within-host competition model, number of dual-carriage steps
PARAMETER ( bool,           whc_neutral,        true );     // for within-host competition model, whether Sr (Rs) compartments transmits only S (R)
PARAMETER ( double,         whc_iota,           0.001 );    // for within-host competition model, (potentially untransmitted) amount of minor carriage in Sr / Rs individuals
PARAMETER ( int,            D_n,                1 );        // for D-type model, number of D-types
PARAMETER ( bool,           D_report_all,       false );    // for D-type model, whether to report the frequency of all D-types or just the aggregate S and R frequencies
PARAMETER ( double,         init_freq_s,        0.001 );    // initial proportion of each group infected with susceptible strain
PARAMETER ( double,         init_freq_r,        0.001 );    // initial proportion of each group infected with resistant strain
PARAMETER ( vector<double>, consumption,        {} );       // vector of consumption rates, potentially to be transformed into treatment rates
PARAMETER ( double,         treatment_course,   -1  );      // for consumption specified in DDDs per thousand people per day and a base time unit of one month, number of DDDs comprising one treatment course; if negative, leave tau alone
PARAMETER ( int,            treatment_regions,  -1 );       // if positive, distribute treatment to M regions modelled as a gamma distribution
PARAMETER ( string,         treatment_shift,    "none" );   // none (no shifting), all (all taus shifted by the same amount), or each (taus for each country shifted by separate amount)
PARAMETER ( vector<double>, population,         {} );       // vector of population sizes for each group, with either no entries, [groups] entries, or [groups / treatment_regions] entries

PARAMETER ( string,         beta,               "2" );      // transmission rate
PARAMETER ( string,         c,                  "0.05" );   // cost of resistance
PARAMETER ( string,         u,                  "1" );      // rate of natural clearance
PARAMETER ( string,         tau,                "0.05" );   // rate of treatment
PARAMETER ( string,         k,                  "1" );      // mixed-carriage model: relative efficiency of co-colonisation
PARAMETER ( string,         kk,                 "0" );      // for single carriage model, relative rate of "superinfection"/"knockout" by sensitive strain of R hosts
PARAMETER ( string,         b,                  "0" );      // mixed-carriage model: rate of ecological succession by sensitive strain
PARAMETER ( string,         b0,                 "0" );      // mixed-carriage model: rate of S_R -> S transition as a multiplier on b
PARAMETER ( string,         a,                  "0" );      // D-type model: strength of balancing selection
PARAMETER ( string,         delta,              "0" );      // D-type model: radius of clearance rate, relative to u and centred on u
PARAMETER ( string,         g,                  "1" );      // assortativity of contact
PARAMETER ( string,         G,                  "-1" );     // assortativity of contact at larger scale, if nonnegative and group_scale > 1; overrides custom waifw matrix
PARAMETER ( string,         waifw,              "-1" );     // custom contact matrix, if first element is nonnegative
PARAMETER ( string,         h,                  "1" );      // proportion of potential hosts in each group
PARAMETER ( string,         psi,                "0" );      // importation of strains; per D-type under model 3
PARAMETER ( string,         psi_r,              "0" );      // fraction of imported strains that are resistant; per D-type under model 3
PARAMETER ( string,         z,                  "" );       // scaling factor for treatment_shift
PARAMETER ( string,         shape,              "1" );      // shape of gamma distribution for use with treatment_regions
PARAMETER ( string,         l_noise,            "0.1" );    // noise parameter for likelihood model
PARAMETER ( string,         v,                  "0" );      // vaccinated fraction
PARAMETER ( string,         dv,                 "1" );      // factor by which clearance is multiplied in vaccinees (typically, enhanced clearance)
PARAMETER ( string,         sv,                 "1" );      // factor by which susceptibility is multiplied in vaccinees (typically, reduced susceptibility)
PARAMETER ( string,         tv,                 "1" );      // factor by which the treatment rate is multiplied in vaccinees (typically, reduced treatment)
PARAMETER ( string,         x,                  "" );       // wildcard parameters
PARAMETER ( Func<string>,   xfunc,              "" );       // wildcard assignment

PARAMETER ( bool,           play_changes,       false );    // If true, read "changes" from a table in a file. Other params set as usual. Possibility for interventions
PARAMETER ( int,            changes_limit,      -1 );       // If positive, only read this many lines from changes file (approximate)
PARAMETER ( string,         changes_file,       "null" );   // File to read changes from
PARAMETER ( vector<string>, interv,             {} );       // For play_changes, intervention to apply: none, vaccine, carriage, treatment, hightreatment
PARAMETER ( vector<double>, interv_coverage,    {} );       // Coverage of intervention; 0 is no coverage, 1 is full coverage
PARAMETER ( vector<double>, interv_strength,    {} );       // Strength of intervention; 0 is no intervention, 1 is full intervention
PARAMETER ( vector<double>, interv_time,        {} );       // Evaluate intervention after this many time units; -1 is until equilibrium
PARAMETER ( bool,           interv_import,      false );    // If true, adjust importation for vaccine and carriage interventions
PARAMETER ( double,         interv_psi,         0 );        // Slope of psi adjustment if interv_import is true (e.g. 2 means psi is 0 when interv_strength is 0.5)

PARAMETER ( bool,           mcmc,           false );        // Run MCMC
PARAMETER ( int,            mcmc_algorithm, 1 );            // Which algorithm to use: DE-MCMC (1), DREAM (2), optimize (3), optimize with small initial step size (4), single-parameter parameter sweep (-N_STEPS)
PARAMETER ( int,            mcmc_burn_in,   1000 );         // Number of burn-in iterations
PARAMETER ( int,            mcmc_thinning,  1 );            // Thinning factor on chains
PARAMETER ( int,            mcmc_R_skip,    -1 );           // If positive, after this many iterations of burn-in, skip to main iterations if all Gelman-Rubin-Brooks R statistics < 1.05.
PARAMETER ( string,         mcmc_resume,    "" );           // Resume MCMC from this output file; all other parameters should be the same except (optionally) trials
PARAMETER ( vector<int>,    mcmc_data_r,    {0} );          // Samples in each country testing positive for resistance
PARAMETER ( vector<int>,    mcmc_data_n,    {1} );          // Number of samples in each country
PARAMETER ( string,         l_model,        "normal_trunc" ); // Likelihood model ("beta" = beta model, "error" = error model, "normal_clamp" = clamped normal,
                                                                // "normal_trunc" = truncated normal, "binomial" = just binomial, "mean" = normal centered on mean)
PARAMETER ( string,         l_y,            "U 0 1" );      // likelihood penalty applied for average carriage prevalence
PARAMETER ( string,         l_d,            "U 0 1" );      // likelihood penalty applied for average dual-carrier prevalence
PARAMETER ( string,         l_rf,           "U 0 1" );      // likelihood penalty applied for average resistance frequency
PARAMETER ( string,         l_yi,           "U 0 1" );      // likelihood penalty applied for per-span-group carriage prevalence
PARAMETER ( bool,           l_weight,       false );        // If true, weight likelihood components by the population of the respective country

PARAMETER ( bool,           find,           false );        // if true, find parameters that satisfy a given equilibrium outcome
PARAMETER ( string,         find_param,     "beta_s" );     // the parameter to focus on
PARAMETER ( double,         find_x0,        0.0 );          // param min
PARAMETER ( double,         find_x1,        1.0 );          // param max
PARAMETER ( int,            find_n,         10 );           // number of tests within range (x0, x1)
PARAMETER ( string,         find_goal,      "Y" );          // either Y (total infected), Yi (infected in ith group), Rfrac (resistance prevalence), or Rfraci (res prev in ith group)
PARAMETER ( int,            find_i,         0 );            // for Yi or Rfraci above, which group to evaluate
PARAMETER ( double,         find_value,     0.5 );          // what the find_goal should equal

PARAMETER ( bool,           autofind,       false );        // whether to run autofind 
PARAMETER ( double,         autofind_y,     0.5 );          // carriage to aim for
PARAMETER ( int,            autofind_i,     -1 );           // if non-negative, this group is where we aim to pass through Rfrac = 50%
PARAMETER ( double,         autofind_tau,   0.125 );        // if autofind_i is non-negative, this treatment rate is where we aim to pass through Rfrac = 50%
PARAMETER ( string,         autofind_cost,  "beta_ratio" ); // which cost to adjust (beta_ratio or xi)

PARAMETER ( int,            trials,         1 );            // Number of trials
PARAMETER ( bool,           root,           false );        // Solve for root of ODEs rather than integrate ODEs
PARAMETER ( bool,           start_previous, false );        // When possible, start solution at previous country equilibrium
PARAMETER ( double,         epsilon,        1e-6 );         // Error tolerance
PARAMETER ( int,            max_epochs,     1000 );         // Maximum epochs (ODE integration) or iterations (ODE root finding)
PARAMETER ( int,            min_epochs,     0 );            // Minimum epochs
PARAMETER ( double,         t_epoch,        100 );          // Time in one epoch
PARAMETER ( double,         t_step,         0.1 );          // Initial time step
PARAMETER ( bool,           euler,          false );        // Use Euler integration -- only for play_changes when interv_time is positive

PARAMETER ( string,         fileout,        "./out.txt" );  // Output file
PARAMETER ( int,            threads,        1 );            // Multithreading
PARAMETER ( bool,           verbose,        true );         // Be wordy

