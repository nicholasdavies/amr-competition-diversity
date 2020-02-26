#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <deque>
#include <string>
#include <stdexcept>
#include <nlopt.hpp>
#include <boost/math/distributions/beta.hpp>
#include <omp.h>
#include "Randomizer/randomizer.h"
#include "Randomizer/distribution.h"

#ifndef MCMC_H
#define MCMC_H

template <typename Observables, typename LikelihoodFunc, typename ReportFunc>
void RandomSample(Randomizer& R, LikelihoodFunc likelihood, ReportFunc report,
    int iterations, std::vector<Distribution>& priors, bool verbose = true)
{
    using namespace std;

// #ifdef _OPENMP
//     if (n_threads > 0)
//         omp_set_num_threads(n_threads);
// #endif

    // Parameter set
    int n_theta = priors.size();
    vector<double> theta(n_theta, 0.0);

    // Assemble target func
    auto target = [&](vector<double>& theta, Observables& obs, double& l)
    {
        double p = 0;
        for (int d = 0; d < n_theta; ++d)
        {
            double pd = priors[d].LogProbability(theta[d]);
            if (pd == -std::numeric_limits<double>::infinity())
            {
                l = -std::numeric_limits<double>::infinity();
                return -std::numeric_limits<double>::infinity();
            }
            p += pd;
        }

        l = likelihood(theta, obs);

        return l + p;
    };

    // Do iterations
    if (verbose)
        cout << "\nSampling...\n";
    for (int i = 0; i < iterations; ++i)
    {
        // Sample from prior
        for (int d = 0; d < n_theta; ++d)
            theta[d] = priors[d].Random(R);

        // Evaluate target function
        Observables obs;
        double l = 0, p = target(theta, obs, l);

        // Report results of this iteration
        report(i, p, 0, l, theta, obs);

        // Print progress
        if (verbose && i % 100 == 0)
            cout << "." << flush;
    }

    if (verbose)
        cout << "\n";
}

template <typename Observables, typename LikelihoodFunc, typename ReportFunc>
void MCMC_Priors(Randomizer& R, LikelihoodFunc likelihood, ReportFunc report,
    int burn_in, int iterations, int n_chains, std::vector<Distribution>& priors,
    bool verbose = true, std::vector<std::string> param_names = std::vector<std::string>(),
    bool in_parallel = false, int n_threads = -1, 
    std::vector<std::vector<double>> init = std::vector<std::vector<double>>(), int init_iter = 1)
{
    using namespace std;

#ifdef _OPENMP
    if (in_parallel && n_threads > 0)
        omp_set_num_threads(n_threads);
#endif

    if (n_chains < 3)
        throw runtime_error("Cannot use DE-MCMC with fewer than 3 chains.");

    // Store acceptance rates
    unsigned int ar_size = 1000;
    unsigned int ar_i = 0;
    bool show_ar = false;
    vector<bool> acceptances(ar_size, false);

    // Storage for chains and settings
    int n_theta = priors.size();
    vector<vector<double>> chains(n_chains, vector<double>(n_theta, 0.0));
    vector<Observables> obs(n_chains, Observables());
    vector<double> p(n_chains, 0.0);    // log probability for each chain
    vector<double> l(n_chains, 0.0);    // log likelihood for each chain

    // Storage for calls to Randomizer - to make thread-safe
    vector<vector<double>> random_perturbations(n_chains, vector<double>(n_theta, 0.0));
    vector<double> random_tests(n_chains, 0.0);

    // Assemble target func
    auto target = [&](vector<double>& theta, Observables& obs, double& l)
    {
        double p = 0;
        for (int d = 0; d < n_theta; ++d)
        {
            double pd = priors[d].LogProbability(theta[d]);
            if (pd == -std::numeric_limits<double>::infinity())
            {
                l = -std::numeric_limits<double>::infinity();
                return -std::numeric_limits<double>::infinity();
            }
            p += pd;
        }

        l = likelihood(theta, obs);

        return l + p;
    };

    // Initialize chains . . .
    if (init.empty())
    {
        // . . . from prior
        for (int c = 0; c < n_chains; ++c)
            for (int d = 0; d < n_theta; ++d)
                chains[c][d] = priors[d].RandomInit(R);
    }
    else
    {
        // . . . from initial values supplied
        if ((int)init.size() != n_chains || (int)init[0].size() != n_theta)
            throw runtime_error("init vector supplied is not the right size.");
        chains = init;
    }

    // Set initial probabilities and observables
    if (verbose)
        cout << "Initializing chains...\n";
    #pragma omp parallel for if(in_parallel) schedule(dynamic)
    for (int c = 0; c < n_chains; ++c)
    {
        p[c] = target(chains[c], obs[c], l[c]);
        if (verbose)
            cout << "." << flush;
    }

    // Do iterations
    if (verbose)
        cout << "\nIterating...\n";
    for (int i = init_iter; i < burn_in + iterations; ++i)
    {
        // Prepare storage and random variates
        for (int c = 0; c < n_chains; ++c)
        {
            for (int d = 0; d < n_theta; ++d)
                random_perturbations[c][d] = R.Normal(0.0, 0.1);
            random_tests[c] = R.Uniform();
        }
        auto saved_chains = chains;
        vector<int> accept(n_chains, 0);

        #pragma omp parallel for if(in_parallel) schedule(dynamic)
        for (int c = 0; c < n_chains; ++c)
        {
            vector<double> theta_p = chains[c];
            int c_from = c;

            // Generate proposal by random jump
            for (int d = 0; d < n_theta; ++d)
                theta_p[d] += random_perturbations[c][d];

            // Calculate log-probability and accept or reject
            Observables obs_p;
            double l_p = 0;
            double p_p = target(theta_p, obs_p, l_p);
            if ( (p_p == -std::numeric_limits<double>::infinity() && p[c_from] == -std::numeric_limits<double>::infinity() && random_tests[c] < 0.5)
                || (p_p > -std::numeric_limits<double>::infinity() && random_tests[c] < exp(p_p - p[c_from])) )
            {
                accept[c_from] = 1;
                chains[c_from] = theta_p;
                obs[c_from] = obs_p;
                p[c_from] = p_p;
                l[c_from] = l_p;
            }
        }

        // Update acceptances
        for (int c = 0; c < n_chains; ++c)
        {
            if (ar_i == ar_size - 1) show_ar = true;
            acceptances[ar_i] = accept[c];
            ar_i = (ar_i + 1) % ar_size;
        }

        // Report results of this iteration
        for (int c = 0; c < n_chains; ++c)
            report(i - burn_in, p[c], c, l[c], chains[c], obs[c]);

        // Print progress
        if (verbose)
        {
            cout << "." << flush;
            if (i % 100 == 0)
            {
                cout << "\n" << (i < burn_in ? "burn-in" : "main") << " iteration " << i - burn_in << ":";
                if (!param_names.empty())
                {
                    cout << "\n         " << setw(12) << right << "log (P)" << setw(12) << right << "log (L)";
                    for (auto n : param_names)
                        cout << setw(12) << right << n;
                }
                for (int c = 0; c < n_chains; ++c)
                {
                    cout << "\nchain" << setw(4) << right << c << setw(12) << right << p[c] << setw(12) << right << l[c];
                    for (int d = 0; d < n_theta; ++d)
                        cout << setw(12) << right << chains[c][d];
                }

                double acceptance_rate = show_ar ? (double)count(acceptances.begin(), acceptances.end(), true) / ar_size : -1;
                cout << "\nacceptance rate: " << acceptance_rate << "\n\n" << flush;
            }
        }
    }
    if (verbose)
        cout << "\n";
}

template <typename Observables, typename LikelihoodFunc, typename ReportFunc>
void DEMCMC_Priors(Randomizer& R, LikelihoodFunc likelihood, ReportFunc report,
    int burn_in, int iterations, int n_chains, std::vector<Distribution>& priors,
    bool verbose = true, std::vector<std::string> param_names = std::vector<std::string>(),
    bool in_parallel = false, int n_threads = -1, bool classic_gamma = false, int R_skip = -1,
    std::vector<std::vector<double>> init = std::vector<std::vector<double>>(), int init_iter = 1)
{
    using namespace std;

#ifdef _OPENMP
    if (in_parallel && n_threads > 0)
        omp_set_num_threads(n_threads);
#endif

    if (n_chains < 3)
        throw runtime_error("Cannot use DE-MCMC with fewer than 3 chains.");

    // Store acceptance rates
    unsigned int ar_size = 1000;
    unsigned int ar_i = 0;
    bool show_ar = false;
    vector<bool> acceptances(ar_size, false);

    // Storage for chains and settings
    int n_theta = priors.size();
    vector<vector<double>> chains(n_chains, vector<double>(n_theta, 0.0));
    vector<Observables> obs(n_chains, Observables());
    vector<double> p(n_chains, 0.0);    // log probability for each chain
    vector<double> l(n_chains, 0.0);    // log likelihood for each chain
    double b = 0.001;

    // Storage for calls to Randomizer - to make thread-safe
    vector<vector<double>> random_perturbations(n_chains, vector<double>(n_theta, 0.0));
    vector<double> random_tests(n_chains, 0.0);
    vector<double> random_gammas(n_chains, 0.0);
    vector<vector<int>> random_chains(n_chains, vector<int>(2, 0));

    // Storage for calculation of Gelman-Rubin-Brooks R diagnostic
    int running_count = 0;
    vector<vector<double>> running_mean(n_chains, vector<double>(n_theta, 0.0));
    vector<vector<double>> running_M2(n_chains, vector<double>(n_theta, 0.0));
    vector<vector<double>> running_s2(n_chains, vector<double>(n_theta, 0.0));
    vector<double> R_hat(n_theta, 0.0);
    const int n_between_R = 100;

    // Assemble target func
    auto target = [&](vector<double>& theta, Observables& obs, double& l)
    {
        double p = 0;
        for (int d = 0; d < n_theta; ++d)
        {
            double pd = priors[d].LogProbability(theta[d]);
            if (pd == -std::numeric_limits<double>::infinity())
            {
                l = -std::numeric_limits<double>::infinity();
                return -std::numeric_limits<double>::infinity();
            }
            p += pd;
        }

        l = likelihood(theta, obs);

        return l + p;
    };

    // Initialize chains . . .
    if (init.empty())
    {
        // . . . from prior
        for (int c = 0; c < n_chains; ++c)
            for (int d = 0; d < n_theta; ++d)
                chains[c][d] = priors[d].RandomInit(R);
    }
    else
    {
        // . . . from initial values supplied
        if ((int)init.size() != n_chains || (int)init[0].size() != n_theta)
            throw runtime_error("init vector supplied is not the right size.");
        chains = init;
    }

    // Set initial probabilities and observables
    if (verbose)
        cout << "Initializing chains...\n";
    #pragma omp parallel for if(in_parallel) schedule(dynamic)
    for (int c = 0; c < n_chains; ++c)
    {
        p[c] = target(chains[c], obs[c], l[c]);
        if (verbose)
            cout << "." << flush;
    }

    // Do iterations
    if (verbose)
        cout << "\nIterating...\n";
    for (int i = init_iter; i < burn_in + iterations; ++i)
    {
        // Prepare storage and random variates
        bool migration = i < burn_in * 0.5 ? R.Bernoulli(0.05) : false;
        vector<int> migration_indices(n_chains, 0);
        if (migration)
        {
            for (int c = 0; c < n_chains; ++c)
                migration_indices[c] = c;
            R.Shuffle(migration_indices.begin(), migration_indices.end());
        }
        for (int c = 0; c < n_chains; ++c)
        {
            for (int d = 0; d < n_theta; ++d)
                random_perturbations[c][d] = R.Uniform(-b, b);
            random_tests[c] = R.Uniform();
            if (!migration)
            {
                if (classic_gamma)
                    random_gammas[c] = (i % 10 == 0 ? 1.0 : 2.38 / sqrt(2 * n_theta));
                else
                    random_gammas[c] = R.Uniform(0.5, 1.0);
                do random_chains[c][0] = R.Discrete(n_chains); while (random_chains[c][0] == c);
                do random_chains[c][1] = R.Discrete(n_chains); while (random_chains[c][1] == c || random_chains[c][1] == random_chains[c][0]);
            }
        }
        auto saved_chains = chains;
        vector<int> accept(n_chains, 0);

        #pragma omp parallel for if(in_parallel) schedule(dynamic)
        for (int c = 0; c < n_chains; ++c)
        {
            vector<double> theta_p = chains[c];
            int c_from = c;

            // Generate proposal, either by migration...
            if (migration)
            {
                c_from = migration_indices[c];
                theta_p = saved_chains[migration_indices[(c + 1) % n_chains]];
                for (int d = 0; d < n_theta; ++d)
                    theta_p[d] += random_perturbations[c][d];
            }
            else // ... or by directed mutation
            {
                for (int d = 0; d < n_theta; ++d)
                    theta_p[d] += random_gammas[c] * (saved_chains[random_chains[c][1]][d] - saved_chains[random_chains[c][0]][d]) + random_perturbations[c][d];
            }

            // Calculate log-probability and accept or reject
            Observables obs_p;
            double l_p = 0;
            double p_p = target(theta_p, obs_p, l_p);
            if ( (p_p == -std::numeric_limits<double>::infinity() && p[c_from] == -std::numeric_limits<double>::infinity() && random_tests[c] < 0.5)
                || (p_p > -std::numeric_limits<double>::infinity() && random_tests[c] < exp(p_p - p[c_from])) )
            {
                accept[c_from] = 1;
                chains[c_from] = theta_p;
                obs[c_from] = obs_p;
                p[c_from] = p_p;
                l[c_from] = l_p;
            }
        }

        // Update acceptances
        for (int c = 0; c < n_chains; ++c)
        {
            if (ar_i == ar_size - 1) show_ar = true;
            acceptances[ar_i] = accept[c];
            ar_i = (ar_i + 1) % ar_size;
        }

        // Update Gelman-Rubin-Brooks R
        bool R_all_ok = true;
        if (R_skip > 0)
        {
            ++running_count;
            for (int c = 0; c < n_chains; ++c)
            {
                for (int d = 0; d < n_theta; ++d)
                {
                    double delta = chains[c][d] - running_mean[c][d];
                    running_mean[c][d] += delta / running_count;
                    double delta2 = chains[c][d] - running_mean[c][d];
                    running_M2[c][d] += delta * delta2;
                }
            }

            // Calculate R every n_between_R generations
            if (i % n_between_R == 0)
            {
                // Finalise running mean and variance
                for (int c = 0; c < n_chains; ++c)
                    for (int d = 0; d < n_theta; ++d)
                        running_s2[c][d] = running_M2[c][d] / (running_count - 1);

                // Calculate statistic for each parameter
                for (int d = 0; d < n_theta; ++d)
                {
                    double M = n_chains;
                    double N = running_count;
                    double W = 0, X = 0, B = 0;
                    for (int c = 0; c < n_chains; ++c)
                    {
                        W += running_s2[c][d];
                        X += running_mean[c][d];
                    }
                    W /= M;
                    X /= M;

                    for (int c = 0; c < n_chains; ++c)
                        B += (running_mean[c][d] - X) * (running_mean[c][d] - X);
                    B *= N / (M - 1);

                    double var = ((N - 1) / N) * W + B / N;
                    R_hat[d] = std::sqrt(var / W);

                    if (R_hat[d] > 1.05)
                        R_all_ok = false;
                }
            }
        }

        // Report results of this iteration
        for (int c = 0; c < n_chains; ++c)
            report(i - burn_in, p[c], c, l[c], chains[c], obs[c]);

        // Print progress
        if (verbose)
        {
            cout << "." << flush;
            if (i % 100 == 0)
            {
                cout << "\n" << (i < burn_in ? "burn-in" : "main") << " iteration " << i - burn_in << ":";
                if (!param_names.empty())
                {
                    cout << "\n         " << setw(12) << right << "log (P)" << setw(12) << right << "log (L)";
                    for (auto n : param_names)
                        cout << setw(12) << right << n;
                }
                for (int c = 0; c < n_chains; ++c)
                {
                    cout << "\nchain" << setw(4) << right << c << setw(12) << right << p[c] << setw(12) << right << l[c];
                    for (int d = 0; d < n_theta; ++d)
                        cout << setw(12) << right << chains[c][d];
                }
                if (R_skip > 0)
                {
                    cout << "\nGelman-Rubin-Brooks diagnostic R ";
                    for (int d = 0; d < n_theta; ++d)
                        cout << setw(12) << right << R_hat[d];
                }

                double acceptance_rate = show_ar ? (double)count(acceptances.begin(), acceptances.end(), true) / ar_size : -1;
                cout << "\nacceptance rate: " << acceptance_rate << "\n\n" << flush;
            }
        }

        // Skip past burn-in if R OK
        if (R_skip > 0 && i > R_skip && i < burn_in && R_all_ok)
        {
            if (verbose)
                cout << "\n\nSkipping to iterations (R < 1.05 reached).\n\n";
            i = burn_in - 1;
        }

    }
    if (verbose)
        cout << "\n";
}

template <typename Observables, typename TargetFunc, typename ReportFunc>
void DEMCMC(Randomizer& R, TargetFunc target, ReportFunc report,
    int burn_in, int iterations, int n_chains,
    std::vector<double> lb, std::vector<double> ub,
    bool verbose = true, std::vector<std::string> param_names = std::vector<std::string>(),
    bool in_parallel = false, int n_threads = -1, bool classic_gamma = false, int R_skip = -1,
    std::vector<std::vector<double>> init = std::vector<std::vector<double>>(), int init_iter = 1)
{
    std::vector<Distribution> priors;
    for (unsigned int b = 0; b < lb.size(); ++b)
        priors.push_back({ Distribution::Uniform, { lb[b], ub[b] } });
    DEMCMC_Priors<Observables, TargetFunc, ReportFunc>
        (R, target, report, burn_in, iterations, n_chains, priors, verbose, param_names, in_parallel, n_threads, classic_gamma, R_skip, init, init_iter);
}

struct DREAM_Options
{
    int burn_in = 1000;
    int iterations = 10000;
    int n_chains = 10;
    int n_threads = -1;
    bool verbose = true;
    std::vector<std::string> param_names = std::vector<std::string>();
    std::vector<std::vector<double>> init_theta = std::vector<std::vector<double>>();
    int init_iter = 1;
};

// TO DO
// just clean up the function below - in particular, jump_thresholds and jump_tests are not needed (just make jump, vector of vector of bools)
// - remove TEMP bit
// - refactor generally
// - rename this something other than DREAM, as it's not exactly the same. Daydream perhaps.
// - add interquartile range stuff...

template <typename Observables, typename LikelihoodFunc, typename ReportFunc>
void DREAM(Randomizer& R, LikelihoodFunc likelihood, ReportFunc report, std::vector<Distribution>& priors, DREAM_Options options = DREAM_Options())
{
    using namespace std;

    // TEMP - to facilitate conversion from DEMCMC to DREAM
    int burn_in = options.burn_in;
    int iterations = options.iterations;
    int n_chains = options.n_chains;
    int n_threads = options.n_threads;
    bool verbose = options.verbose;
    vector<string> param_names = options.param_names;
    vector<vector<double>> init = options.init_theta;
    int init_iter = options.init_iter;


#ifdef _OPENMP
    if (n_threads > 1)
        omp_set_num_threads(n_threads);
#endif

    if (n_chains < 3)
        throw runtime_error("Cannot use DREAM with fewer than 3 chains.");

    // Store acceptance rates
    unsigned int ar_size = 8192;
    unsigned int ar_i = 0;
    bool show_ar = false;
    vector<bool> acceptances(ar_size, false);

    // Storage for chains and settings
    int n_theta = priors.size();
    vector<vector<double>> chains(n_chains, vector<double>(n_theta, 0.0));
    vector<Observables> obs(n_chains, Observables());
    vector<double> p(n_chains, 0.0);    // log probability for each chain
    vector<double> l(n_chains, 0.0);    // log likelihood for each chain
    double b = 0.05;
    double bstar = 1e-6;
    vector<double> p50(n_chains, 0.0);

    // Storage for calls to Randomizer - to make thread-safe
    vector<vector<double>> random_perturbations(n_chains, vector<double>(n_theta, 0.0));
    vector<double> random_extensions(n_chains, 0.0);
    vector<double> random_tests(n_chains, 0.0);
    vector<vector<int>> random_chains(n_chains, vector<int>(2, 0));
    vector<double> random_jumpthresholds(n_chains, 0.0);
    vector<vector<double>> random_jumptests(n_chains, vector<double>(n_theta, 0.0));

    // Assemble target func
    auto target = [&](vector<double>& theta, Observables& obs, double& l)
    {
        double p = 0;
        for (int d = 0; d < n_theta; ++d)
        {
            double pd = priors[d].LogProbability(theta[d]);
            if (pd == -std::numeric_limits<double>::infinity())
            {
                l = -std::numeric_limits<double>::infinity();
                return -std::numeric_limits<double>::infinity();
            }
            p += pd;
        }

        l = likelihood(theta, obs);

        return l + p;
    };

    // Initialize chains . . .
    if (init.empty())
    {
        // . . . from prior
        for (int c = 0; c < n_chains; ++c)
            for (int d = 0; d < n_theta; ++d)
                chains[c][d] = priors[d].RandomInit(R);
    }
    else
    {
        // . . . from initial values supplied
        if ((int)init.size() != n_chains || (int)init[0].size() != n_theta)
            throw runtime_error("init vector supplied is not the right size.");
        chains = init;
    }

    // Set initial probabilities and observables
    if (verbose)
        cout << "Initializing chains...\n";
    #pragma omp parallel for if(n_threads > 1) schedule(dynamic)
    for (int c = 0; c < n_chains; ++c)
    {
        p[c] = target(chains[c], obs[c], l[c]);
        if (verbose)
            cout << "." << flush;
    }

    // Do iterations
    if (verbose)
        cout << "\nIterating...\n";
    for (int i = init_iter; i < burn_in + iterations; ++i)
    {
        // Prepare storage and random variates
        bool migration = i < burn_in * 0.5 ? R.Bernoulli(0.05) : false;
        vector<int> migration_indices(n_chains, 0);
        if (migration)
        {
            for (int c = 0; c < n_chains; ++c)
                migration_indices[c] = R.Uniform(n_chains);
            /// for (int c = 0; c < n_chains; ++c)
            ///     migration_indices[c] = c;
            /// R.Shuffle(migration_indices.begin(), migration_indices.end());
        }
        for (int c = 0; c < n_chains; ++c)
        {
            for (int d = 0; d < n_theta; ++d)
                random_perturbations[c][d] = R.Normal(0, bstar);
            random_extensions[c] = R.Uniform(-b, b);
            random_tests[c] = R.Uniform();
            do random_chains[c][0] = R.Discrete(n_chains); while (random_chains[c][0] == c);
            do random_chains[c][1] = R.Discrete(n_chains); while (random_chains[c][1] == c || random_chains[c][1] == random_chains[c][0]);
            random_jumpthresholds[c] = (1. + R.Discrete(3)) / 3.;
            for (int d = 0; d < n_theta; ++d)
                random_jumptests[c][d] = R.Uniform();
        }
        auto saved_chains = chains;
        vector<int> accept(n_chains, 0);

        #pragma omp parallel for if(n_threads > 1) schedule(dynamic)
        for (int c = 0; c < n_chains; ++c)
        {
            vector<double> theta_p = chains[c];
            int c_from = c;

            // Generate proposal, either by migration...
            if (migration)
            {
                theta_p = saved_chains[migration_indices[c]];
                for (int d = 0; d < n_theta; ++d)
                    theta_p[d] += random_perturbations[c][d];

                /// c_from = migration_indices[c];
                /// theta_p = saved_chains[migration_indices[(c + 1) % n_chains]];
                /// for (int d = 0; d < n_theta; ++d)
                ///     theta_p[d] += random_perturbations[c][d];
            }
            else // ... or by directed mutation
            {
                int n_theta_changing = 0;
                for (int d = 0; d < n_theta; ++d)
                    if (random_jumptests[c][d] < random_jumpthresholds[c])
                        ++n_theta_changing;
                n_theta_changing = max(1, n_theta_changing);

                double gamma = (1 + random_extensions[c]) * (i % 5 == 0 ? 1.0 : 2.38 / sqrt(2 * n_theta_changing));
                for (int d = 0; d < n_theta; ++d)
                {
                    if (random_jumptests[c][d] < random_jumpthresholds[c])
                    {
                        theta_p[d] += gamma * (saved_chains[random_chains[c][1]][d] - saved_chains[random_chains[c][0]][d])
                            + random_perturbations[c][d];
                    }
                }
            }

            // Calculate log-probability and accept or reject
            Observables obs_p;
            double l_p = 0;
            double p_p = target(theta_p, obs_p, l_p);
            if ( (p_p == -std::numeric_limits<double>::infinity() && p[c_from] == -std::numeric_limits<double>::infinity() && random_tests[c] < 0.5)
                || (p_p > -std::numeric_limits<double>::infinity() && random_tests[c] < exp(p_p - p[c_from])) )
            {
                accept[c_from] = 1;
                chains[c_from] = theta_p;
                obs[c_from] = obs_p;
                p[c_from] = p_p;
                l[c_from] = l_p;
            }
        }

        // Update acceptances
        for (int c = 0; c < n_chains; ++c)
        {
            if (ar_i == ar_size - 1) show_ar = true;
            acceptances[ar_i] = accept[c];
            ar_i = (ar_i + 1) % ar_size;
        }

        // Report results of this iteration
        for (int c = 0; c < n_chains; ++c)
            report(i - burn_in, p[c], c, l[c], chains[c], obs[c]);

        // Print progress
        if (verbose)
        {
            cout << "." << flush;
            if (i % 100 == 0)
            {
                cout << "\n" << (i < burn_in ? "burn-in" : "main") << " iteration " << i - burn_in << ":";
                if (!param_names.empty())
                {
                    cout << "\n         " << setw(12) << right << "log (P)" << setw(12) << right << "log (L)";
                    for (auto n : param_names)
                        cout << setw(12) << right << n;
                }
                for (int c = 0; c < n_chains; ++c)
                {
                    cout << "\nchain" << setw(4) << right << c << setw(12) << right << p[c] << setw(12) << right << l[c];
                    for (int d = 0; d < n_theta; ++d)
                        cout << setw(12) << right << chains[c][d];
                }

                double acceptance_rate = show_ar ? (double)count(acceptances.begin(), acceptances.end(), true) / ar_size : -1;
                cout << "\nacceptance rate: " << acceptance_rate << "\n\n" << flush;
            }
        }
    }
    if (verbose)
        cout << "\n";
}

template <typename Observables, typename PFunc>
struct OptimizeData
{
    PFunc target;
    Observables obs;
    std::vector<Distribution> priors;
    bool verbose;
    unsigned int count;

    static double Test(const std::vector<double>& theta, std::vector<double>& grad, void* f_data)
    {
        (void) grad;
        auto& THIS = *(OptimizeData<Observables, PFunc>*)(f_data);
        ++THIS.count;

        // Assemble prior probability
        double p = 0;
        for (unsigned int d = 0; d < theta.size(); ++d)
        {
            double pd = THIS.priors[d].LogProbability(theta[d]);
            if (pd == std::numeric_limits<double>::lowest())
                return std::numeric_limits<double>::lowest();
            p += pd;
        }

        // Get target probability (i.e. likelihood)
        double t = THIS.target(theta, THIS.obs);

        if (THIS.verbose)
        {
            if (THIS.count % 100 == 0)
            {
                std::cout << "Optimize step " << THIS.count << ": ";
                for (auto th : theta)
                    std::cout << " " << th;
                std::cout << ", ll = " << t << ", lp = " << t + p << ".\n";
            }
        }

        return t + p;
    }
};


template <typename Observables, typename LikelihoodFunc, typename ReportFunc>
void Optimize_Priors(Randomizer* R, LikelihoodFunc likelihood, ReportFunc report, std::vector<Distribution>& priors,
    int trials = 1, double tolerance = 1e-6, bool global = true, bool verbose = false, std::vector<double>* initial = 0, double stepsize = -1)
{
    using namespace std;

    int n_theta = priors.size();
    vector<double> lb, ub;
    for (unsigned int i = 0; i < priors.size(); ++i)
    {
        lb.push_back(priors[i].LowerBound());
        ub.push_back(priors[i].UpperBound());
    }

    for (int i = 0; i < trials; ++i)
    {
        OptimizeData<Observables, LikelihoodFunc> f_data;
        f_data.target = likelihood;
        f_data.priors = priors;
        f_data.verbose = verbose;
        f_data.count = 0;
        vector<double> theta(n_theta, 0.0);

        if (initial)
            theta = *initial;
        else for (int d = 0; d < n_theta; ++d)
        {
            if (R)
                theta[d] = priors[d].RandomInit(*R);
            else
                theta[d] = (lb[d] + ub[d]) * 0.5;
        }

        double p;
        int error = 0;

        if (global)
        {
            if (verbose)
                cout << "Global optimization:\n";
            nlopt::opt go(nlopt::GN_ISRES, n_theta);
            go.set_max_objective(OptimizeData<Observables, LikelihoodFunc>::Test, &f_data);
            go.set_lower_bounds(lb);
            go.set_upper_bounds(ub);

            go.set_maxeval(100);
            go.set_ftol_abs(tolerance);

            try {
                go.optimize(theta, p);
            }
            catch (std::exception& e) {
                cout << "NLOPT exception: " << e.what() << "\n";
                error = -1;
            }
        }

        if (verbose)
            cout << "Local refinement:\n";

        nlopt::opt lo(nlopt::LN_COBYLA, n_theta);
        lo.set_max_objective(OptimizeData<Observables, LikelihoodFunc>::Test, &f_data);
        lo.set_lower_bounds(lb);
        lo.set_upper_bounds(ub);

        lo.set_ftol_abs(tolerance);

        if (stepsize > 0)
            lo.set_initial_step(stepsize);

        do {
            error = 0;
            try {
                lo.optimize(theta, p);
            }
            catch (std::exception& e) {
                cout << "NLOPT exception: " << e.what() << "\n";
                if (R != 0)
                {
                    error = -1;
                    cout << "Trying different step size.\n";
                    lo.set_initial_step(0.001 * R->Exponential(1.0));
                }
            }
        } while (error < 0);

        // error = 0;
        // try {
        //     lo.optimize(theta, p);
        // }
        // catch (std::exception& e) {
        //     cout << "NLOPT exception: " << e.what() << "\n";
        //     error = -1;

        //     std::cout << "Error: at theta";
        //     for (auto th : theta)
        //         std::cout << " " << th;
        //     std::cout << ".\n";

        //     std::cout << "Attempting smaller step size.\n";
        //     lo.set_initial_step(stepsize > 0 ? stepsize / 1000 : 0.001);

        //     error = 0;
        //     try {
        //         lo.optimize(theta, p);
        //     }
        //     catch (std::exception& e) {
        //         cout << "NLOPT exception: " << e.what() << "\n";
        //         error = -1;

        //         std::cout << "Error: at theta";
        //         for (auto th : theta)
        //             std::cout << " " << th;
        //         std::cout << ".\n";
        //     }
        // }

        double lprior = 0;
        for (unsigned int d = 0; d < theta.size(); ++d)
            lprior += priors[d].LogProbability(theta[d]);

        report(f_data.count, p, error, p - lprior, theta, f_data.obs);
    }
}

template <typename Observables, typename TargetFunc, typename ReportFunc>
void Optimize(Randomizer* R, TargetFunc target, ReportFunc report,
    std::vector<double> lb, std::vector<double> ub, int trials = 1, bool verbose = false)
{
    std::vector<Distribution> priors;
    for (unsigned int b = 0; b < lb.size(); ++b)
        priors.push_back({ Distribution::Uniform, { lb[b], ub[b] } });
    Optimize_Priors<Observables, TargetFunc, ReportFunc>
        (R, target, report, priors, trials, verbose);
}

#endif
