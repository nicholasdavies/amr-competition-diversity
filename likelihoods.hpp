// likelihoods.hpp
// Likelihood- and MCMC-related functions.

#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <iostream>
#include <omp.h>
#include "MCMC/mcmc.h"

// Helper functions: probability mass functions.
// Calculate the probability mass at k for Binomial(n, p).
double BinomialPMF(int n, double p, int k)
{
    if (k > n || k < 0) return 0.0;
    return boost::math::pdf(boost::math::binomial(n, p), k);
}

// Helper to calculate the log binomial coefficient for N choose n; inspired by Stan's implementation (mc-stan.org)
double LogBinomialCoefficient(double N, double n)
{
    using boost::math::lgamma;
    const double Threshold = 1000;
    if (N - n < Threshold)
        return lgamma(N + 1) - lgamma(n + 1) - lgamma(N + 1 - n);
    return n * log(N - n) + (N + 0.5) * log(N / (N - n)) + 1. / (12. * N) - n - 1. / (12. * (N - n)) - lgamma(n + 1);
}

// Calculate the log probability mass at k for Binomial(n, p).
double LogBinomialPMF(int n, double p, int k)
{
    if (k > n || k < 0)
        return -std::numeric_limits<double>::infinity();

    return LogBinomialCoefficient(n, k) + k * log(p) + (n - k) * log(1 - p);
}

// Calculate the probability density at x for Beta(a, b)
double BetaPDF(double a, double b, double x)
{
    return boost::math::pdf(boost::math::beta_distribution<double>(a, b), x);
}

// Calculate the log probability density at x for Beta(a, b)
double LogBetaPDF(double a, double b, double x)
{
    using boost::math::lgamma;
    if (x < 0 || x > 1)     return -std::numeric_limits<double>::infinity();
    if (a == 1 && b == 1)   return 0;
    return lgamma(a + b) - lgamma(a) - lgamma(b) + (a - 1) * log(x) + (b - 1) * log(1 - x);
}

// Calculate the probability density at x for Normal(mu, sigma)
double NormPDF(double mu, double sigma, double x)
{
    return boost::math::pdf(boost::math::normal_distribution<double>(mu, sigma), x);
}

// Calculate the log probability density at x for Normal(mu, sigma)
double LogNormPDF(double mu, double sigma, double x)
{
    return -(x - mu) * (x - mu) / (2 * sigma * sigma) - log(sigma) - 0.918938533204673;
}

// Calculate the cumulative density at x for Normal(mu, sigma)
double NormCDF(double mu, double sigma, double x)
{
    return boost::math::cdf(boost::math::normal_distribution<double>(mu, sigma), x);
}

// Calculate the inverse of the CDF for Normal(mu, sigma)
double InvNormCDF(double mu, double sigma, double q)
{
    return boost::math::quantile(boost::math::normal_distribution<double>(mu, sigma), q);
}

// Clamp the argument x between 0 and 1.
double Clamp(double x)
{
    return max(0.0, min(1.0, x));
}

// Helper for beta model
extern "C" double beta_f(double x, void * params)
{
    double alpha = ((double*)params)[0];
    double beta = ((double*)params)[1];
    int data_r = ((double*)params)[2];
    int data_n = ((double*)params)[3];
    return BetaPDF(alpha, beta, x) * BinomialPMF(data_n, x, data_r);
}

// Helper for normal model
extern "C" double normal_f(double x, void * params)
{
    double mu = ((double*)params)[0];
    double sigma = ((double*)params)[1];
    int data_r = ((double*)params)[2];
    int data_n = ((double*)params)[3];

    return NormPDF(mu, sigma, x) * BinomialPMF(data_n, x, data_r);
}

// Model predicted resistance prevalence is model_R;
// Empirical resistance prevalence is data_r resistant isolates out of data_n total isolates.
double DoLogLikelihood(double model_R, int data_r, int data_n, double l_noise)
{
    if (*GSL_Workspace == 0)
        *GSL_Workspace = gsl_integration_workspace_alloc(1000);
  
    double likelihood = 0;

    if (P.l_model == "binomial")
    {
        return LogBinomialPMF(data_n, max(1e-6, min(1 - 1e-6, model_R)), data_r);
    }
    else if (P.l_model == "error")
    {
        for (int n_in_error = 0; n_in_error <= data_n; ++n_in_error)
        {
            int n_not_in_error = data_n - n_in_error;
            double P_n_in_error = BinomialPMF(data_n, l_noise, n_in_error);
            for (int r_in_error = 0; r_in_error <= data_r; ++r_in_error)
            {
                int r_not_in_error = data_r - r_in_error;
                likelihood += P_n_in_error * BinomialPMF(n_in_error, 0.5, r_in_error) * BinomialPMF(n_not_in_error, model_R, r_not_in_error);
            }
        }
    }
    else if (P.l_model == "beta")
    {
        // Parameters of beta distribution
        double a = model_R * (l_noise - 2) + 1;
        double b = (1 - model_R) * (l_noise - 2) + 1;

        double params[4] = { a, b, (double)data_r, (double)data_n };
        gsl_function F;
        F.function = &beta_f;
        F.params = &params;
        double error;

        if (int ret = gsl_integration_qags(&F, 0, 1, 0, 1e-7, 1000, *GSL_Workspace, &likelihood, &error))
            cout << "\nWARNING: GSL integrator returned error in integration of beta noise.\n"
                 << "params " << a << " " << b << " " << data_r << " " << data_n << " " << likelihood << " " << error << " error code " << ret << "\n";
    }
    else if (P.l_model == "normal_clamp" || P.l_model == "normal_trunc")
    {
        double params[4] = { model_R, l_noise, (double)data_r, (double)data_n };
        gsl_function F;
        F.function = &normal_f;
        F.params = &params;
        double error;

        if (int ret = gsl_integration_qags(&F, 0, 1, 0, 1e-7, 1000, *GSL_Workspace, &likelihood, &error))
            cout << "\nWARNING: GSL integrator returned error in integration of normal noise.\n"
                 << "params " << model_R << " " << l_noise << " " << data_r << " " << data_n << " " << likelihood << " " << error << " error code " << ret << "\n";

        if (P.l_model == "normal_clamp")
            likelihood += NormCDF(model_R, l_noise, 0.) * BinomialPMF(data_n, 0., data_r)
                 +  (1. - NormCDF(model_R, l_noise, 1.)) * BinomialPMF(data_n, 1., data_r);
        else // i.e. normal_trunc
            likelihood /= (NormCDF(model_R, l_noise, 1) - NormCDF(model_R, l_noise, 0));
    }
    else if (P.l_model == "mean")
    {
        return LogNormPDF((data_r + 1.0) / (data_n + 2.0), l_noise, model_R);
    }
    else if (P.l_model == "none")
    {
        return 0;
    }
    else
    {
        throw runtime_error("Invalid likelihood model specified.");
    }

    return log(likelihood);
}

double GetLogLikelihoodFromResults(Results& res)
{
    if (res.zero_likelihood)
        return -std::numeric_limits<double>::infinity();

    // Calculate model-predicted resistant carriers & total carriers
    vector<double> r(res.Rs.size() / P.group_span, 0.0);
    vector<double> inf(r.size(), 0.0);
    vector<double> h(r.size(), 0.0);

    for (unsigned int k = 0; k < r.size(); ++k)
    {
        for (int g = 0; g < P.group_span; ++g)
        {
            r[k]     += res.Rs   [k * P.group_span + g] * V->h[k * P.group_span + g];
            inf[k]   += res.Infs [k * P.group_span + g] * V->h[k * P.group_span + g];
            h[k]     +=                                   V->h[k * P.group_span + g];
        }
    }

    // Calculate log-likelihood of each data point
    double ll = 0;
    double avg_population = accumulate(P.population.begin(), P.population.end(), 0.0) / P.population.size();
    for (unsigned int k = 0; k < P.mcmc_data_r.size(); ++k)
    {
        // Resistance prevalence - fit to model prediction plus dispersion from "l_noise" model
        double modelR = (inf[k] == 0 ? 0 : Clamp(r[k] / inf[k]));
        double pll = DoLogLikelihood(modelR, P.mcmc_data_r[k], P.mcmc_data_n[k], V->l_noise);
        if (P.l_weight)
            ll += pll * P.population[k] / avg_population;
        else
            ll += pll;

        if (!P.l_yi.empty())
        {
            Distribution distYi(P.l_yi);
            ll += distYi.LogProbability(Clamp(inf[k] / h[k]));
        }
    }

    // Fit average carriage and average resistance prevalence to given distribution
    if (!P.l_y.empty())
    {
        Distribution distY(P.l_y);
        ll += distY.LogProbability(res.InfTotal);
    }
    if (!P.l_d.empty())
    {
        Distribution distD(P.l_d);
        ll += distD.LogProbability(res.SR / res.InfTotal);
    }
    if (!P.l_rf.empty())
    {
        Distribution distRf(P.l_rf);
        ll += distRf.LogProbability(res.RFrac);
    }
    
    return ll;
}

void MCMCSetParameter(const std::vector<double>& theta, int& index, std::string& P_param, double& v_param)
{
    if (P_param.size() > 0 && isalpha(P_param[0]))
        v_param = theta[index++];
};

void MCMCSetParameter(const std::vector<double>& theta, int& index, std::string& P_param, vector<double>& v_param)
{
    vector<string> parts = Unserialize(P_param, ",;");

    if (parts.size() == 1)
    {
        if (parts[0].size() > 0 && isalpha(parts[0][0]))
        {
            double x = theta[index++];
            for (auto& v : v_param)
                v = x;
        }
    }
    else if (parts.size() == v_param.size())
    {
        for (unsigned int p = 0; p < parts.size(); ++p)
            if (parts[p].size() > 0 && isalpha(parts[p][0]))
                v_param[p] = theta[index++];
    }
    else
    {
        throw runtime_error("Size mismatch for parameter " + P_param);
    }
};

void MCMCSetParameters(const std::vector<double>& theta)
{
    int index = 0;

    #define ODE_PARAM(type, name) MCMCSetParameter(theta, index, P.name, V->name);
    DO_ODE_PARAMS;
    #undef ODE_PARAM
}

double MCMCTarget(const std::vector<double>& theta, Results& observables)
{
    MCMCSetParameters(theta);

    observables = Trial(-1);

    double ll = GetLogLikelihoodFromResults(observables);

    return ll;
}

void DoMCMC(std::ostream& out)
{
    SetParameters(0);
    vector<Distribution> priors;
    vector<string> param_names;

    auto add_parameter = [&](string& P_param, string name)
    {
        vector<string> parts = Unserialize(P_param, ",;");
        for (unsigned int p = 0; p < parts.size(); ++p)
        {
            vector<string> subcode = Unserialize(parts[p]);
            if (isalpha(subcode[0][0]))
            {
                priors.push_back(Distribution(subcode));
                param_names.push_back(parts.size() == 1 ? name : name + to_string(p));
            }
        }
    };

    #define ODE_PARAM(type, name) add_parameter(P.name, #name);
    DO_ODE_PARAMS;
    #undef ODE_PARAM

    auto report = [&](int trial, double lp, int chain, double ll, std::vector<double>& params, Results& obs)
    {
        if (trial % P.mcmc_thinning == 0)
            ReportResults(out, obs, true, trial, lp, chain, ll, params);
    };

    int n_chains = priors.size() * 2;
    if (P.mcmc_algorithm >= 3)
        n_chains = 1;

    // Optionally, resume MCMC from specified output file
    vector<vector<double>> init;
    int init_iter = 1;
    if (!P.mcmc_resume.empty())
    {
        istringstream in(Exec("tail -n " + to_string(n_chains) + " " + P.mcmc_resume + " | cut -f 1,5-" + to_string(5 + param_names.size() - 1)));
        for (int c = 0; c < n_chains; ++c)
        {
            init.push_back(vector<double>());
            in >> init_iter;
            for (unsigned int d = 0; d < param_names.size(); ++d)
            {
                double x; in >> x;
                init.back().push_back(x);
            }
        }
        init_iter = init_iter + P.mcmc_burn_in + 1;
    }

    if (P.mcmc_algorithm == 1)
    {
        DEMCMC_Priors<Results>(R, MCMCTarget, report, P.mcmc_burn_in, (P.trials + n_chains - 1) / n_chains,
            n_chains, priors, true, param_names, P.threads > 1, P.threads, true, P.mcmc_R_skip, init, init_iter);
    }
    else if (P.mcmc_algorithm == 2)
    {
        DREAM_Options options;
        options.burn_in = P.mcmc_burn_in;
        options.n_chains = priors.size();
        options.iterations = (P.trials + options.n_chains - 1) / options.n_chains;
        options.verbose = true;
        options.param_names = param_names;
        options.n_threads = P.threads;
        options.init_theta = init;
        options.init_iter = init_iter;
        DREAM<Results>(R, MCMCTarget, report, priors, options);
    }
    else if (P.mcmc_algorithm < 0)
    {
        auto par = find(param_names.begin(), param_names.end(), P.find_param);
        if (par == param_names.end())
            throw logic_error("Unrecognised param name.");
        unsigned int par_index = par - param_names.begin();

        double xmin = priors[par_index].LowerBound();
        double xmax = priors[par_index].UpperBound();

        for (int i = 0; i < -P.mcmc_algorithm; ++i)
        {
            auto priors2 = priors;
            double x = xmin + ((xmax - xmin) * i) / (-P.mcmc_algorithm - 1.0);
            vector<double> params {x, x + 1e-6};
            priors2[par_index] = Distribution(Distribution::Uniform, params);
            Optimize_Priors<Results>(&R, MCMCTarget, report, priors2, 1, P.epsilon, false, P.verbose);
        }
    }
    else if (P.mcmc_algorithm >= 3)
    {
        Optimize_Priors<Results>(&R, MCMCTarget, report, priors, 1, P.epsilon, false, P.verbose, init.empty() ? 0 : &init[0], P.mcmc_algorithm == 4 ? 0.001 : -1);
    }
}