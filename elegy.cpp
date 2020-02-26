#include "Config/config.h"
#include "Randomizer/randomizer.h"
#include "String/string.hpp"
#include "Boundary/boundary.h"
#include "Thread/thread.h"
#include <boost/numeric/odeint.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <memory>
#include <array>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <unistd.h>
#include <cstdio>
using namespace boost::numeric::odeint;
using namespace boost::math;
using namespace std;

//// REFACTORING TODO / WISH LIST
//// - For play_changes mode, save "original" trial number as well as intervention name, strength, & time. This will help with matching interventions. This will require changing elegy.interventions.R.
//// - For play_changes mode, perhaps separately save trial number and intervention details? (less important, but will simplify code in elegy.interventions.R)
//// - Save actual index of e.g. x (& other parameters) rather than 0 for first, 1 for second. This will also require changing config files.

// ODE parameters
#define DO_ODE_PARAMS \
    ODE_PARAM(vector<double>,  beta)            \
    ODE_PARAM(vector<double>,  c)               \
    ODE_PARAM(vector<double>,  u)               \
    ODE_PARAM(vector<double>,  tau)             \
    ODE_PARAM(vector<double>,  k)               \
    ODE_PARAM(vector<double>,  kk)              \
    ODE_PARAM(vector<double>,  b)               \
    ODE_PARAM(vector<double>,  b0)              \
    ODE_PARAM(vector<double>,  a)               \
    ODE_PARAM(vector<double>,  delta)           \
    ODE_PARAM(vector<double>,  g)               \
    ODE_PARAM(vector<double>,  G)               \
    ODE_PARAM(vector<double>,  waifw)           \
    ODE_PARAM(vector<double>,  h)               \
    ODE_PARAM(vector<double>,  psi)             \
    ODE_PARAM(vector<double>,  psi_r)           \
    ODE_PARAM(vector<double>,  z)               \
    ODE_PARAM(vector<double>,  shape)           \
    ODE_PARAM(double,          l_noise)         \
    ODE_PARAM(double,          v)               \
    ODE_PARAM(double,          dv)              \
    ODE_PARAM(double,          sv)              \
    ODE_PARAM(vector<double>,  tv)              \
    ODE_PARAM(vector<double>,  x)

struct ODEParams
{
    #define ODE_PARAM(type, name) type name;
    DO_ODE_PARAMS;
    #undef ODE_PARAM
};

void PrintParam(ostream& out, double v)
{
    out << v << "\n";
}

void PrintParam(ostream& out, vector<double>& v)
{
    for (auto& x : v)
        out << x << " ";
    out << "\n";
}

void PrintParams(ostream& out, ODEParams& v)
{
    #define ODE_PARAM(type, name) { out << #name << " = "; PrintParam(out, v.name); }
    DO_ODE_PARAMS;
    #undef ODE_PARAM
}

// Global objects
Parameters P;
Randomizer R;
OMPWrapper<gsl_integration_workspace*> GSL_Workspace(0);
OMPWrapper<ODEParams> V;

vector<string> Unserialize(string s, string sep = " \t\n")
{
    vector<string> sub;
    string::size_type start = 0, end = 0;
    while (end != string::npos)
    {
        end = s.find_first_of(sep, start);
        if (end == start)
        {
            ++start;
            continue;
        }
        sub.push_back(trim(s.substr(start, end - start)));
        start = end + 1;
    }
    return sub;
}

string Exec(string cmd)
{
    array<char, 128> buffer;
    string result;
    unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);
    if (!pipe) throw runtime_error("popen() failed!");
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
        result += buffer.data();
    return result;
}

void SetParameter(double& param, string name, string code, int trial, int subscript = -1);

void SetParameter(vector<double>& param, string name, string code, int trial, int subscript = -1)
{
    if (code.empty())
        return;

    // Pipes separate trials; semicolons/commas separate groups within a trial.
    if (code.find('|') != string::npos)
    {
        vector<string> trial_parts = Unserialize(code, "|");
        if (trial_parts.size() != (unsigned int)P.trials)
            throw runtime_error("Mismatch with number of trials for parameter " + name + ".");
        code = trial_parts[trial];
    }

    if (code.find(',') != string::npos || code.find(';') != string::npos)
    {
        vector<string> parts = Unserialize(code, ",;");
        if (parts.size() != param.size())
            throw runtime_error("Mismatch with number of elements for parameter " + name + ".");
        for (unsigned int i = 0; i < param.size(); ++i)
            SetParameter(param[i], name, parts[i], trial);
        return;
    }

    try
    {
        if (subscript < 0)
            fill(param.begin(), param.end(), stod(code));
        else if (subscript < (int)param.size())
            param[subscript] = stod(code);
    }
    catch (invalid_argument&)
    {
        if (!isalpha(code[0]))
            throw runtime_error("Could not interpret code " + code + " for parameter " + name + ".");
    }
}

void SetParameter(double& param, string name, string code, int trial, int subscript)
{
    if (subscript >= 0)
        throw runtime_error("Attempting to index-assign code " + code + " to atomic variable " + name + ".");
    vector<double> v(1, param);
    SetParameter(v, name, code, trial, subscript);
    param = v[0];
}

void LowLevelSet(double& a, int s, double v)
{
    (void)s;
    a = v;
}

void LowLevelSet(vector<double>& a, int s, double v)
{
    if (s < 0) fill(a.begin(), a.end(), v);
    else if (s >= (int)a.size()) throw runtime_error("Subscript " + to_string(s) + " exceeds bounds of parameter.");
    else a[s] = v;
}

bool DoXFunc(int trial = 0)
{
    #pragma omp critical
    {
        for (size_t i = 0; i < V->x.size(); ++i)
            P.xfunc << V->x[i];
        P.xfunc << trial;

        int n_ret = P.xfunc.CallMultRet();

        if (n_ret > 1 && n_ret % 3 != 0)
        {
            throw std::runtime_error("P.xfunc must return a multiple of 3 number of values.");
        }

        if (n_ret > 1)
        {
            for (int i = 0; i < n_ret; i += 3)
            {
                string p_name;
                int p_subscript;
                double p_value;
                P.xfunc >> p_name >> p_subscript >> p_value;
                #define ODE_PARAM(type, name) else if (p_name == #name) LowLevelSet(V->name, p_subscript, p_value);
                if (false) ;
                DO_ODE_PARAMS
                else throw std::runtime_error("P.xfunc attempting to set non-existent ODE parameter " + p_name + ".");
                #undef ODE_PARAM
            }
        }
        P.xfunc.FlushRet();
    }

    // Construct compartment sizes from populations
    if (!P.population.empty())
    {
        double total_population = accumulate(P.population.begin(), P.population.end(), 0.0);
        if ((int)P.population.size() == P.groups)
            for (int g = 0; g < P.groups; ++g)
                V->h[g] = P.population[g] / total_population;
        else if ((int)P.population.size() * P.treatment_regions == P.groups)
            for (int g = 0; g < P.groups; ++g)
                V->h[g] = P.population[g / P.treatment_regions] / (total_population * P.treatment_regions);
        else
            throw std::runtime_error("Population sizes not compatible with number of groups and treatment regions.");
    }

    // Construct treatment rates from consumption rates (in DDD)
    vector<double> treatment;
    if (P.treatment_course >= 0)
    {
        treatment = P.consumption;
        transform(treatment.begin(), treatment.end(), treatment.begin(), [](double x) { return x * 365.25 / (12000 * P.treatment_course); });
    }
    else
    {
        treatment = V->tau;
    }

    // Shift treatment rates with z
    if (P.treatment_shift == "none" && V->z.size() == 0)
        ;
    else if (P.treatment_shift == "all" && V->z.size() == 1)
        transform(treatment.begin(), treatment.end(), treatment.begin(), [](double x) { return x * V->z[0]; });
    else if (P.treatment_shift == "each" && V->z.size() == treatment.size())
        transform(treatment.begin(), treatment.end(), V->z.begin(), treatment.begin(), [](double x, double z) { return x * z; });
    else
        throw std::runtime_error("treatment_shift must be either none, all, or each, and the size of z must correspond.");

    // Assign to sub-national regions
    if (P.treatment_regions > 0)
    {
        if (P.groups % P.treatment_regions != 0)
            throw runtime_error("number of groups must be a multiple of the number of treatment regions.");

        int n_countries = P.groups / P.treatment_regions;

        if (V->shape[0] <= 0)
            return false;

        auto GammaQuantile = [](double q, double shape, double scale) {
            return boost::math::quantile(boost::math::gamma_distribution<double>(shape, scale), q);
        };

        // Assign treatment rates to tau
        for (int c = 0; c < n_countries; ++c)
        {
            for (int r = 0; r < P.treatment_regions; ++r)
            {
                double X0 = boost::math::gamma_q(V->shape[c] + 1, GammaQuantile(double(r) / P.treatment_regions, V->shape[c], 1.)), X1 = 0;
                if (r < P.treatment_regions - 1)
                    X1 = boost::math::gamma_q(V->shape[c] + 1, GammaQuantile(double(r + 1) / P.treatment_regions, V->shape[c], 1.));
                V->tau[c * P.treatment_regions + r] = treatment[c] * P.treatment_regions * (X0 - X1);
            }
        }
    }
    else
    {
        V->tau = treatment;
    }

    return true;
}

bool SetParameters(int trial)
{
    // Individual parameters
    #define ODE_PARAM(type, name) SetParameter(V->name, #name, P.name, trial);
    DO_ODE_PARAMS;
    #undef ODE_PARAM

    return DoXFunc(trial);
}

// Results of a run
struct Results
{
    Results()
     : S(0), R(0), SR(0), RTotal(0), InfTotal(0), RFrac(0), max_epochs(-1), abc_pass(false), zero_likelihood(false) { }

    double S, R, SR, RTotal, InfTotal, RFrac;
    int max_epochs;
    vector<double> compartments;
    vector<double> Rs, Duals, Infs;
    bool abc_pass;
    bool zero_likelihood;

    vector<double> saved_Y;
};

#include "metapopulation.hpp"
#include "subpopulations.hpp"
#include "results.hpp"

Results FunctionalFit(int trial);

Results Trial(int trial = 0, vector<double>* iY = 0, string in = "", double ic = 0, double is = 0, double it0 = 0, double it1 = 0)
{
    switch (P.model)
    {
        case 0: { return FunctionalFit(trial); }
        case 1: { Metapopulation<PopulationSingle> mp;      return mp.Trial(trial, iY, in, ic, is, it0, it1); }
        case 2: { Metapopulation<PopulationWithinHost> mp;  return mp.Trial(trial, iY, in, ic, is, it0, it1); }
        case 3: { Metapopulation<PopulationDType> mp;       return mp.Trial(trial, iY, in, ic, is, it0, it1); }

        default:
            throw runtime_error("Unsupported model.");
            break;
    }
}

#include "likelihoods.hpp"

class FindTest
{
public:
    string param, goal;
    double value;
    int index;

    FindTest(string p, string g, double v, int i = 0)
     : param(p), goal(g), value(v), index(i) { }

    bool operator()(double x)
    {
        P.Set(param, to_string(x));
        SetParameters(0);
        auto res = Trial(0);
        if (goal == "Y")
            return res.InfTotal > value;
        else if (goal == "Yi")
            return res.Infs[index] > value;
        else if (goal == "Rfrac")
            return res.RFrac > value;
        else if (goal == "Rfraci")
            return res.Rs[index] / res.Infs[index] > value;
        else
            throw std::runtime_error("Unrecognised find goal.");
    }
};

Results FunctionalFit(int trial)
{
    SetParameters(trial);
    Results res;
    for (int g = 0; g < P.groups; ++g)
    {
        double RFrac = 0;
        for (int j = 0; j < (int)V->x.size(); ++j)
            RFrac += V->x[j] * pow(V->tau[g], j);
        RFrac = Clamp(RFrac);
        res.compartments.push_back(RFrac);
        res.Rs.push_back(RFrac);
        res.Duals.push_back(0);
        res.Infs.push_back(1);
    }
    return res;
}

int main(int argc, char* argv[])
{
    // Read parameters
    P.Read(argc, argv);

    // Initialisation
    gsl_set_error_handler_off();
    PopulationWithinHost::Setup();
    auto elements = [](string& V) -> int { if (V.empty()) return 0; else return count(V.begin(), V.end(), ',') + 1; };
    for (auto& v : V)
    {
        int ng = P.groups, nG = 1;
        if (P.treatment_regions > 0) 
        { 
            ng = P.groups / P.treatment_regions;
            nG = P.groups / P.treatment_regions; 
        }

        v.beta.resize(P.groups, 0);
        v.c.resize(P.groups, 0);
        v.u.resize(P.groups, 0);
        v.tau.resize(P.groups, 0);
        v.k.resize(P.groups, 0);
        v.kk.resize(P.groups, 0);
        v.b.resize(P.groups, 0);
        v.b0.resize(P.groups, 0);
        v.a.resize(P.groups, 0);
        v.delta.resize(P.groups, 0);
        v.g.resize(ng, 0);
        v.G.resize(nG, 0);
        v.waifw.resize(elements(P.waifw), 0);
        v.h.resize(P.groups, 0);
        v.tv.resize(P.groups, 0);
        v.psi.resize(elements(P.psi), 0);
        v.psi_r.resize(elements(P.psi_r), 0);
        v.z.resize(elements(P.z), 0);
        v.shape.resize(max(1, P.groups / P.treatment_regions), 0);
        v.x.resize(elements(P.x), 0);
    }

    if (!P.find)
        P.Write(cout);
    
    if (P.play_changes)
    {
        if (!(P.interv.size() == P.interv_coverage.size() && P.interv.size() == P.interv_strength.size() && P.interv.size() == P.interv_time.size()))
            throw runtime_error("interv, interv_coverage, interv_strength and interv_time must all have the same number of entries.");

        omp_set_num_threads(P.threads);

        // Preprocess file
        string sorted_file = string(P_tmpdir) + "/elegy.XXXXXXXX";
        mktemp(&*sorted_file.begin());
        cout << "Processing changes file to " << sorted_file << "...\n";

        // write header to sorted_file
        system(("grep -m 1 -h ^[^#] " + P.changes_file + " > " + sorted_file).c_str());

        // select, sort, and optionally sample lines to sorted_file
        if (P.changes_limit > 0)
        {
            double span = stod(Exec("grep ^[0-9] " + P.changes_file + " | wc -l"));
            if (P.changes_limit > 1)
                span = (span - 1) / (P.changes_limit - 1);
            
            system(("grep -h ^[0-9] " + P.changes_file +
                    " | sort -n -k3,3 -k1,1 -S128M" +
                    " | awk 'BEGIN {t = 1; s = " + to_string(span) + ";} NR >= t - 1e-6 { print $0; t += s; }'" +
                    " >> " + sorted_file).c_str());
        }
        else
        {
            system(("grep -h ^[0-9] " + P.changes_file +
                    " | sort -n -k3,3 -k1,1 -S128M >> " + sorted_file).c_str());
        }

        cout << "Running changes...\n";

        // Load changes
        ifstream changes(sorted_file);
        string line;
        vector<string> columns, lines, outlines;
        vector<int> counts;
        int tot_lines = 0, col0 = 0, col1 = 0;

        ofstream fout(P.fileout);
        Header(fout, true, true);

        // Evaluate each line of changes file
        while (changes.good())
        {
            // Read some lines
            do
            {
                if (line.empty() || line[0] == '#' || line[0] == '-')
                {
                    continue;
                }
                else if (columns.empty())
                {
                    columns = unserialize(line, "\t");
                    auto c0 = find(columns.begin(), columns.end(), "ll");
                    auto c1 = find(columns.begin(), columns.end(), "S");
                    if (c0 == columns.end() || c1 == columns.end() || c0 >= c1 - 1)
                    {
                        cout << "Supplied changes file not formatted for processing (params must be between ll and S columns).\n";
                        remove(sorted_file.c_str());
                        exit(1);
                    }
                    col0 = c0 - columns.begin() + 1;
                    col1 = c1 - columns.begin();
                    columns = vector<string>(columns.begin() + col0, columns.begin() + col1);
                }
                else
                {
                    vector<string> s = unserialize(line, "\t");
                    s = vector<string>(s.begin() + col0, s.begin() + col1);
                    string pline = serialize(s, "\t");

                    if (!lines.empty() && pline == lines.back())
                    {
                        ++counts.back();
                    }
                    else
                    {
                        if (lines.size() >= 100)
                            break;
                        lines.push_back(pline);
                        counts.push_back(1);
                    }
                }
            } while (getline(changes, line));

            // Process the lines
            outlines.resize(lines.size());
            #pragma omp parallel for schedule(dynamic)
            for (unsigned int l = 0; l < lines.size(); ++l)
            {
                // Set parameters for run
                vector<string> run_data = unserialize(lines[l], "\t");
                for (unsigned int p = 0; p < columns.size(); ++p)
                {
                    string alpha = columns[p];
                    int num = -1;
                    auto sep = columns[p].find_first_of("01234567890");
                    if (sep != string::npos)
                    {
                        alpha = columns[p].substr(0, sep);
                        num = stoi(columns[p].substr(sep));
                    }

                    #define ODE_PARAM(type, name) if (alpha == #name) SetParameter(V->name, #name, run_data[p], tot_lines + l, num);
                    DO_ODE_PARAMS;
                    #undef ODE_PARAM
                }

                // Run all interventions
                ostringstream sout;
                Results res_pre, res;
                bool run_pre = any_of(P.interv_time.begin(), P.interv_time.end(), [&](double t) { return t >= 0; });
                if (run_pre) res_pre = Trial(tot_lines + l);

                for (unsigned int i = 0; i < P.interv.size(); ++i)
                {
                    // Do run
                    ODEParams V0 = *V;

                    // If running to equilibrium, start from previous run if we can (as initial conditions shouldn't matter), for efficiency.
                    if (P.interv_time[i] < 0)
                        res = Trial(tot_lines + l, res.saved_Y.size() > 0 ? &res.saved_Y : 0, P.interv[i], P.interv_coverage[i], P.interv_strength[i], 0.0, P.interv_time[i]);
                    // If running same intervention to a later time step, start from previous run if we can.
                    else if (i > 0 && P.interv[i] == P.interv[i - 1] && P.interv_coverage[i] == P.interv_coverage[i - 1]
                             && P.interv_strength[i] == P.interv_strength[i - 1] && P.interv_time[i] >= P.interv_time[i - 1])
                        res = Trial(tot_lines + l, &res.saved_Y, P.interv[i], P.interv_coverage[i], P.interv_strength[i], P.interv_time[i - 1], P.interv_time[i]);
                    // Otherwise, start from saved run above (which is pre-intervention scenario).
                    else
                        res = Trial(tot_lines + l, &res_pre.saved_Y, P.interv[i], P.interv_coverage[i], P.interv_strength[i], 0.0, P.interv_time[i]);

                    ostringstream interv_name;
                    interv_name << P.interv[i] << ";" << P.interv_coverage[i] << ";" << P.interv_strength[i] << ";" << P.interv_time[i];
                    ReportResults(sout, res, true, -99999, 1, -1, 1, vector<double>(), interv_name.str());
                    *V = V0;
                }
                outlines[l] = sout.str();
            }

            // Output the results
            for (int l = 0; l < (int)outlines.size(); ++l)
                for (int c = 0; c < counts[l]; ++c)
                    fout << outlines[l];

            // Reset to read more lines
            cout << "\n" << flush;
            tot_lines += lines.size();
            outlines.clear();
            lines.clear();
            counts.clear();
        }
        cout << "\n";

        // Remove temporary file
        remove(sorted_file.c_str());
    }
    else if (P.mcmc)
    {
        if (P.fileout == P.mcmc_resume)
            throw runtime_error("Attempting to resume from and write to the same file.");
        ofstream fout(P.fileout);
        Header(fout, true, true);
        DoMCMC(fout);
    }
    else if (P.find)
    {
        cout << FindBoundary(FindTest(P.find_param, P.find_goal, P.find_value, P.find_i), P.find_x0, P.find_x1, P.find_n, true, 1);
    }
    else if (P.autofind)
    {
        // Part I: find
        string suffix = (P.autofind_i >= 0 ? "i" : "");
        string saved_tau = P.tau;
        cout << "Finding beta_s...\n";
        FindBoundary(FindTest("beta_s", "Y" + suffix, P.autofind_y, 0), 0.0, 10.0, 50, true, 1, true);
        cout << "Finding " << P.autofind_cost << "...\n";
        if (P.autofind_i < 0)
            P.Set("tau", to_string(P.autofind_tau));
        FindBoundary(FindTest(P.autofind_cost, "Rfrac" + suffix, 0.5, P.autofind_i), 0.0, 10.0, 50, true, 0, true);
        P.Set("tau", saved_tau);

        // Part II: run
        P.fileout += ".autofind.txt";
        ofstream fout(P.fileout);
        Header(fout, true, true);
        for (int tr = 0; tr < P.trials; ++tr)
        {
            // Do a trial
            Results res = Trial(tr);

            // Save results and update progress
            ReportResults(fout, res, true, tr);
            if (tr % 100 == 0)
                cout << "." << flush;
        }
        cout << "\n";
    }

    else
    {
        ofstream fout(P.fileout);
        Header(fout, true, false);

        for (int tr = 0; tr < P.trials; ++tr)
        {
            // Do a trial
            Results res = Trial(tr);

            // Save results and update progress
            ReportResults(fout, res, false, tr);
            if (tr % 100 == 0)
                cout << "." << flush;
        }
        cout << "\n";
    }

    // Cleanup
    for (auto& w : GSL_Workspace)
        gsl_integration_workspace_free(w);

    return 0;
}