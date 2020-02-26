// metapopulation.hpp
// General implementation of a metapopulation with various network structures and plugin dynamics for subpopulations.

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <nlopt.hpp>

double AbsoluteDifference(vector<double>& y1, vector<double>& y2)
{
    double s = 0.0;
    for (size_t i = 0; i < y1.size(); ++i)
        s += fabs(max(0.0, y1[i]) - max(0.0, y2[i]));
    return s;
}

// class Metapopulation<typename Subpopulation>
// Subpopulation must have four static members:
//  ::Elements(), the number of explicit elements/compartments in a subpopulation
//  ::InitialState(), returning a vector<double> with the initial state for a given subpopulation
//  ::Gradient(const int g, const vector<double>& Y, vector<double>& dYdt, const vector<double>& contact, const double t), updating dYdt for subpopulation g according to current state Y, contacts, and time t
//  vector<string> ::Compartments(), returning the names of the compartments that Report gives
//  vector<double> ::Report(const vector<double>& Y), returning the frequencies of the compartments
template <typename Subpopulation>
class Metapopulation
{
public:
    Metapopulation() { }

    Results Trial(int trial, vector<double>* init_Y = 0, string intervention = "", double intervention_coverage = 0,
        double intervention_strength = 0, double intervention_t0 = 0, double intervention_t1 = 0)
    {
        Results res;

        if (!SetParameters(trial))
        {
            res.zero_likelihood = true;
            auto fake_Y = InitialState();
            Outcome(fake_Y, res, 0);
            return res;
        }
        MakeContactMatrix();
        SingleRun(res, init_Y, intervention, intervention_coverage, intervention_strength, intervention_t0, intervention_t1);
        
        // Calculate overall resistance fraction
        res.RFrac = res.RTotal / res.InfTotal;

        return res;
    }

    void SingleRun(Results& res, vector<double>* init_Y, string intervention, double intervention_coverage, double intervention_strength, double intervention_t0, double intervention_t1)
    {
        Intervene(intervention, intervention_coverage, intervention_strength);

        vector<double> Y, Y0;
        if (init_Y == 0)
            Y0 = Y = InitialState();
        else
            Y0 = Y = *init_Y;

        int most_epochs = 0;

        // Iterate over each separated block of groups
        int group_leap = P.separate_groups ? P.group_connectivity : P.groups;
        for (g0 = 0; g0 < P.groups; g0 += group_leap)
        {
            // Extract block for integration
            g1 = g0 + group_leap;
            vector<double> y(Y.begin() + g0 * Subpopulation::Elements(), Y.begin() + g1 * Subpopulation::Elements());
            vector<double> y0(y);

            // Solve for equilibrium of y, or run y until equilibrium or max_epochs, or from intervention_t0 to intervention_t1 if running intervention
            int epochs = 1;
            if (P.root && intervention_t1 <= 0)
            {
                if (g0 > 0 && P.start_previous && !init_Y)
                    y = vector<double>(Y.begin() + (g0 - group_leap) * Subpopulation::Elements(), Y.begin() + g0 * Subpopulation::Elements());

                epochs = FindRoot(y);

                if (epochs < 0) // If root finding was unsuccessful, alternate ODE integration and more root finding
                {
                    y = y0;
                    for (epochs = 0; epochs < P.max_epochs; ++epochs)
                    {
                        integrate(*this, y, 0.0, P.t_epoch, P.t_step);
                        if (AbsoluteDifference(y, y0) < P.epsilon)
                            break;
                        y0 = y;
                        if (FindRoot(y) > 0)
                            break;
                        y = y0;
                    }
                }
            }
            else if (!intervention.empty() && intervention_t1 >= intervention_t0)
            {
                if (P.euler)    // testing...
                {
                    auto dydt = y;
                    int counter = 0;
                    for (double t = intervention_t0; t < intervention_t1; t += P.t_step, ++counter)
                    {
                        (*this)(y, dydt, t);
                        transform(y.begin(), y.end(), dydt.begin(), y.begin(), [&](double x, double dx) { return x + P.t_step * dx; });
                        if (counter % 100 == 0)
                            cout << t << " " << flush;
                    }
                    cout << "\n";
                }
                else
                    integrate(*this, y, intervention_t0, intervention_t1, P.t_step);
            }
            else if (P.euler) // integrate with above!
            {
                auto dydt = y;
                for (double t = intervention_t0; t < P.max_epochs * P.t_epoch; t += P.t_step)
                {
                    (*this)(y, dydt, t);
                    transform(y.begin(), y.end(), dydt.begin(), y.begin(), [&](double x, double dx) { return x + P.t_step * dx; });
                }
            }
            else for (epochs = 0; epochs < P.max_epochs; ++epochs)
            {
                integrate(*this, y, epochs * P.t_epoch, (epochs + 1) * P.t_epoch, P.t_step);
                double diff = AbsoluteDifference(y, y0);
                y0 = y;
                if (diff < P.epsilon && epochs >= P.min_epochs)
                    break;
            }

            // Save results of block run
            copy(y.begin(), y.end(), Y.begin() + g0 * Subpopulation::Elements());

            if (epochs >= P.max_epochs && P.verbose)
                cout << "max epochs reached.\n";

            most_epochs = max(most_epochs, epochs);
        }

        Outcome(Y, res, most_epochs);
    }

    void operator() (const vector<double>& Y, vector<double>& dYdt, const double t)
    {
        if (g1 == g0 + 1 && P.skip_contact)
        {
            Subpopulation::Gradient(g0, Y, dYdt, Y, t);
        }
        else
        {
            for (int c = 0; c < g1 - g0; ++c)
            {
                // Calculate total contact with individuals of each class for this group
                contact.assign(Subpopulation::Elements(), 0.0);
                for (int d = 0; d < g1 - g0; ++d)
                    for (int e = 0; e < Subpopulation::Elements(); ++e)
                        contact[e] += p[c + g0][d + g0] * Y[Subpopulation::Elements() * d + e];

                // Calculate gradients within subpopulation
                vector<double> subY(Y.begin() + c * Subpopulation::Elements(), Y.begin() + (c + 1) * Subpopulation::Elements());
                vector<double> subdYdt(dYdt.begin() + c * Subpopulation::Elements(), dYdt.begin() + (c + 1) * Subpopulation::Elements());
                Subpopulation::Gradient(c + g0, subY, subdYdt, contact, t);
                copy(subdYdt.begin(), subdYdt.end(), dYdt.begin() + c * Subpopulation::Elements());
            }
        }
    }

private:
    vector<double> contact;         // Assembled contact rates
    vector<vector<double>> p;       // Contact matrix: rate of group-i individuals contacting group-j individuals.
    int g0, g1;                     // Span of groups to integrate over in operator()

    void MakeContactMatrix()
    {
        if (!P.waifw.empty() && (int)P.waifw.size() != P.group_connectivity * P.group_connectivity)
            throw std::runtime_error("Waifw matrix size does not match connectivity of groups.");

        if (P.groups % P.group_connectivity != 0)
            throw std::runtime_error("Number of groups must be a multiple of group connectivity.");

        p = vector<vector<double>>(P.groups, vector<double>(P.groups, 0.0)); // contact matrix

        // Two-scale network (e.g. within and between country)
        if (P.group_span > 1 && V->G[0] >= 0)
        {
            for (g0 = 0; g0 < P.groups; g0 += P.group_connectivity)
            {
                g1 = g0 + P.group_connectivity;
                double SH = accumulate(V->h.begin() + g0, V->h.begin() + g1, 0.0); // total size of set of countries

                // G is "within country" assortativity, while g is "within province" assortativity.
                double H = 1, h = 1;
                for (int d = g0; d < g1; ++d)
                {
                    H = V->h[d] / SH;              // rel. size of contacted province to set of countries
                    h = V->h[d] / accumulate(V->h.begin() +  (d / P.group_span) * P.group_span,
                                             V->h.begin() + ((d / P.group_span) + 1) * P.group_span,
                                             0.0); // rel. size of contacted province to its country

                    for (int c = g0; c < g1; ++c)
                    {
                        // G is "within country" assortativity, while g is "within province" assortativity.
                        double G = V->G[c / P.group_span], g = V->g[c / P.group_span];

                        if (c == d) // same province
                            p[c][d] = G * (g + (1 - g) * h) + (1 - G) * H;
                        else if (c / P.group_span == d / P.group_span)  // provinces in same country
                            p[c][d] = G * (1 - g) * h + (1 - G) * H;
                        else // provinces in different countries
                            p[c][d] = (1 - G) * H;
                    }
                }
            }

        }
        // One-scale network or WAIFW matrix
        else
        {
            double H; // total subnetwork size
            auto waifw = [&](int c, int d, int g0) -> double { // if no waifw matrix is provided, random mixing is assumed
                if (P.waifw.size() > 0)
                    return P.waifw[(c - g0) * P.group_connectivity + d - g0];
                return V->h[d] / H;
            };

            for (g0 = 0; g0 < P.groups; g0 += P.group_connectivity)
            {
                g1 = g0 + P.group_connectivity;
                H = accumulate(V->h.begin() + g0, V->h.begin() + g1, 0.0); // total subnetwork size

                // The parameter g controls the deviation in assortativity.
                // This is a piecewise linear interpolation where -1 is random mixing, 0 is the WAIFW matrix, and 1 is assortative mixing.
                // This means that when no WAIFW matrix is provided, all g in [-1,0] are equivalent.
                for (int c = g0; c < g1; ++c)
                {
                    for (int d = g0; d < g1; ++d)
                    {
                        if (V->g[c] <= 0)
                            p[c][d] = (1 + V->g[c]) * waifw(c, d, g0) - V->g[c] * V->h[d] / H;
                        else
                            p[c][d] = (1 - V->g[c]) * waifw(c, d, g0) + V->g[c] * (c == d ? 1.0 : 0.0);
                    }
                }
            }
        }
    }

    static vector<double> InitialState()
    {
        vector<double> state;
        state.reserve(Subpopulation::Elements() * P.groups);
        vector<double> y0 = Subpopulation::InitialState();
        for (int g = 0; g < P.groups; ++g)
            state.insert(state.end(), y0.begin(), y0.end());
        return state;
    }

    void Outcome(const vector<double>& Y, Results& res, int epochs)
    {
        res.max_epochs = max(epochs, res.max_epochs);
        res.saved_Y = Y;

        double total_size = accumulate(V->h.begin(), V->h.end(), 0.0);

        // Calculate average properties and collate group contents
        for (int c = 0; c < P.groups; ++c)
        {
            // Add compartments, plus calculated fraction of disease which is resistant, to output vector
            double so, ro, b, r;
            vector<double> subY(Y.begin() + c * Subpopulation::Elements(), Y.begin() + (c + 1) * Subpopulation::Elements());
            auto group = Subpopulation::Report(subY, so, ro, b, r);

            res.S += so * V->h[c] / total_size;
            res.R += ro * V->h[c] / total_size;
            res.SR += b * V->h[c] / total_size;
            res.RTotal += r * V->h[c] / total_size;
            res.InfTotal += (so + ro + b) * V->h[c] / total_size;
            double RFrac = r / (so + ro + b);

            res.compartments.insert(res.compartments.end(), group.begin(), group.end());
            res.compartments.push_back(RFrac);

            res.Rs.push_back(r);
            res.Duals.push_back(b);
            res.Infs.push_back(so + ro + b);
        }
    }

    void Intervene(string intervention, double intervention_coverage, double intervention_strength)
    {
        if (intervention == "vaccine")
        {
            V->v = intervention_coverage;
            V->sv = 1. - intervention_strength;
            if (P.interv_import)
                for (auto& ps : V->psi)
                    ps *= max(0.0, 1. - P.interv_psi * intervention_strength);
        }
        else if (intervention == "clearance")
        {
            V->v = intervention_coverage;
            if (intervention_strength == 1)
                throw runtime_error("Cannot implement a full-strength clearance-hastening vaccine, because the clearance rate is undefined.");
            V->dv = 1. / (1. - intervention_strength);
            if (P.interv_import)
                for (auto& ps : V->psi)
                    ps *= max(0.0, 1. - P.interv_psi * intervention_strength);
        }
        else if (intervention == "treatment")
        {
            V->v = intervention_coverage;
            for (auto& tv : V->tv)
                tv = 1. - intervention_strength;
        }
        else if (intervention == "hightreatment")
        {
            V->v = intervention_coverage;
            // weighted average of treatment rates we want to reach in covered individuals
            double target = inner_product(V->h.begin(), V->h.end(), V->tau.begin(), 0.0) * (1 - intervention_strength);
            auto func = [&](double ceiling) { return inner_product(V->h.begin(), V->h.end(), V->tau.begin(), 0.0,
                                                                   [](double x, double y) { return x + y; },
                                                                   [&](double h, double tau) { return h * min(ceiling, tau); }) > target; };
            double ceiling = FindBoundary(func, *min_element(V->tau.begin(), V->tau.end()), *max_element(V->tau.begin(), V->tau.end()), 20, false, 1, false);
            for (int g = 0; g < P.groups; ++g)
                V->tv[g] = min(ceiling, V->tau[g]) / V->tau[g];
        }
        else if (intervention == "hightreatment_between")
        {
            V->v = intervention_coverage;
            // weighted average of treatment rates we want to reach in covered individuals
            double target = inner_product(V->h.begin(), V->h.end(), V->tau.begin(), 0.0) * (1 - intervention_strength);
            // maximum treatment rate in each country
            vector<double> maxes(P.groups / P.group_span, 0.0);
            for (int c = 0; c < (int)maxes.size(); ++c)
                maxes[c] = *max_element(V->tau.begin() + c * P.group_span, V->tau.begin() + (c + 1) * P.group_span);
            // find ceiling multiplier
            auto func = [&](double ceilingm) { double tau = 0.0;
                                                for (int g = 0; g < P.groups; ++g)
                                                    tau += V->h[g] * V->tau[g] * min(1.0, ceilingm / maxes[g / P.group_span]);
                                                return tau > target; };
            double ceilingm = FindBoundary(func, *min_element(V->tau.begin(), V->tau.end()), *max_element(V->tau.begin(), V->tau.end()), 20, false, 1, false);
            // apply intervention
            for (int g = 0; g < P.groups; ++g)
                V->tv[g] = min(1.0, ceilingm / maxes[g / P.group_span]);
        }
        else if (intervention == "hightreatment_within")
        {
            V->v = intervention_coverage;
            // weighted average of treatment rates we want to reach in covered individuals
            double target = inner_product(V->h.begin(), V->h.end(), V->tau.begin(), 0.0) * (1 - intervention_strength);
            // maximum treatment rate in each country
            vector<double> maxes(P.groups / P.group_span, 0.0);
            for (int c = 0; c < (int)maxes.size(); ++c)
                maxes[c] = *max_element(V->tau.begin() + c * P.group_span, V->tau.begin() + (c + 1) * P.group_span);
            // find ceiling multiplier
            auto func = [&](double ceilingm) { double tau = 0.0;
                                                for (int g = 0; g < P.groups; ++g)
                                                    tau += V->h[g] * min(ceilingm * maxes[g / P.group_span], V->tau[g]);
                                                return tau > target; };
            double ceilingm = FindBoundary(func, 0, 1, 20, false, 1, false);
            // apply intervention
            for (int g = 0; g < P.groups; ++g)
                V->tv[g] = min(ceilingm * maxes[g / P.group_span], V->tau[g]) / V->tau[g];
        }
        else if (intervention != "none" && intervention != "")
        {
            throw runtime_error("Unrecognised intervention type " + intervention);
        }
    }

    int FindRoot(vector<double>& Y)
    {
        // Set up root solver
        gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrid, Y.size());
        gsl_multiroot_function f = { &RootGradient, Y.size(), this };
        gsl_vector *x = gsl_vector_alloc(Y.size());
        for (size_t i = 0; i < Y.size(); ++i)
            gsl_vector_set(x, i, Y[i]);
        gsl_multiroot_fsolver_set(s, &f, x);

        // Iterate
        int status = 0, iter = 0;
        do
        {
            ++iter;
            status = gsl_multiroot_fsolver_iterate(s);

            if (status)   // check if solver is stuck
                break;

            status = gsl_multiroot_test_residual(s->f, P.epsilon);
        } while (status == GSL_CONTINUE && iter < P.max_epochs);

        // Get results
        bool negative = false, exceeds_simplex = false;
        double running_sum = 0.0;
        size_t gs = P.vaccinate ? Subpopulation::Elements() / 2 : Subpopulation::Elements();
        for (size_t i = 0; i < Y.size(); ++i)
        {
            Y[i] = gsl_vector_get(s->x, i);

            // Clean up small negative values
            if (-P.epsilon < Y[i] && Y[i] < 0)
                Y[i] = 0;

            // Zero compartments where there are no hosts
            if (P.vaccinate)
                if ((V->v == 1 && (i / gs) % 2 == 0)
                 || (V->v == 0 && (i / gs) % 2 == 1))
                    Y[i] = 0;

            // Check for anomalies
            if (Y[i] < 0)
                negative = true;
            running_sum += max(0.0, Y[i]);
            if (i % gs == gs - 1)
            {
                double limit = P.vaccinate ? ((i / gs) % 2 ? V->v : 1 - V->v) : 1.0;
                if (running_sum > limit)
                    exceeds_simplex = true;
                running_sum = 0;
            }
        }

        // Check validity of results
        if (status) // GSL could not find the root
        {
            cout << "*";
            iter = -1;
        }
        else if (negative) // A root component is negative
        {
            cout << "-";
            iter = -2;
        }
        else if (exceeds_simplex) // A population group has more than 100% of individuals as carriers
        {
            cout << "!";
            iter = -3;
        }

        // Clean up and exit
        gsl_multiroot_fsolver_free(s);
        gsl_vector_free(x);
        return iter;
    }

    static int RootGradient(const gsl_vector* x, void* params, gsl_vector* f)
    {
        auto parent = static_cast<Metapopulation<Subpopulation>*>(params);
        vector<double> Y(x->size, 0.0), dYdt(x->size, 0.0);

        // Get vector and calculate gradient
        for (size_t i = 0; i < x->size; ++i)
            Y[i] = gsl_vector_get(x, i);
        (*parent)(Y, dYdt, 0.);

        for (size_t i = 0; i < x->size; ++i)
            gsl_vector_set(f, i, dYdt[i]);

        return GSL_SUCCESS;
    }
};
