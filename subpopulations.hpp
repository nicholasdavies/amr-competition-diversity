// subpopulations.hpp
// dynamics of subpopulations within a metapopulation.
// In all these, a model component in [brackets] indicates assumed based on a density constraint

void VaxAdjust(vector<double>& y0, bool guess)
{
    auto size = y0.size();
    if (P.vaccinate)
    {
        y0.reserve(size * 2);
        y0.insert(y0.end(), y0.begin(), y0.end());
        transform(y0.begin(), y0.begin() + size, y0.begin(), [](double y) { return y * (1 - V->v); });
        transform(y0.begin() + size, y0.end(), y0.begin() + size, [=](double y) { return y * (guess ? V->v * V->sv / V->dv : V->v); });
    }
}

// Single-infection
// [X], IS, IR
class PopulationSingle
{
public:
    static int Elements() { return 2 * (P.vaccinate ? 2 : 1); }
    enum Components { I_S, I_R, I_Sv, I_Rv };

    static vector<double> InitialState()
    {
        vector<double> y0(2, 0.0);

        y0[0] = P.init_freq_s;
        y0[1] = P.init_freq_r;
        VaxAdjust(y0, false);

        return y0;
    }

    static void Gradient(const int g, const vector<double>& Y, vector<double>& dYdt, const vector<double>& contact, const double t)
    {
        (void) t;
        double S_tot = max(0.0, contact[I_S] + (P.vaccinate ? contact[I_Sv] : 0.0));
        double R_tot = max(0.0, contact[I_R] + (P.vaccinate ? contact[I_Rv] : 0.0));
        double lambda_S = V->beta[g]                 * S_tot + V->psi[0] * (1 - V->psi_r[0]);
        double lambda_R = V->beta[g] * (1 - V->c[g]) * R_tot + V->psi[0] *      V->psi_r[0] ;

        Gradient2(I_S, Y, dYdt, 1 - V->v, lambda_S, lambda_R, V->tau[g], V->u[g], V->kk[g]);
        if (P.vaccinate)
            Gradient2(I_Sv, Y, dYdt, V->v, lambda_S * V->sv, lambda_R * V->sv, V->tau[g] * V->tv[g], V->u[g] * V->dv, V->kk[g]);
    }

    static void Gradient2(int j, const vector<double>& Y, vector<double>& dYdt, double size, double lambda_S, double lambda_R, double tau, double u, double kk)
    {
        double S = max(0.0, Y[j + I_S]);
        double R = max(0.0, Y[j + I_R]);
        double X = max(0.0, size - S - R);

        dYdt[j + I_S] = lambda_S * X        // infection of susceptibles with S
                    - (u + tau) * S         // reversion to susceptibility because of treatment or natural clearance
                    + kk * lambda_S * R;    // superinfection of R hosts by sensitive strain

        dYdt[j + I_R] = lambda_R * X        // infection of susceptibles with R
                    - u * R                 // reversion to susceptibility because of natural clearance
                    - kk * lambda_S * R;    // superinfection of R hosts by sensitive strain
    }

    static vector<double> Report(const vector<double>& Y, double& S_only, double& R_only, double& both, double& R_total)
    {
        S_only = Y[I_S] + (P.vaccinate ? Y[I_Sv] : 0.0);
        R_only = Y[I_R] + (P.vaccinate ? Y[I_Rv] : 0.0);
        both = 0;
        R_total = R_only;
        return { 1 - S_only - R_only, S_only, R_only };
    }

    static vector<string> Compartments()
    {
        return { "X", "S", "R" };
    }
};

// Within-Host
// [X], S, Sr, ..., Rs, R
class PopulationWithinHost
{
public:
    static int Elements() { return (2 + P.whc_steps) * (P.vaccinate ? 2 : 1); }
    static vector<vector<double>> F;

    // Initialise F vector, if needed
    static void Setup()
    {
        if (P.whc_steps < 2)
            throw runtime_error("whc_steps must be 2 or greater.");

        F.resize(P.whc_steps + 1);

        double a = log(P.whc_iota) - log(1 - P.whc_iota);
        double as = 2 * a / (P.whc_steps - 1);

        // Fill F vector with 1, 0.999, ..., 0.001, 0
        F[P.whc_steps].push_back(1.0);
        for (int i = 0; i < P.whc_steps; ++i, a -= as)
            F[P.whc_steps].push_back(1. / (1. + exp(a)));
        F[P.whc_steps].push_back(0.0);

        if (P.whc_neutral)
        {
            F[P.whc_steps][1] = 1.0;
            F[P.whc_steps][P.whc_steps] = 0.0;
        }
    }

    static vector<double> InitialState()
    {
        vector<double> y0(2 + P.whc_steps, 0.0);

        y0[0] = P.init_freq_s;
        y0[1 + P.whc_steps] = P.init_freq_r;
        VaxAdjust(y0, false);

        return y0;
    }

    static void Gradient(const int g, const vector<double>& Y, vector<double>& dYdt, const vector<double>& contact, const double t)
    {
        (void) t;
        double lambda_S = V->psi[0] * (1 - V->psi_r[0]);
        double lambda_R = V->psi[0] * V->psi_r[0];
        double duals = 0, dualsv = 0;
        const int vj = 2 + P.whc_steps;

        for (int i = 0; i < vj; ++i)
        {
            lambda_S += V->beta[g]                 *      F[P.whc_steps][i]  * (contact[i] + (P.vaccinate ? contact[vj + i] : 0.0));
            lambda_R += V->beta[g] * (1 - V->c[g]) * (1 - F[P.whc_steps][i]) * (contact[i] + (P.vaccinate ? contact[vj + i] : 0.0));
            if (i > 0 && i <= P.whc_steps)
            {
                duals += Y[i];
                if (P.vaccinate)
                    dualsv += Y[vj + i];
            }
        }

        Gradient2(g, 0, Y, dYdt, (1 - V->v) - Y[0] - duals - Y[1 + P.whc_steps], duals, lambda_S, lambda_R, V->tau[g], V->u[g]);
        if (P.vaccinate)
            Gradient2(g, vj, Y, dYdt, V->v - Y[vj] - dualsv - Y[vj + 1 + P.whc_steps], dualsv, lambda_S * V->sv, lambda_R * V->sv, V->tau[g] * V->tv[g], V->u[g] * V->dv);
    }

    static void Gradient2(const int g, const int j, const vector<double>& Y, vector<double>& dYdt, double X, double duals, double lambda_S, double lambda_R, double tau, double u)
    {
        const int S = j, Sr = j + 1, Rs = j + P.whc_steps, R = j + 1 + P.whc_steps;

        dYdt[S] = lambda_S * X
            - (u + tau) * Y[S]
            - V->k[g] * lambda_R * Y[S]
            + V->b0[g] * V->b[g] * Y[Sr];
        dYdt[Sr] = V->k[g] * lambda_R * Y[S]
            - (u + tau) * Y[Sr]
            + V->b[g] * Y[Sr + 1]
            - V->b0[g] * V->b[g] * Y[Sr];

        for (int d = Sr + 1; d < Rs; ++d)
            dYdt[d] = - (u + tau) * Y[d]
                + V->b[g] * (Y[d + 1] - Y[d]);

        dYdt[Rs] = V->k[g] * lambda_S * Y[R]
            - (u + tau) * Y[Rs]
            - V->b[g] * Y[Rs];
        dYdt[R] = lambda_R * X
            - u * Y[R]
            + tau * duals
            - V->k[g] * lambda_S * Y[R];
    }

    static vector<double> Report(const vector<double>& Y, double& S_only, double& R_only, double& both, double& R_total)
    {
        const int S = 0, R = 1 + P.whc_steps;
        const int vj = 2 + P.whc_steps;

        S_only = Y[S] + (P.vaccinate ? Y[vj + S] : 0.0);
        R_only = Y[R] + (P.vaccinate ? Y[vj + R] : 0.0);
        both = 0;
        R_total = 0;

        vector<double> ret(3 + P.whc_steps, 0.0);
        for (int i = 0; i < 2 + P.whc_steps; ++i)
        {
            double yi = Y[i] + (P.vaccinate ? Y[vj + i] : 0.0);
            R_total += (1 - F[P.whc_steps][i]) * yi;
            if (i > 0 && i <= P.whc_steps)
                both += yi;
            ret[i + 1] = yi;
        }
        ret[0] = 1 - S_only - R_only - both;

        return ret;
    }

    static vector<string> Compartments()
    {
        int sub = P.whc_steps / 2;
        int mid = P.whc_steps % 2;

        vector<string> comp = { "X" };
        for (int i = 0; i <= sub; ++i)
            comp.push_back("S" + string(i, 'r'));
        if (mid)
            comp.push_back("SR");
        for (int i = 0; i <= sub; ++i)
            comp.push_back("R" + string(sub - i, 's'));
        return comp;
    }
};

vector<vector<double>> PopulationWithinHost::F;

// D-type
// [X], IS0, IR0, IS1, IR1, ..., ISN-1, IRN-1
class PopulationDType
{
public:
    static int Elements() { return (2 * P.D_n) * (P.vaccinate ? 2 : 1); }

    static vector<double> InitialState()
    {
        vector<double> y0;
        for (int d = 0; d < P.D_n; ++d)
            y0.insert(y0.end(), { P.init_freq_s / P.D_n, P.init_freq_r / P.D_n });
        VaxAdjust(y0, false);

        return y0;
    }

    static void Gradient(const int g, const vector<double>& Y, vector<double>& dYdt, const vector<double>& contact, const double t)
    {
        (void) t;
        const int vj = 2 * P.D_n;
        double Ysumn = accumulate(Y.begin(), Y.begin() + vj, 0.), Ysumv = 0;
        if (P.vaccinate) Ysumv = accumulate(Y.begin() + vj, Y.end(), 0.);
        double X = (1 - V->v) - Ysumn, Xv = V->v - Ysumv;

        for (int d = 0; d < P.D_n; ++d)
        {
            double yd = Y[2 * d + 0] + Y[2 * d + 1] + (P.vaccinate ? Y[vj + 2 * d + 0] + Y[vj + 2 * d + 1] : 0.0);
            double b = pow(1. - yd / (Ysumn + Ysumv) + 1. / P.D_n, V->a[g]); // balancing selection on d-type

            double S_tot = max(0.0, contact[2 * d + 0] + (P.vaccinate ? contact[vj + 2 * d + 0] : 0.0));
            double R_tot = max(0.0, contact[2 * d + 1] + (P.vaccinate ? contact[vj + 2 * d + 1] : 0.0));
            double lambda_S = V->beta[g]                 * S_tot + V->psi[d % V->psi.size()] * (1 - V->psi_r[d % V->psi_r.size()]) / P.D_n;
            double lambda_R = V->beta[g] * (1 - V->c[g]) * R_tot + V->psi[d % V->psi.size()] *      V->psi_r[d % V->psi_r.size()]  / P.D_n;
            double u = V->u[g] * (1 + V->delta[g] * (2 * d / (P.D_n - 1.) - 1));

            dYdt[2 * d + 0] = lambda_S * b * X - (u + V->tau[g]) * Y[2 * d + 0];
            dYdt[2 * d + 1] = lambda_R * b * X -               u * Y[2 * d + 1];

            if (P.vaccinate)
            {
                dYdt[vj + 2 * d + 0] = lambda_S * V->sv * b * Xv - (u * V->dv + V->tau[g] * V->tv[g]) * Y[vj + 2 * d + 0];
                dYdt[vj + 2 * d + 1] = lambda_R * V->sv * b * Xv -  u * V->dv                         * Y[vj + 2 * d + 1];
            }
        }
    }

    static vector<double> Report(const vector<double>& Y, double& S_only, double& R_only, double& both, double& R_total)
    {
        const int vj = 2 * P.D_n;
        S_only = 0; R_only = 0; both = 0;
        for (int d = 0; d < P.D_n; ++d)
        {
            S_only += Y[2 * d + 0] + (P.vaccinate ? Y[vj + 2 * d + 0] : 0.0);
            R_only += Y[2 * d + 1] + (P.vaccinate ? Y[vj + 2 * d + 1] : 0.0);
        }
        R_total = R_only;

        if (P.D_report_all)
        {
            vector<double> ret = { 1 - S_only - R_only };
            for (int d = 0; d < P.D_n; ++d)
            {
                ret.push_back(Y[2 * d + 0] + (P.vaccinate ? Y[vj + 2 * d + 0] : 0.0));
                ret.push_back(Y[2 * d + 1] + (P.vaccinate ? Y[vj + 2 * d + 1] : 0.0));
            }
            return ret;
        }

        return { 1 - S_only - R_only, S_only, R_only };
    }

    static vector<string> Compartments()
    {
        if (P.D_report_all)
        {
            vector<string> comp = { "X" };
            for (int d = 0; d < P.D_n; ++d)
            {
                comp.push_back("S" + to_string(d) + "d");
                comp.push_back("R" + to_string(d) + "d");
            }
            return comp;
        }

        return { "X", "S", "R" };
    }
};
