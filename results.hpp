// results.hpp

int NVaryingParams(std::string& P_param)
{
    int n = 0;
    vector<string> parts = Unserialize(P_param, ";,");
    for (auto& p : parts)
        if (p.size() > 0 && isalpha(p[0])) ++n;
    return n;
};

void Header(ostream& out, bool write_params, bool condensed)
{
    if (write_params)
        P.Write(out, "# ");

    auto param_header = [&](std::string& P_param, std::string name, bool condensed)
    {
        int n = NVaryingParams(P_param);
        if (n == 1 || !condensed)
            out << "\t" << name;
        else if (n > 1)
            for (int i = 0; i < n; ++i)
                out << "\t" << name << i;
    };

    out << "trial";
    if (P.mcmc)
        out << "\tlp\tchain\tll";

    #define ODE_PARAM(type, name) param_header(P.name, #name, condensed);
    DO_ODE_PARAMS;
    #undef ODE_PARAM

    out << "\tS\tR\tSR\tRfrac\tepochs";

    vector<string> compartments;
    switch (P.model)
    {
        case 0: break;
        case 1: compartments = PopulationSingle::Compartments(); break;
        case 2: compartments = PopulationWithinHost::Compartments(); break;
        case 3: compartments = PopulationDType::Compartments(); break;

        default:
            throw runtime_error("Unsupported model.");
            break;
    }
    compartments.push_back("Rfrac");

    for (int g = 0; g < P.groups; ++g)
        for (unsigned int c = 0; c < compartments.size(); ++c)
            out << "\t" << compartments[c] << g;

    out << "\n";
}


void ParamReport(ostream& out, string& P_params, double v_param, bool condensed)
{
    int n = NVaryingParams(P_params);

    if (n == 1 || !condensed)
        out << "\t" << v_param;
}

void ParamReport(ostream& out, string& P_params, vector<double>& v_param, bool condensed)
{
    if (!condensed)
    {
        out << "\t";
        for (unsigned int p = 0; p < v_param.size(); ++p)
            out << v_param[p] << (p == v_param.size() - 1 ? "" : ";");
    }
    else
    {
        vector<string> parts = Unserialize(P_params, ";,");
        for (unsigned int p = 0; p < parts.size(); ++p)
            if (parts[p].size() > 0 && isalpha(parts[p][0]))
                out << "\t" << v_param[p];
    }
}

void ReportResults(ostream& out, Results& r, bool condensed, int trial, double lp = 1, int chain = -1, double ll = 1, vector<double> special_params = vector<double>(), string trial_name = "")
{
    if (!trial_name.empty())
        out << trial_name;
    else
        out << trial;
    if (P.mcmc)
        out << "\t" << lp << "\t" << chain << "\t" << ll;

    if (special_params.empty() || !condensed)
    {
        #define ODE_PARAM(type, name) ParamReport(out, P.name, V->name, condensed);
        DO_ODE_PARAMS;
        #undef ODE_PARAM
    }
    else for (double& p : special_params)
        out << "\t" << p;

    out << "\t" << r.S << "\t" << r.R << "\t" << r.SR << "\t" << r.RFrac << "\t" << r.max_epochs;
    for (double c : r.compartments)
        out << "\t" << c;
    out << "\n";
}