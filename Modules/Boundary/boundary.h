// boundary.h

#include <stdexcept>
#include <string>

template <typename IndicatorFunc>
double FindBoundary(IndicatorFunc f, double x0, double x1, unsigned int n, bool check_bounds = true, int sign = 0, bool verbose = false)
{
    using namespace std;

    auto message = [verbose](string message) { if (verbose) cout << message << flush; };

    if (!check_bounds && sign == 0)
        throw runtime_error("FindBoundary: cannot guess sign of indicator function without checking bounds.");

    bool b = sign > 0;
    message("Boundary: ");
    if (check_bounds)
    {
        message("[bounds");
        bool f0 = f(x0);
        bool f1 = f(x1);
        message("] ");
        if (f0 == f1)
            throw runtime_error(string("FindBoundary: indicator function evaluates to ") + (f0 ? "true" : "false") + string(" at both ends of range."));
        if (sign == 0)
            b = f1;
        else if (b != f1)
            throw runtime_error("FindBoundary: indicator function inconsistent with sign given.");
    }

    message("[bifurcations");
    double xp = (x0 + x1) / 2;
    for (unsigned int i = 0; i < n; ++i)
    {
        if (f(xp) == b)
            x1 = xp;
        else
            x0 = xp;
        xp = (x0 + x1) / 2;
        message(".");
    }
    message("] ");
    message(" boundary identified at " + to_string(xp) + ".\n");

    return xp;
}