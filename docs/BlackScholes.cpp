#include <Rcpp.h>
#include <math.h>
#include <vector>

using namespace Rcpp;

double normCDF(double x) {
    return erfc( -x * M_SQRT1_2 ) / 2;
}

// [[Rcpp::export]]
std::vector<double> BlackScholes_Cpp(
        const double S0, 
        const double K, 
        const double t, 
        const double r, 
        const double v ) {
    
    const double d1 = (log(S0/K) + (r + v*v/2) * t) / (v * sqrt(t));
    const double d2 = d1 - v * sqrt(t);
    
    const double C = S0 * normCDF(d1) - K * exp( -r*t ) * normCDF(d2);
    const double P = -S0 * normCDF(-d1) + K * exp( -r*t ) * normCDF(-d2);
    
    // std::cout << "Call Price: \t" << C << std::endl;
    // std::cout << "Put Price: \t" << P << std::endl;
    
    std::vector<double> prices{ C, P };
    
    return prices;
}

// END OF FILE
