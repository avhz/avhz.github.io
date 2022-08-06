
// Header files and namespaces -------------------------------------------------

#include <Rcpp.h>
#include <math.h>
#include <numeric>
#include <random>
#include <vector>
#include <functional> 

using namespace Rcpp;

// Non-exported functions ------------------------------------------------------

std::vector<double> linspace(double a, double b, int num) {
    std::vector<double> v(num);
    for (int i = 0; i < num; i++) { v[i] = a + i * ( (b - a) / num ); }
    return v;
}

std::vector<double> cumulative_sum(std::vector<double> x){
    std::vector<double> v(x.size());
    std::partial_sum(x.begin(), x.end(), v.begin(), std::plus<double>());
    return v;
}

std::vector<double> standard_normal(int n) {
    // random device class instance, source of 'true' randomness for initializing random seed
    std::random_device rd;
    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 gen(rd());

    float sample;
    std::vector<double> v(n);
    for(int i = 0; i < n; i++) {
        // instance of class std::normal_distribution with specific mean and stddev
        std::normal_distribution<float> d(0, 1);
        // get random number with normal distribution using gen as random source
        sample = d(gen);
        // add to vector
        v[i] = sample;
    }
    return v;
}

// Exported functions ----------------------------------------------------------

// [[Rcpp::export]]
std::vector<double> GBM_Cpp(
    double tau,         // time to expiry
    double r,           // risk free rate
    double sigma,       // volatility
    double S0,          // initial stock price
    int N               // number of steps between 0 and T
) {
    
    // length of each time step
    double dt = tau / N;
    
    // Vector for GBM paths
    std::vector<double> St(N+1);
    St[0] = S0;
    // fill(St.begin(), St.end(), S0);
    
    // vector of time points
    std::vector<double> time = linspace(0.0, tau, N);
    
    // standard normal sample of N elements
    std::vector<double> Z = standard_normal(N);
    
    // Brownian motion increments
    std::vector<double> dW(Z.size());

    for (int i = 0; i < (int)Z.size(); i++) {
        dW[i] = Z[i] * std::sqrt(dt);
    }
    
    // Brownian motion at each time (N+1 elements)
    std::vector<double> W{0};     
    std::vector<double> cumsum_dW = cumulative_sum(dW);
    W.insert(W.end(), cumsum_dW.begin(), cumsum_dW.end());
    
    for (int i = 1; i < (int)St.size(); i++) {
        St[i] = S0 * exp((r - sigma * sigma / 2) * time[i] + sigma * W[i]);
    }
    
    St.pop_back();

    return( St );
}

// END OF FILE.