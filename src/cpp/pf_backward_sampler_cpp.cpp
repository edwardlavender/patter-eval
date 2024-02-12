#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector dtruncgamma_cpp(NumericVector x, double shape = 15, double scale = 15, double mobility = 500) {
  int n = x.size();
  NumericVector tt(n);

  for (int i = 0; i < n; ++i) {
    if (x[i] <= mobility) {
      tt[i] = R::dgamma(x[i], shape, scale, 0) / (R::pgamma(mobility, shape, scale, 1, 0));
    } else {
      x[i] = 0;
    }
  }

  return tt;
}

/*** R
x <- c(0, 50, 100, 400, 500, 1000)
(a <- dtruncgamma(x))
(b <- dtruncgamma_cpp(x))
all.equal(a, b)
*/

//' @title Density function
//' @description Evaluate f(s_{t} | s_{t - 1})
//' @param ppast A multi-row matrix of x,y coordinates for previous particles
//' @param pnow A one-row matrix of coordinates for the current particle
// [[Rcpp::export]]

NumericVector dpropose_cpp(NumericMatrix ppast,
                           NumericVector pnow,
                           double shape,
                           double scale,
                           double mobility) {

  // Calculate distances between current particle and all previous particles
  NumericVector x_diff = ppast(_, 0) - pnow[0];
  NumericVector y_diff = ppast(_, 1) - pnow[1];
  NumericVector distances = sqrt(pow(x_diff, 2) + pow(y_diff, 2));

  // Compute densities
  NumericVector dens = dtruncgamma_cpp(distances, shape, scale, mobility);

  return dens;
}

//' @title Sampling function
//' @description Sample s_{t-1} using f(s_{t} | s_{t - 1})
//' @param ppast A matrix of x, y coordinates for previous particles
//' @param dens A vector of movement densities into ppast locations
// [[Rcpp::export]]

NumericVector sample_cpp(NumericMatrix ppast, NumericVector dens) {
  int np = ppast.nrow();

  // Validation
  double dens_sum = sum(dens);
  if (dens_sum == 0) {
    stop("All densities are zero.");
  }

  // Normalise densities
  NumericVector dens_norm = dens / dens_sum;

  // Sample
  int index = 0;
  double u = R::runif(0, 1);
  double dens_cum = dens_norm[0];
  while (u > dens_cum && index < np - 1) {
    index++;
    dens_cum += dens_norm[index];
  }

  return ppast(index, _);
}

//' @title Backward sampling algorithm
//' @value The function returns a numeric matrix of matrix, with one row for each path and one column for each time step.
// [[Rcpp::export]]

NumericMatrix pf_backward_sampler_cpp(List particles, double shape, double scale, double mobility) {

  // Identify final particle samples (coordinates)
  // A. Define number time steps
  int nt = particles.size();
  // B. Identify particles for final time step
  NumericMatrix pfinal = particles[nt - 1];
  // C. Define number of finishing particles
  int np = pfinal.nrow();

  // Define output matrix of paths
  // * One row for each path
  // * Two columns per time step (for the x and y coordinates)
  int nc = nt * 2;
  NumericMatrix paths(np, nc);

  // Loop over each particle
  for (int i = 0; i < np; ++i) {

    // Initiate loop with final coordinates
    paths(i, nc - 2) = pfinal(i, 0);
    paths(i, nc - 1) = pfinal(i, 1);

    // Loop backwards over time steps
    for (int j = nt; j >= 2; --j) {

      Rcpp::Rcout << "Particle " << i + 1 << " (" << "timestep " << j << ")" << std::endl;

      // Define matrix indices
      int iy_now = j * 2 - 1;
      int ix_now = j * 2 - 2;
      int iy_past = j * 2 - 3;
      int ix_past = j * 2 - 4;

      // Identify current particle coordinates
      NumericVector pnow = NumericVector::create(paths(i, ix_now), paths(i, iy_now));
      Rcpp::Rcout << pnow << std::endl;

      // Identify previous particles coordinates
      // j is time step of current element; j - 1 is corresponding index; previous particles are in j - 2
      NumericMatrix ppast = particles[j - 2];

      // Calculate movement densities
      NumericVector dens = dpropose_cpp(ppast, pnow, shape, scale, mobility);

      // Sample a previous particle & record
      NumericVector pxy = sample_cpp(ppast, dens);
      paths(i, ix_past) = pxy[0];
      paths(i, iy_past) = pxy[1];
    }
  }

  return paths;
}
