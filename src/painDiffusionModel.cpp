// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

namespace pgu{
  arma::mat init_kernel()
  {
    int width = 3;
    arma::mat kernel(width, width);
    for(int x = 0; x < width; x++) {
      for(int y = 0; y < width; y++){
        kernel(x,y) = 1.0;
      }
    }
    kernel(0,0) = 0.0;
    kernel(0,2) = 0.0;
    kernel(2,2) = 0.0;
    kernel(2,0) = 0.0;
    kernel(1,1) = -4.0;
    return kernel;
  }
  
  int toroidal_edge_correction(int max_idx, int idx)
  {
    if (idx<0)
      return idx + max_idx;
    if(idx >= max_idx)
      return idx - max_idx;
    return idx;
  }
  
  arma::mat convolution(arma::mat matrix, arma::mat kernel)
  {
    double integrated_intensity;
    int x_corr, y_corr;
    int m = matrix.n_rows;
    int n = matrix.n_cols;
    arma::mat matrix_new(m, n);
    for(int x = 0; x < m; x++) {
      for(int y = 0; y < n; y++){
        integrated_intensity = 0.0;
        for(int i = -1; i <= 1; i++){
          for(int j = -1; j <= 1; j++){
            x_corr = pgu::toroidal_edge_correction(m, x - i);
            y_corr = pgu::toroidal_edge_correction(n, y - j);
            integrated_intensity += kernel(i+1, j+1) * matrix(x_corr, y_corr);
          }
        }
        matrix_new(x,y) = integrated_intensity;
      }
    }
    return matrix_new;
  }
  
  double decay_reaction(double x, double k)
  {
    if (x>0)
    {
      return k * x;
    }
    return 0;
  }
  
  int heaviside(double x, double omega)
  {
    if(x<omega)
    {
      return(0);
    }
    return(1);
  }
  
  double hill(double x, double ec50, double n)
  {
    return(1.0 / (1.0 + std::pow(ec50/x, n)));
  }
  
  arma::cube update_system(arma::cube compartments, arma::mat kernel, double q, double ec50, double n, double kin,  double diff_coeff, double tau, double omega, double delta_t)
  {
    int u = compartments.n_rows;
    int v = compartments.n_cols;
    double delta_clearance;
    double delta_effect;
    double delta_decay;
    double delta_diffusion;
    arma::cube compartments_new(u, v, 4);
    arma::mat convolved_compartment = pgu::convolution(compartments.slice(2), kernel);
    for(int x = 0; x < u; x++) {
      for(int y = 0; y < v; y++){
        delta_clearance = pgu::decay_reaction(compartments(x, y, 0), q);
        delta_effect = pgu::decay_reaction(pgu::hill(compartments(x, y, 0), ec50, n), kin);
        delta_decay = pgu::decay_reaction(compartments(x, y, 2), 1.0/tau);
        delta_diffusion = diff_coeff * convolved_compartment(x,y);
        if(compartments(x, y, 2) <= 0)
        {
          delta_decay = 0;
        }
        compartments_new(x, y, 0) = compartments(x, y, 0) - delta_t * delta_clearance;
        compartments_new(x, y, 1) = compartments(x, y, 1) + delta_t * delta_clearance;
        compartments_new(x, y, 2) = compartments(x, y, 2) + delta_t * (delta_diffusion + delta_effect - delta_decay);
        compartments_new(x, y, 3) = pgu::heaviside(compartments_new(x, y, 2), omega);
      }
    }
    return compartments_new;
  }
}

//' @title simulate_system
//' @description Simulate the time evaluation of the diffusion of intracutaneously applied capsaicin-evoked pain mediated by the spinal cord.
//' @details The model is build from 4 compartments. 1 The peripheral compartment; 2 The central compartment; 3 The transit compartment; 4 The effect compartment
//' @param compartments Tow-dimensional representation of model compartments.
//' @param q Clearance of the capsaicin conentration in the peripheral compartment of the subcutaneous fat tissue.
//' @param ec50 Effect activity triggered by capsaicin concentration of the peripheral compartment in the transit compartment.
//' @param n The Hill coefficient as a measure of ultrasensitivity.
//' @param kin Apparent effect production rate in the transit compartment.
//' @param diff_coeff Lateral diffusion coefficient for the pain stimulus diffusion within the transit compartment.
//' @param tau Decay time that describes the decay of the effect in the transit compartment.
//' @param omega Threshold of severity perception of a pain stimulus.
//' @param t Duration of the simulation.
//' @param inc Number of simulated increments.
//' @return Tow-dimensional representation of model compartments. 4 Compartments.
//' @export
// [[Rcpp::export]]
arma::cube simulate_system(arma::cube compartments, double q, double ec50, double n, double kin, double diff_coeff, double tau, double omega, double t, double inc)
{
  double delta_t = t/inc;
  arma::mat kernel = pgu::init_kernel();
  for(int i = 0; i < inc; i++){
    compartments = pgu::update_system(compartments, kernel, q, ec50, n, kin, diff_coeff, tau, omega, delta_t); 
  }
  return compartments;
}
