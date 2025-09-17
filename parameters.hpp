#ifndef INCLUDED_parameters_hpp_
#define INCLUDED_parameters_hpp_

// Parameters of STOLAS
const double sigma = pow(2.,-3.); // ksigma = 2pi sigma exp(N) / L, nsigma = sigma exp(N)
const double dn = 1; // Thickness of nsigma sphere shell
const int NLnoise = pow(2,6); // Box size L
const int totalnoiseNo = pow(2,4); // The number of chunks
const int NL = pow(NLnoise,3)/totalnoiseNo; // Box size L for each noisemap
// const int divnumber = pow(2,5); // The number of chunks
const double dN = 0.01/log2(exp(1)); // e-folds step
const double Nprec = 1e-7; // Precision of e-foldings
const double dlogn = 0.1; // Width of bin in power spectrum

int internumber = 0;
double nsigmareset = 8.;
int totalstep = ceil(log((NLnoise/2-1)/sigma)/dN); // Total number of time step with noise
int firststep = ceil(log(nsigmareset/sigma)/dN);
int itpstep = ceil((log(nsigmareset/sigma)-log(nsigmareset/sigma/2.))/dN);
int aninum = 10; // dividing number for animation
const bool EoI_noise = false; // Adding the noise until EoN
int MeanNumber = 20;

// Outputs
const bool szeta = true; // Output the zeta
const bool sfield = false; // Output the animation
const bool strajectory = false; // Output the trajectory
const bool spower = true; // Output the power spectrum
const bool sweight = false; // Output the weight
const bool scompaction = false; // Output the compaction
bool sanimation = false; // Output the animation

// Directory name of saved data, you can change after "make clean" in your terminal.
const std::string sdatadir = "data";
const std::string sourcedir = "../source";
const std::string noisefiledir = "noisedata/noiselist_";
const std::string noisefilenamediv = "/part_";

#endif
