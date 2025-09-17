#ifndef INCLUDED_model_hpp_
#define INCLUDED_model_hpp_

#define NFIELDS 1 // Number of waterfall
#define NFLOOP for (int nf = 0; nf < NFIELDS; nf++)

extern inline double pw2(double x) {return x*x;}
extern inline double pw3(double x) {return x*x*x;}
extern inline double pw4(double x) {return x*x*x*x;}

// type
typedef std::vector<double> state_type;

// Model parameters
const std::string model = "test"; // Name of the model

const double mm = 1.e-5;

// Initial value of fields
const double PSI_INIT = 13.0;;
const double DPSI_INIT = 0.;

const double EoNvalue = 1.5;

// Compute psi_r
double Psir(const state_type PHI) {
  double wfterm = 0.;
  NFLOOP wfterm += pw2(PHI[2*nf]);
  return wfterm;
}

// Potential
double VV(const state_type &phi) {
  double wfterm = 0;
  NFLOOP{
    wfterm += pw2(phi[2*nf]);
  }
  return 0.5*pw2(mm)*wfterm;
}

// Derivative
double Vpsi(const state_type &phi, int NField) {
  double wfterm = phi[2*NField];
	return pw2(mm)*wfterm;
}

double hubble(const state_type &phi) {
  double wfkterm = 0;
  NFLOOP{
    wfkterm += pw2(phi[2*nf+1])/2.;
  }
  return sqrt((phi[1]*phi[1]/2. + wfkterm + VV(phi))/3.);
}

double ep(const state_type &phi) {
  double HH = hubble(phi);
  double wfkterm = 0;
  NFLOOP{
    wfkterm += pw2(phi[2*nf+1])/2.;
  }
  return wfkterm/HH/HH;
}

// pow spectrum of fields
double calPphi(const state_type &phi) {
  double HH = hubble(phi);
  return pw2(HH/2./M_PI);
}

double calPpsi(const state_type &phi) {
  double HH = hubble(phi);
  return pw2(HH/2./M_PI);
}

// The condition at the end of inflation
double EoI(const state_type &phi) {
  return 0.8 - ep(phi);
}

double EoN(const state_type &phi) {
  return EoNvalue + pw2(mm)/VV(phi);
}

// EoM
void dphidN(const state_type &phi, state_type &dxdt, const double t) {
  double HH = hubble(phi);

  NFLOOP{
    double p2 = phi[2*nf+1]; // pi
    dxdt[2*nf] = p2/HH;
    dxdt[2*nf+1] = -3.*p2 - Vpsi(phi,nf)/HH;
  }
  
}

#endif