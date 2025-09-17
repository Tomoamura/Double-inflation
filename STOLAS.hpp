#ifndef INCLUDED_STOLAS_
#define INCLUDED_STOLAS_

#define _USR_MATH_DEFINES
#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
#include <complex>
#include <sys/time.h>

#include <boost/numeric/odeint.hpp>

#define euler_gamma 0.57721566490153286061

#include "parameters.hpp"
#include "model.hpp"
#include "src/vec_op.hpp"
#include "src/fft.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


#include <random>
std::mt19937 engine(1); // Fixed seed
// std::random_device seed; // random seed
// std::mt19937 engine(seed());
std::normal_distribution<> dist(0., 1.);


// useful macro
#define LOOP for(int i = 0; i < NL; i++) for(int j = 0; j < NL; j++) for(int k = 0; k < NL; k++)
#define LOOPLONG for(int i = 0; i < NLnoise; i++) for(int j = 0; j < NLnoise; j++) for(int k = 0; k < NLnoise; k++)


const std::string Nfileprefix = sdatadir + "/" + model + "/Nmap_";
const std::string fieldfileprefix = sdatadir + "/" + model + "/field_";
const std::string trajectoryfileprefix = sdatadir + "/" + model + "/trajectory_";
const std::string animationfileprefix = sdatadir + "/" + model + "/animation/animation_";
const std::string powfileprefix = sdatadir + "/" + model + "/power_";
const std::string powsfileprefix = sdatadir + "/" + model + "/powers";
bool noisefilefail, Nfilefail, sanisuperH = false;


int Nn, noisefiledirNo, noisefileNo;
std::ofstream Nfile, fieldfile, fieldfileA, trajectoryfile, powfile, powsfile;
std::vector<std::vector<double>> noisedataTr;
std::vector<std::vector<std::vector<double>>> noisedata(NFIELDS);
std::vector<double> Ndata;
std::vector<std::vector<std::vector<std::complex<double>>>> Nmap3D;

state_type phii(2*NFIELDS,0.);
std::vector<std::vector<double>> Phidata(NLnoise*NLnoise*NLnoise, std::vector<double>(2*NFIELDS,0.));
std::vector<std::vector<double>> phievol(NLnoise*NLnoise*NLnoise, std::vector<double>(2*NFIELDS,0.));
std::vector<std::vector<std::vector<double>>> phianimation;

std::vector<double> Nnoise(NLnoise*NLnoise*NLnoise, 0);
std::vector<double> Naverage(NLnoise*NLnoise*NLnoise, 0);
std::vector<double> Ntotal(NLnoise*NLnoise*NLnoise, 0);


bool checknoisefile();
bool checkNfilefail();
 
void dNmap(int noisefileNo);
void animation(int NoisefileDirNo, int noisefileNo);
void spectrum(std::vector<double> Ndata, int noisefiledirNo);

std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix);
bool compareByFirstColumn(const std::vector<double>& a, const std::vector<double>& b);

void STOLAS(int NoisefileDirNo, int NoisefileNo, int InterpolatingNo);


// -- functions -----------------------
void initialize(){
  Ndata = std::vector<double>(NLnoise*NLnoise*NLnoise,0);
  NFLOOP{
    phii[2*nf] = PSI_INIT;
    phii[2*nf+1] = DPSI_INIT;
  }

  LOOPLONG{
    phievol[NLnoise*NLnoise*i + NLnoise*j + k] = phii;
  }

  if (sanimation) {
    phianimation=std::vector<std::vector<std::vector<double>>> (int(totalstep/aninum),std::vector<std::vector<double>>(NLnoise*NLnoise*NLnoise, std::vector<double>(2*NFIELDS,0.)));
  }
}


void STOLAS(int NoisefileDirNo, int NoisefileNo, int InterpolatingNo) {
  std::ifstream noisefile;
  
  NFLOOP{
    noisefile.open(noisefiledir + std::to_string(NoisefileDirNo+nf) + noisefilenamediv + std::to_string(NoisefileNo) + "_" + std::to_string(InterpolatingNo) + std::string(".bin"), std::ios::binary);
    noisefilefail = noisefile.fail();
    
    if (!noisefile.fail()) {
      while (true) {
        std::vector<double> vv(NL);
        noisefile.read(reinterpret_cast<char*>(vv.data()), NL * sizeof(double));

        if (noisefile.gcount() != NL * sizeof(double)) {
          break;
        }
        noisedataTr.push_back(vv);
      }
      noisedata[nf] = transpose(noisedataTr);
      noisedataTr.clear();
      noisefile.close();
    }
  }

}


void dNmap(int NoiseNo, int InterpolatingNo) {
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NL; i++) {
    int numstep = 0;
    double N = dN*(totalstep + InterpolatingNo*itpstep) + Nnoise[i + NL*NoiseNo];
    state_type phi = phievol[i + NL*NoiseNo];

    // Find zero crossing time
    typedef boost::numeric::odeint::runge_kutta_dopri5<state_type>base_stepper_type;
    auto stepper = make_dense_output(1.0e-15, 1.0e-15, base_stepper_type());

    stepper.initialize(phi, N, dN);
    state_type phil(2*NFIELDS);
    state_type phir(2*NFIELDS);
    state_type phim(2*NFIELDS);

    while (true){
      stepper.do_step(dphidN);
      phi = stepper.current_state();
      N = stepper.current_time();
      
      if (EoI(phi)<0){
        double precphi = 1.e+2;
        double Nl = N;
        double Nr = N - stepper.current_time_step();
        double Nmid = 0;

        while (precphi>Nprec){
          stepper.calc_state(Nl, phil);
          stepper.calc_state(Nr, phir);
          Nmid = (Nl*EoI(phir) - Nr*EoI(phil)) / (EoI(phir)-EoI(phil));
          stepper.calc_state(Nmid, phim);
          if (EoI(phim)>0){
            Nr = Nmid;
          }
          else {
            Nl = Nmid;
          }
          precphi = std::fabs(EoI(phim));
        }

        N = Nmid;
        break;
      }

      if(strajectory && i==0 && NoiseNo==0){
        trajectoryfile << N << ' ' << phi[0] << ' ';
        double PSISUM = Psir(phi);
        trajectoryfile << PSISUM << ' ' << EoI(phi) << ' ' << hubble(phi) << std::endl;
      }
    }

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      Ndata[i + NL*NoiseNo] = N;
      phievol[i + NL*NoiseNo] = phi;
    }
  }
}


void evolution(int NoiseNo) {
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NL; i++) {
    double N = 0;
    int numstep = 0;
    int animationcount = 0; // for animation
    double sqrt_dN = std::sqrt(dN);
    state_type phi = phievol[i + NL*NoiseNo];

    // stepper
    boost::numeric::odeint::runge_kutta4<state_type> stepper_noise;

    for (size_t n=0; n<noisedata[0][0].size(); n++) {
      double psiamp = sqrt(calPpsi(phi));
      stepper_noise.do_step(dphidN, phi, N, dN);
      N += dN;
      double dw = 0.;
      NFLOOP{
        dw = noisedata[nf][i][n]; //dist(engine);//
        phi[2*nf] += psiamp * dw * sqrt_dN;
      }

      if(strajectory && i==0 && NoiseNo==0){
        double Nd = N;
        if(sanisuperH) Nd += dN*firststep;
        trajectoryfile << Nd << ' ' << phi[0] << ' ';
        double PSISUM = Psir(phi);
        trajectoryfile << PSISUM << ' ' << EoI(phi) << ' ' << hubble(phi) << std::endl;
      }

      animationcount++;
      if(sanimation && animationcount>aninum-1){
        int ef = int(n/aninum);
        if(sanisuperH) ef += int(firststep/aninum);
        phianimation[ef][i + NL*NoiseNo] = phi;
        animationcount=0;
      }
    }


#ifdef _OPENMP
#pragma omp critical
#endif
    {
      phievol[i + NL*NoiseNo] = phi;
    }
  }
}


void evolutionNoise(int NoiseNo) {
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NL; i++) {
    double N = 0;
    double sqrt_dN = std::sqrt(dN);
    state_type phi = phievol[i + NL*NoiseNo];

    // stepper
    boost::numeric::odeint::runge_kutta4<state_type> stepper_noise;

    while (EoN(phi)>0) {
      double psiamp = sqrt(calPpsi(phi));
      stepper_noise.do_step(dphidN, phi, N, dN);
      N += dN;
      NFLOOP{
        double dw = dist(engine);
        phi[2*nf] += psiamp * dw * sqrt_dN;
      }

      if(strajectory && i==0 && NoiseNo==0){
        trajectoryfile << N+dN*totalstep << ' ' << phi[0] << ' ';
        double PSISUM = Psir(phi);
        trajectoryfile << PSISUM << ' ' << EoI(phi) << ' ' << hubble(phi) << std::endl;
      }
    }

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      phievol[i + NL*NoiseNo] = phi;
      Nnoise[i + NL*NoiseNo] = N;
    }
  }
}


void OpenFiles(int NoisefiledirNo, int Interpolatingnumber){
  std::string NLfilename = std::to_string(NLnoise) + std::string("_") + std::to_string(NFIELDS) + std::string("_") + std::to_string(NoisefiledirNo);
  std::string InterFileName = NLfilename + std::string("_") + std::to_string(Interpolatingnumber) + std::string(".dat");

  Nfile.open(Nfileprefix + InterFileName);
  Nfile << std::setprecision(10);
  Nfilefail = Nfile.fail();

  if (sfield) {
    fieldfile.open(fieldfileprefix + InterFileName);
    fieldfile << std::setprecision(10);
  }
  
  if (strajectory) {
    trajectoryfile.open(trajectoryfileprefix + InterFileName);
    trajectoryfile << std::setprecision(10);
  }
  
}


void save_zeta(int NoiseNo, int Interpolatingnumber){
  // double noisestepN = dN*(totalstep + Interpolatingnumber*itpstep);
  for (int i=0; i<NL; i++){
    double Nsave = 0.;
    if(EoI_noise) {
      Nsave = Naverage[i + NL*NoiseNo];// + noisestepN;
    }
    else {
      Nsave = Ndata[i + NL*NoiseNo];// + Nnoise[i + NL*NoiseNo] + noisestepN;
    }
    // double Nsave = Ndata[i + NL*NoiseNo] + noisestepN;
    // if(EoI_noise) Nsave += Nnoise[i + NL*NoiseNo];
    Nfile << i + NL*NoiseNo << ' ' << Nsave << std::endl;
  }
}


void save_field(int NoiseNo){
  for (int i=0; i<NL; i++) {
    fieldfile << i + NL*NoiseNo << ' ';
    NFLOOP{
      fieldfile << phievol[i + NL*NoiseNo][2*nf] << ' ';
    }
    fieldfile << std::endl;
  }
}


inline int index(int i, int j, int k) {
  return i * NLnoise * NLnoise + j * NLnoise + k;
}

inline int indexsmall(int i, int j, int k) {
  return i * NLnoise * NLnoise/4 + j * NLnoise/2 + k;
}

inline int INCREMENT(int i) {
  // return (i == NLnoise - 1) ? 0 : i + 1;
  return (i == NLnoise - 1) ? i : i + 1;
}

inline int DECREMENT(int i) {
  // return (i == 0) ? NLnoise - 1 : i - 1;
  return i - 1;
}


double LinearInterplation(std::vector<std::vector<std::vector<double>>> LatticeField, int x, int y, int z){
  int xp = INCREMENT(x);
  int xm = DECREMENT(x);
  int yp = INCREMENT(y);
  int ym = DECREMENT(y);
  int zp = INCREMENT(z);
  int zm = DECREMENT(z);
  int oddMask = (x % 2) * 4 + (y % 2) * 2 + (z % 2);

  switch (oddMask) {
    case 0:
      return LatticeField[x/2][y/2][z/2];
    
    case 4: // on the x-axis
      return 0.5*(LatticeField[xp/2][y/2][z/2] + LatticeField[xm/2][y/2][z/2]);
      
    case 2: // on the y-axis
      return 0.5*(LatticeField[x/2][yp/2][z/2] + LatticeField[x/2][ym/2][z/2]);
    
    case 1: // on the z-axis
      return 0.5*(LatticeField[x/2][y/2][zp/2] + LatticeField[x/2][y/2][zm/2]);

    case 6: // on the xy-plane
      return 0.25*( LatticeField[xp/2][yp/2][z/2] + LatticeField[xp/2][ym/2][z/2]
                  + LatticeField[xm/2][yp/2][z/2] + LatticeField[xm/2][ym/2][z/2]);

    case 5: // on the xz-plane
      return 0.25*( LatticeField[xp/2][y/2][zp/2] + LatticeField[xp/2][y/2][zm/2]
                  + LatticeField[xm/2][y/2][zp/2] + LatticeField[xm/2][y/2][zm/2]);

    case 3: // on the yz-plane
      return 0.25*( LatticeField[x/2][yp/2][zp/2] + LatticeField[x/2][yp/2][zm/2]
                   + LatticeField[x/2][ym/2][zp/2] + LatticeField[x/2][ym/2][zm/2]);
    
    case 7: // on the mid-point
      return 0.125*(LatticeField[xp/2][yp/2][zp/2] + LatticeField[xp/2][yp/2][zm/2]
                  + LatticeField[xp/2][ym/2][zp/2] + LatticeField[xp/2][ym/2][zm/2]
                  + LatticeField[xm/2][yp/2][zp/2] + LatticeField[xm/2][yp/2][zm/2]
                  + LatticeField[xm/2][ym/2][zp/2] + LatticeField[xm/2][ym/2][zm/2]);

    default:
      std::cout << "Bad interpolation" << std::endl;
      return -1;
  }
}

int findMaxZeta(){
  std::vector<double> Nblock(8,0);

  LOOPLONG{
    int idxsum = (i < NLnoise/2 ? 0 : 1) + (j < NLnoise/2 ? 0 : 2) + (k < NLnoise/2 ? 0 : 4);
    int idxf = index(i,j,k);
    NFLOOP{
      Nblock[idxsum] += Ndata[idxf];
    }
  }
  
  return std::distance(Nblock.begin(), std::max_element(Nblock.begin(), Nblock.end()));
}


std::vector<int> findNMaxBox(std::vector<double> NData){
  int maxNpoint = std::distance(NData.begin(), std::max_element(NData.begin(), NData.end()));

  int x = maxNpoint/NLnoise/NLnoise, y = (maxNpoint%(NLnoise*NLnoise))/NLnoise, z = maxNpoint%NLnoise;
  int xb = x - NLnoise/4, yb = y - NLnoise/4, zb = z - NLnoise/4;

  if(x > NLnoise/2) xb = NLnoise/2;
  if(xb<0) xb = 0;
  if(y > NLnoise/2) yb = NLnoise/2;
  if(yb<0) yb = 0;
  if(z > NLnoise/2) zb = NLnoise/2;
  if(zb<0) zb = 0;

  std::vector<int> ShiftVector{xb,yb,zb};

  return ShiftVector;
}


std::vector<std::vector<double>> InterpolatingPhi(std::vector<int> Shift){
  int NLsmall = NLnoise/2;

  std::vector<std::vector<std::vector<std::vector<double>>>> LatticeField(2*NFIELDS+2,std::vector<std::vector<std::vector<double>>>(NLsmall,std::vector<std::vector<double>>(NLsmall,std::vector<double>(NLsmall,0.))));
  std::vector<std::vector<double>> phinew(NLnoise*NLnoise*NLnoise, std::vector<double>(2*NFIELDS,0));

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int x = 0; x < NLsmall; x++) {
    for (int y = 0; y < NLsmall; y++) {
      for (int z = 0; z < NLsmall; z++) {
        int idx = index(x+Shift[0],y+Shift[1],z+Shift[2]);
        for (int nf = 0; nf < 2*NFIELDS; nf++) {
          LatticeField[nf][x][y][z] = phievol[idx][nf];
        }
      }
    }
  }

  // interpolation
  int complete=0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int x = 0; x < NLnoise; x++) {
    for (int y = 0; y < NLnoise; y++) {
      for (int z = 0; z < NLnoise; z++) {
        int idx = index(x,y,z);
        for (int nf = 0; nf < 2*NFIELDS; nf++) {
          phinew[idx][nf] = LinearInterplation(LatticeField[nf],x,y,z);
        }
      }
    }
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      complete++;
      std::cout << "\rInterpolation       : " << std::setw(3) << 100*complete/NLnoise << "%" << std::flush;
    }
  }
  std::cout << std::endl;

  return phinew;
}


std::ifstream fieldinitfile;
std::vector<std::vector<double>> fieldindata;
std::vector<std::vector<double>> getphifile(const int order, int NoisefiledirNo){
  std::vector<std::vector<double>> FieldinData;
  int NLsmall = NLnoise/2;
  std::string FieldFileName = fieldfileprefix + std::to_string(NLnoise) + "_" + std::to_string(NoisefiledirNo) + "_" + std::to_string(order-1) + std::string(".dat");
  std::cout << FieldFileName << std::endl;

  fieldinitfile.open(FieldFileName);
  if (fieldinitfile.fail()) std::cout << "field init file fail" << std::endl;

  double f1, f2, f3, f4, f5;
  while (fieldinitfile >> f1 >> f2 >> f3 >> f4 >> f5) {
    FieldinData.push_back({f1, f2, f3, f4, f5});
  }
  fieldinitfile.close();

  std::sort(FieldinData.begin(), FieldinData.end(), compareByFirstColumn);

  return FieldinData;
}

bool checknoisefile() {
  return noisefilefail;
}

bool checkNfilefail() {
  return Nfilefail;
}


// Export animation
void animation(int NoisefiledirNo, int EachField){
  std::cout << "ExportAnimation" << std::endl;

  if (EachField==0){
    for (int n = 0; n < phianimation.size(); n++) {
      std::string fFileName = std::to_string(NLnoise) + std::string("_") + std::to_string(NFIELDS) + std::string("_") + std::to_string(NoisefiledirNo) + std::string("_") + std::to_string(n);
      fieldfileA.open(animationfileprefix + fFileName + std::string(".dat"));
      fieldfileA << std::setprecision(10);

      for (int i=0; i<phianimation[n].size(); i++) {
        double PSISUM = 0;
        NFLOOP{
          PSISUM += pw2(phianimation[n][i][2*nf]);
        }
        fieldfileA << sqrt(PSISUM) << ' ';
      }
      fieldfileA << std::endl;
      fieldfileA.close();
    }
  }
  else{ // export each field
    NFLOOP{
      std::string fFileName = std::to_string(NLnoise) + std::string("_") + std::to_string(NFIELDS) + std::string("_") + std::to_string(noisefiledirNo);
      fieldfileA.open(animationfileprefix + fFileName + std::string("_") + std::to_string(nf) + std::string(".dat"));
      fieldfileA << std::setprecision(10);

      for (int n = 0; n < phianimation.size(); n++) {
        for (int i=0; i<phianimation[n].size(); i++) {
          fieldfileA << phianimation[n][i][2*nf] << ' ';
        }
        fieldfileA << std::endl;
      }
      fieldfileA.close();
    }
  }

  // std::string fFileName = std::to_string(NLnoise) + std::string("_") + std::to_string(NFIELDS) + std::string("_") + std::to_string(NoisefiledirNo);
  // fieldfileA.open(animationfileprefix + std::string("phi_") + fFileName + std::string(".dat"));

  // fieldfileA << std::setprecision(16);
  // for (int n = 0; n < phianimation.size(); n++) {
  //   for (int i=0; i<phianimation[n].size(); i++) {
  //     fieldfileA << phianimation[n][i][0] << ' ';
  //   }
  //   fieldfileA << std::endl;
  // }
  // fieldfileA.close();
}

std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix) {
  if (matrix.empty() || matrix[0].empty()) return {}; // 空の場合の処理

  size_t rows = matrix.size();
  size_t cols = matrix[0].size();
  std::vector<std::vector<double>> result(cols, std::vector<double>(rows));

  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      result[j][i] = matrix[i][j];
    }
  }
  return result;
}

bool compareByFirstColumn(const std::vector<double>& a, const std::vector<double>& b) {
  return a[0] < b[0];
}

// Calculate power spectrum
void spectrum(std::vector<double> Ndata, int noisefiledirNo) {
  // powfile.open(powfileprefix + std::to_string(NLnoise) + std::string("_") + std::to_string(noisefiledirNo) + std::string(".dat"));
  // powfile << std::setprecision(10);
  std::vector<std::vector<std::vector<std::complex<double>>>> Nmap3D = std::vector<std::vector<std::vector<std::complex<double>>>>(NLnoise, std::vector<std::vector<std::complex<double>>>(NLnoise, std::vector<std::complex<double>>(NLnoise, 0)));

  for (int i=0; i<NLnoise*NLnoise*NLnoise; i++) {
    int x=i/NLnoise/NLnoise ,y=(i%(NLnoise*NLnoise))/NLnoise, z=i%NLnoise;
    if(EoI_noise) {
      Nmap3D[x][y][z] = Ndata[i];
    }
    else {
      Nmap3D[x][y][z] = Ndata[i];
    }
  }
  std::vector<std::vector<std::vector<std::complex<double>>>> Nk=fft_fftw(Nmap3D);
  
  powsfile.open(powsfileprefix + std::string(".dat"), std::ios::app);
  powsfile << std::setprecision(10);
  int imax = ceil(log(NLnoise/2)/dlogn);
  std::vector<double> disc_power(imax, 0);

  LOOPLONG{
    int nxt, nyt, nzt; // shifted index
    if (i<=NLnoise/2) {
      nxt = i;
    } else {
      nxt = i-NLnoise;
    }

    if (j<=NLnoise/2) {
      nyt = j;
    } else {
      nyt = j-NLnoise;
    }

    if (k<=NLnoise/2) {
      nzt = k;
    } else {
      nzt = k-NLnoise;
    }
    
    double rk=nxt*nxt+nyt*nyt+nzt*nzt;
    // powfile << sqrt(rk) << " " << norm(Nk[i][j][k])/NLnoise/NLnoise/NLnoise/NLnoise/NLnoise/NLnoise << std::endl;

    double LogNk = log(sqrt(rk));
    double calPk = norm(Nk[i][j][k])/NLnoise/NLnoise/NLnoise/NLnoise/NLnoise/NLnoise;
    for (size_t ii = 0; ii < imax; ii++) {
      if ((dlogn*ii-LogNk)<dlogn/2. && -dlogn/2.<=(dlogn*ii-LogNk) && rk!=0) {
        disc_power[ii] += calPk/dlogn;
        break;
      }
    }
  }
  powsfile << noisefiledirNo << " ";
  for (size_t ii = 0; ii < imax; ii++) {
    powsfile << disc_power[ii] << " " ;
  }
  powsfile << std::endl;
  powsfile.close();
  // powfile.close();
  std::cout << "ExportPowerSpectrum" << std::endl;
}


#endif
