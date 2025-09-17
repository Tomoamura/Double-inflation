#include "STOLAS.hpp"

// main function
int main(int argc, char* argv[])
{
  if (argc!=2) {
    std::cout << "Specify the noise file number correctly." << std::endl;
    return -1;
  }
  
  // ---------- start timer ----------
  struct timeval Nv;
  struct timezone Nz;
  double before, after;
  
  gettimeofday(&Nv, &Nz);
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // ---------------------------------

#ifdef _OPENMP
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
#endif

  int noisefiledirNo = atoi(argv[1]);
  // internumber = 10;//atoi(argv[1]);

  std::cout << "Noise file No. : " << noisefiledirNo << std::endl;
  std::cout << "model : " << model << std::endl;

  // initialize fields
  initialize();

  // Solve fields
  for (size_t interpolatingnumber = 0; interpolatingnumber < internumber+1; interpolatingnumber++) {
    OpenFiles(noisefiledirNo, interpolatingnumber);
    
    if (checkNfilefail()) {
      std::cout << "The export file couldn't be opened. 'mkdir data'" << std::endl;
      return -1;
    }

    // fields evolution
    for (int noiseNo = 0; noiseNo < totalnoiseNo; noiseNo++) {
      STOLAS(noisefiledirNo,noiseNo,interpolatingnumber);
      if (checknoisefile()) {
        std::cout << "The noise file couldn't be opened." << std::endl;
        return -1;
      }

      evolution(noiseNo);

      std::cout << "\rLatticeSimulation   : " << std::setw(3) << 100*noiseNo/(totalnoiseNo-1) << "%" << std::flush;
    }
    std::cout << std::endl;


    // keep fields map
    Phidata = phievol;

    // Solve untill superH
    sanisuperH = true;
    // if (interpolatingnumber==internumber) {
      for (int noiseNo = 0; noiseNo < totalnoiseNo; noiseNo++) {
        STOLAS(noisefiledirNo,noiseNo,999);
        if (checknoisefile()) {
          std::cout << "The noise file couldn't be opened." << std::endl;
          return -1;
        }
    
        evolution(noiseNo);

        // Take an average for each patches
        if(EoI_noise) {
          std::vector<std::vector<double>> PhidataAv = phievol;
          for (int av = 0; av < MeanNumber; av++) {
            evolutionNoise(noiseNo); // Adding the noise until EoN()
            dNmap(noiseNo,interpolatingnumber);
            Ntotal += Ndata;
            Ndata.assign(Ndata.size(), 0.0); // reset vector
            phievol = PhidataAv;
          }
          for (int i=0; i<NL; i++) Naverage[i + NL*noiseNo] = Ntotal[i + NL*noiseNo]/(double)MeanNumber;
        }
        else {
          // evolutionNoise(noiseNo); // Adding the noise until EoN()
          dNmap(noiseNo,interpolatingnumber);
          Naverage = Ndata;
        }

        save_zeta(noiseNo,interpolatingnumber); // svave delta N map
        if(sfield) save_field(noiseNo);

        std::cout << "\rCompute delta N map : " << std::setw(3) << 100*noiseNo/(totalnoiseNo-1) << "%" << std::flush;
      }
      std::cout << std::endl;
      if(spower) spectrum(Naverage,interpolatingnumber);
    // }
    Nfile.close();
    fieldfile.close();

    // ---------- stop timer ----------
    gettimeofday(&Nv, &Nz);
    after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
    std::cout << after - before << " sec." << std::endl;
    before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
    // --------------------------------

    // Interpolation
    if (interpolatingnumber!=internumber) {
      // int maxNpoint = std::distance(Ndata.begin(), std::max_element(Ndata.begin(), Ndata.end()));

      // int x = maxNpoint/NLnoise/NLnoise, y = (maxNpoint%(NLnoise*NLnoise))/NLnoise, z = maxNpoint%NLnoise;
      // int xb = x - NLnoise/4, yb = y - NLnoise/4, zb = z - NLnoise/4;

      // if(x > NLnoise/2) xb = NLnoise/2;
      // if(xb<0) xb = 0;
      // if(y > NLnoise/2) yb = NLnoise/2;
      // if(yb<0) yb = 0;
      // if(z > NLnoise/2) zb = NLnoise/2;
      // if(zb<0) zb = 0;

      // std::vector<int> shift = findNMaxBox(Ndata);
      std::vector<int> shift{NLnoise/4,NLnoise/4,NLnoise/4};

      phievol = Phidata; // reset field values
      Ntotal.assign(Ntotal.size(), 0.0); // reset vector
      phievol = InterpolatingPhi(shift);
      std::cout << "Number of interpolation: " << interpolatingnumber+1 << std::endl;
    }

    // ---------- stop timer ----------
    gettimeofday(&Nv, &Nz);
    after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
    std::cout << after - before << " sec." << std::endl;
    before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
    // --------------------------------
  }

  if(sanimation) animation(noisefiledirNo,0);


  // ---------- stop timer ----------
  // gettimeofday(&Nv, &Nz);
  // after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // std::cout << after - before << " sec." << std::endl;
  // std::cout << std::endl;
  // --------------------------------
}
