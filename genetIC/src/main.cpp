#include <stdlib.h>
#include <sstream>
#include <ctime>
#include <iostream>

#include "bindings.hpp"
#include "tools/logging.hpp"
#include "ic.hpp"

int runGenetic(string paramFilename);

#ifdef USE_COMPLEX_FFT
using ICf = ICGenerator<complex<float>>;
using ICd = ICGenerator<complex<double>>;
#else
using ICf = ICGenerator<float>;
using ICd = ICGenerator<double>;
#endif

void header(ostream &outf, std::string prefix="") {
  outf << prefix << "genetIC v1.5.0, compiled: " << __DATE__ << ", " << __TIME__ << endl;
  time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%b %d %Y, %X", &tstruct);
  outf << prefix << "                 runtime: " << buf << endl << endl;
#ifdef GIT_VERSION
  if (sizeof(GIT_VERSION) > 1)
    outf << prefix << "git HEAD:" << GIT_VERSION << endl;
  if (sizeof(GIT_MODIFIED) > 1) {
    outf << prefix << "At the time of compilation, the following files had been modified:" << endl;
    std::istringstream modifFiles(GIT_MODIFIED);
    while(!modifFiles.eof()) {
      std::string thisFile;
      modifFiles >> thisFile;
      outf << prefix << "  " << thisFile << endl;
    }
  }
#endif


}


template<typename IC>
int runGenetic(string paramFilename) {
  IC generator;

  ofstream outf;
  outf.open("IC_output.params");

  if (!outf.is_open()) {
    logging::entry() << "Error: could not open output file" << endl;
    exit(1);
  }

  header(outf, "# ");
  header(cerr);

  // Initialising FFTW is not strictly required (it will be done later if missed here) but, since it might print
  // some output, it's neater to do it straight away:
  tools::numerics::fourier::initialise();

  runInterpreter(generator, paramFilename, outf);

  return 0;
}


void usageMessage() {
  using namespace std;
  cout << "Usage: genetIC paramfile [-f] [-c <cache-size>]" << endl << endl
       << " The paramfile is a text file of commands (see example provided with genetIC distribution)." << endl << endl
       << " If option -f is specified, genetIC uses float (32-bit) instead of double (64-bit) internally." << endl
       << " The output format is unaffected by the internal bit depth." << endl << endl
       << " If option -c <cache-size> is specified, the transfer function cache is limited to the specified"
          " number of entries. This can be used to reduce memory usage, but may slow down the code." << endl;
}

int main(int argc, char *argv[]) {
  using namespace std;

  bool useFloat = false;

  if (argc<2) {
    usageMessage();
    return -1;
  }


  string fname(argv[1]);

  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "-f") == 0) {
      useFloat = true;
    } else if (strcmp(argv[i], "-c") == 0) {
      if (i + 1 >= argc) {
        cerr << "Error: -c option requires an argument" << endl;
        return -1;
      }
      cosmology::lru_cache_size = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
      usageMessage();
      return 0;
    } else {
      if (i != argc - 1) {
        cerr << "Error: unknown option " << argv[i] << endl;
        return -1;
      }
    }
  }


  if(useFloat) {
    return runGenetic<ICf>(fname);
  } else {
    return runGenetic<ICd>(fname);
  }


}
