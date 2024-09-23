#include "bindings.hpp"
#include "ic.hpp"

void usageMessage() {
  using namespace std;
  cout << "Usage: genetIC_mapper <paramfile1> <paramfile2> <id-file-input> <id-file-output>" << endl << endl
       << " Processes the geometry in paramfile1 and paramfile2, then takes a file of particle IDs from" << endl
        << " id-file-input and maps them from the geometry in paramfile1 to the geometry in paramfile2." << endl
        << " The mapped IDs are written to id-file-output." << endl;
}

int main(int argc, char *argv[]) {
  using namespace std;

  bool useFloat = false;

  if (argc<4) {
    usageMessage();
    return -1;
  }

  std::string fname1(argv[1]);
  std::string fname2(argv[2]);
  std::string idFileInput(argv[3]);
  std::string idFileOutput(argv[4]);

  logging::entry() << "Processing geometry in " << fname1 << std::endl;
  dummyic::DummyICGenerator<double> generator1;
  {
    logging::IndentWhileInScope temporaryIndent;
    runInterpreter<ICGenerator<double>>(generator1, fname1);
  }

  dummyic::DummyICGenerator<double> generator2(&generator1);
  logging::entry() << "Processing geometry in " << fname2 << std::endl;
  {
    logging::IndentWhileInScope temporaryIndent;
    runInterpreter<ICGenerator<double>>(generator2, fname2);
  }



  logging::entry() << "Loading IDs from " << idFileInput << std::endl;
  generator1.loadID(idFileInput);

  logging::entry() << "Writing IDs to " << idFileOutput << std::endl;

  generator2.ICGenerator<double>::dumpID(idFileOutput);

  logging::entry() << "Done." << std::endl;


}
