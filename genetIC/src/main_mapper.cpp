#include "bindings.hpp"
#include "ic.hpp"

void usageMessage() {
  using namespace std;
  cout << "Usage(1): genetIC_mapper <paramfile1> <paramfile2> <id-file-input> <id-file-output>" << endl << endl
       << " Processes the geometry in paramfile1 and paramfile2, then takes a file of particle IDs from" << endl
        << " id-file-input and maps them from the geometry in paramfile1 to the geometry in paramfile2." << endl
        << " The mapped IDs are written to id-file-output." << endl;
  cout << "Usage(2): genetIC_mapper <paramfile> <id-file-output>" << endl << endl
       << " Processes the geometry in paramfile and writes flagged IDs to id-file-output." << endl;
}

int main(int argc, char *argv[]) {
  using namespace std;

  bool useFloat = false;

  if (argc!=5 && argc!=3) {
    logging::entry() << "argc = " << argc << std::endl;
    usageMessage();
    return -1;
  }

  if(argc == 5) {

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
    generator1.clearCellFlags();
    generator2.clearCellFlags();

    // the input IDs are relative to the output for generator1, not any input mapper that may be active:
    generator1.clearInputMapper();

    generator1.loadID(idFileInput);

    // generator2.propagateFlagsToRefinedCells(generator1.getMultiLevelContext());
    generator1.propagateFlagsToRefinedCells(generator2.getMultiLevelContext());



    logging::entry() << "Writing IDs to " << idFileOutput << std::endl;

    generator2.ICGenerator<double>::dumpID(idFileOutput);

    logging::entry() << "Done." << std::endl;
  } else {
    std::string fname(argv[1]);
    std::string idFileOutput(argv[2]);

    logging::entry() << "Processing geometry in " << fname << std::endl;
    dummyic::DummyICGenerator<double> generator;
    {
      logging::IndentWhileInScope temporaryIndent;
      runInterpreter<ICGenerator<double>>(generator, fname);
    }

    logging::entry() << "Writing flagged IDs to " << idFileOutput << std::endl;
    generator.ICGenerator<double>::dumpID(idFileOutput);

    logging::entry() << "Done." << std::endl;
  }


}
