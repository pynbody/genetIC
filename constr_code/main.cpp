#include "parser.hpp"
#include "ICgauss_deriv_clean.cc"
#include "ic.cpp"

typedef IC<MyFloat> ICf;

void setup_parser(ClassDispatch<ICf,void> &dispatch) {
    dispatch.add_class_route("Om",&ICf::setOmegaM0);
    dispatch.add_class_route("Ol",&ICf::setOmegaLambda0);
    dispatch.add_class_route("s8",&ICf::setSigma8);
    dispatch.add_class_route("Boxl",&ICf::setBoxLen);
    dispatch.add_class_route("zin",&ICf::setZ0);
    dispatch.add_class_route("n",&ICf::setn);
    dispatch.add_class_route("output",&ICf::setOutputMode);
    dispatch.add_class_route("seed",&ICf::setSeed);
    dispatch.add_class_route("camb",&ICf::setCambDat);
    dispatch.add_class_route("outdir",&ICf::setOutDir);
    dispatch.add_class_route("gadgetformat",&ICf::setGadgetFormat);

/*
    dispatch.add_class_route("IDfile",&ICf::loadID);
    dispatch.add_class_route("append_IDfile",&ICf::appendID);
    dispatch.add_class_route("select_sphere",&ICf::selectSphere);
    dispatch.add_class_route("centre_max",&ICf::centreDenmax);
    dispatch.add_class_route("centre_on",&ICf::centreParticle);
    dispatch.add_class_route("order",&ICf::reorderBuffer);
    dispatch.add_class_route("truncate",&ICf::truncateBuffer);
    dispatch.add_class_route("calculate",&ICf::calculate);
    dispatch.add_class_route("constrain",&ICf::constrain);
    */

}


int main(int argc, char *argv[]) {
    using namespace std;

    if(argc!=2)
    {
        cerr << "Usage: ./galaxylab paramfile" << endl;
        return -1;
    }

    ifstream inf;
    inf.open(argv[1]);

    if(!inf.is_open()) {
        cerr << "Error: could not open parameter file " << argv[1] << endl;
        exit(1);
    }

    // The main program is contained in this class:
    ICf generator;

    // Set up the command interpreter to issue commands to main_generator
    ClassDispatch<ICf,void> dispatch(generator);

    setup_parser(dispatch);

    // Read and act on the commands
    dispatch.run_loop(inf);

    // All done - write out
    generator.write();

    // Finished
    return 0;
}
