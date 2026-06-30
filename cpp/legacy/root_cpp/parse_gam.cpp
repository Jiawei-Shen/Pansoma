#include <iostream>
#include <fstream>
#include <vector>
#include "vg.pb.h"  // Protobuf-generated header file

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <gam_file>" << endl;
        return 1;
    }

    string gam_filename = argv[1];
    ifstream gam_file(gam_filename, ios::binary);

    if (!gam_file) {
        cerr << "Error: Cannot open GAM file: " << gam_filename << endl;
        return 1;
    }

    // Read and parse each alignment in the GAM file
    vg::Alignment alignment;
    while (alignment.ParseFromIstream(&gam_file)) {
        cout << "Read Name: " << alignment.name() << endl;
        cout << "Sequence: " << alignment.sequence() << endl;
        cout << "Mapping Quality: " << alignment.mapping_quality() << endl;
        cout << "--------------------------" << endl;
    }

    gam_file.close();
    return 0;
}
