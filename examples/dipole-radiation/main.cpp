
#include <iostream>
#include <string>

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

extern void read_nfff_grid();
extern void process_nfff_data();
extern void get_radiation_pattern();

int main (int, char *[])
{
    // Try block to detect exceptions raised by any of the calls inside it
    try {
        // Turn off the auto-printing when failure occurs so that we can
        // handle the errors appropriately
        Exception::dontPrint();

        read_nfff_grid();

        process_nfff_data();

        get_radiation_pattern();

    }  // end of try block

    // catch failure caused by the H5File operations
    catch (FileIException error) {
        error.printError();
        return -1;
    }

    // catch failure caused by the DataSet operations
    catch(DataSetIException error) {
        error.printError();
        return -1;
    }

    return 0;  // successfully terminated
}
