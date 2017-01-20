
#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkRectilinearGrid.h>
#include <vtkSmartPointer.h>
#include <pugixml.hpp>

#define RESULT_FILE "result.pvd"
#define NFFF_FILE "nfff.xml"

extern "C" double getEx(int i, int j, int k);
extern "C" double getEy(int i, int j, int k);
extern "C" double getEz(int i, int j, int k);
extern "C" double getJx(int i, int j, int k);
extern "C" double getJy(int i, int j, int k);
extern "C" double getJz(int i, int j, int k);
extern "C" double getHx(int i, int j, int k);
extern "C" double getHy(int i, int j, int k);
extern "C" double getHz(int i, int j, int k);

extern "C" double getE(int i, int j, int k);
extern "C" double getJ(int i, int j, int k);
extern "C" double getH(int i, int j, int k);

extern "C" void init_data_pvd()
{
    // получили документ
    pugi::xml_document doc;

    // добавили "VTKFile" узел
    pugi::xml_node file = doc.append_child("VTKFile");

    //задали атрибуты
    file.append_attribute("type") = "Collection";
    file.append_attribute("version") = "0.1";
    file.append_attribute("byte_order") = "LittleEndian";
    file.append_attribute("compressor") = "vtkZLibDataCompressor";

    //добавили "Collection"
    pugi::xml_node coll = file.append_child("Collection");

    //сохранили
    doc.save_file(RESULT_FILE);
}

static void _append_data_set_section(const char *file_name, const double time)
{
    int state;
    pugi::xml_document doc;
    char time_str[32];

    pugi::xml_parse_result result = doc.load_file(RESULT_FILE);

    pugi::xml_node file = doc.child("VTKFile");
    pugi::xml_node coll = file.child("Collection");

    //добавили "DataSet"
    pugi::xml_node data = coll.append_child("DataSet");

    sprintf(time_str, "%le", time);
    data.append_attribute("timestep") = time_str;
    data.append_attribute("group") = "";
    data.append_attribute("part") = "0";
    data.append_attribute("file") = file_name;

    //сохранили
    state = doc.save_file(RESULT_FILE);
}

extern "C" void save_data_vtr(const char *file_name, const double time,
                              const int nx, const double *xi,
                              const int ny, const double *yi,
                              const int nz, const double *zi)
{
    int i, j, k;

    // Create a rectilinear grid by defining three arrays specifying the
    // coordinates in the x-y-z directions.
    vtkFloatArray *xCoords = vtkFloatArray::New();
    xCoords->SetName("X");
    for (i = 0; i <= nx; ++i) {
        xCoords->InsertNextValue(xi[i]);
    }

    vtkFloatArray *yCoords = vtkFloatArray::New();
    yCoords->SetName("Y");
    for (i = 0; i <= ny; ++i) {
        yCoords->InsertNextValue(yi[i]);
    }

    vtkFloatArray *zCoords = vtkFloatArray::New();
    zCoords->SetName("Z");
    for (i = 0; i <= nz; ++i) {
        zCoords->InsertNextValue(zi[i]);
    }

    // Create a grid
    vtkSmartPointer<vtkRectilinearGrid> rectilinearGrid =
        vtkSmartPointer<vtkRectilinearGrid>::New();

    rectilinearGrid->SetDimensions(nx+1, ny+1, nz+1);
    rectilinearGrid->SetXCoordinates(xCoords);
    rectilinearGrid->SetYCoordinates(yCoords);
    rectilinearGrid->SetZCoordinates(zCoords);

    float pnts[3];

    vtkFloatArray *vectorE = vtkFloatArray::New();
    vectorE->SetName("E");
    vectorE->SetNumberOfComponents(3);
    for (k = 1; k <= nz; ++k) {
        for (j = 1; j <= ny; ++j) {
            for (i = 1; i <= nx; ++i) {
                pnts[0] = getEx(i, j, k);
                pnts[1] = getEy(i, j, k);
                pnts[2] = getEz(i, j, k);
                vectorE->InsertNextTuple(pnts);
            }
        }
    }
    rectilinearGrid->GetCellData()->AddArray(vectorE);

    vtkFloatArray *vectorH = vtkFloatArray::New();
    vectorH->SetName("H");
    vectorH->SetNumberOfComponents(3);
    for (k = 1; k <= nz; ++k) {
        for (j = 1; j <= ny; ++j) {
            for (i = 1; i <= nx; ++i) {
                pnts[0] = getHx(i, j, k);
                pnts[1] = getHy(i, j, k);
                pnts[2] = getHz(i, j, k);
                vectorH->InsertNextTuple(pnts);
            }
        }
    }
    rectilinearGrid->GetCellData()->AddArray(vectorH);

    // Write file
    vtkSmartPointer<vtkXMLRectilinearGridWriter> writer =
        vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
    writer->SetFileName(file_name);
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(structuredGrid);
#else
    writer->SetInputData(rectilinearGrid);
#endif
    writer->Write();

    _append_data_set_section(file_name, time);

    xCoords->Delete();
    yCoords->Delete();
    zCoords->Delete();
    vectorE->Delete();
    vectorH->Delete();
}

extern "C" void init_data_nfff(const double time)
{
    // получили документ
    pugi::xml_document doc;

    // добавили "NFFF" узел
    pugi::xml_node nfff = doc.append_child("NFFF");

    //задали атрибуты
    nfff.append_attribute("time") = time;

    //добавили "HDF5"
    pugi::xml_node hdf5 = nfff.append_child("HDF5");

    //сохранили
    doc.save_file(NFFF_FILE);
}

extern "C" void append_nfff_data(const char *file_name, const double time)
{
    int state;
    pugi::xml_document doc;

    pugi::xml_parse_result result = doc.load_file(NFFF_FILE);

    pugi::xml_node nfff = doc.child("NFFF");
    pugi::xml_node hdf5 = nfff.child("HDF5");

    //добавили "data"
    pugi::xml_node data = hdf5.append_child("data");

    data.append_attribute("timestep") = time;
    data.append_attribute("file") = file_name;

    //сохранили
    state = doc.save_file(NFFF_FILE);
}
