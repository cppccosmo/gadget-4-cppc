#include <mpi.h>
#include <stdio.h>
#include <hdf5.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "gadgetconfig.h"
#include "data/constants.h"
#include "data/dtypes.h"
#include "data/macros.h"
#include "io/io.h"
hid_t my_H5Acreate(hid_t loc_id, const char *attr_name, hid_t type_id, hid_t space_id, hid_t acpl_id);
hid_t my_H5Screate(H5S_class_t type);
herr_t my_H5Aclose(hid_t attr_id, const char *attr_name);
herr_t my_H5Awrite(hid_t attr_id, hid_t mem_type_id, const void *buf, const char *attr_name);
herr_t my_H5Sclose(hid_t dataspace_id, H5S_class_t type);
herr_t my_H5Tclose(hid_t type_id);

void IO_Def::write_compile_time_options_in_hdf5(hid_t handle)
{
hid_t hdf5_dataspace, hdf5_attribute;
double val;
hid_t atype = H5Tcopy(H5T_C_S1);
H5Tset_size(atype, 1);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);hdf5_attribute = my_H5Acreate(handle, "ADDITIONAL_GRID" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ADDITIONAL_GRID");
my_H5Aclose(hdf5_attribute, "ADDITIONAL_GRID");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);hdf5_attribute = my_H5Acreate(handle, "CB_PHASE" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "CB_PHASE");
my_H5Aclose(hdf5_attribute, "CB_PHASE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);hdf5_attribute = my_H5Acreate(handle, "DOUBLEPRECISION" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 2;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "DOUBLEPRECISION");
my_H5Aclose(hdf5_attribute, "DOUBLEPRECISION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);hdf5_attribute = my_H5Acreate(handle, "MFLR_RST" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "MFLR_RST");
my_H5Aclose(hdf5_attribute, "MFLR_RST");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);hdf5_attribute = my_H5Acreate(handle, "MULTIPOLE_ORDER" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "MULTIPOLE_ORDER");
my_H5Aclose(hdf5_attribute, "MULTIPOLE_ORDER");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);hdf5_attribute = my_H5Acreate(handle, "NGENIC" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1024;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "NGENIC");
my_H5Aclose(hdf5_attribute, "NGENIC");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);hdf5_attribute = my_H5Acreate(handle, "NSOFTCLASSES" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 2;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "NSOFTCLASSES");
my_H5Aclose(hdf5_attribute, "NSOFTCLASSES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);hdf5_attribute = my_H5Acreate(handle, "NTYPES" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 2;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "NTYPES");
my_H5Aclose(hdf5_attribute, "NTYPES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);hdf5_attribute = my_H5Acreate(handle, "PERIODIC" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "PERIODIC");
my_H5Aclose(hdf5_attribute, "PERIODIC");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);hdf5_attribute = my_H5Acreate(handle, "PMGRID" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1024;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "PMGRID");
my_H5Aclose(hdf5_attribute, "PMGRID");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);hdf5_attribute = my_H5Acreate(handle, "POSITIONS_IN_32BIT" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "POSITIONS_IN_32BIT");
my_H5Aclose(hdf5_attribute, "POSITIONS_IN_32BIT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);hdf5_attribute = my_H5Acreate(handle, "POWERSPEC_ON_OUTPUT" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "POWERSPEC_ON_OUTPUT");
my_H5Aclose(hdf5_attribute, "POWERSPEC_ON_OUTPUT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);hdf5_attribute = my_H5Acreate(handle, "SELFGRAVITY" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SELFGRAVITY");
my_H5Aclose(hdf5_attribute, "SELFGRAVITY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);hdf5_attribute = my_H5Acreate(handle, "THERMAL_VEL_IC" , atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "THERMAL_VEL_IC");
my_H5Aclose(hdf5_attribute, "THERMAL_VEL_IC");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

my_H5Tclose(atype);
}

