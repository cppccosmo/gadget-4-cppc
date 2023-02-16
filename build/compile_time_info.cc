#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gadgetconfig.h"
#include "data/dtypes.h"
#include "data/allvars.h"
#include "main/main.h"
void output_compile_time_options(void)
{
printf(
"    ADDITIONAL_GRID\n"
"    CB_PHASE\n"
"    DOUBLEPRECISION=2\n"
"    MFLR_RST\n"
"    MULTIPOLE_ORDER=1\n"
"    NGENIC=1024\n"
"    NSOFTCLASSES=2\n"
"    NTYPES=2\n"
"    PERIODIC\n"
"    PMGRID=1024\n"
"    POSITIONS_IN_32BIT\n"
"    POWERSPEC_ON_OUTPUT\n"
"    SELFGRAVITY\n"
"    THERMAL_VEL_IC\n"
);
}