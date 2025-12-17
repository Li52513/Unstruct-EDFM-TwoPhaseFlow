#include   "SinglePhase_CO2_TH_withoutWell.h"
#include  "SinglePhase_CO2_TH_withWell.h"
#include    "TestCapillaryTerm_v2.h"
#include    "run_IMPES_Iteration_TimeTerm_AnalyticalTest.h"
#include    "run_IMPES_Iteration_TwoPhase_BL_Numerical.h"
#include "run_IMPES_Iteration_TwoPhase_withWell.h"
#include "SinglePhase_TH_CO2_withwell_revised.h"



int main()
{
    return  run_IMPES_Iteration_TwoPhase_WellCase();
}

