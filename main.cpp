#include   "SinglePhase_CO2_TH_withoutWell.h"
#include  "SinglePhase_CO2_TH_withWell.h"
#include    "TestCapillaryTerm_v2.h"
#include    "run_IMPES_Iteration_TimeTerm_AnalyticalTest.h"
#include    "run_IMPES_Iteration_TwoPhase_BL_Numerical.h"
#include    "run_IMPES_Iteration_TwoPhase_withWell.h"
#include    "run_FC_P_IMPES_I_TwoPhase_withWell.h"
#include    "SinglePhase_TH_CO2_withwell_revised.h"
#include    "EDFM_withFracture_Geomtrycreate.h"
#include    "2D-EDFM-SinglePhase-CO2-HT.h"



int main()
{
//return  run_IMPES_Iteration_TwoPhase_WellCase();
    return run_FC_P_IMPES_I_TwoPhase_WellCase();
}
