/**
 * @file Test_ValidationSuite.h
 * @brief Systematic validation test suite: 8 physics scenarios x 6 topology/property variants.
 *
 * Naming convention:
 *   Val_T{N}_{Const|Var}_{NoFrac|SingleFrac|CrossFrac}()
 *
 *   T1 = no-well single-phase pressure diffusion
 *   T2 = no-well single-phase thermal diffusion
 *   T3 = no-well two-phase Buckley-Leverett
 *   T4 = no-well two-phase + heat (BL + thermal front)
 *   T5 = with-well single-phase Theis
 *   T6 = with-well single-phase thermal front
 *   T7 = with-well two-phase radial BL
 *   T8 = with-well two-phase radial BL + heat
 *
 * Const = ConstantWater / Corey  ->  analytical solution comparison
 * Var   = Water/CO2 EOS          ->  cross-compare with Const result
 * NoFrac / SingleFrac / CrossFrac ->  matrix-only / one diagonal frac / X fractures
 */
#pragma once

namespace ValidationSuite {

// === T1: no-well single-phase pressure diffusion ===
void Val_T1_Const_NoFrac();
void Val_T1_Const_SingleFrac();
void Val_T1_Const_CrossFrac();
void Val_T1_Var_NoFrac();
void Val_T1_Var_SingleFrac();
void Val_T1_Var_CrossFrac();

// === T2: no-well single-phase thermal diffusion ===
void Val_T2_Const_NoFrac();
void Val_T2_Const_SingleFrac();
void Val_T2_Const_CrossFrac();
void Val_T2_Var_NoFrac();
void Val_T2_Var_SingleFrac();
void Val_T2_Var_CrossFrac();

// === T3: no-well two-phase Buckley-Leverett ===
void Val_T3_Const_NoFrac();
void Val_T3_Const_SingleFrac();
void Val_T3_Const_CrossFrac();
void Val_T3_Var_NoFrac();
void Val_T3_Var_SingleFrac();
void Val_T3_Var_CrossFrac();

// === T4: no-well two-phase + heat ===
void Val_T4_Const_NoFrac();
void Val_T4_Const_SingleFrac();
void Val_T4_Const_CrossFrac();
void Val_T4_Var_NoFrac();
void Val_T4_Var_SingleFrac();
void Val_T4_Var_CrossFrac();

// === T5: with-well single-phase Theis ===
void Val_T5_Const_NoFrac();
void Val_T5_Const_SingleFrac();
void Val_T5_Const_CrossFrac();
void Val_T5_Var_NoFrac();
void Val_T5_Var_SingleFrac();
void Val_T5_Var_CrossFrac();

// === T6: with-well single-phase thermal front ===
void Val_T6_Const_NoFrac();
void Val_T6_Const_SingleFrac();
void Val_T6_Const_CrossFrac();
void Val_T6_Var_NoFrac();
void Val_T6_Var_SingleFrac();
void Val_T6_Var_CrossFrac();

// === T7: with-well two-phase radial BL ===
void Val_T7_Const_NoFrac();
void Val_T7_Const_SingleFrac();
void Val_T7_Const_CrossFrac();
void Val_T7_Var_NoFrac();
void Val_T7_Var_SingleFrac();
void Val_T7_Var_CrossFrac();

// === T8: with-well two-phase radial BL + heat ===
void Val_T8_Const_NoFrac();
void Val_T8_Const_SingleFrac();
void Val_T8_Const_CrossFrac();
void Val_T8_Var_NoFrac();
void Val_T8_Var_SingleFrac();
void Val_T8_Var_CrossFrac();

// === Grouped runners ===
void Run_All_T1();
void Run_All_T2();
void Run_All_T3();
void Run_All_T4();
void Run_All_T5();
void Run_All_T6();
void Run_All_T7();
void Run_All_T8();

void Run_ValidationSuite_All();

} // namespace ValidationSuite
