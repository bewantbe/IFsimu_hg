#include "stdafx.h"
#include "loop.h"
#include "neuron.h"
#include "raster.h"
#include "datainput.h"
#include "datahandling.h"
#include "poisson_input.h"

// Force Eigen3 to use MKL when icc is detected
#if defined(__ICC) || defined(__INTEL_COMPILER)
#  define EIGEN_USE_MKL_ALL
#endif
#define NDEBUG  // Turn off debug in Eigen3
// Seems that icc is not good at compiling template library (such as Eigen3)
// To have a good performance, you need to have Eigen version >= 3.2
#include <Eigen/Dense>

void (*ode_solver)(neuron*, double, double);

#if POISSON_INPUT_USE
#if SMOOTH_CONDUCTANCE_USE
// use smooth conductance HE and HI
void force_inputR( neuron *tempneu, int index_neuron)
{
  // excitatory conductance
//  tempneu[index_neuron].value[Stepsmooth_Con] += Strength_Exinput;
  tempneu[index_neuron].value[Stepsmooth_Con] += g_arr_poisson_strength_E[index_neuron];
  // inhibitory conductance
//  tempneu[index_neuron].value[2*Stepsmooth_Con] += Strength_Ininput;
  tempneu[index_neuron].value[2*Stepsmooth_Con] += g_arr_poisson_strength_I[index_neuron];
}



#else  // of SMOOTH_CONDUCTANCE_USE
// only use excitatory conductance and inhibitory conductance: GE and GI
void force_input( neuron *tempneu, int index_neuron)
{
  // excitatory conductance
//  tempneu[index_neuron].value[1] += Strength_Exinput*2.0/(Time_ExCon);
  tempneu[index_neuron].value[1] += g_arr_poisson_strength_E[index_neuron];
  // inhibitory conductance
//  tempneu[index_neuron].value[Stepsmooth_Con+1] += Strength_Ininput*5.0/(Time_InCon);
  tempneu[index_neuron].value[Stepsmooth_Con+1] += g_arr_poisson_strength_I[index_neuron];
}
#endif  // SMOOTH_CONDUCTANCE_USE

void voltage_dt(int index_neuron,double t,double m,double h,double n,double gE,double gI,double v,double &dv)
{
# if EXPONENTIAL_IF_USE
  dv = - Con_Leakage*(v - Vot_Leakage) - Con_sodium*m*m*m*h*(v-Vot_sodium) - Con_potassium*n*n*n*n*(v - Vot_potassium)
          - gE*(v - Vot_Excitatory) - gI*(v - Vot_Inhibitory) + Con_Leakage*VOT_DELTAT*exp((v-VOT_TAKEOFF)/VOT_DELTAT);
# else
  dv = - Con_Leakage*(v - Vot_Leakage) - Con_sodium*m*m*m*h*(v-Vot_sodium) - Con_potassium*n*n*n*n*(v - Vot_potassium)
          - gE*(v - Vot_Excitatory) - gI*(v - Vot_Inhibitory);
# endif
}

void whole_dt(const neuron * const neu_val, neuron *neu_dt, double *volt, int index_neuron, double t)
{
  //dv_dt
  const double * const &neu_i_val = neu_val[index_neuron].value;
  double *&neu_i_dt = neu_dt[index_neuron].value;
  double v  = neu_i_val[0];
  double gE = neu_i_val[1];
  double gI = neu_i_val[Stepsmooth_Con+1];
  double m  = neu_i_val[2*Stepsmooth_Con+1];
  double h  = neu_i_val[2*Stepsmooth_Con+2];
  double n  = neu_i_val[2*Stepsmooth_Con+3];

# if EXPONENTIAL_IF_USE
  neu_i_dt[0] = - Con_Leakage*(v - Vot_Leakage) - Con_sodium*m*m*m*h*(v-Vot_sodium) - Con_potassium*n*n*n*n*(v - Vot_potassium)
          - gE*(v - Vot_Excitatory) - gI*(v - Vot_Inhibitory) + Con_Leakage*VOT_DELTAT*exp((v-VOT_TAKEOFF)/VOT_DELTAT);
# else
  neu_i_dt[0] = - Con_Leakage*(v - Vot_Leakage) - Con_sodium*m*m*m*h*(v-Vot_sodium) - Con_potassium*n*n*n*n*(v - Vot_potassium)
          - gE*(v - Vot_Excitatory) - gI*(v - Vot_Inhibitory);
# endif

//mhn_dt
  {
    double e1 = exp(-v+2.5);
    double e2 = sqrt(e1*exp(-2.5));
    neu_i_dt[2*Stepsmooth_Con+1] = (v-2.5)/(1-e1)*(1-m) - 4*exp(-v/1.8)*m;
    neu_i_dt[2*Stepsmooth_Con+2] = 0.07*e2*(1-h) - 1/(1+e1*exp(0.5))*h;
    neu_i_dt[2*Stepsmooth_Con+3] = 0.1*(v-1.0)/(1-e1*exp(-1.5))*(1-n) - 0.125*sqrt(sqrt((e2)))*n;
  }

//conductance_dt
int index_firing = -1;
# if SMOOTH_CONDUCTANCE_USE
//smooth_dt
  for (int con_index = 1; con_index<Stepsmooth_Con; con_index++){
    neu_i_dt[con_index] = -neu_i_val[con_index]/Time_ExCon + neu_i_val[con_index+1];
    neu_i_dt[Stepsmooth_Con+con_index] = -neu_i_val[Stepsmooth_Con+con_index]/Time_InCon + neu_i_val[Stepsmooth_Con+con_index+1];
  }
  neu_i_dt[Stepsmooth_Con] = -neu_i_val[Stepsmooth_Con]/Time_ExConR ;
  neu_i_dt[2*Stepsmooth_Con] = -neu_i_val[2*Stepsmooth_Con]/Time_InConR ;

#if CORTICAL_STRENGTH_NONHOMO
  for (index_firing = 0; index_firing<g_num_neu; index_firing++) {
    if (index_firing != index_neuron) {
      if(volt[index_firing]== 0) continue;
      neu_i_dt[Stepsmooth_Con] +=
        (neu_val[index_neuron].type * Strength_CorEE
         + (1-neu_val[index_neuron].type) * Strength_CorIE)
        * neu_val[index_firing].type
        * cortical_matrix[index_neuron][index_firing]
        * volt[index_firing];

      neu_i_dt[2*Stepsmooth_Con] +=
        (neu_val[index_neuron].type * Strength_CorEI
         + (1-neu_val[index_neuron].type) * Strength_CorII)
        * (1-neu_val[index_firing].type)
        * cortical_matrix[index_neuron][index_firing]
        *volt[index_firing];
    }
  }
#else
  for (index_firing = 0; index_firing<g_num_neu; index_firing++) {
    if (index_firing != index_neuron) {
      if(volt[index_firing]== 0) continue;
      neu_i_dt[Stepsmooth_Con] +=
        neu_val[index_firing].type * neu_val[index_neuron].type    * Strength_CorEE
        *volt[index_firing]
        + neu_val[index_firing].type * (1-neu_val[index_neuron].type)* Strength_CorIE
        *volt[index_firing];

      neu_i_dt[2*Stepsmooth_Con] +=
        (1-neu_val[index_firing].type) * neu_val[index_neuron].type    * Strength_CorEI
        *volt[index_firing]
        + (1-neu_val[index_firing].type) * (1-neu_val[index_neuron].type)* Strength_CorII
        *volt[index_firing];
    }
  }
#endif

#else//non_smooth conductance
  neu_i_dt[1] = -neu_i_val[1]/Time_ExConR;
  neu_i_dt[Stepsmooth_Con+1] = -neu_i_val[Stepsmooth_Con+1]/Time_InConR;

#if CORTICAL_STRENGTH_NONHOMO
  for (index_firing = 0;index_firing<g_num_neu;index_firing++) {
    if (index_firing != index_neuron) {
      if(volt[index_firing]== 0) continue;
      neu_i_dt[1] +=
        neu_val[index_firing].type    * neu_val[index_neuron].type
        *Strength_CorEE * cortical_matrix[index_neuron][index_firing]
        *volt[index_firing]
        + neu_val[index_firing].type * (1-neu_val[index_neuron].type)
        *Strength_CorIE * cortical_matrix[index_neuron][index_firing]
        *volt[index_firing];

      neu_i_dt[Stepsmooth_Con+1] +=
        (1-neu_val[index_firing].type)    * neu_val[index_neuron].type
        *Strength_CorEI * cortical_matrix[index_neuron][index_firing]
        *volt[index_firing]
        + (1-neu_val[index_firing].type) * (1-neu_val[index_neuron].type)
        *Strength_CorII * cortical_matrix[index_neuron][index_firing]
        *volt[index_firing];
    }
  }
#else
for (index_firing = 0;index_firing<g_num_neu;index_firing++) {
  if (index_firing != index_neuron) {
    if(volt[index_firing]== 0) continue;
    neu_i_dt[1] +=
      neu_val[index_firing].type * neu_val[index_neuron].type    * Strength_CorEE
      *volt[index_firing]
      + neu_val[index_firing].type * (1-neu_val[index_neuron].type)* Strength_CorIE
      *volt[index_firing];

    neu_i_dt[Stepsmooth_Con+1] +=
      (1-neu_val[index_firing].type) * neu_val[index_neuron].type    * Strength_CorEI
      *volt[index_firing]
      + (1-neu_val[index_firing].type) * (1-neu_val[index_neuron].type)* Strength_CorII
      *volt[index_firing];
  }
}
#endif

#endif//nonsmooth end
}

#if SMOOTH_CONDUCTANCE_USE && CORTICAL_STRENGTH_NONHOMO && !EXPONENTIAL_IF_USE
void whole_dt_single(const double * const &neu_i_val, double *&neu_i_dy, double t)
{
  //dv_dt
  double v  = neu_i_val[0];
  double gE = neu_i_val[1];
  double gI = neu_i_val[Stepsmooth_Con+1];
  double m  = neu_i_val[2*Stepsmooth_Con+1];
  double h  = neu_i_val[2*Stepsmooth_Con+2];
  double n  = neu_i_val[2*Stepsmooth_Con+3];

  neu_i_dy[0] = - Con_Leakage*(v-Vot_Leakage) - Con_sodium*m*m*m*h*(v-Vot_sodium) - Con_potassium*n*n*n*n*(v-Vot_potassium)
                - gE*(v-Vot_Excitatory) - gI*(v-Vot_Inhibitory);

//mhn_dt
  {
    double e1 = exp(-v+2.5);
    double e2 = sqrt(e1*exp(-2.5));
    neu_i_dy[2*Stepsmooth_Con+1] = (v-2.5)/(1-e1)*(1-m) - 4*exp(-v/1.8)*m;
    neu_i_dy[2*Stepsmooth_Con+2] = 0.07*e2*(1-h) - h/(1+e1*exp(0.5));
    neu_i_dy[2*Stepsmooth_Con+3] = 0.1*(v-1.0)/(1-e1*exp(-1.5))*(1-n) - 0.125*sqrt(sqrt((e2)))*n;
  }

  //conductance_dt
  for (int con_index = 1; con_index<Stepsmooth_Con; con_index++){
    neu_i_dy[               con_index] = -neu_i_val[               con_index]*(1.0/Time_ExCon) + neu_i_val[con_index+1];
    neu_i_dy[Stepsmooth_Con+con_index] = -neu_i_val[Stepsmooth_Con+con_index]*(1.0/Time_InCon) + neu_i_val[Stepsmooth_Con+con_index+1];
  }
  neu_i_dy[  Stepsmooth_Con] = -neu_i_val[  Stepsmooth_Con]*(1.0/Time_ExConR);
  neu_i_dy[2*Stepsmooth_Con] = -neu_i_val[2*Stepsmooth_Con]*(1.0/Time_InConR);
}

void whole_dt_vector(const neuron * const neu_val, neuron *neu_dy, double *volt, double t)
{
  for (int i = 0; i < g_num_neu; i++) {
    volt[i] = (neu_val[i].value[0]>4 ? 1/(1+exp(-5*(neu_val[i].value[0] - 8.5))) : 0.0);
  }
  for (int i = 0; i < g_num_neu; i++) {
    whole_dt_single(neu_val[i].value, neu_dy[i].value, t);
  }

  // Compute influence from the network (H^E and H^I)
  for (int i = 0; i < g_num_neu_ex; i++) {
    if(volt[i] == 0) continue;
    for (int j = 0; j < g_num_neu_ex; j++) {
      neu_dy[j].value[Stepsmooth_Con] +=
        Strength_CorEE * cortical_matrix[j][i] * volt[i];
    }
    for (int j = g_num_neu_ex; j < g_num_neu; j++) {
      neu_dy[j].value[Stepsmooth_Con] +=
        Strength_CorIE * cortical_matrix[j][i] * volt[i];
    }
  }
  for (int i = g_num_neu_ex; i < g_num_neu; i++) {
    if(volt[i] == 0) continue;
    for (int j = 0; j < g_num_neu_ex; j++) {
      neu_dy[j].value[2*Stepsmooth_Con] +=
        Strength_CorEI * cortical_matrix[j][i] * volt[i];
    }
    for (int j = g_num_neu_ex; j < g_num_neu; j++) {
      neu_dy[j].value[2*Stepsmooth_Con] +=
        Strength_CorII * cortical_matrix[j][i] * volt[i];
    }
  }
}

void get_dx(Eigen::ArrayXXd &dx, const Eigen::ArrayXXd &xm, double t)
{
  const Eigen::ArrayXd &v  = xm.col(0);
  const Eigen::ArrayXd &gE = xm.col(1);
  const Eigen::ArrayXd &hE = xm.col(2);
  const Eigen::ArrayXd &gI = xm.col(3);
  const Eigen::ArrayXd &hI = xm.col(4);
  const Eigen::ArrayXd &m  = xm.col(5);
  const Eigen::ArrayXd &h  = xm.col(6);
  const Eigen::ArrayXd &n  = xm.col(7);

  // cost 12%
  dx.col(0) = -Con_Leakage * (v - Vot_Leakage)
              -Con_sodium * m*m*m*h * (v-Vot_sodium)
              -Con_potassium * n*n*n*n * (v-Vot_potassium)
              -gE * (v-Vot_Excitatory) - gI * (v-Vot_Inhibitory);

  // cost 42%
  // mhn
  Eigen::ArrayXd e1(g_num_neu);
  Eigen::ArrayXd e2(g_num_neu);
  e1 = 2.5-v;
  e1 = exp(e1);
  e2 = v*(-1/1.8);
  e2 = exp(e2);
  dx.col(5) = (v-2.5)/(1-e1)*(1-m) - 4*e2*m;
  e2 = e1*exp(-2.5);
  e2 = sqrt(e2);
  dx.col(6) = 0.07*e2*(1-h) - h/(1+e1*exp(0.5));
  e2 = sqrt(e2);
  e2 = sqrt(e2);
  dx.col(7) = 0.1*(v-1.0)/(1-e1*exp(-1.5))*(1-n) - 0.125*e2*n;

  //conductance_dt
  dx.col(1) = -gE * (1.0/Time_ExCon) + hE;
  dx.col(3) = -gI * (1.0/Time_InCon) + hI;
  dx.col(2) = -hE * (1.0/Time_ExConR);
  dx.col(4) = -hI * (1.0/Time_InConR);

  // cost 7%
  // Compute influence from the network (H^E and H^I)
  double synaptic_volt;
  for (int i = 0; i < g_num_neu_ex; i++) {
    synaptic_volt = xm(i,0)>4 ? 1/(1+exp(-5*(xm(i,0) - 8.5))) : 0.0;
    if(synaptic_volt == 0) continue;
    for (int j = 0; j < g_num_neu_ex; j++) {
      dx(j, 2) += Strength_CorEE * cortical_matrix[j][i] * synaptic_volt;
    }
    for (int j = g_num_neu_ex; j < g_num_neu; j++) {
      dx(j, 2) += Strength_CorIE * cortical_matrix[j][i] * synaptic_volt;
    }
  }
  for (int i = g_num_neu_ex; i < g_num_neu; i++) {
    synaptic_volt = xm(i,0)>4 ? 1/(1+exp(-5*(xm(i,0) - 8.5))) : 0.0;
    if(synaptic_volt == 0) continue;
    for (int j = 0; j < g_num_neu_ex; j++) {
      dx(j, 4) += Strength_CorEI * cortical_matrix[j][i] * synaptic_volt;
    }
    for (int j = g_num_neu_ex; j < g_num_neu; j++) {
      dx(j, 4) += Strength_CorII * cortical_matrix[j][i] * synaptic_volt;
    }
  }
}
#endif // if SMOOTH_CONDUCTANCE_USE && CORTICAL_STRENGTH_NONHOMO && !EXPONENTIAL_IF_USE

#else   // not POISSON_INPUT_USE

double external_current(int index_neuron, double t)
{
  return Current_0;//(Current_0 + Current_1*cos(2*M_PI*Rate_input*t+phase[index_neuron]));
}

void voltage_dt(int index_neuron,double t,double m,double h,double n,double gE,double gI,double v,double &dv)
{
# if EXPONENTIAL_IF_USE
  dv = - Con_Leakage*(v - Vot_Leakage) - gE*(v - Vot_Excitatory) - gI*(v - Vot_Inhibitory) + Con_Leakage*VOT_DELTAT*exp((v-VOT_TAKEOFF)/VOT_DELTAT)
          + external_current(index_neuron, t);
# else
  dv = - Con_Leakage*(v - Vot_Leakage) - Con_sodium*m*m*m*h*(v-Vot_sodium) - Con_potassium*n*n*n*n*(v - Vot_potassium)
       - gE*(v - Vot_Excitatory) - gI*(v - Vot_Inhibitory) + external_current(index_neuron, t);
#endif
}

void whole_dt(const neuron *neu_val, neuron *neu_dt, double *volt, int index_neuron, double t)
{
//dv_dt
double v = neu_val[index_neuron].value[0];
double sca_v = neu_val[index_neuron].value[0]*10.0-65.0 ;//rescale the voltage
double gE = neu_val[index_neuron].value[1];
double gI = neu_val[index_neuron].value[Stepsmooth_Con+1];
double m = neu_val[index_neuron].value[2*Stepsmooth_Con+1];
double h = neu_val[index_neuron].value[2*Stepsmooth_Con+2];
double n = neu_val[index_neuron].value[2*Stepsmooth_Con+3];
#if EXPONENTIAL_IF_USE
  neu_dt[index_neuron].value[0] = - Con_Leakage*(v - Vot_Leakage) - gE*(v - Vot_Excitatory) - gI*(v - Vot_Inhibitory) + Con_Leakage*VOT_DELTAT*exp((v-VOT_TAKEOFF)/VOT_DELTAT)
          + external_current(index_neuron, t);
#else
  neu_dt[index_neuron].value[0] = - Con_Leakage*(v - Vot_Leakage) - Con_sodium*m*m*m*h*(v-Vot_sodium) - Con_potassium*n*n*n*n*(v - Vot_potassium)
        - gE*(v - Vot_Excitatory) - gI*(v - Vot_Inhibitory) + external_current(index_neuron, t);
#endif

//mhn_dt
  neu_dt[index_neuron].value[2*Stepsmooth_Con+1] = 0.1*(sca_v+40.0)/(1-exp(-(sca_v+40.0)/10.0))*(1-m) - 4*exp(-(sca_v+65.0)/18.0)*m;

  neu_dt[index_neuron].value[2*Stepsmooth_Con+2] = 0.07*exp(-(sca_v+65.0)/20.0)*(1-h) - exp((sca_v+35.0)/10)/(1+exp((sca_v+35.0)/10))*h;

  neu_dt[index_neuron].value[2*Stepsmooth_Con+3] = 0.01*(sca_v+55.0)/(1-exp(-(sca_v+55.0)/10.0))*(1-n) - 0.125*exp(-(sca_v+65.0)/80.0)*n;

//conductance_dt
int index_firing = -1;
# if SMOOTH_CONDUCTANCE_USE
//smooth_dt
 for (int con_index = 1;con_index<Stepsmooth_Con;con_index++){
        neu_dt[index_neuron].value[con_index] = -neu_val[index_neuron].value[con_index]/Time_ExCon + neu_val[index_neuron].value[con_index+1];
        neu_dt[index_neuron].value[Stepsmooth_Con+con_index] = -neu_val[index_neuron].value[Stepsmooth_Con+con_index]/Time_InCon + neu_val[index_neuron].value[Stepsmooth_Con+con_index+1];
    }
  neu_dt[index_neuron].value[Stepsmooth_Con] = -neu_val[index_neuron].value[Stepsmooth_Con]/Time_ExConR + external_current(index_neuron, t);
  neu_dt[index_neuron].value[2*Stepsmooth_Con] = -neu_val[index_neuron].value[2*Stepsmooth_Con]/Time_InConR ;

#if CORTICAL_STRENGTH_NONHOMO
  for (index_firing = 0; index_firing<g_num_neu; index_firing++) {
    if (index_firing != index_neuron) {
      if(volt[index_firing]== 0) continue;
      neu_dt[index_neuron].value[Stepsmooth_Con] +=
        neuRK[index_firing].type    * neuRK[index_neuron].type
        *Strength_CorEE * cortical_matrix[index_neuron][index_firing]
        *volt[index_firing]
        + neuRK[index_firing].type * (1-neuRK[index_neuron].type)
        *Strength_CorIE * cortical_matrix[index_neuron][index_firing]
        *volt[index_firing];

      neu_dt[index_neuron].value[2*Stepsmooth_Con] +=
        (1-neuRK[index_firing].type)    * neuRK[index_neuron].type
        *Strength_CorEI * cortical_matrix[index_neuron][index_firing]
        *volt[index_firing]
        + (1-neuRK[index_firing].type) * (1-neuRK[index_neuron].type)
        *Strength_CorII * cortical_matrix[index_neuron][index_firing]
        *volt[index_firing];
    }
  }
#else
  for (index_firing = 0; index_firing<g_num_neu; index_firing++) {
    if (index_firing != index_neuron) {
      if(volt[index_firing]== 0) continue;
      neu_dt[index_neuron].value[Stepsmooth_Con] +=
        neuRK[index_firing].type * neuRK[index_neuron].type    * Strength_CorEE
        *volt[index_firing]
        + neuRK[index_firing].type * (1-neuRK[index_neuron].type)* Strength_CorIE
        *volt[index_firing];

      neu_dt[index_neuron].value[2*Stepsmooth_Con] +=
        (1-neuRK[index_firing].type) * neuRK[index_neuron].type    * Strength_CorEI
        *volt[index_firing]
        + (1-neuRK[index_firing].type) * (1-neuRK[index_neuron].type)* Strength_CorII
        *volt[index_firing];
    }
  }
#endif

#else//non_smooth conductance
  neu_dt[index_neuron].value[1] = -neu_val[index_neuron].value[1]/Time_ExConR+external_current(index_neuron, t);
  neu_dt[index_neuron].value[Stepsmooth_Con+1] = -neu_val[index_neuron].value[Stepsmooth_Con+1]/Time_InConR;

#if CORTICAL_STRENGTH_NONHOMO
  for (index_firing = 0;index_firing<g_num_neu;index_firing++) {
    if (index_firing != index_neuron) {
      if(volt[index_firing]== 0) continue;
      neu_dt[index_neuron].value[1] +=
        neu_val[index_firing].type    * neu_val[index_neuron].type
        *Strength_CorEE * cortical_matrix[index_neuron][index_firing]
        *volt[index_firing]
        + neu_val[index_firing].type * (1-neu_val[index_neuron].type)
        *Strength_CorIE * cortical_matrix[index_neuron][index_firing]
        *volt[index_firing];

      neu_dt[index_neuron].value[Stepsmooth_Con+1] +=
        (1-neu_val[index_firing].type)    * neu_val[index_neuron].type
        *Strength_CorEI * cortical_matrix[index_neuron][index_firing]
        *volt[index_firing]
        + (1-neu_val[index_firing].type) * (1-neu_val[index_neuron].type)
        *Strength_CorII * cortical_matrix[index_neuron][index_firing]
        *volt[index_firing];
    }
  }
#else
for (index_firing = 0;index_firing<g_num_neu;index_firing++) {
  if (index_firing != index_neuron) {
    if(volt[index_firing]== 0) continue;
    neu_dt[index_neuron].value[1] +=
      neu_val[index_firing].type * neu_val[index_neuron].type    * Strength_CorEE
      *volt[index_firing]
      + neu_val[index_firing].type * (1-neu_val[index_neuron].type)* Strength_CorIE
      *volt[index_firing];

    neu_dt[index_neuron].value[Stepsmooth_Con+1] +=
      (1-neu_val[index_firing].type) * neu_val[index_neuron].type    * Strength_CorEI
      *volt[index_firing]
      + (1-neu_val[index_firing].type) * (1-neu_val[index_neuron].type)* Strength_CorII
      *volt[index_firing];
  }
}
#endif

#endif//nonsmooth ends
}

#endif
// deal with the function 1/(1+exp(-(volt1+volt2)/deno))

///note by jyl: This function is supposed to increase the speed of simulation in the exponential part
///using the hermit interpolation, though not applied. If interested, have a try! :)
double hermexp(double volt1, double volt2, double deno)
{
    double f = 1/(1+exp(-(volt1+volt2)/deno));
    /*double der0 =0.25;
    double g0 = 0.5;
    double v = (volt1+volt2)/deno;
    double f;
    if (v<xleft)
    f = maximum((v-xleft)*dgdt+gleft,0.0);
    else{
        if (v>-xleft)
        f = minimum((v+xleft)*dgdt+1-gleft,1.0);
        else{
                    // use gleft, gright, dleft, dright to construct a cubic polynomial
                    // basic function f1 satisfies: f1(left)=1, f1(0)=0, f1'(left)=0, f1'(0)=0
                    double f1 = gleft*(2*v+0-3*xleft)*(v-0)*(v-0)/(0-xleft)/(0-xleft)/(0-xleft);
                    // basic function f2 satisfies: f2(left)=0, f2(0)=1, f2'(left)=0, f2'(0)=0
                    double f2 = g0*(3*0-2*v-xleft)*(v-xleft)*(v-xleft)/(0-xleft)/(0-xleft)/(0-xleft);
                    // basic function f3 satisfies: f3(a)=0, f3(b)=0, f3'(a)=1, f3'(b)=0
                    double f3 = dgdt*(v-xleft)*(v-0)*(v-0)/(0-xleft)/(0-xleft);
                    // basic function f4 satisfies: f4(a)=0, f4(b)=0, f4'(a)=0, f4'(b)=1
                    double f4 = der0*(v-xleft)*(v-xleft)*(v-0)/(0-xleft)/(0-xleft);

                    f=f1+f2+f3+f4;
                if (v > 0)
                    f = 1-f;

        }
    }*/
    return f;
}




#if runge_kutta == runge_kutta2
void runge_kutta2(neuron *tempneu, double subTstep, double t_evolution)
{
        // dy/dt = f(t,y) t = t(n), y = y(n)
  double ini_time = t_evolution;
  double mid_time = t_evolution + subTstep/2;

  int i = -1;
  int j = -1;
  //initialize the information of the network
  for (i = 0; i<g_num_neu; i++) {
    for (j = 0; j<2*Stepsmooth_Con+4; j++) {
      neuRK[i].value[j]  = tempneu[i].value[j];
    }
    vol[i]=(neuRK[i].value[0]>4?(exp((tempneu[i].value[0]-8.5)*5)/(1+exp((tempneu[i].value[0]-8.5)*5))):0);
  } ///jyl: I set this value(4) to increase the computing speed. Once the voltage is larger than this value, the simulation program
  ///would take into account. Note that the error caused by this technique is smaller than the accuracy of RK method

  //compute the dy1 = f(t(n), y(n)), here y,f,dy1 are 2*stepsmooth+4-dimension vectors
  for (i = 0; i<g_num_neu; i++) {
    whole_dt(neuRK, neu_d1, vol, i, ini_time);
    //update the information of the network at midtime, which will be used in the following RK2 method
    for (j = 0; j<2*Stepsmooth_Con+4; j++)
      neuRK[i].value[j]  = tempneu[i].value[j]  + neu_d1[i].value[j]  * subTstep/2;
  }
  for (i = 0; i<g_num_neu; i++)
    vol[i]=(neuRK[i].value[0]>4?exp((tempneu[i].value[0] + neu_d1[i].value[0]  * subTstep/2-8.5)*5)/(1+exp((tempneu[i].value[0] + neu_d1[i].value[0]  * subTstep/2-8.5)*5)):0);


  // dy2 = f(t(n)+h/3, y(n)+dy1*h/2)
  for (i = 0; i<g_num_neu; i++)
    whole_dt(neuRK, neu_d2, vol, i, mid_time);


//update the final data
//dy = dy2
  for (i = 0; i<g_num_neu; i++) {
    for (j=0; j<2*Stepsmooth_Con+4; j++) {
      tempneu[i].value[j] = tempneu[i].value[j] + neu_d2[i].value[j] *subTstep;
    }
  }
}
#endif

#if runge_kutta == runge_kutta3
void runge_kutta3(neuron *tempneu, double subTstep, double t_evolution)
{
    // dy/dt = f(t,y) t = t(n), y = y(n)
  double ini_time = t_evolution;
  double one_third_time = t_evolution + subTstep/3;
  double two_third_time = t_evolution + 2*subTstep/3;

  int i = -1;
  int j = -1;
  //initialize the information of the network
  for (i = 0; i<g_num_neu; i++) {
    for (j = 0; j<2*Stepsmooth_Con+4; j++) {
      neuRK[i].value[j]  = tempneu[i].value[j];
    }
    vol[i]=(neuRK[i].value[0]>4?(exp((tempneu[i].value[0]-8.5)*5)/(1+exp((tempneu[i].value[0]-8.5)*5))):0);
  }

  //compute the dy1 = f(t(n), y(n)), here y,f,dy1 are 2*stepsmooth+4-dimension vectors
  for (i = 0; i<g_num_neu; i++) {
    whole_dt(neuRK, neu_d1, vol, i, ini_time);
    //update the information of the network at midtime, which will be used in the following RK2 method
    for (j = 0; j<2*Stepsmooth_Con+4; j++)
      neuRK[i].value[j]  = tempneu[i].value[j]  + neu_d1[i].value[j]  * subTstep/3;
  }
  for (i = 0; i<g_num_neu; i++)
    vol[i]=(neuRK[i].value[0]>4?exp((tempneu[i].value[0] + neu_d1[i].value[0]  * subTstep/2-8.5)*5)/(1+exp((tempneu[i].value[0] + neu_d1[i].value[0]  * subTstep/2-8.5)*5)):0);


  // dy2 = f(t(n)+h/3, y(n)+dy1*h/3)
  for (i = 0; i<g_num_neu; i++) {
    whole_dt(neuRK, neu_d2, vol, i, one_third_time);
    //update the information of the network at midtime, which will be used in the following RK4 method
    for (j = 0; j<2*Stepsmooth_Con+4; j++)
      neuRK[i].value[j]  = tempneu[i].value[j]  + neu_d2[i].value[j]  * subTstep*2/3;
  }
  for (i = 0; i<g_num_neu; i++){
    vol[i]=(neuRK[i].value[0]>4?exp((tempneu[i].value[0] + neu_d2[i].value[0]  * subTstep/2-8.5)*5)/(1+exp((tempneu[i].value[0] + neu_d2[i].value[0]  * subTstep/2-8.5)*5)):0);
  }


//dy3 = f(t(n)+h/2, y(n)+dy2*h*2/3)
  for (i = 0; i<g_num_neu; i++)
    whole_dt(neuRK, neu_d3, vol, i, two_third_time);


//update the final data
//dy = (dy1+3*dy3)/4
  for (i = 0; i<g_num_neu; i++) {
    for (j=0; j<2*Stepsmooth_Con+4; j++) {
      tempneu[i].value[j] = tempneu[i].value[j] + (neu_d1[i].value[j] + 3*neu_d3[i].value[j]) *subTstep/4;
    }
  }
}
#endif

void runge_kutta4(neuron *tempneu, double subTstep, double t_evolution)
{
  int i, j;

  // dy1 = f(t(n), y(n)), here y,f,dy1 are 2*stepsmooth+4-dimension vectors
  for (i = 0; i<g_num_neu; i++) {
    for (j = 0; j<2*Stepsmooth_Con+4; j++) {
      neuRK[i].value[j]  = tempneu[i].value[j];
    }
  }
  whole_dt_vector(neuRK, neu_d1, vol, t_evolution);

  // dy2 = f(t(n)+h/2, y(n)+dy1*h/2)
  for (i = 0; i<g_num_neu; i++) {
    for (j = 0; j<2*Stepsmooth_Con+4; j++) {
      neuRK[i].value[j]  = tempneu[i].value[j] + neu_d1[i].value[j] * (subTstep/2);
    }
  }
  whole_dt_vector(neuRK, neu_d2, vol, t_evolution + subTstep/2);

  // dy3 = f(t(n)+h/2, y(n)+dy2*h/2)
  for (i = 0; i<g_num_neu; i++){
    for (j = 0; j<2*Stepsmooth_Con+4; j++) {
      neuRK[i].value[j]  = tempneu[i].value[j] + neu_d2[i].value[j] * (subTstep/2);
    }
  }
  whole_dt_vector(neuRK, neu_d3, vol, t_evolution + subTstep/2);

  // dy4 = f(t(n)+h, y(n)+dy3*h)
  for (i = 0; i<g_num_neu; i++) {
    for (j = 0; j<2*Stepsmooth_Con+4; j++) {
      neuRK[i].value[j]  = tempneu[i].value[j] + neu_d3[i].value[j] * subTstep;
    }
  }
  whole_dt_vector(neuRK, neu_d4, vol, t_evolution + subTstep);

  // update the final data
  // dy = (dy1+2*dy2+2*dy3+dy4)/6
  for (i = 0; i<g_num_neu; i++) {
    for (j=0; j<2*Stepsmooth_Con+4; j++) {
      tempneu[i].value[j] = tempneu[i].value[j] + (neu_d1[i].value[j] + 2*neu_d2[i].value[j] + 2*neu_d3[i].value[j] + neu_d4[i].value[j]) *subTstep/6;
    }
  }
}

void runge_kutta4_vec(neuron *tempneu, double subTstep, double t_evolution)
{
  int i, j;
  Eigen::ArrayXXd xt(g_num_neu, 2*Stepsmooth_Con+4);  // all variables of all neurons
  Eigen::ArrayXXd k1(g_num_neu, 2*Stepsmooth_Con+4);
  Eigen::ArrayXXd k2(g_num_neu, 2*Stepsmooth_Con+4);
  Eigen::ArrayXXd k3(g_num_neu, 2*Stepsmooth_Con+4);
  Eigen::ArrayXXd k4(g_num_neu, 2*Stepsmooth_Con+4);

  for (i = 0; i<g_num_neu; i++) {
    for (j = 0; j<2*Stepsmooth_Con+4; j++) {
      xt(i,j) = tempneu[i].value[j];
    }
  }

  get_dx(k1, xt, t_evolution);

  get_dx(k2, xt + k1 * (subTstep/2), t_evolution + subTstep/2);

  get_dx(k3, xt + k2 * (subTstep/2), t_evolution + subTstep/2);

  get_dx(k4, xt + k3 * subTstep, t_evolution + subTstep);

  xt += (k1 + 2*k2 + 2*k3 + k4) * (subTstep/6);

  for (i = 0; i<g_num_neu; i++) {
    for (j = 0; j<2*Stepsmooth_Con+4; j++) {
      tempneu[i].value[j] = xt(i,j);
    }
  }
}


///modified 2013/1/13 12:51

void hermit(double a, double b, double va, double vb,
            double dva, double dvb, double x, double &fx)
{
  // use v(a), v(b), dv(a)/dt, dv(b)/dt to construct a cubic polynomial
  // basic function f1 satisfies: f1(a)=1, f1(b)=0, f1'(a)=0, f1'(b)=0
  double f1 = va*(2*x+b-3*a)*(x-b)*(x-b)/(b-a)/(b-a)/(b-a);
  // basic function f2 satisfies: f2(a)=0, f2(b)=1, f2'(a)=0, f2'(b)=0
  double f2 = vb*(3*b-2*x-a)*(x-a)*(x-a)/(b-a)/(b-a)/(b-a);
  // basic function f3 satisfies: f3(a)=0, f3(b)=0, f3'(a)=1, f3'(b)=0
  double f3 = dva*(x-a)*(x-b)*(x-b)/(b-a)/(b-a);
  // basic function f4 satisfies: f4(a)=0, f4(b)=0, f4'(a)=0, f4'(b)=1
  double f4 = dvb*(x-a)*(x-a)*(x-b)/(b-a)/(b-a);

  // the polynomial of v(x) - Vot_Threshold
  fx = f1 + f2 + f3 + f4 - Vot_Threshold;
}

double root_search(void (*func)(double a, double b, double va, double vb,
                                double dva, double dvb, double x, double &fx),
                   double x1, double x2,
                   double fx1, double fx2,
                   double dfx1, double dfx2, double xacc)
{
  int j;
  double tempx1,tempx2,dx,f,fmid,xmid,root;

  // for firing time case, fx1<0, fmid>0
  (*func)(x1,x2,fx1,fx2,dfx1,dfx2,x1,f);
  (*func)(x1,x2,fx1,fx2,dfx1,dfx2,x2,fmid);
  /**************************************
    if (fabs(x2-x1)<xacc)
    {
      return x1;
    }
  ***************************************/
  if (f*fmid > 0) {
    P_ERR("The neuron does not fire at the end point?!\n");
    fprintf(stderr, "voltage difference at the beginning time and ending time: %g;%g\n",f,fmid);
    return x1;
  }

  tempx1 = x1;
  tempx2 = x2;
  for (j=0; j<Maxnum_search; j++) {
    dx = tempx2 - tempx1;
    xmid = tempx1 + dx/2;
    (*func)(x1,x2,fx1,fx2,dfx1,dfx2,xmid,fmid);
    if (fmid <= 0.0)
      tempx1 = xmid;
    else
      tempx2 = xmid;
    // the interval is small enough or already find the root
    if (fabs(fmid)<xacc) {
      root = xmid;
      return root;
    }
  }
  static int cnt_bad;
  cnt_bad++;
  P_ERR("Too many bisections in root searching!\n");
  // at most ONE error per 1000ms per neuron in total
  if (cnt_bad > 0.001*g_comp_time*g_num_neu) {
    fprintf(stderr,"There are over 1 warning per second per neuron (too many bisections in root searching)!\n");
    fprintf(stderr,"Program terminated.\n");
    exit(3);
  }
  return xmid;
}

void spiking_time(neuron *bef_neu, neuron *cur_neu, int index_neuron, double t_evolution, double Ta,
                  double Tb, double &firing_time)
{

  double va  = bef_neu[index_neuron].value[0];
  double dva = 0; // initialization
  double gEa = bef_neu[index_neuron].value[1];
  double gIa = bef_neu[index_neuron].value[Stepsmooth_Con+1];
  double ma  = bef_neu[index_neuron].value[2*Stepsmooth_Con+1];
  double ha  = bef_neu[index_neuron].value[2*Stepsmooth_Con+2];
  double na  = bef_neu[index_neuron].value[2*Stepsmooth_Con+3];

  voltage_dt(index_neuron,t_evolution,ma,ha,na,gEa,gIa,va,dva);


  double vb  = cur_neu[index_neuron].value[0];
  double dvb = 0; // initialization
  double gEb = cur_neu[index_neuron].value[1];
  double gIb = cur_neu[index_neuron].value[Stepsmooth_Con+1];
  double mb  = cur_neu[index_neuron].value[2*Stepsmooth_Con+1];
  double hb  = cur_neu[index_neuron].value[2*Stepsmooth_Con+2];
  double nb  = cur_neu[index_neuron].value[2*Stepsmooth_Con+3];

  voltage_dt(index_neuron,(t_evolution+Tb-Ta),mb,hb,nb,gEb,gIb,vb,dvb);

  firing_time = root_search(hermit, Ta, Tb, va, vb, dva, dvb, root_acc);
}


#if POISSON_INPUT_USE
void next_poisson_spike( double begin_time,
                         int *tempbegin_poisson_index, double &end_time, int &next_neu_index)
{
  double temp_end_time = 2*Tstep;//initialization
  double com_poisson_time = 0.0;
  int poisson_index = -1;
  int temp_next_neu_index = -1;

  for (int i=0; i < g_num_neu; i++)
  {
    poisson_index = tempbegin_poisson_index[i];
    if (poisson_index == -1) continue;//this means the current calculating neuron has no poisson spikes in the Tstep any more
    com_poisson_time = poisson_input[i].vect_value[poisson_index];
    if (temp_end_time >= com_poisson_time) {
      temp_end_time = com_poisson_time;
      temp_next_neu_index = i;
    }
  }
  end_time = temp_end_time;
  next_neu_index = temp_next_neu_index;

  return;
}
#endif


// this function is used to compute all the information of the whole network in a smooth subTstep [begin_time,end_time]
void sub_network_evolve(neuron *tempneu,  double begin_time, double end_time)
{
  //copy the information of all the neurons
  int index_neuron;
  for (index_neuron = 0; index_neuron<g_num_neu; index_neuron++) {
    neuron_copy(former_neu[index_neuron],tempneu[index_neuron]);
  }
  double subTstep = end_time - begin_time;
  double t_evolution = time_evolution + begin_time;
  double firing_time = -1.0;
//evovle the system from t_evolution+begin
  ode_solver(tempneu, subTstep, t_evolution);

//search the spiking time of the index_neuron
  for (index_neuron = 0; index_neuron<g_num_neu; index_neuron++) {
    if(former_neu[index_neuron].value[0]<Vot_Threshold && tempneu[index_neuron].value[0] >= Vot_Threshold) {
      spiking_time(former_neu,tempneu, index_neuron, time_evolution, begin_time,
                   end_time, firing_time);
      tempneu[index_neuron].state_neuron = STATE_SPIKING;
      //insert the spike
      raster_insert(spike_list,index_neuron,firing_time+time_evolution);
    }
    if(former_neu[index_neuron].value[0] >= Vot_Threshold && tempneu[index_neuron].value[0]<Vot_Threshold) {
      tempneu[index_neuron].state_neuron = STATE_NORMAL;
    }
  }
}

void network_initialization()
{
  // initialize the spike_list of the current computing time step
  raster_allocate(spike_list, RASTER_SPIKE_SIZE);

  //******************** poisson input spikes generation
#if POISSON_INPUT_USE
  int i;
  for (i=0; i<g_num_neu; i++) {
    poisson_generator(i, poisson_input[i]);
    if (poisson_input[i].vect_size<0)
      g_begin_poisson_index[i] = -1;
    else
      g_begin_poisson_index[i] = 0;
  }

#endif
}

// evolve the whole neuron system from t(0) to t(0)+Tstep
void network_evolution()
{
  int j;
  double begin_time = 0;
  double end_time = 0;


#if POISSON_INPUT_USE

  int next_neu_index = -1;//initialization
  int *&begin_poisson_index = g_begin_poisson_index;

  do {
    next_poisson_spike(begin_time,
                       begin_poisson_index, end_time, next_neu_index);

    if (g_fp_save_poisson_events!=NULL && next_neu_index != -1) {
      fprintf(g_fp_save_poisson_events, "%d\t%.17e\n", next_neu_index, time_evolution+end_time);
    }

    if(end_time>Tstep) end_time = Tstep;

    //renew the begin index of the poisson spike
    if (next_neu_index !=-1) {
      j = begin_poisson_index[next_neu_index];
      if (j!=-1 && j<poisson_input[next_neu_index].vect_size-1) {
        begin_poisson_index[next_neu_index]++;
      } else {
        begin_poisson_index[next_neu_index] = -1;
      }
    }
    if (begin_time<end_time)
      // this judgement is used for the case that "end_time = begin_time = 0" which means the next poisson spike
      // just occur right at the begining of the interval.
      sub_network_evolve(neu, begin_time, end_time);

    //add the force of poisson spike
    if (next_neu_index >= 0) {
#if SMOOTH_CONDUCTANCE_USE
      force_inputR(neu,next_neu_index);
#else
      force_input(neu,next_neu_index);
#endif
    }
    //calculating the next subTstep
    begin_time = end_time;
  } while(end_time<Tstep  && next_neu_index != -1 );
  //***********finish the evolution from t(0) to t(0)+Tstep

  //copy the spike information of the time interval
  raster_quicksort(spike_list);
  if (spike_list.ras_index>0) {
    for (j=0; j<spike_list.ras_index; j++) {
      raster_insert(RAS,spike_list.array_index[j],spike_list.array_firingtime[j]);
    }
  }
#else//network with poisson spikes ends
  begin_time = 0;
  end_time = Tstep;

  //we evolve the whole system directly from t(0) to t(0) + Tstep
  //since this time interval is already smooth
  sub_network_evolve(neu, begin_time, end_time);
  //***********finish the evolution from t(0) to t(0)+Tstep

  //copy the spike information of the time interval
  raster_quicksort(spike_list);
  if (spike_list.ras_index>0) {
    for (j=0; j<spike_list.ras_index; j++) {
      raster_insert(RAS,spike_list.array_index[j],spike_list.array_firingtime[j]);
    }
  }
#endif//network with no poisson spikes ends

  raster_destroy(spike_list);
///2013/1/24 16:06

}

void compute_perstep()
{
  int i;
  if (time_evolution >= last_time + g_comp_time) {
    RUN_DONE = 1;                   // program stop signal
    return ;
  }
  if (g_num_neu >= 12) {  // choose one faster
    ode_solver = runge_kutta4_vec;
  } else {
    ode_solver = runge_kutta4;
  }

  network_initialization();
  network_evolution();

  time_evolution += Tstep;
  // only update the data that needed
  int imax, neuron_index = 0, var_index = 0;
  if (g_no_graphic) {
    if (g_cond_out)
      imax = g_num_neu*(2*Stepsmooth_Con+4);
    else
      imax = g_num_neu;
  } else {
    imax = g_num_neu*size_neuronvar;
  }
///jyl strobe_veri is written to check the accuracy
  /*
  for (i = 0; i<imax;i++){
      var_index = i/g_num_neu;
      neuron_index = i%g_num_neu;
      double v = neu[neuron_index].value[var_index];
      strobe_veri(GLOBAL_STRA[i],time_evolution,v);
  }
*/
  if (g_b_RC_filter) {
    for (i=0; i<g_num_neu; i++) {

      double v = neu[i].value[0];

      strobeupdateRCFilter(GLOBAL_STRA[i], time_evolution, Tstep, v);
    }

    for (i=g_num_neu; i<imax; i++) {
      neuron_index = i%g_num_neu;
      var_index = i/g_num_neu;
      strobeupdateRCFilter(GLOBAL_STRA[i], time_evolution, Tstep,
                           neu[neuron_index].value[var_index]);
    }
//    orignal version
//    for (i=0; i<imax; i++) {
//      neuron_index = i%g_num_neu;
//      var_index = i/g_num_neu;
//      strobeupdateRCFilter(GLOBAL_STRA[i], time_evolution, Tstep,
//                           neu[neuron_index].value[var_index]);
//    }
/*
  } else{
      for (i=0;i<g_num_neu;i++){
        strobeupdate(GLOBAL_STRA[i], time_evolution, Tstep,
                     neu[i].value[2*Stepsmooth_Con+3]);
      }
  }
  */
    }
    else {
    for (i=0; i<imax; i++) {
      neuron_index = i%g_num_neu;
      var_index = i/g_num_neu;
      strobeupdate(GLOBAL_STRA[i], time_evolution, Tstep,
                   neu[neuron_index].value[var_index]);
    }
  }



  ///XYY: output neuron data to text file every SLIGHT_BIN
  /**edited by jyl
    structure of GLOBAL_STRA[] (same as neu)
    without smooth jump
      g_num_neu*[0,1) V   voltage
      g_num_neu*[1,2) GE  Ex_conductance
      g_num_neu*[2,3) GI  In_conductance
      g_num_neu*[3,4) m   ion channel conductance
      g_num_neu*[4,5) h   ion channel conductance
      g_num_neu*[5,6) n   ion channel conductance
    with smooth jump e.g. Stepsmooth_Com = 4
      g_num_neu*[0,1) V   voltage
      g_num_neu*[1,2) GE  Ex_conductance
      g_num_neu*[2,3) HE  multi-step smoothness
      g_num_neu*[3,4) IE  multi-step smoothness
      g_num_neu*[4,5) JE  multi-step smoothness
      g_num_neu*[5,6) KE  multi-step smoothness
      g_num_neu*[6,7) GI
      g_num_neu*[7,8) HI  multi-step smoothness
      g_num_neu*[8,9) II  multi-step smoothness
      g_num_neu*[9,10) JI  multi-step smoothness
      g_num_neu*[10,11) KI  multi-step smoothness
      g_num_neu*[11,12) m   ion channel conductance
      g_num_neu*[12,13) h   ion channel conductance
      g_num_neu*[13,14) n   ion channel conductance
  **/

  if (g_b_save_while_cal) {
    static int oldtab = 0;
    if (GLOBAL_STRA[0]->tab != oldtab) {          // if data is updated
      // loop for each neuron, see the structure of GLOBAL_STRA
      if (g_b_save_use_binary) {                  // save volt in binary format
        for (i = 0; i < g_num_neu; i++)
          fwrite(&(GLOBAL_STRA[i]->data[oldtab]), sizeof(double), 1, g_fout);
      } else {
        for (i = 0; i < g_num_neu; i++)
          fprintf(g_fout, "%19.16f ", GLOBAL_STRA[i]->data[oldtab]);
        fprintf(g_fout,"\n");
      }
      if (g_cond_out) {                           // save conductance
       if (g_b_save_use_binary) {                  // save volt in binary format
        for (i = g_num_neu; i < 2*g_num_neu; i++)
          fwrite(&(GLOBAL_STRA[i]->data[oldtab]), sizeof(double), 1, g_cond_out);
        for (i = (Stepsmooth_Con+1)*g_num_neu; i < (Stepsmooth_Con+2)*g_num_neu; i++)
          fwrite(&(GLOBAL_STRA[i]->data[oldtab]), sizeof(double), 1, g_cond_out);
      } else {
        for (i = g_num_neu; i < 2*g_num_neu; i++)
          fprintf(g_cond_out, "%19.16f ", GLOBAL_STRA[i]->data[oldtab]);
        for (i = (Stepsmooth_Con+1)*g_num_neu; i < (Stepsmooth_Con+2)*g_num_neu; i++)
          fprintf(g_cond_out, "%19.16f ", GLOBAL_STRA[i]->data[oldtab]);
        fprintf(g_cond_out,"\n");
      }
      }
      oldtab = GLOBAL_STRA[0]->tab;               // position to the last data
    }
  }
}
/// Save data at the tail of program
void LastRun()
{
  if (!g_b_quiet && !g_no_graphic)
    printf("t_end: %g ms\n", time_evolution);
  // save voltages
  if (!g_b_save_while_cal) {
    for (int k = 0; k < GLOBAL_STRA[0]->tab; k++) {
      // loop for each neuron, see the structure of GLOBAL_STRA
      for (int i = 0*g_num_neu; i < 1*g_num_neu; i++) {
        fprintf(g_fout, "%19.16f ", GLOBAL_STRA[i]->data[k]);
      }
      fprintf(g_fout,"\n");
    }
    fclose(g_fout);
  } else {
    fclose(g_fout);
  }
  if (g_cond_out) {
    fclose(g_cond_out);
  }
  if (g_fp_save_poisson_events!=NULL) {
    fclose(g_fp_save_poisson_events);
  }

  /// record firing time
  // neuron index starting from one
  if (g_ras_path[0]) {
    FILE *frasout = fopen(g_ras_path, "w");
    if (frasout==NULL) {
      printf("\nError: Fail to open \"%s\" for spike time file output!\n", g_ras_path);
    } else {
      for (int k=0; k<RAS.ras_index; k++) {
        fprintf(frasout, "%d\t%f\n", RAS.array_index[k]+1, RAS.array_firingtime[k]);
      }
      fclose(frasout);
    }

    /// count firing number (of each neuron) in each interval outTstep
//    const double outTstep = 32*Tstep;
//    FILE *fout = fopen("./data/z.txt", "w");
//    int tsc = 0;
//    int *z = (int*)calloc(g_num_neu, sizeof(int));
//    for (int i=0; i<RAS.ras_index; i++) {
//      if (RAS.array_firingtime[i] < outTstep*tsc) {
//        z[RAS.array_index[i]]++;
//      } else {
//        for (int j=0; j<g_num_neu; j++) {
//          if (z[j])
//            fprintf(fout,"%3d", z[j]);
//          else
//            fprintf(fout,"   ");
//        }
//        fprintf(fout,"\n");
//        memset(z, 0, g_num_neu*sizeof(z[0]));
//        tsc++;
//        i--;
//      }
//    }
//    for (int j=0; j<g_num_neu; j++)
//      fprintf(fout,"%3d", z[j]);
//    fprintf(fout,"\n");
//    fclose(fout);
//    free(z);
  }

  // give some basic statistic quantities
  int *frct = (int*)calloc(g_num_neu, sizeof(int));
  if (frct == NULL) {
    fprintf(stderr, "Warning: No enough memory for counting firing numbers.\n");
    fprintf(stderr, "Skipped this step.\n\n");
    return ;
  }
  for (int k=0; k<RAS.ras_index; k++)
    frct[RAS.array_index[k]]++;
  if (!g_b_quiet) {
    printf("\nFiring count:\n");
    for (int i=0; i<g_num_neu; i++)
      printf("%7d ", frct[i]);
    printf("\nAverage firing interval(ms):\n");
    for (int i=0; i<g_num_neu; i++)
      printf("%7.2f ", time_evolution/frct[i]);
    printf("\n");
  }/// print average spiking time
  if (g_spike_interval_path[0]) {
    FILE *fout = fopen(g_spike_interval_path, "w");
    if (fout==NULL) {
      printf("Error: Fail to open \"%s\" for average firing interval output!\n",
             g_spike_interval_path);
    } else {
      for (int i=0; i<g_num_neu; i++)
        fprintf(fout, "%.15g ", time_evolution/frct[i]);
      fclose(fout);
    }
  }
  free(frct);
}
