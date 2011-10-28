#include "stdafx.h"
#include "loop.h"
#include "neuron.h"
#include "raster.h"
#include "datainput.h"
#include "datahandling.h"
#include "poisson_input.h"

#if POISSON_INPUT_USE
#if SMOOTH_CONDUCTANCE_USE
// use smooth conductance HE and HI
inline void force_inputR( neuron *tempneu, int index_neuron)
{
  // excitatory conductance
//  tempneu[index_neuron].value[Stepsmooth_Con] += Strength_Exinput;
  tempneu[index_neuron].value[Stepsmooth_Con] += g_arr_poisson_strength_E[index_neuron];
  // inhibitory conductance
//  tempneu[index_neuron].value[2*Stepsmooth_Con] += Strength_Ininput;
  tempneu[index_neuron].value[2*Stepsmooth_Con] += g_arr_poisson_strength_I[index_neuron];
}
#else
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
#endif
#endif

// only use excitatory conductance and inhibitory conductance: GE and GI
void force_cortic( neuron *tempneu, int index_neuron, int index_spikingneuron)
{
  if (index_neuron != index_spikingneuron) {
#if CORTICAL_STRENGTH_NONHOMO
    tempneu[index_neuron].value[1] +=
       tempneu[index_spikingneuron].type *    tempneu[index_neuron].type
       * Strength_CorEE * cortical_matrix[index_neuron][index_spikingneuron]
     + tempneu[index_spikingneuron].type * (1-tempneu[index_neuron].type)
       * Strength_CorIE * cortical_matrix[index_neuron][index_spikingneuron];

    tempneu[index_neuron].value[Stepsmooth_Con+1] +=
       (1-tempneu[index_spikingneuron].type) *    tempneu[index_neuron].type
       * Strength_CorEI * cortical_matrix[index_neuron][index_spikingneuron]
     + (1-tempneu[index_spikingneuron].type) * (1-tempneu[index_neuron].type)
       * Strength_CorII * cortical_matrix[index_neuron][index_spikingneuron];
#else
    tempneu[index_neuron].value[1] +=
       tempneu[index_spikingneuron].type *    tempneu[index_neuron].type  * Strength_CorEE
     + tempneu[index_spikingneuron].type * (1-tempneu[index_neuron].type) * Strength_CorIE;

    tempneu[index_neuron].value[Stepsmooth_Con+1] +=
       (1-tempneu[index_spikingneuron].type) *    tempneu[index_neuron].type  * Strength_CorEI
     + (1-tempneu[index_spikingneuron].type) * (1-tempneu[index_neuron].type) * Strength_CorII;
#endif
  }
}

// use smooth conductance HE and HI
void force_corticR(neuron *tempneu, int index_neuron, int index_spikingneuron)
{
  if (index_neuron != index_spikingneuron) {
#if CORTICAL_STRENGTH_NONHOMO
    tempneu[index_neuron].value[Stepsmooth_Con] +=
       tempneu[index_spikingneuron].type *    tempneu[index_neuron].type
       * Strength_CorEE * cortical_matrix[index_neuron][index_spikingneuron]
     + tempneu[index_spikingneuron].type * (1-tempneu[index_neuron].type)
       * Strength_CorIE * cortical_matrix[index_neuron][index_spikingneuron];

    tempneu[index_neuron].value[2*Stepsmooth_Con] +=
       (1-tempneu[index_spikingneuron].type) *    tempneu[index_neuron].type
       * Strength_CorEI * cortical_matrix[index_neuron][index_spikingneuron]
     + (1-tempneu[index_spikingneuron].type) * (1-tempneu[index_neuron].type)
       * Strength_CorII * cortical_matrix[index_neuron][index_spikingneuron];
#else
    tempneu[index_neuron].value[Stepsmooth_Con] +=
       tempneu[index_spikingneuron].type *    tempneu[index_neuron].type  * Strength_CorEE
     + tempneu[index_spikingneuron].type * (1-tempneu[index_neuron].type) * Strength_CorIE;

    tempneu[index_neuron].value[2*Stepsmooth_Con] +=
       (1-tempneu[index_spikingneuron].type) *    tempneu[index_neuron].type  * Strength_CorEI
     + (1-tempneu[index_spikingneuron].type) * (1-tempneu[index_neuron].type) * Strength_CorII;
#endif
  }
}

// use non-smooth conductance
void conductance_decay( neuron *tempneu, int index_neuron, double subTstep)
{
// renew the Excitatory conductance
  tempneu[index_neuron].value[1] =
    tempneu[index_neuron].value[1]*exp(-subTstep/Time_ExCon);

// renew the Inhibitory conductance
  tempneu[index_neuron].value[Stepsmooth_Con+1] =
    tempneu[index_neuron].value[Stepsmooth_Con+1]*exp(-subTstep/Time_InCon);
}

// use smooth conductance
void conductance_evolve( neuron *tempneu, int index_neuron, double subTstep)
{
//******************************************************************************
//********************** renew the Excitatory conductance
  // renew the ex-conductance with highest smoothness
  tempneu[index_neuron].value[1] =
    tempneu[index_neuron].value[1]*exp(-subTstep/Time_ExCon)
   +1.0*Time_ExCon*Time_ExConR/(Time_ExConR-Time_ExCon)
    *(exp(-subTstep/Time_ExConR) - exp(-subTstep/Time_ExCon))
    *tempneu[index_neuron].value[Stepsmooth_Con];

  tempneu[index_neuron].value[Stepsmooth_Con] =
    tempneu[index_neuron].value[Stepsmooth_Con]*exp(-subTstep/Time_ExConR);

//******************************************************************************
//********************* renew the Inhibitory conductance
  // renew the in-conductance with highest smoothness
  tempneu[index_neuron].value[Stepsmooth_Con+1] =
    tempneu[index_neuron].value[Stepsmooth_Con+1]*exp(-subTstep/Time_InCon)
   +1.0*Time_InCon*Time_InConR/(Time_InConR-Time_InCon)
    *(exp(-subTstep/Time_InConR) - exp(-subTstep/Time_InCon))
    *tempneu[index_neuron].value[2*Stepsmooth_Con];

  tempneu[index_neuron].value[2*Stepsmooth_Con] =
    tempneu[index_neuron].value[2*Stepsmooth_Con]*exp(-subTstep/Time_InConR);
}

double external_current(int index_neuron, double t)
{
  return (Current_0 + Current_1*cos(2*PI*Rate_input*t+phase[index_neuron]));
}

// cortical interaction by conductance, external input either by current or by conductance
void voltage_dt(int index_neuron, double t, double gE, double gI, double v, double &dv_dt)
{
#if POISSON_INPUT_USE
  dv_dt = - Con_Leakage*(v - Vot_Leakage) - gE*(v - Vot_Excitatory)
          - gI*(v - Vot_Inhibitory);
#else
  dv_dt = - Con_Leakage*(v - Vot_Leakage) - gE*(v - Vot_Excitatory)
          - gI*(v - Vot_Inhibitory) + external_current(index_neuron, t);
#endif
}

void runge_kutta2(neuron *tempneu, int index_neuron, double subTstep,
                  double t_evolution, double &temp_vot,
                  void (*dvdt)(int index_neuron, double t, double gE, double gI,
                               double v, double &dv_dt))
{
  // dy/dt = f(t,y) t = t(n), y = y(n)
  double m1, m2;
  double newstep = subTstep; // time step needs to be changed for awaken neuron

  // the neuron is in refractory period
  tempneu[index_neuron].value[size_neuronvar-1] += ( 1 - tempneu[index_neuron].state_neuron )*newstep;

  // the neuron is in refractory period but will be awaken in this time!
  if (tempneu[index_neuron].value[size_neuronvar-1] >= TIME_REFRACTORY) {
    newstep = tempneu[index_neuron].value[size_neuronvar-1] - TIME_REFRACTORY;
    tempneu[index_neuron].state_neuron = STATE_ACTIVE;
    tempneu[index_neuron].value[size_neuronvar-1] = 0;
  }

  if (tempneu[index_neuron].state_neuron == STATE_ACTIVE) {
    double v0 = tempneu[index_neuron].value[0];
// if the neuron is awaken in this time step, the value of conductance should be used
// with the value of the awaking point******************************************
#if SMOOTH_CONDUCTANCE_USE
    double g_E0 = tempneu[index_neuron].value[1]*exp((newstep-subTstep)/Time_ExCon)
                  + 1.0*Time_ExCon*Time_ExConR/(Time_ExConR-Time_ExCon)
                  *(exp((newstep-subTstep)/Time_ExConR) - exp((newstep-subTstep)/Time_ExCon))
                  *tempneu[index_neuron].value[Stepsmooth_Con];
    double H_E0 = tempneu[index_neuron].value[Stepsmooth_Con]*exp((newstep-subTstep)/Time_ExConR);

    double g_I0 = tempneu[index_neuron].value[Stepsmooth_Con+1]*exp((newstep-subTstep)/Time_InCon)
                  + 1.0*Time_InCon*Time_InConR/(Time_InConR-Time_InCon)
                  *(exp((newstep-subTstep)/Time_InConR) - exp((newstep-subTstep)/Time_InCon))
                  *tempneu[index_neuron].value[2*Stepsmooth_Con];
    double H_I0 = tempneu[index_neuron].value[2*Stepsmooth_Con]*exp((newstep-subTstep)/Time_InConR);

    double g_E1 = g_E0*exp(-newstep/2/Time_ExCon)
                  + 1.0*Time_ExCon*Time_ExConR/(Time_ExConR-Time_ExCon)
                  *(exp(-newstep/2/Time_ExConR) - exp(-newstep/2/Time_ExCon))*H_E0;
    double g_I1 = g_I0*exp(-newstep/2/Time_InCon)
                  + 1.0*Time_InCon*Time_InConR/(Time_InConR-Time_InCon)
                  *(exp(-newstep/2/Time_InConR) - exp(-newstep/2/Time_InCon))*H_I0;
#else
    double g_E0 = tempneu[index_neuron].value[1]*exp((newstep-subTstep)/Time_ExCon);
    double g_I0 = tempneu[index_neuron].value[Stepsmooth_Con+1]*exp((newstep-subTstep)/Time_InCon);
    double g_E1 = g_E0*exp(-newstep/2/Time_ExCon);
    double g_I1 = g_I0*exp(-newstep/2/Time_InCon);
#endif

    double ini_time = t_evolution + (subTstep - newstep);
    double mid_time = t_evolution + (subTstep - newstep) + newstep/2;

    // m1 = f(t(n), y(n))
    (*dvdt)(index_neuron, ini_time, g_E0, g_I0, v0, m1);

    // m2 = f(t(n)+h/2, y(n)+m1*h/2)
    (*dvdt)(index_neuron, mid_time, g_E1, g_I1, v0+m1/2*newstep, m2);

    temp_vot = v0 + m2*newstep;
  }
}

void runge_kutta3(neuron *tempneu, int index_neuron, double subTstep,
                  double t_evolution, double &temp_vot,
                  void (*dvdt)(int index_neuron, double t, double gE, double gI,
                               double v, double &dv_dt))
{
  // dy/dt = f(t,y) t = t(n), y = y(n)
  double m1, m2, m3;
  double newstep = subTstep; // time step needs to be changed for awaken neuron

  // the neuron is in refractory period
  tempneu[index_neuron].value[size_neuronvar-1] +=
    (1 - tempneu[index_neuron].state_neuron)*newstep;

  // the neuron is in refractory period but will be awaken in this time!
  if (tempneu[index_neuron].value[size_neuronvar-1] >= TIME_REFRACTORY) {
    newstep = tempneu[index_neuron].value[size_neuronvar-1] - TIME_REFRACTORY;
    tempneu[index_neuron].state_neuron = STATE_ACTIVE;
    tempneu[index_neuron].value[size_neuronvar-1] = 0;
  }

  if (tempneu[index_neuron].state_neuron == STATE_ACTIVE) {
    double v0 = tempneu[index_neuron].value[0];
// if the neuron is awaken in this time step, the value of conductance should be used
// with the value of the awaking point******************************************
#if SMOOTH_CONDUCTANCE_USE
    double g_E0 = tempneu[index_neuron].value[1]*exp((newstep-subTstep)/Time_ExCon)
                  + 1.0*Time_ExCon*Time_ExConR/(Time_ExConR-Time_ExCon)
                  *(exp((newstep-subTstep)/Time_ExConR) - exp((newstep-subTstep)/Time_ExCon))
                  *tempneu[index_neuron].value[Stepsmooth_Con];
    double H_E0 = tempneu[index_neuron].value[Stepsmooth_Con]*exp((newstep-subTstep)/Time_ExConR);

    double g_I0 = tempneu[index_neuron].value[Stepsmooth_Con+1]*exp((newstep-subTstep)/Time_InCon)
                  + 1.0*Time_InCon*Time_InConR/(Time_InConR-Time_InCon)
                  *(exp((newstep-subTstep)/Time_InConR) - exp((newstep-subTstep)/Time_InCon))
                  *tempneu[index_neuron].value[2*Stepsmooth_Con];
    double H_I0 = tempneu[index_neuron].value[2*Stepsmooth_Con]*exp((newstep-subTstep)/Time_InConR);

    double g_E1 = g_E0*exp(-newstep/3/Time_ExCon)
                  + 1.0*Time_ExCon*Time_ExConR/(Time_ExConR-Time_ExCon)
                  *(exp(-newstep/3/Time_ExConR) - exp(-newstep/3/Time_ExCon))*H_E0;
    double g_I1 = g_I0*exp(-newstep/3/Time_InCon)
                  + 1.0*Time_InCon*Time_InConR/(Time_InConR-Time_InCon)
                  *(exp(-newstep/3/Time_InConR) - exp(-newstep/3/Time_InCon))*H_I0;

    double g_E2 = g_E0*exp(-2*newstep/3/Time_ExCon)
                  + 1.0*Time_ExCon*Time_ExConR/(Time_ExConR-Time_ExCon)
                  *(exp(-2*newstep/3/Time_ExConR) - exp(-2*newstep/3/Time_ExCon))*H_E0;
    double g_I2 = g_I0*exp(-2*newstep/3/Time_InCon)
                  + 1.0*Time_InCon*Time_InConR/(Time_InConR-Time_InCon)
                  *(exp(-2*newstep/3/Time_InConR) - exp(-2*newstep/3/Time_InCon))*H_I0;
#else
    double g_E0 = tempneu[index_neuron].value[1]*exp((newstep-subTstep)/Time_ExCon);
    double g_I0 = tempneu[index_neuron].value[Stepsmooth_Con+1]*exp((newstep-subTstep)/Time_InCon);
    double g_E1 = g_E0*exp(-newstep/3/Time_ExCon);
    double g_I1 = g_I0*exp(-newstep/3/Time_InCon);
    double g_E2 = g_E0*exp(-2*newstep/3/Time_ExCon);
    double g_I2 = g_I0*exp(-2*newstep/3/Time_InCon);
#endif

    double ini_time = t_evolution + (subTstep - newstep);
    double one_third_time = t_evolution + (subTstep - newstep) + newstep/3;
    double two_third_time = t_evolution + (subTstep - newstep) + 2*newstep/3;

    // m1 = f(t(n), y(n))
    (*dvdt)(index_neuron, ini_time, g_E0, g_I0, v0, m1);

    // m2 = f(t(n)+h/3, y(n)+m1*h/3)
    (*dvdt)(index_neuron, one_third_time, g_E1, g_I1, v0+m1/3*newstep, m2);

    // m3 = f(t(n)+2*h/3, y(n)+2*m2*h/3)
    (*dvdt)(index_neuron, two_third_time, g_E2, g_I2, v0+2*m2/3*newstep, m3);

    temp_vot = v0 + (m1 + 3*m3)*newstep/4;
  }
}

void runge_kutta4(neuron *tempneu, int index_neuron, double subTstep,
                  double t_evolution, double &temp_vot,
                  void (*dvdt)(int index_neuron, double t, double gE, double gI,
                               double v, double &dv_dt))
{
  // dy/dt = f(t,y) t = t(n), y = y(n)
  double m1, m2, m3, m4;
  double newstep = subTstep; // time step needs to be changed for awaken neuron

  // the neuron is in refractory period
  tempneu[index_neuron].value[size_neuronvar-1] +=
    (1 - tempneu[index_neuron].state_neuron)*newstep;

  // the neuron is in refractory period but will be awaken in this time!
  if (tempneu[index_neuron].value[size_neuronvar-1] >= TIME_REFRACTORY) {
    newstep = tempneu[index_neuron].value[size_neuronvar-1] - TIME_REFRACTORY;
    tempneu[index_neuron].state_neuron = STATE_ACTIVE;
    tempneu[index_neuron].value[size_neuronvar-1] = 0;

  }

  if (tempneu[index_neuron].state_neuron == STATE_ACTIVE) {
    double v0 = tempneu[index_neuron].value[0];
// if the neuron is awaken in this time step, the value of conductance should be used
// with the value of the awaking point******************************************
#if SMOOTH_CONDUCTANCE_USE
    double g_E0 = tempneu[index_neuron].value[1]*exp((newstep-subTstep)/Time_ExCon)
                  + 1.0*Time_ExCon*Time_ExConR/(Time_ExConR-Time_ExCon)
                  *(exp((newstep-subTstep)/Time_ExConR) - exp((newstep-subTstep)/Time_ExCon))
                  *tempneu[index_neuron].value[Stepsmooth_Con];
    double H_E0 = tempneu[index_neuron].value[Stepsmooth_Con]*exp((newstep-subTstep)/Time_ExConR);

    double g_I0 = tempneu[index_neuron].value[Stepsmooth_Con+1]*exp((newstep-subTstep)/Time_InCon)
                  + 1.0*Time_InCon*Time_InConR/(Time_InConR-Time_InCon)
                  *(exp((newstep-subTstep)/Time_InConR) - exp((newstep-subTstep)/Time_InCon))
                  *tempneu[index_neuron].value[2*Stepsmooth_Con];
    double H_I0 = tempneu[index_neuron].value[2*Stepsmooth_Con]*exp((newstep-subTstep)/Time_InConR);

    double g_E1 = g_E0*exp(-newstep/2/Time_ExCon)
                  + 1.0*Time_ExCon*Time_ExConR/(Time_ExConR-Time_ExCon)
                  *(exp(-newstep/2/Time_ExConR) - exp(-newstep/2/Time_ExCon))*H_E0;
    double g_I1 = g_I0*exp(-newstep/2/Time_InCon)
                  + 1.0*Time_InCon*Time_InConR/(Time_InConR-Time_InCon)
                  *(exp(-newstep/2/Time_InConR) - exp(-newstep/2/Time_InCon))*H_I0;

    double g_E2 = g_E0*exp(-newstep/Time_ExCon)
                  + 1.0*Time_ExCon*Time_ExConR/(Time_ExConR-Time_ExCon)
                  *(exp(-newstep/Time_ExConR) - exp(-newstep/Time_ExCon))*H_E0;
    double g_I2 = g_I0*exp(-newstep/Time_InCon)
                  + 1.0*Time_InCon*Time_InConR/(Time_InConR-Time_InCon)
                  *(exp(-newstep/Time_InConR) - exp(-newstep/Time_InCon))*H_I0;
#else
    double g_E0 = tempneu[index_neuron].value[1]*exp((newstep-subTstep)/Time_ExCon);
    double g_I0 = tempneu[index_neuron].value[Stepsmooth_Con+1]*exp((newstep-subTstep)/Time_InCon);
    double g_E1 = g_E0*exp(-newstep/2/Time_ExCon);
    double g_I1 = g_I0*exp(-newstep/2/Time_InCon);
    double g_E2 = g_E0*exp(-newstep/Time_ExCon);
    double g_I2 = g_I0*exp(-newstep/Time_InCon);
#endif

    double ini_time = t_evolution + (subTstep - newstep);
    double two_fourth_time = t_evolution + (subTstep - newstep) + newstep/2;
    double end_time = t_evolution + (subTstep - newstep) + newstep;

    // m1 = f(t(n), y(n))
    (*dvdt)(index_neuron, ini_time, g_E0, g_I0, v0, m1);

    // m2 = f(t(n)+h/2, y(n)+m1*h/2)
    (*dvdt)(index_neuron, two_fourth_time, g_E1, g_I1, v0+m1/2*newstep, m2);

    // m3 = f(t(n)+h/2, y(n)+m2*h/2)
    (*dvdt)(index_neuron, two_fourth_time, g_E1, g_I1, v0+m2/2*newstep, m3);

    // m4 = f(t(n)+h, y(n)+m3*h)
    (*dvdt)(index_neuron, end_time, g_E2, g_I2, v0+m3*newstep, m4);

    temp_vot = v0 + (m1 + 2*m2 + 2*m3 + m4)*newstep/6;

  } else { // neuron is in the refractory period!
    temp_vot = VOT_RESET;
  }
}

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
  P_ERR("Too many bisections in root searching!\n");
  return xmid;
}

void spiking_time(neuron *tempneu, int index_neuron, double t_evolution, double Ta,
                  double Tb, double vot_cross, double &firing_time)
{
  double subTstep = Tb - Ta;

  double va = tempneu[index_neuron].value[0];
  double dva = 0; // initialization
  double gEa = tempneu[index_neuron].value[1];
  double gIa = tempneu[index_neuron].value[Stepsmooth_Con+1];

  voltage_dt(index_neuron,t_evolution,gEa,gIa,va,dva);

  double vb = vot_cross;
  double dvb = 0; // initialization

#if SMOOTH_CONDUCTANCE_USE
  double gEb = tempneu[index_neuron].value[1]*exp(-subTstep/Time_ExCon)
               + 1.0*Time_ExCon*Time_ExConR/(Time_ExConR-Time_ExCon)
               *(exp(-subTstep/Time_ExConR) - exp(-subTstep/Time_ExCon))
               *tempneu[index_neuron].value[Stepsmooth_Con];
  double gIb = tempneu[index_neuron].value[Stepsmooth_Con+1]*exp(-subTstep/Time_InCon)
               + 1.0*Time_InCon*Time_InConR/(Time_InConR-Time_InCon)
               *(exp(-subTstep/Time_InConR) - exp(-subTstep/Time_InCon))
               *tempneu[index_neuron].value[2*Stepsmooth_Con];
#else
  double gEb = tempneu[index_neuron].value[1]*exp(-subTstep/Time_ExCon);
  double gIb = tempneu[index_neuron].value[Stepsmooth_Con+1]*exp(-subTstep/Time_InCon);
#endif

  voltage_dt(index_neuron,(t_evolution+Tb-Ta),gEb,gIb,vb,dvb);

  firing_time = root_search(hermit, Ta, Tb, va, vb, dva, dvb, root_acc);
}

#if POISSON_INPUT_USE
// given two cortical spike time T1 and T2 (T1<T2), evolve each neuron from
// T1 to T2, thus the whole system evolves from T1 to T2
// if there is some neuron firing between [T1, T2], record the nearest time TT
// evolve the whole system from T1 to TT!
//***************** renewing conductance process is done outside this program!
void single_neuron_test(neuron *tempneu, int index_neuron, int firing_neuron,
                        double begin_time, double end_time,
                        int *begin_poisson_index, struct raster &temp_spike_list, int &jump)
{
  // because the system may need to be pulled back, therefore some variables are
  // just temporal such as: *tempneu, *begin_poisson_index, end_time

  // consider if there is a neuron jump from A(quiet) group to A(spike) group!
  jump = 0;
  //************** cortical force to conductance is done outside this function!

  // record the starting point of poisson input in this interval!
  int j = begin_poisson_index[index_neuron];
  // for each neuron i, [T1, T2] is divided by poisson input spikes t(i1),t(i2),...
  // evolve interval [T1, t(i1)] and then [t(i2), t(i3)], so on and so fourth
  double subbegin_time = begin_time;
  double subend_time = end_time; // record the middle point of poisson input!
  double sub_tstep = end_time - begin_time;

  // do hermit interpolation if a neuron fires
  double firing_time = 0;

  // record the voltage once it fires, if neuron is during refractory period
  // and will not wake up, then runge_kutta method just do nothing!
  // temp_vot will be the same value as the voltage of tempneu[i]
  double temp_vot = tempneu[index_neuron].value[0];
//  double temp_threshold = 0; // record the value of threshold once it fires

  // if there is an external input current!
  double t_evolution = time_evolution+begin_time;
  //****************************************************************************
  //**************  no poisson input in the whole interval!
  if (j >= poisson_input[index_neuron].vect_size) {
    if (index_neuron==firing_neuron) {
#if SMOOTH_CONDUCTANCE_USE
      conductance_evolve(tempneu, index_neuron, sub_tstep);
#else
      conductance_decay(tempneu, index_neuron, sub_tstep);
#endif
      tempneu[index_neuron].value[0] = VOT_RESET;
      tempneu[index_neuron].state_neuron = STATE_REFRACTORY;
      return;
    }

#if RK2_RUN
    runge_kutta2(tempneu, index_neuron, sub_tstep, t_evolution, temp_vot,
                 voltage_dt);
#endif
#if RK3_RUN
    runge_kutta3(tempneu, index_neuron, sub_tstep, t_evolution, temp_vot,
                 voltage_dt);
#endif
#if RK4_RUN
    runge_kutta4(tempneu, index_neuron, sub_tstep, t_evolution, temp_vot,
                 voltage_dt);
#endif

    // neuron is firing without poisson input spikes
    if(temp_vot >= Vot_Threshold) {
      spiking_time(tempneu, index_neuron, t_evolution, subbegin_time,
                   subend_time, temp_vot, firing_time);
      raster_insert(temp_spike_list, index_neuron, firing_time);
      jump++;
      return;
    }
    tempneu[index_neuron].value[0] = temp_vot;

#if SMOOTH_CONDUCTANCE_USE
    conductance_evolve(tempneu, index_neuron, sub_tstep);
#else
    conductance_decay(tempneu, index_neuron, sub_tstep);
#endif

    return;
  } // end if there is no poisson spikes
//******************************* 2007.10.19 1:50AM ****************************
  while ((j<poisson_input[index_neuron].vect_size)
         && (poisson_input[index_neuron].vect_value[j]<=end_time)) {
    subend_time = poisson_input[index_neuron].vect_value[j];
    sub_tstep = subend_time - subbegin_time;
    if (index_neuron == firing_neuron) {
#if SMOOTH_CONDUCTANCE_USE
      conductance_evolve(tempneu, index_neuron, sub_tstep);
      force_inputR(tempneu, index_neuron);
#else
      conductance_decay(tempneu, index_neuron, sub_tstep);
      force_input(tempneu, index_neuron);
#endif
      subbegin_time = subend_time;
      j++;
      continue;
    }

#if RK2_RUN
    runge_kutta2(tempneu, index_neuron, sub_tstep, t_evolution, temp_vot,
                 voltage_dt);
#endif
#if RK3_RUN
    runge_kutta3(tempneu, index_neuron, sub_tstep, t_evolution, temp_vot,
                 voltage_dt);
#endif
#if RK4_RUN
    runge_kutta4(tempneu, index_neuron, sub_tstep, t_evolution, temp_vot,
                 voltage_dt);
#endif

//******************************************************************************
    if(temp_vot >= Vot_Threshold) {
      spiking_time(tempneu, index_neuron, t_evolution, subbegin_time, subend_time,
                   temp_vot, firing_time);
      raster_insert(temp_spike_list, index_neuron, firing_time);
      jump++;
      return;
    }

//******************************************************************************
    // neuron is not firing at this sub-interval,
    // renew the voltage and conductance, then go ahead to the next sub-interval
    tempneu[index_neuron].value[0] = temp_vot;

#if SMOOTH_CONDUCTANCE_USE
    conductance_evolve(tempneu, index_neuron, sub_tstep);
    force_inputR(tempneu, index_neuron);
#else
    conductance_decay(tempneu, index_neuron, sub_tstep);
    force_input(tempneu, index_neuron);
#endif
    subbegin_time = subend_time;
    t_evolution += sub_tstep;
    j++;
  } // end of while (for poisson spikes)

  // finish the last sub-interval!
  //*******************************************************************
  sub_tstep = end_time - subbegin_time;

  if (index_neuron == firing_neuron) {

#if SMOOTH_CONDUCTANCE_USE
    conductance_evolve(tempneu, index_neuron, sub_tstep);
#else
    conductance_decay(tempneu, index_neuron, sub_tstep);
#endif

//********************************* 2007.10.19 18:39PM *************************
    tempneu[index_neuron].value[0] = VOT_RESET;
    tempneu[index_neuron].state_neuron = STATE_REFRACTORY;

    begin_poisson_index[index_neuron] = j;
    return;
  }

#if RK2_RUN
  runge_kutta2(tempneu, index_neuron, sub_tstep, t_evolution, temp_vot,
               voltage_dt);
#endif
#if RK3_RUN
  runge_kutta3(tempneu, index_neuron, sub_tstep, t_evolution, temp_vot,
               voltage_dt);
#endif
#if RK4_RUN
  runge_kutta4(tempneu, index_neuron, sub_tstep, t_evolution, temp_vot,
               voltage_dt);
#endif

//******************************************************************************
  if(temp_vot >= Vot_Threshold) {
// this is just for finishing the last step, so the ending time should be "end_time"
// instead of "subend_time"*****************************************************
    spiking_time(tempneu, index_neuron, t_evolution, subbegin_time, end_time,
                 temp_vot, firing_time);
    raster_insert(temp_spike_list, index_neuron, firing_time);
    jump++;
    return;

  }

//******************************************************************************
  // neuron is not firing at this sub-interval,
  // renew the voltage and conductance, then go ahead to the next sub-interval
  tempneu[index_neuron].value[0] = temp_vot;

#if SMOOTH_CONDUCTANCE_USE
  conductance_evolve(tempneu, index_neuron, sub_tstep);
#else
  conductance_decay(tempneu, index_neuron, sub_tstep);
#endif

  // cortical force to conductance is done outside this function!
  begin_poisson_index[index_neuron] = j;
}
#endif

// This is for non-poisson input only current input case!!!
void single_neuron_test(neuron *tempneu, int index_neuron, int firing_neuron,
                        double begin_time, double end_time,
                        struct raster &temp_spike_list, int &jump)
{
  // because the system may need to be pulled back, therefore some variables are
  // just temporal such as: *tempneu, end_time

  // consider if there is a neuron jump from A(quiet) group to A(spike) group!
  jump = 0;
  //************** cortical force to conductance is done outside this function!

  // do hermit interpolation if a neuron fires
  double firing_time = 0;
  double temp_vot = tempneu[index_neuron].value[0]; // record the voltage once it fires

  // if there is an external input current!
  double sub_tstep = end_time - begin_time;

  double t_evolution = time_evolution + begin_time;
  //****************************************************************************
  //**************  no poisson input in the whole interval!
  if (index_neuron==firing_neuron) {
#if SMOOTH_CONDUCTANCE_USE
    conductance_evolve(tempneu, index_neuron, sub_tstep);
#else
    conductance_decay(tempneu, index_neuron, sub_tstep);
#endif
    tempneu[index_neuron].value[0] = VOT_RESET;
    tempneu[index_neuron].state_neuron = STATE_REFRACTORY;

    return;
  }
#if RK2_RUN
  runge_kutta2(tempneu, index_neuron, sub_tstep, t_evolution, temp_vot,
               voltage_dt);
#endif

#if RK3_RUN
  runge_kutta3(tempneu, index_neuron, sub_tstep, t_evolution, temp_vot,
               voltage_dt);
#endif

#if RK4_RUN
  runge_kutta4(tempneu, index_neuron, sub_tstep, t_evolution, temp_vot,
               voltage_dt);
#endif

//******************************************************************************
  // neuron is firing without poisson input spikes
  if(temp_vot >= Vot_Threshold) {
    spiking_time(tempneu, index_neuron, t_evolution, begin_time,
                 end_time, temp_vot, firing_time);
    raster_insert(temp_spike_list, index_neuron, firing_time);

    jump++;
    return;
  }

//******************************************************************************
  tempneu[index_neuron].value[0] = temp_vot;

#if SMOOTH_CONDUCTANCE_USE
  conductance_evolve(tempneu, index_neuron, sub_tstep);
#else
  conductance_decay(tempneu, index_neuron, sub_tstep);
#endif

}

void next_cortical_spike(neuron *tempneu, double begin_time,
                         int *tempbegin_poisson_index)
{
  struct raster new_spike_list;
  raster_initialize(new_spike_list);
  raster_allocate(new_spike_list, spike_list.ras_size);

  for (int i=0; i<g_num_neu; i++) {
#if POISSON_INPUT_USE
    tmp_tempbegin_poisson_index[i] = tempbegin_poisson_index[i];
#endif
    neuron_copy_2(tmp_tempneu[i], tempneu[i]);
  }
  // search the A(spike) list tempneurons first to get the most probable next spiking time!
  for (int i=0; i<spike_list.ras_index; i++) {
    double end_time = Tstep;
    int num_jump = 0;
    int index_tempneu = spike_list.array_index[i];
    // This is for determining the next cortical spike by searching
#if POISSON_INPUT_USE
    // choose begin_time as the last firing time, and end_time as the whole timestep!
    single_neuron_test(tmp_tempneu, index_tempneu, -1, begin_time, end_time,
                       tmp_tempbegin_poisson_index, new_spike_list, num_jump);
#else
    single_neuron_test(tmp_tempneu, index_tempneu, -1, begin_time, end_time,
                       new_spike_list, num_jump);
#endif
    // the tempneuron will not be firing in this time interval!
    if(num_jump == 0) {
      // still remain this tempneuron in the new list
      // but put the firing time of this tempneuron to later time
      raster_insert(new_spike_list, index_tempneu, 2*Tstep);
    }
  }
  // get a sequential increment spiking time list
  raster_destroy(spike_list);
  raster_copy(spike_list, new_spike_list);
  raster_destroy(new_spike_list);

  raster_quicksort(spike_list);  // sort the new firing events!

  return;
}

void network_initialization()
{
  // record if there is some neuron firing in the middle of the time interval!
  int i;
  int jump;

#if POISSON_INPUT_USE
  int *&tempbegin_poisson_index = g_tempbegin_poisson_index;
#endif

  double begin_time = 0;
  double end_time = Tstep;
  raster_allocate(spike_list, RASTER_SPIKE_SIZE);

  //******************** poisson input spikes generation
#if POISSON_INPUT_USE
  // last_input is dynamically distributed space outside this function!
  for (i=0; i<g_num_neu; i++) {
    poisson_generator(i, poisson_input[i]);
  }
#endif

  //****************************************************************************
  // construct a temporal neuron variable in case the system needs to be pulled back!
  neuron *&tempneu = tmp_tempneu2;        // use global temporary variable

  for (i=0; i<g_num_neu; i++) {
    neuron_copy_raw_static(tempneu[i], neu[i]);
#if POISSON_INPUT_USE
    tempbegin_poisson_index[i] = 0;
#endif
  }

  //****************************************************************************
  // get the initial A(spike) list!
  int sum_spikes = 0;

  for (i=0; i<g_num_neu; i++) {
    jump = 0;
#if POISSON_INPUT_USE // start of using poisson input
    single_neuron_test(tempneu, i, -1, begin_time, end_time,
                       tempbegin_poisson_index, spike_list, jump);
#else // not using poisson input
    single_neuron_test(tempneu, i, -1, begin_time, end_time,
                       spike_list, jump);
#endif
    sum_spikes += jump;
  }
  //****************************************************************************
  // there is no spike in this time interval!
  if (sum_spikes == 0) {
    for (i=0; i<g_num_neu; i++) {
      neuron_copy_raw_static(neu[i],tempneu[i]);
    }
    return;
  } // end for no neuron spiking case!
  //****************************** release the temporal variables
  // sort the spike list!
  raster_quicksort(spike_list);
}

// evolve the whole neuron system from t(0) to t(0)+Tstep
void network_evolution()
{
  int i;
  int index_firingneuron;
  int jump; // record if there is some neuron firing in the middle of the time interval!
  double begin_time = 0;
  double end_time = 0;

#if POISSON_INPUT_USE
  int *&begin_poisson_index = g_begin_poisson_index;
  int *&tempbegin_poisson_index = g_tempbegin_poisson_index;
#else
  int *tempbegin_poisson_index = NULL;
#endif

  //****************************************************************************
  // there is no spike in this time step, evolve system to the next time step!
  if (spike_list.ras_index <= 0) {
    raster_destroy(spike_list);
    return;
  }

  // there is some spikes in this time step, set up a temporal neuron variable
  // to evolve system to a sub-timestep to see if there is any spikes in this sub-timestep
  neuron *&tempneu = tmp_tempneu2;

  for (i=0; i<g_num_neu; i++) {
#if POISSON_INPUT_USE
    begin_poisson_index[i] = 0;
    tempbegin_poisson_index[i] = 0;
#endif
//    neuron_copy_2(tempneu[i], neu[i]);
    neuron_copy_raw_static(tempneu[i], neu[i]);
  }

  //****************************************************************************
  // get the earliest firing time from the A(spike) list!
  int sum_spikes = 0;
  do { // big circle!
    //*******************************************************************
    // evolve the system to a sub-timestep!
    if (spike_list.ras_index <=0) {
      end_time = Tstep;
      index_firingneuron = -1;
    } else {
      end_time = spike_list.array_firingtime[0];
      index_firingneuron = spike_list.array_index[0];
    }
//*******************************************
    if (end_time > Tstep) {
      end_time = Tstep;
      index_firingneuron = -1;
    }
//*********************************************
    if (index_firingneuron>=0 &&
        tempneu[index_firingneuron].state_neuron == STATE_REFRACTORY) {
      tempneu[index_firingneuron].value[size_neuronvar-1] = 0;  // fix the Bug 010.
    }
    sum_spikes = 0; // record if there is a spike occuring before end_time
    for (i=0; i<g_num_neu; i++) {
      jump = 0;
#if POISSON_INPUT_USE
      single_neuron_test(tempneu, i, index_firingneuron, begin_time, end_time,
                         tempbegin_poisson_index, spike_list, jump);
#else
      single_neuron_test(tempneu, i, index_firingneuron, begin_time, end_time,
                         spike_list, jump);
#endif
      sum_spikes += jump;
    }

    //*******************************************************************
    // if there is no spike in this sub-timestep! go to the next sub-timestep!
    if (sum_spikes == 0) {
      if(end_time<Tstep || (index_firingneuron!=-1)) {
        raster_insert(RAS, index_firingneuron, (end_time+time_evolution));
        //*****************************************************************
        // release the firing neuron from the A(spike) list
        raster_erase(spike_list, 0);
        for (i=0; i<g_num_neu; i++) {
          // renew the conductance at the end of sub-timestep
#if SMOOTH_CONDUCTANCE_USE
          force_corticR(tempneu, i, index_firingneuron);
#else
          force_cortic(tempneu, i, index_firingneuron);
#endif
        }
        // obtain the next sub-time step
        begin_time = end_time;
        next_cortical_spike(tempneu, begin_time, tempbegin_poisson_index);
      }
      for (i=0; i<g_num_neu; i++) {
        neuron_copy_raw_static(neu[i],tempneu[i]);
#if POISSON_INPUT_USE
        // renew the start index of poisson input
        begin_poisson_index[i] = tempbegin_poisson_index[i];
#endif
      }
    } else { // there is some guy firing before the supposed first firing guy!
      for (i=0; i<g_num_neu; i++) {
        // revert the value of neuron
        neuron_copy_raw_static(tempneu[i], neu[i]);
#if POISSON_INPUT_USE
        // revert the start index of poisson input
        tempbegin_poisson_index[i] = begin_poisson_index[i];
#endif
      }
      raster_quicksort(spike_list);
    }
  } while (sum_spikes!=0 || end_time<Tstep);

  raster_destroy(spike_list);
}

void compute_perstep()
{
  int i;

  if (time_evolution >= last_time + g_comp_time) {
    RUN_DONE = 1;                   // program stop signal
    return ;
  }

  network_initialization();
  network_evolution();

  time_evolution += Tstep;

  // only update the data that needed
  int imax, neuron_index = 0, var_index = 0;
  if (g_no_graphic) {
    if (g_cond_out)
      imax = g_num_neu*2;
    else
      imax = g_num_neu;
  } else {
    imax = g_num_neu*size_neuronvar;
  }
  if (g_b_RC_filter) {
    for (i=0; i<imax; i++) {
      neuron_index = i%g_num_neu;
      var_index = i/g_num_neu;
      strobeupdateRCFilter(GLOBAL_STRA[i], time_evolution, Tstep,
                           neu[neuron_index].value[var_index]);
    }
  } else {
    for (i=0; i<imax; i++) {
      neuron_index = i%g_num_neu;
      var_index = i/g_num_neu;
      strobeupdate(GLOBAL_STRA[i], time_evolution, Tstep,
                   neu[neuron_index].value[var_index]);
    }
  }

  ///XYY: output neuron data to text file every SLIGHT_BIN
  /**
    structure of GLOBAL_STRA[] (same as neu)
    without smooth jump
      g_num_neu*[0,1) V   voltage
      g_num_neu*[1,2) GE  Ex_conductance
      g_num_neu*[2,3) GI  In_conductance
      g_num_neu*[3,4) ref_time
    with smooth jump
      g_num_neu*[0,1) voltage data
      g_num_neu*[1,2) GE
      g_num_neu*[2,3) HE  multi-step smoothness
      g_num_neu*[3,4) GI
      g_num_neu*[4,5) HI  multi-step smoothness
      g_num_neu*[5,6) ref_time
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
        for (i = g_num_neu; i < 2*g_num_neu; i++)
          fprintf(g_cond_out, "%19.16f ", GLOBAL_STRA[i]->data[oldtab]);
        fprintf(g_cond_out,"\n");
      }
      oldtab = GLOBAL_STRA[0]->tab;               // position to the last data
    }
  }
}

/// Save data at the tail of program
void LastRun()
{
  if (!g_b_quiet)
    printf("t_end: %g ms\n", time_evolution);
  // save voltages
  if (!g_b_save_while_cal) {
    FILE *fout = fopen(g_staffsave_path, "w");
    for (int k = 0; k < GLOBAL_STRA[0]->tab; k++) {
      // loop for each neuron, see the structure of GLOBAL_STRA
      for (int i = 0*g_num_neu; i < 1*g_num_neu; i++) {
        fprintf(fout, "%19.16f ", GLOBAL_STRA[i]->data[k]);
      }
      fprintf(fout,"\n");
    }
    fclose(fout);
  } else {
    fclose(g_fout);
  }
  if (g_cond_out) {
    fclose(g_cond_out);
  }

  /// record firing time
  // neuron index starting from one
  if (g_ras_path[0]) {
    FILE *frasout = fopen(g_ras_path, "w");
    for (int k=0; k<RAS.ras_index; k++) {
      fprintf(frasout, "%d\t%f\n", RAS.array_index[k]+1, RAS.array_firingtime[k]);
    }
    fclose(frasout);

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
  }
  if (g_spike_interval_path[0]) {
    FILE *fout = fopen(g_spike_interval_path, "w");
    for (int i=0; i<g_num_neu; i++)
      fprintf(fout, "%g ", time_evolution/frct[i]);
    fclose(fout);
  }
  free(frct);
}
