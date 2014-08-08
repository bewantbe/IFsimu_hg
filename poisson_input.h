///mark: modified on 2013/1/20
#ifndef _POISSON_INPUT_H_
#define _POISSON_INPUT_H_

#if SMOOTH_CONDUCTANCE_USE

//this computes the effect of the poisson spike force towards the target neuron in the condition of smooth conductance
//specifically, the poisson force is targeted on KE & KI
void force_inputR( neuron *tempneu, int index_neuron);

#else
//this computes the effect of the poisson spike force towards the target neuron in the condition of nonsmooth conductance
//specifically, the poisson force is targeted on GE & GI
void force_input( neuron *tempneu, int index_neuron);
#endif

// this function is used to compute the derivative of all the variables in neuron
void whole_dt(const neuron * const neu_val, neuron *neu_dt, double *volt, int index_neuron, double t);

//this is used to search the immediate next poisson spike of the whole network
void next_poisson_spike(neuron *tempneu, double begin_time,
                        int *tempbegin_poisson_index, double &end_time, int &next_neu_index);

double hermexp(double volt1, double volt2, double deno);

// this function is used in the computation of the spiking time
void voltage_dt(int index_neuron,double t,double m,double h,double n,double gE,double gI,double v,double &dv);
// this computes the external input current at time t (only for current input case)
double external_current(int index_neuron, double t);

// The following are the Runge-Kutta algorithm with different orders of convergence
void runge_kutta2(neuron *tempneu, double subTstep, double t_evolution);

void runge_kutta3(neuron *tempneu, double subTstep, double t_evolution);

void runge_kutta4(neuron *tempneu, double subTstep, double t_evolution);

// this uses four values to form a hermitian polynomial on the given interval [a,b]
// "va", "vb" are the given values of voltage at the point "a" and "b",
// "dva" and "dvb" are the given values of voltage derivative at the point "a" and "b"
// "x" is some point in [a,b], "fx" is the corresponding value of the hermitian polynomial
void hermit(double a, double b, double va, double vb, double dva, double dvb,
            double x, double &fx);

// search the root of function "func" at the given interval [x1,x2] under
// the accuracy of "xacc"
double root_search(void (*func)(double a, double b, double va, double vb,
                                double dva, double dvb, double x, double &fx), double x1, double x2,
                   double fx1, double fx2, double dfx1, double dfx2, double xacc);

// this determines the spike time of neuron at the given interval
void spiking_time(neuron *bef_neu, neuron *cur_neu, int index_neuron, double t_evolution, double Ta,
                  double Tb, double &firing_time);

// this function is used to compute all the information of the whole network in a smooth subTstep [begin_time,end_time]
void sub_network_evolve(neuron *tempneu, double begin_time, double end_time);

// initial setup of the whole neuron system (for reference trajectory)
void network_initialization();

// evolve the whole neuron system (for reference trajectory) from t(0) to t(0)+Tstep
void network_evolution();

// this just evolves the whole system for one time step ahead
void compute_perstep();

void LastRun();

#endif
