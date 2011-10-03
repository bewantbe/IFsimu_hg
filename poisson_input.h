#ifndef _POISSON_INPUT_H_
#define _POISSON_INPUT_H_

#if SMOOTH_CONDUCTANCE_USE
// same as the above but for smoothed conductance case, it describes the jump
// of H by the external Poisson spike train
inline void force_inputR( neuron *tempneu, int index_neuron);
#else
// this is used for non-smooth conductance case, it describes the jump of
// conductance G by the external Poisson spike train
void force_input( neuron *tempneu, int index_neuron);
#endif

// this is used for non-smooth conductance case, it describes the jump of
// conductance G by other neurons' spikes
void force_cortic( neuron *tempneu, int index_neuron, int index_spikingneuron);

// same as the above but for smoothed conductance case, it describes the jump
// of H by other neurons' spikes where H is non-smoothed part
void force_corticR( neuron *tempneu, int index_neuron, int index_spikingneuron);

// this is used for non-smooth conductance case, it describes the evolution of
// conductance between two adjacent cortical spikes
void conductance_decay( neuron *tempneu, int index_neuron, double subTstep);

// same as the above but for smoothed conductance case
void conductance_evolve( neuron *tempneu, int index_neuron, double subTstep);

// this computes the external input current at time t (only for current input case)
double external_current(int index_neuron, double t);

// this computes the derivative of voltage at time t
void voltage_dt(int index_neuron, double t, double gE, double gI, double v, double &dv_dt);

// The following are the Runge-Kutta algorithm with different orders of convergence
void runge_kutta2(neuron *tempneu, int index_neuron, double subTstep, double t_evolution,
                  double &temp_vot,
                  void (*dvdt)(int index_neuron, double t, double gE, double gI, double v,
                               double &dv_dt));

void runge_kutta3(neuron *tempneu, int index_neuron, double subTstep, double t_evolution,
                  double &temp_vot,
                  void (*dvdt)(int index_neuron, double t, double gE, double gI, double v,
                               double &dv_dt));

void runge_kutta4(neuron *tempneu, int index_neuron, double subTstep, double t_evolution,
                  double &temp_vot,
                  void (*dvdt)(int index_neuron, double t, double gE, double gI, double v,
                               double &dv_dt));

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
void spiking_time(neuron *tempneu, int index_neuron, double t_evolution, double Ta,
                  double Tb, double vot_cross, double &firing_time);

// This is for poisson input case!!!
void single_neuron_test(neuron *tempneu, int index_neuron,
                        int firing_neuron, double begin_time, double end_time,
                        int *begin_poisson_index, struct raster &temp_spike_list, int &jump);

// This is for current input case!!!
void single_neuron_test(neuron *tempneu, int index_neuron,
                        int firing_neuron, double begin_time, double end_time,
                        struct raster &temp_spike_list, int &jump);

// this determines the next cortical spike time at the given interval
void next_cortical_spike(neuron *tempneu, double begin_time,
                         int *tempbegin_poisson_index);

// initial setup of the whole neuron system (for reference trajectory)
void network_initialization();

// evolve the whole neuron system (for reference trajectory) from t(0) to t(0)+Tstep
void network_evolution();

// this just evolves the whole system for one time step ahead
void compute_perstep();

void LastRun();

#endif
