#ifndef _NEURON_H_
#define _NEURON_H_

// this is structue to save neuron's information
struct neuron {
  int type;   // 1 for Excitatory and 0 for Inhibitory
  int state_neuron; // 1 for active and 0 for refractory

  /* this indicates whether use smoothed conductance or not, if it equals 1, then
   the conductance is described only by one ODE with infinite rising time scale
   and a finite decay time scale; if it equals 2, then the conductance is
   described by two ODEs with both finite rising and decay time scale. */
  int step_conductance; // step of smoothness for conductance

  // bool variable to indicate that whether this data structure has
  // been initialized or not
  bool is_allocated;

  // voltage, Ex_conductance with multi-step smoothness
  // In_conductance with multi-step smoothness, ref_time
  // size equals (2*step_conductance+2)
  double * value;
};

// destroy the data structure
void neuron_destroy (neuron &neu);

// initial setup of data structure
void neuron_initialize(neuron& neu);

// set up data value in the structure
void neuron_set_value(neuron& neu, int typ, int state_neu, int step_con);

// copy data from another data structure
void neuron_copy(neuron &neu, const neuron &src);
void neuron_copy_2(neuron &neu, const neuron &src);
void neuron_copy_raw_static(neuron &neu, const neuron &src);

#endif
