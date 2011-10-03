#include "stdafx.h"
#include "neuron.h"

void neuron_destroy (neuron &neu)
{
  if(neu.is_allocated == true) {
    free(neu.value);
    neu.value = NULL;
    neu.is_allocated = false;
  }
}

void neuron_initialize(neuron &neu)
{
  neu.value = NULL;
  neu.is_allocated = false;
}

void neuron_set_value(neuron &neu, int typ, int state_neu, int step_con)
{
  neuron_destroy(neu);

  neu.type = typ;
  neu.state_neuron = state_neu;
  neu.step_conductance = step_con;

  neu.value = (double *)calloc(2*step_con+2, sizeof(double));
  P_NULL_ERR(neu.value, "Error: neuron_set_value: Allocation failed.");

  neu.is_allocated = true;
}

void neuron_copy(neuron &neu, const neuron &src)
{
  // free the allocated memory for neu first
  if(neu.is_allocated == true) {
    neuron_destroy(neu);
  }
  neu.is_allocated = src.is_allocated;
  if(neu.is_allocated == false) {
    fprintf(stderr, "the source neuron is null!\n");
    return;
  } else {
    neu.type = src.type;
    neu.state_neuron = src.state_neuron;
    neu.step_conductance = src.step_conductance;
    neu.value = (double*)malloc((2*src.step_conductance+2)*sizeof(double));
    assert(neu.value);
    memcpy(neu.value, src.value, (2*src.step_conductance+2)*sizeof(double));
  }
}

void neuron_copy_2(neuron &neu, const neuron &src)
{
  if(src.is_allocated == false) {
    fprintf(stderr, "the source neuron is null!\n");
    return;
  }
  if(neu.is_allocated == false) {
    neu.value = (double*)malloc((2*src.step_conductance+2)*sizeof(double));
    assert(neu.value);
    neu.is_allocated = true;
  }
  neu.type = src.type;
  neu.state_neuron = src.state_neuron;
  neu.step_conductance = src.step_conductance;
  memcpy(neu.value, src.value, (2*src.step_conductance+2)*sizeof(double));
}

// copy data, assum neu and src have the same type
void neuron_copy_raw_static(neuron &neu, const neuron &src)
{
//  neu.type = src.type;
  neu.state_neuron = src.state_neuron;
//  neu.step_conductance = Stepsmooth_Con;
  memcpy(neu.value, src.value, (2*Stepsmooth_Con+2)*sizeof(double));
}
