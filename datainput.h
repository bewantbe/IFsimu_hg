#ifndef _DATAINPUT_H_
#define _DATAINPUT_H_

// read data from input file to set up some global variables
void readinput(char *filename);

int read_cortical_matrix(const char *filepath, double **cor_mat);
int Str2Arr(const char *c_str, double *a);
int ReadPR(const char *filename, double *r);
int ReadPS(const char *filename, double *ae, double *ai);

// initial setup of global variables, such as assigning space ...
void setglobals();

// initial setup for the external drive, the neurons' information, ...
void input_initialization();

#if POISSON_INPUT_USE
// generate poisson input spike train
void poisson_generator(int index_neuron, vector &Timing_input);
#endif

// release all data structure
void data_dump();

#endif
