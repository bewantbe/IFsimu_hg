#include "stdafx.h"
#include "random.h"
#include "datainput.h"
#include "datahandling.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

using std::cout;
using std::endl;

int readinput(char *filename)
{
  int verbose = g_b_verbose_debug;
  char varname[128];
  double temp = 0.0;

  if (filename==NULL || filename[0]=='\0')  // no path?
    return 1;

  std::ifstream data_readin(filename);

  if (!data_readin.good())  // bad file name?
    return -1;

  while (!data_readin.eof()) {
    data_readin>>varname;
    data_readin.ignore(2);

    if (strcmp(varname,"Tstep")==0) {
      data_readin>>Tstep;
      if (verbose) {
        cout<<"Time step read to be: "<<Tstep<<endl;
      }
      continue;
    }
    if (strcmp(varname,"Rate_input")==0) {
      data_readin>>Rate_input;
      if (verbose) {
        cout<<"Poisson input rate read to be: "<<Rate_input<<endl;
      }
      continue;
    }
#if POISSON_INPUT_USE
    if (strcmp(varname,"Strength_Exinput")==0) {
      data_readin>>Strength_Exinput;
      if (verbose) {
        cout<<"Poisson input strength to Ex. neurons read to be: "
            <<Strength_Exinput<<endl;
      }
      continue;
    }
    if (strcmp(varname,"Strength_Ininput")==0) {
      data_readin>>Strength_Ininput;
      if (verbose) {
        cout<<"Poisson input strength to In. neurons read to be: "
            <<Strength_Ininput<<endl;
      }
      continue;
    }
#else
    if (strcmp(varname,"average_current")==0) {
      data_readin>>Current_0;
      if (verbose) {
        cout<<"average external current read to be: "
            <<Current_0<<endl;
      }
      continue;
    }
    if (strcmp(varname,"amplitude_current")==0) {
      data_readin>>Current_1;
      if (verbose) {
        cout<<"average external current read to be: "
            <<Current_1<<endl;
      }
      continue;
    }
#endif
    if (strcmp(varname,"Strength_CorEE")==0) {
      data_readin>>Strength_CorEE;
      if (verbose) {
        cout<<"Cortical strength of Ex. to Ex. neurons read to be: "
            <<Strength_CorEE<<endl;
      }
      continue;
    }
    if (strcmp(varname,"Strength_CorIE")==0) {
      data_readin>>Strength_CorIE;
      if (verbose) {
        cout<<"Cortical strength of Ex. to In. neurons read to be: "
            <<Strength_CorIE<<endl;
      }
      continue;
    }
    if (strcmp(varname,"Strength_CorII")==0) {
      data_readin>>Strength_CorII;
      if (verbose) {
        cout<<"Cortical strength of In. to In. neurons read to be: "
            <<Strength_CorII<<endl;
      }
      continue;
    }
    if (strcmp(varname,"Strength_CorEI")==0) {
      data_readin>>Strength_CorEI;
      if (verbose) {
        cout<<"Cortical strength of In. to Ex. neurons read to be: "
            <<Strength_CorEI<<endl;
      }
      continue;
    }
    if (strcmp(varname,"initial_pertub_Vot")==0) {
      data_readin>>initial_pertub_Vot;
      if (verbose) {
        cout<<"initial random seed for voltage read to be: "
            <<initial_pertub_Vot<<endl;
      }
      continue;
    }
    if (strcmp(varname,"initial_pertub_Ex")==0) {
      data_readin>>initial_pertub_Ex;
      if (verbose) {
        cout<<"initial random seed for excitatory conductance read to be: "
            <<initial_pertub_Ex<<endl;
      }
      continue;
    }
#if SMOOTH_CONDUCTANCE_USE
    if (strcmp(varname,"initial_pertub_Ex_H")==0) {
      data_readin>>initial_pertub_Ex_H;
      if (verbose) {
        cout<<"initial random seed for subexcitatory conductance read to be: "
            <<initial_pertub_Ex_H<<endl;
      }
      continue;
    }
#endif
    if (strcmp(varname,"initial_pertub_In")==0) {
      data_readin>>initial_pertub_In;
      if (verbose) {
        cout<<"initial random seed for inhibitory conductance read to be: "
            <<initial_pertub_In<<endl;
      }
      continue;
    }
#if SMOOTH_CONDUCTANCE_USE
    if (strcmp(varname,"initial_pertub_In_H")==0) {
      data_readin>>initial_pertub_In_H;
      if (verbose) {
        cout<<"initial random seed for subinhibitory conductance read to be: "
            <<initial_pertub_In_H<<endl;
      }
      continue;
    }
#endif
    if (strcmp(varname,"initial_seed")==0) {
      data_readin>>initial_seed;
      if (verbose) {
        cout<<"initial seed for poisson input generator read to be: "
            <<initial_seed<<endl;
      }
      continue;
    }
    if (strcmp(varname,"last_time")==0) {
      data_readin>>last_time;
      if (verbose) {
        cout<<"ending time of the last run read to be: "
            <<last_time<<endl;
      }
      continue;
    }
    // if there is no corresponding variable, just go to the next line to read!!
    data_readin>>temp;
    if (verbose) {
      cout<<"no variables to match!"<<endl;
    }
  }
  data_readin.close();

  return 0;
}

#if CORTICAL_STRENGTH_NONHOMO
int read_cortical_matrix(const char *filepath, double **cor_mat)
{
  // all set to NaN, for error checking below
  for (int i=0; i<g_num_neu; i++)
    memset(cor_mat[i], 255, g_num_neu*sizeof(double));

  std::ifstream fin(filepath);
  std::string str;
  int i;
  for (i=0; i<g_num_neu; i++) {       // read the matrix
    getline(fin, str);
    std::istringstream sin(str);
    int j;
    for (j=0; j<g_num_neu; j++) {
      sin>>cor_mat[i][j];
      if (isnan(cor_mat[i][j])) break;
    }
    if (j!=g_num_neu) break;
  }
  if (!fin.good() || i != g_num_neu) {
    fin.close();
    if (filepath && filepath[0]=='-' && filepath[1]=='\0') {
      for (int i=0; i<g_num_neu; i++)
        for (int j=0; j<g_num_neu; j++)
          cor_mat[i][j] = 1;
      return 1;
    } else {
      return -1;
    }
  }
  fin.close();

  return 0;
}
#endif

bool IsEmptyStr(const char *st)
{
  int j=0;
  char c;
  while ((c=st[j])) {
    if (c!=' ' && c!='\t' && c!='\n') break;
    j++;
  }
  return c=='\0';
}

// conver string to array, won't change the value of the one not given
// return the number of numbers in st
// No error checking, use with careful
int Str2Arr(const char *c_str, double *a, int size)
{
  std::istringstream sin(c_str);
  std::string sst;
  int cnt=0;
  while (sin>>sst) {
    std::string::size_type idx = sst.find('@');
    if (idx == std::string::npos) {              // not found '@'
      std::istringstream numin(sst);
      if (cnt>=size) return -cnt-1;    // out of range
      numin>>a[cnt];
    } else {
      sst.at(idx) = ' ';
      std::istringstream numin(sst);
      double r=0;
      int j=0;
      numin>>r>>j;
      if (j>size || j<1) return -cnt-1;       // out of range
      a[j-1] = r;
    }
    cnt++;
  }
  return cnt;
}

int ReadPR(const char *filename, double *r, int size)
{
  if (filename==NULL || filename[0]=='\0') {
//    for (int j=0; j<g_num_neu; j++) r[j] = 1;
    return 1;
  }
  std::ifstream fin(filename);
  std::string str;
  getline(fin, str);
  return (Str2Arr(str.c_str(), r, size)==g_num_neu) ? 0 : -1;
}

// Read poisson input strength from file filename
int ReadPS(const char *filename, double *ae, double *ai)
{
  if (filename==NULL || filename[0]=='\0') {
//    for (int j=0; j<g_num_neu; j++) {
//      ae[j] = 1;
//      ai[j] = 1;
//    }
    return 1;
  }

  std::ifstream fin(filename);
  std::string str;
  getline(fin, str);
  std::istringstream sin(str);
  int j=0;                      // reads ae
  while (sin>>ae[j]) {
    j++;
    if (j>g_num_neu) break;
  }
  if (j!=g_num_neu) return -1;

  str.clear();                  // is there data for ai?
  getline(fin, str);
  if (IsEmptyStr(str.c_str())) {
    for (j=0; j<g_num_neu; j++)
      ai[j] = 1;
    return 0;
  }

  j=0;                          // reads ai
  while (sin>>ai[j]) {
    j++;
    if (j>g_num_neu) break;
  }
  if (j!=g_num_neu) return -1;
  return 0;
}

// temp space used in next_cortical_spike @ poisson_input.cpp
neuron *tmp_tempneu = NULL;
// temp space used in network_initialization, network_evolution @ poisson_input.cpp
neuron *tmp_tempneu2 = NULL;
//
int *tmp_tempbegin_poisson_index = NULL;
// temporal vector to instore the poisson input timing in this interval
//vector tmp_vec1_poisson_generator;
// uesd in network_initialization() and network_evolution()
int *g_tempbegin_poisson_index = NULL;
int *g_begin_poisson_index = NULL;

#define ALLOCRT(rt,v) rt=(int)(v!=NULL)&(rt)?-1:0;

//********************************** 2007.10.22 3:15AM *************************
int setglobals()
{
//  const char *err_st_mem = "Error: setglobals(): Fail to allocate memory.";
  int i, rt=-1;

  g_num_neu = g_num_neu_ex + g_num_neu_in;
  int RASTER_SIZE = g_num_neu*5;

  raster_initialize(spike_list);
  raster_initialize(RAS);
  raster_allocate(RAS, RASTER_SIZE);

#if CORTICAL_STRENGTH_NONHOMO
  cortical_matrix = (double**)malloc(g_num_neu*sizeof(double*));
  ALLOCRT(rt, cortical_matrix);
//  P_NULL_ERR(cortical_matrix, err_st_mem);
  for (i=0; i<g_num_neu; i++) {
    cortical_matrix[i] = (double*)calloc(g_num_neu, sizeof(double));
    ALLOCRT(rt, cortical_matrix[i]);
//    P_NULL_ERR(cortical_matrix[i], err_st_mem);
  }
#endif

  neu          = (neuron *)malloc(sizeof(neuron)*g_num_neu);
  tmp_tempneu  = (neuron *)malloc(sizeof(neuron)*g_num_neu);
  tmp_tempneu2 = (neuron *)malloc(sizeof(neuron)*g_num_neu);
  ALLOCRT(rt, neu);
  ALLOCRT(rt, tmp_tempneu);
  ALLOCRT(rt, tmp_tempneu2);
//  P_NULL_ERR(neu,          err_st_mem);
//  P_NULL_ERR(tmp_tempneu,  err_st_mem);
//  P_NULL_ERR(tmp_tempneu2, err_st_mem);
  for (i=0; i<g_num_neu; i++) {
    neuron_initialize(neu[i]);
    neuron_initialize(tmp_tempneu[i]);
    neuron_initialize(tmp_tempneu2[i]);
  }

  // intialize the voltage, Exconductance and Inconductance
  for( i = 0; i < g_num_neu; i++ ) {
    if ( i < g_num_neu_ex ) {
      neuron_set_value(neu[i],          Type_Exneuron, STATE_ACTIVE, Stepsmooth_Con);
      neuron_set_value(tmp_tempneu[i],  Type_Exneuron, STATE_ACTIVE, Stepsmooth_Con);
      neuron_set_value(tmp_tempneu2[i], Type_Exneuron, STATE_ACTIVE, Stepsmooth_Con);
    } else {
      neuron_set_value(neu[i],          Type_Inneuron, STATE_ACTIVE, Stepsmooth_Con);
      neuron_set_value(tmp_tempneu[i],  Type_Inneuron, STATE_ACTIVE, Stepsmooth_Con);
      neuron_set_value(tmp_tempneu2[i], Type_Inneuron, STATE_ACTIVE, Stepsmooth_Con);
    }
    ALLOCRT(rt, neu[i].value);
    ALLOCRT(rt, tmp_tempneu[i].value);
    ALLOCRT(rt, tmp_tempneu2[i].value);
//    P_NULL_ERR(neu[i].value,          err_st_mem);
//    P_NULL_ERR(tmp_tempneu[i].value,  err_st_mem);
//    P_NULL_ERR(tmp_tempneu2[i].value, err_st_mem);
  }

#if POISSON_INPUT_USE
  // poission_input[i] records the timing of poisson input spikes
  // in some time interval for i'th neuron!
  poisson_input = (vector *)malloc(sizeof(vector)*g_num_neu);
  ALLOCRT(rt, poisson_input);
//  P_NULL_ERR(poisson_input, err_st_mem);
  for (i=0; i<g_num_neu; i++) {
    vector_initialize(poisson_input[i]);
    vector_set_value(poisson_input[i], Maxnum_input, 0.0);
    ALLOCRT(rt, poisson_input[i].vect_value);
//    P_NULL_ERR(poisson_input[i].vect_value, err_st_mem);
  }

  // record the initial random seed for each neuron!
  initialseed_neuron = (long *)malloc(sizeof(long)*g_num_neu);
  last_input = (double *)malloc(sizeof(double)*g_num_neu);
  ALLOCRT(rt, last_input);
  ALLOCRT(rt, initialseed_neuron);
//  P_NULL_ERR(last_input, err_st_mem);
//  P_NULL_ERR(initialseed_neuron, err_st_mem);

  ran_iy = (long *)malloc(sizeof(long)*g_num_neu);
  ran_iv = (long **)malloc(sizeof(long*)*g_num_neu);
  ALLOCRT(rt, ran_iy);
  ALLOCRT(rt, ran_iv);
//  P_NULL_ERR(ran_iy, err_st_mem);
//  P_NULL_ERR(ran_iv, err_st_mem);

  for (i=0; i<g_num_neu; i++) {
    ran_iy[i] = 0;
    ran_iv[i] = (long*)malloc(sizeof(long)*NTAB);
  }

  tmp_tempbegin_poisson_index = (int *)malloc(sizeof(int)*g_num_neu);
  g_tempbegin_poisson_index = (int *)malloc(sizeof(int)*g_num_neu);
  g_begin_poisson_index = (int *)malloc(sizeof(int)*g_num_neu);

  g_arr_poisson_rate = (double *)malloc(g_num_neu*sizeof(double));
  g_arr_poisson_strength_E = (double *)malloc(g_num_neu*sizeof(double));
  g_arr_poisson_strength_I = (double *)malloc(g_num_neu*sizeof(double));

  ALLOCRT(rt, (void*)(
      (long)tmp_tempbegin_poisson_index  & (long)g_tempbegin_poisson_index
    & (long)g_begin_poisson_index        & (long)g_arr_poisson_rate
    & (long)g_arr_poisson_strength_E     & (long)g_arr_poisson_strength_I
    ));

  for (int j=0; j<g_num_neu; j++) {
    g_arr_poisson_rate[j] = 1;
    g_arr_poisson_strength_E[j] = 1;
    g_arr_poisson_strength_I[j] = 1;
  }
#else
  // each neuron has different phase for external input current!
  phase = (double *)malloc(g_num_neu*sizeof(double));
  ALLOCRT(rt, phase);
//  P_NULL_ERR(phase, err_st_mem);
  for (i=0; i<g_num_neu; i++) {
    phase[i] = 2*M_PI*i/g_num_neu;
  }
#endif

  GLOBAL_STRA = (struct strobe **)tcalloc(
    g_num_neu*size_neuronvar/*number of variables*/, sizeof(struct strobe *));

  return -(1+rt);   // -1 when fail, 0 when success
}

void input_initialization()
{
  int i;
  int NONE_INHIBITORY;

  if (g_num_neu_in) NONE_INHIBITORY = 0; else NONE_INHIBITORY = 1;

  for( i = 0; i < g_num_neu; i++ ) {
    neu[i].value[0] = ran0(&initial_pertub_Vot);
    neu[i].value[1] = ran0(&initial_pertub_Ex);
    neu[i].value[Stepsmooth_Con+1] = (1 - NONE_INHIBITORY)*ran0(&initial_pertub_In);
#if SMOOTH_CONDUCTANCE_USE
    neu[i].value[Stepsmooth_Con] = ran0(&initial_pertub_Ex_H);
    neu[i].value[2*Stepsmooth_Con] = (1 - NONE_INHIBITORY)*ran0(&initial_pertub_In_H);
#endif
  }
#if POISSON_INPUT_USE
  srand(initial_seed);
  for( i = 0; i < g_num_neu; i++ ) {
    initialseed_neuron[i] = (long) -rand();
//    last_input[i] = -log(RANDOM(initialseed_neuron+i,ran_iy,ran_iv,i))/Rate_input;
    last_input[i] = -log(RANDOM(initialseed_neuron+i,ran_iy,ran_iv,i))/g_arr_poisson_rate[i];
  }
#endif

  for (i=0; i<g_num_neu*size_neuronvar; i++) {
    GLOBAL_STRA[i]=strobemake(WINDOW_BIN_LENGTH/*ms*/,
      SLIGHT_BIN/*ms*/,0/*we don't want to maintain cyclical average*/);
  }

  time_evolution = last_time;
}

#if POISSON_INPUT_USE
// generate poisson input spikes for each neuron in each time interval!
// must make sure Timing_input has at least Maxnum_input space
void poisson_generator(int index_neuron, vector &Timing_input)
{
  // " >= " means dealing with the input spike at the starting point of the time interval!
  if (last_input[index_neuron] >= Tstep) {
    // there is no poisson spike input in this timestep interval!
    last_input[index_neuron] -= Tstep; // continue to deal with the next interval!
    // set the first value less than zero: no input spikes in this interval
    Timing_input.vect_size = -1;
    Timing_input.vect_index = 0;
    return;
  }

  int size = 0;
  double spiketiming_wait = last_input[index_neuron];
  do {
    Timing_input.vect_value[size] = spiketiming_wait;
    spiketiming_wait += -log(RANDOM(initialseed_neuron+index_neuron,ran_iy,ran_iv,index_neuron))/g_arr_poisson_rate[index_neuron];
    size++;
    if(size >= Maxnum_input) {
      std::cerr<<"Error: There are more spikes exceeding the capacity of vector!"<<endl;
      std::cerr<<"Please enlarge the Maximum size of vector!"<<endl;
      exit(1);
    }
  } while (spiketiming_wait < Tstep);
  // Timing_input contains all poisson input spikes in this time interval!
  // the first value is recorded as last_input[index_neuron]
  // the last value is recorded as the (size-1)th value of vect1
  Timing_input.vect_size = size;
  Timing_input.vect_index = 0;

  // last_input records the nearest poisson input spike waiting time beyond this time interval!
  last_input[index_neuron] = spiketiming_wait - Tstep;
}
#endif

void data_dump()
{
  for(int i = 0; i < g_num_neu; i++ ) {
    neuron_destroy(tmp_tempneu[i]);
    neuron_destroy(tmp_tempneu2[i]);
    neuron_destroy(neu[i]);
#if CORTICAL_STRENGTH_NONHOMO
    free(cortical_matrix[i]);
#endif
#if POISSON_INPUT_USE
    free(ran_iv[i]);
    vector_destroy(poisson_input[i]);
#endif
  }
  free(tmp_tempneu);                  tmp_tempneu = NULL;
  free(tmp_tempneu2);                 tmp_tempneu2 = NULL;
  raster_destroy(RAS);
  raster_destroy(spike_list);
  free(neu);                          neu = NULL;
#if POISSON_INPUT_USE
  free(poisson_input);                poisson_input = NULL;
  free(g_arr_poisson_rate);           g_arr_poisson_rate = NULL;
  free(g_arr_poisson_strength_E);     g_arr_poisson_strength_E = NULL;
  free(g_arr_poisson_strength_I);     g_arr_poisson_strength_I = NULL;
  free(g_tempbegin_poisson_index);    g_tempbegin_poisson_index = NULL;
  free(g_begin_poisson_index);        g_begin_poisson_index = NULL;
  free(tmp_tempbegin_poisson_index);  tmp_tempbegin_poisson_index = NULL;
  free(ran_iy);                       ran_iy = NULL;
  free(ran_iv);                       ran_iv = NULL;
  free(last_input);                   last_input = NULL;
  free(initialseed_neuron);           initialseed_neuron = NULL;
#else
  free(phase);                        phase = NULL;
#endif
#if CORTICAL_STRENGTH_NONHOMO
  free(cortical_matrix);              cortical_matrix = NULL;
#endif
}
