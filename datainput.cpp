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
    ///jyl:m,h,n
    if (strcmp(varname,"initial_pertub_m")==0) {
      data_readin>>initial_pertub_m;
      if (verbose) {
        cout<<"initial random seed for fast activation sodium ion channel variable m read to be: "
            <<initial_pertub_m<<endl;
      }
      continue;
    }
    if (strcmp(varname,"initial_pertub_h")==0) {
      data_readin>>initial_pertub_h;
      if (verbose) {
        cout<<"initial random seed for slow inactivation sodium ion channel variable h read to be: "
            <<initial_pertub_h<<endl;
      }
      continue;
    }
    if (strcmp(varname,"initial_pertub_n")==0) {
      data_readin>>initial_pertub_n;
      if (verbose) {
        cout<<"initial random seed for slow activation potassium ion channel variable n read to be: "
            <<initial_pertub_n<<endl;
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
        cout<<"initial random seed for subexcitatory conductance Ex_H read to be: "
            <<initial_pertub_Ex_H<<endl;
      }
      continue;
    }
    if (strcmp(varname,"initial_pertub_Ex_I")==0) {
      data_readin>>initial_pertub_Ex_I;
      if (verbose) {
        cout<<"initial random seed for subexcitatory conductance Ex_I read to be: "
            <<initial_pertub_Ex_I<<endl;
      }
      continue;
    }
    if (strcmp(varname,"initial_pertub_Ex_J")==0) {
      data_readin>>initial_pertub_Ex_J;
      if (verbose) {
        cout<<"initial random seed for subexcitatory conductance Ex_J read to be: "
            <<initial_pertub_Ex_J<<endl;
      }
      continue;
    }
    if (strcmp(varname,"initial_pertub_Ex_K")==0) {
      data_readin>>initial_pertub_Ex_K;
      if (verbose) {
        cout<<"initial random seed for subexcitatory conductance Ex_K read to be: "
            <<initial_pertub_Ex_K<<endl;
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
        cout<<"initial random seed for subinhibitory conductance In_H read to be: "
            <<initial_pertub_In_H<<endl;
      }
      continue;
    }
    if (strcmp(varname,"initial_pertub_In_I")==0) {
      data_readin>>initial_pertub_In_I;
      if (verbose) {
        cout<<"initial random seed for subinhibitory conductance In_I read to be: "
            <<initial_pertub_In_I<<endl;
      }
      continue;
    }
    if (strcmp(varname,"initial_pertub_In_J")==0) {
      data_readin>>initial_pertub_In_J;
      if (verbose) {
        cout<<"initial random seed for subinhibitory conductance In_J read to be: "
            <<initial_pertub_In_J<<endl;
      }
      continue;
    }
    if (strcmp(varname,"initial_pertub_In_K")==0) {
      data_readin>>initial_pertub_In_K;
      if (verbose) {
        cout<<"initial random seed for subinhibitory conductance In_K read to be: "
            <<initial_pertub_In_K<<endl;
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
    if (g_b_verbose_debug)
      std::cerr<<"str:"<<str<<std::endl;
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

// Convert string to array, won't change the value of the one not given
// c_str format: value1[@position1] ...   (position start from 1)
// Return the number of numbers in c_str, or the negitave position of error
// Only have simple error checking, use with careful
int Str2Arr(const char *c_str, double *a, int size)
{
  std::istringstream sin(c_str);
  std::string sst;
  int cnt=0, j=0;
  while (sin>>sst) {
    std::string::size_type idx = sst.find('@');
    if (idx == std::string::npos) {              // not found '@'
      std::istringstream numin(sst);
      if (++j>size) return -cnt-1;               // out of range
      numin>>a[j-1];
    } else {
      sst.at(idx) = ' ';                         // apart the two numbers
      std::istringstream numin(sst);
      double r=0;
      numin>>r>>j;
      if (j>size || j<1) return -cnt-1;          // out of range
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

#define CHK_MEM_RET(v) {if((v)==NULL)return-1;}

//********************************** 2007.10.22 3:15AM *************************
int setglobals()
{
  int i;

  g_num_neu = g_num_neu_ex + g_num_neu_in;
  int RASTER_SIZE = g_num_neu*5;

  raster_initialize(spike_list);
  raster_initialize(RAS);
  raster_allocate(RAS, RASTER_SIZE);

#if CORTICAL_STRENGTH_NONHOMO
  cortical_matrix = (double**)malloc(g_num_neu*sizeof(double*));
  CHK_MEM_RET(cortical_matrix);
  for (i=0; i<g_num_neu; i++) {
    cortical_matrix[i] = (double*)calloc(g_num_neu, sizeof(double));
    CHK_MEM_RET(cortical_matrix[i]);
  }
#endif
  vol          = (double *)malloc(sizeof(double)*g_num_neu);
  CHK_MEM_RET(vol);
  neu          = (neuron *)malloc(sizeof(neuron)*g_num_neu);
  CHK_MEM_RET(neu);
  former_neu   = (neuron *)malloc(sizeof(neuron)*g_num_neu);
  CHK_MEM_RET(former_neu);
  neuRK        = (neuron *)malloc(sizeof(neuron)*g_num_neu);
  CHK_MEM_RET(neuRK);
  neu_d1       = (neuron *)malloc(sizeof(neuron)*g_num_neu);
  CHK_MEM_RET(neu_d1);
  neu_d2       = (neuron *)malloc(sizeof(neuron)*g_num_neu);
  CHK_MEM_RET(neu_d2);
  neu_d3       = (neuron *)malloc(sizeof(neuron)*g_num_neu);
  CHK_MEM_RET(neu_d3);
  neu_d4       = (neuron *)malloc(sizeof(neuron)*g_num_neu);
  CHK_MEM_RET(neu_d4);
  for (i=0; i<g_num_neu; i++) {
    neuron_initialize(neu[i]);
    neuron_initialize(former_neu[i]);
    neuron_initialize(neuRK[i]);
    neuron_initialize(neu_d1[i]);
    neuron_initialize(neu_d2[i]);
    neuron_initialize(neu_d3[i]);
    neuron_initialize(neu_d4[i]);
  }

  // intialize the voltage, Exconductance and Inconductance
  for( i = 0; i < g_num_neu; i++ ) {
    if ( i < g_num_neu_ex ) {
      neuron_set_value(neu[i],          Type_Exneuron, STATE_NORMAL, Stepsmooth_Con);
      neuron_set_value(former_neu[i],   Type_Exneuron, STATE_NORMAL, Stepsmooth_Con);
      neuron_set_value(neuRK[i],        Type_Exneuron, STATE_NORMAL, Stepsmooth_Con);
      neuron_set_value(neu_d1[i],       Type_Exneuron, STATE_NORMAL, Stepsmooth_Con);
      neuron_set_value(neu_d2[i],       Type_Exneuron, STATE_NORMAL, Stepsmooth_Con);
      neuron_set_value(neu_d3[i],       Type_Exneuron, STATE_NORMAL, Stepsmooth_Con);
      neuron_set_value(neu_d4[i],       Type_Exneuron, STATE_NORMAL, Stepsmooth_Con);

    } else {
      neuron_set_value(neu[i],          Type_Inneuron, STATE_NORMAL, Stepsmooth_Con);
      neuron_set_value(former_neu[i],   Type_Inneuron, STATE_NORMAL, Stepsmooth_Con);
      neuron_set_value(neuRK[i],        Type_Inneuron, STATE_NORMAL, Stepsmooth_Con);
      neuron_set_value(neu_d1[i],       Type_Inneuron, STATE_NORMAL, Stepsmooth_Con);
      neuron_set_value(neu_d2[i],       Type_Inneuron, STATE_NORMAL, Stepsmooth_Con);
      neuron_set_value(neu_d3[i],       Type_Inneuron, STATE_NORMAL, Stepsmooth_Con);
      neuron_set_value(neu_d4[i],       Type_Inneuron, STATE_NORMAL, Stepsmooth_Con);
    }
    CHK_MEM_RET(neu[i].value);
    CHK_MEM_RET(former_neu[i].value);
    CHK_MEM_RET(neuRK[i].value);
    CHK_MEM_RET(neu_d1[i].value);
    CHK_MEM_RET(neu_d2[i].value);
    CHK_MEM_RET(neu_d3[i].value);
    CHK_MEM_RET(neu_d4[i].value);
  }

#if POISSON_INPUT_USE
  // poission_input[i] records the timing of poisson input spikes
  // in some time interval for i'th neuron!
  poisson_input = (vector *)malloc(sizeof(vector)*g_num_neu);
  CHK_MEM_RET(poisson_input);
  for (i=0; i<g_num_neu; i++) {
    vector_initialize(poisson_input[i]);
    vector_set_value(poisson_input[i], Maxnum_input, 0.0);
    CHK_MEM_RET(poisson_input[i].vect_value);
  }

  // record the initial random seed for each neuron!
  initialseed_neuron = (long *)malloc(sizeof(long)*g_num_neu);
  last_input = (double *)malloc(sizeof(double)*g_num_neu);
  CHK_MEM_RET(last_input);
  CHK_MEM_RET(initialseed_neuron);

  ran_iy = (long *)malloc(sizeof(long)*g_num_neu);
  ran_iv = (long **)malloc(sizeof(long*)*g_num_neu);
  CHK_MEM_RET(ran_iy);
  CHK_MEM_RET(ran_iv);

  for (i=0; i<g_num_neu; i++) {
    ran_iy[i] = 0;
    ran_iv[i] = (long*)malloc(sizeof(long)*NTAB);
  }

  g_begin_poisson_index = (int *)malloc(sizeof(int)*g_num_neu);

  g_arr_poisson_rate = (double *)malloc(g_num_neu*sizeof(double));
  g_arr_poisson_strength_E = (double *)malloc(g_num_neu*sizeof(double));
  g_arr_poisson_strength_I = (double *)malloc(g_num_neu*sizeof(double));

  CHK_MEM_RET((int*)(
     g_begin_poisson_index!=NULL       && g_arr_poisson_rate!=NULL
    && g_arr_poisson_strength_E!=NULL    && g_arr_poisson_strength_I!=NULL
    ));

  for (int j=0; j<g_num_neu; j++) {
    g_arr_poisson_rate[j] = 1;
    g_arr_poisson_strength_E[j] = 1;
    g_arr_poisson_strength_I[j] = 1;
  }
#else
  // each neuron has different phase for external input current!
  phase = (double *)malloc(g_num_neu*sizeof(double));
  CHK_MEM_RET(phase);
  for (i=0; i<g_num_neu; i++) {
    phase[i] = 2*M_PI*i/g_num_neu;
  }
#endif

  GLOBAL_STRA = (struct strobe **)tcalloc(
    g_num_neu*size_neuronvar/*number of variables*/, sizeof(struct strobe *));

  return 0;
}

void input_initialization()
{
  int i;
  int NONE_INHIBITORY;

  if (g_num_neu_in) NONE_INHIBITORY = 0; else NONE_INHIBITORY = 1;

  for( i = 0; i < g_num_neu; i++ ) {
#if EXPONENTIAL_IF_USE
    neu[i].value[0] = 0.5;//VOT_RESET+ran0(&initial_pertub_Vot)
                     //*(VOT_TAKEOFF-VOT_RESET);
#else // use I&F model
    neu[i].value[0] = ran0(&initial_pertub_Vot);
#endif
//    neu[i].value[0] = ran0(&initial_pertub_Vot);
    neu[i].value[1] = ran0(&initial_pertub_Ex);
    neu[i].value[Stepsmooth_Con+1] = ran0(&initial_pertub_In);
    neu[i].value[2*Stepsmooth_Con+1] = ran0(&initial_pertub_m);
    neu[i].value[2*Stepsmooth_Con+2] = ran0(&initial_pertub_h);
    neu[i].value[2*Stepsmooth_Con+3] = ran0(&initial_pertub_n);
#if SMOOTH_CONDUCTANCE_USE
    neu[i].value[2] = ran0(&initial_pertub_Ex_H);
    neu[i].value[3] = ran0(&initial_pertub_Ex_I);
    neu[i].value[4] = ran0(&initial_pertub_Ex_J);
    neu[i].value[Stepsmooth_Con] = ran0(&initial_pertub_Ex_K);
    neu[i].value[Stepsmooth_Con+2] = (1 - NONE_INHIBITORY)*ran0(&initial_pertub_In_H);
    neu[i].value[Stepsmooth_Con+3] = (1 - NONE_INHIBITORY)*ran0(&initial_pertub_In_I);
    neu[i].value[Stepsmooth_Con+4] = (1 - NONE_INHIBITORY)*ran0(&initial_pertub_In_J);
    neu[i].value[2*Stepsmooth_Con] = (1 - NONE_INHIBITORY)*ran0(&initial_pertub_In_K);
    //xyy: for test
    neu[i].value[0] = 0.000027756626542950875;
    neu[i].value[1] = 0;  // gE
    neu[i].value[2] = 0;  // hE
    neu[i].value[3] = 0;  // gI
    neu[i].value[4] = 0;  // hI
    neu[i].value[5] = 0.052934217620863984;  // m
    neu[i].value[6] = 0.5961110463468279;  // h
    neu[i].value[7] = 0.31768116757978115;  // n
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
    neuron_destroy(neu[i]);
#if CORTICAL_STRENGTH_NONHOMO
    free(cortical_matrix[i]);
#endif
#if POISSON_INPUT_USE
    free(ran_iv[i]);
    vector_destroy(poisson_input[i]);
#endif
  }
  raster_destroy(RAS);
  raster_destroy(spike_list);
  free(neu);                          neu = NULL;
  free(former_neu);                   former_neu = NULL;
  free(neuRK);                        neuRK = NULL;
  free(neu_d1);                       neu_d1 = NULL;
  free(neu_d2);                       neu_d2 = NULL;
  free(neu_d3);                       neu_d3 = NULL;
  free(neu_d4);                       neu_d4 = NULL;

#if POISSON_INPUT_USE
  free(poisson_input);                poisson_input = NULL;
  free(g_arr_poisson_rate);           g_arr_poisson_rate = NULL;
  free(g_arr_poisson_strength_E);     g_arr_poisson_strength_E = NULL;
  free(g_arr_poisson_strength_I);     g_arr_poisson_strength_I = NULL;
  free(g_begin_poisson_index);        g_begin_poisson_index = NULL;
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
