// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently

#if !defined(AFX_STDAFX_H__A9DB83DB_A9FD_11D0_BFD1_444553540000__INCLUDED_)
#define AFX_STDAFX_H__A9DB83DB_A9FD_11D0_BFD1_444553540000__INCLUDED_

#define WIN32_LEAN_AND_MEAN       // Exclude rarely-used stuff from Windows headers
#define _CRT_SECURE_NO_DEPRECATE  // prevent MSVS say it's unsafe to use scanf etc

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>

#if defined(_WIN32) || defined(__WIN32__)
#define _WINDOWS_USE_
#endif

#ifdef _WINDOWS_USE_
#include <windows.h>
#else  // Linux, no Mac OS support
#include <unistd.h>
#define Sleep(msec) usleep((msec)*1000)  // replace Sleep() in winbase.h (in Windows)
#endif

#ifndef isnan
static inline int isnan(double x) { return x != x; }
#endif

//#ifndef INT_MAX
//#if __WORDSIZE != 16
//#define INT_MAX 2147483647
//#define UINT_MAX 4294967295U
//#else
//#define INT_MAX 32767
//#define INT_MAX 65535U
//#endif
//#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include "raster.h"
#include "neuron.h"
#include "loop.h"

// these are used in open GL
#define ESCAPE_KEY 27
#define SPACE_KEY 32
#define ENTER_KEY 13
#define Q_KEY 113
#define NUM_NEU_MAX (32767)

// this is used to locate the drawn menu position when using open GL
#define menu_top 3

#define PI 3.141592653589793
#define maximum(A,B) ((A) > (B) ? (A) : (B))
#define minimum(A,B) ((A) < (B) ? (A) : (B))
// usually B<C, keep B<=A<=C with up/down limite
#define crop(A,B,C) ((A) < (B) ? (B) : ((A) > (C) ? (C) : (A)))
// usually B<C, keep B<=A< C using periodization, i.e. A mod (C-B) + B, but only one period
#define periodize(A,B,C) ((A) < (B) ? (A) + (C) - (B) : ((A) >= (C) ? (A) - (C) + (B) : (A)))

// whether use poisson input or current input
#define POISSON_INPUT_USE 1              ///XYY: modify this
// whether consider the homogeneous cortical strength (a scalar independent of
// neuron index) or not! It is useful when using a N-by-N matrix to indicate the
// strength between neurons
#define CORTICAL_STRENGTH_NONHOMO 1      ///XYY: modify this
// the number of excitatory neurons
#define Number_Exneuron 2                ///XYY: modify this
// the number of inhibitory neurons
#define Number_Inneuron 0
// the final stopping time of computation in millisecond
#define COMP_TIME  100000.0               /// modify this
// saving time length (the interval is [T-20, T]) for voltage, conductance
// and also new raster events
#define TIME_LENGTH 20.0

// whether use the smoothed conductance (meaning that we use more than one ODE to
// describe the rising and decay time scale of conductance)
#define SMOOTH_CONDUCTANCE_USE 0          /// modify this

#if SMOOTH_CONDUCTANCE_USE
// this means using two ODEs to describe the dynamics of conductance
// therefore, the conductance is smooth with both rising and decay time scale
#define Stepsmooth_Con 2
#else
// this means using only one ODE to describe the dynamics of condcuctance
// therefore, the conductance only has decay time scale and infinitely fast
// rising time scale (instantaneous jump)
#define Stepsmooth_Con 1
#endif
// the length of "data" stored in the neuron structure
// the first is saving voltage, the last is saving refractory residence time
// the middle saves the conductance part, there are maybe three cases:
// (1) Ex_conductance after two step of smoothness, (2) Ex_conductance after one
//     step of smoothness, (3) Ex_conductance without smoothness
// (1) In_conductance after two step of smoothness, (2) In_conductance after one
//     step of smoothness, (3) In_conductance without smoothness
#define size_neuronvar  (2*Stepsmooth_Con+2)

// decay time scale of Ex_conductance
#define Time_ExCon      2.0
// if we use smoothed excitatory conductance instead of jump conductance, we
// should have another ODE. this is the time scale in that ODE
#define Time_ExConR     0.5
// decay time scale of In_conductance
#define Time_InCon      5.0
// the same thing but for inhibitory conductance
#define Time_InConR     0.8
// the length of refractory period
#define TIME_REFRACTORY 2.0
// the reset value for voltage after spike
#define VOT_RESET       0.0
// the threshold value for voltage
#define Vot_Threshold   1.0
// the resting value of leakage conductance
#define Vot_Leakage     0.0
// the resting value of excitatory conductance
#define Vot_Excitatory  (14.0/3.0)
// the resting value of inhibitory conductance
#define Vot_Inhibitory  (-2.0/3.0)
// the leakage conductance "G^L"
#define Con_Leakage     0.05

// The following are global variables used!
// this indicates whether the running program is finished or not
extern int RUN_DONE;
extern bool g_no_graphic;
// these are global variables used in open GL
extern GLfloat xdepth;
extern GLfloat ydepth;
extern GLfloat zdepth;
extern int GLUT_WINDOW_ID;
extern int GLOBAL_WINDOW_WIDTH;
extern int GLOBAL_WINDOW_HEIGHT;
extern int FULL_SCREEN;
// this is the whole time length of data shown in the window screen
// normally, this is choosen as 512 ms.
extern int WINDOW_BIN_LENGTH;
// this is the time length of each data point to show, usually, this
// is chosen as 0.5 or 1 ms
extern double SLIGHT_BIN;
// this number indicates the position of cursor in the menu
extern int FIDDLE_PARAMETER;
// this indicates after how many time steps, for the data to show on the screen
extern int STEPS_PER_DRAW;
extern int STD_VIEW;
// this is the strobe structure used in open GL
extern struct strobe **GLOBAL_STRA;
// this number indicates the corresponding variable shown in the window screen
extern int DRAW_FLAG;

extern int GRAYSCALE;
// the initial time for the run
extern double GLOBAL_TI;

// The following variables are the input reading from input file!
// this is time step used in the run
extern double Tstep;
// i.e. COMP_TIME
extern double g_comp_time;
// this is the rate either in the current input case or poisson input case
extern double Rate_input;
// this is the Poisson input strength to excitatory conductance
// it is not used in current input case
extern double Strength_Exinput;
// the same as above but the Poisson input strength to inhibitory conductance
// normally, it is chosen as zero
extern double Strength_Ininput;
// the strength of excitatory neuron to excitatory neuron
extern double Strength_CorEE;
// the strength of excitatory neuron to inhibitory neuron
extern double Strength_CorIE;
// the strength of inhibitory neuron to inhibitory neuron
extern double Strength_CorII;
// the strength of inhibitory neuron to excitatory neuron
extern double Strength_CorEI;
// in current input case, the sinusoidal drive is given by
// Current_0 + Current_1*sin(2*pi*Rate_input*t+phase)
extern double Current_0;
extern double Current_1;
extern double *phase;

// this saves the spike list in each time step
extern struct raster spike_list;

// this saves all the spike event during the computation
extern struct raster RAS;

#if CORTICAL_STRENGTH_NONHOMO
// this matrix saves the cortical strength between neurons in the network
extern double** cortical_matrix;
// this is the random seed for generating random numbers
extern long CORTICAL_SEED;
#endif

// this saves neuronal network information
extern struct neuron *neu;
// this saves the generated Poisson spike train for each neuron
extern struct vector *poisson_input;
// random seed for generating Poisson process
extern long *initialseed_neuron;
// this saves the first spike time of Poisson input outside the current time step
extern double *last_input;

// these are used for generating random variables
extern long initial_pertub_Vot;
extern long initial_pertub_Ex;
extern long initial_pertub_In;
#if SMOOTH_CONDUCTANCE_USE
extern long initial_pertub_Ex_H;
extern long initial_pertub_In_H;
#endif
extern unsigned int initial_seed;
extern long* ran_iy;
extern long** ran_iv;
extern int iset;
extern double gset;

// this saves the current time of computation, the dimension is "ms"
extern double time_evolution;
// this saves the starting time of the run, normally it is chosen as zero
extern double last_time;

extern bool g_b_verbose;
extern bool g_b_quiet;
extern bool g_b_save_while_cal;
extern bool g_b_save_use_binary;
extern bool g_b_RC_filter;

extern double g_RC_filter_ci;          // filter coefficient
extern double g_RC_filter_co;          // filter coefficient

extern int g_num_neu;
extern int g_num_neu_ex;
extern int g_num_neu_in;

extern double g_strength_corEE;

#if POISSON_INPUT_USE
extern double g_poisson_rate;
extern double g_poisson_strength;
extern double *g_arr_poisson_rate;
extern double *g_arr_poisson_strength_E;
extern double *g_arr_poisson_strength_I;
extern int *tmp_tempbegin_poisson_index;
extern int *g_tempbegin_poisson_index;
extern int *g_begin_poisson_index;
#endif

extern char g_staffsave_path[];
extern char g_ras_path[];
extern char g_conductance_path[];
extern char g_spike_interval_path[];

extern neuron *tmp_tempneu;
extern neuron *tmp_tempneu2;

extern FILE *g_fout;
extern FILE *g_cond_out;

// macro defintion for neurons with excitatory type and inhibitory type
#define Type_Exneuron  1
#define Type_Inneuron  0
// macro definition for neurons' state, whether it is ouside refractory or inside it
#define STATE_ACTIVE  1
#define STATE_REFRACTORY  0

// estimated maximum times for searching a root, it is used in root_search function
#define Maxnum_search  50
// estimated maximum numerber of poisson input in each time step, if the input rate
// is really high, we should choose much larger number!
#define Maxnum_input  200
// the tolerance error for root searching
#define root_acc  (1.0e-12)
// the tolerance error to differ two values
#define LOW_BOUND (1.0e-16)

// These are used for saving raster events!
// estimated data length for raster structure (for saving all spike events)
extern int RASTER_SIZE;
// estimated data length for periodically traced vector structure
//extern int VECTOR_SIZE;
// this number means the least tolerance space in the raster structure
// if the space left in the raster structure is less than this number, the new space
// will be assigned to the raster structure
#define RASTER_ENDURE_SIZE  5
// additional assigned space for raster structure
#define RASTER_ADD_SIZE 50
// additional assigned space for vector structure
#define VECTOR_ADD_SIZE 50
// initial estimate data length for raster structure (for saving spike events in
// each time step)
#define RASTER_SPIKE_SIZE 20

// this uses different orders of Runge-Kutta algorithm
#define RK4_RUN 1
#define RK3_RUN 0
#define RK2_RUN 0

// macro definition for random number generating functions
#define RANDOM ran1
#define NTAB 32 // ran1() using

// #ifdef DEBUG
#define P_NULL_ERR(p,st) {if ((p)==NULL) fprintf(stderr, "%s\n", (st));}
// might use perror(st) instead.
// variadic macros, introduced in C99
//#define P_ERR(...) {fprintf(stderr,__VA_ARGS__);}
//#define P_ERR(args...) {fprintf(stderr,##agrs);}
#define P_ERR(ps) {fprintf(stderr,ps);}
//#define P_ERR(ps,p1) {fprintf(stderr,ps,p1);}
//#define P_ERR(ps,p1,p2) {fprintf(stderr,ps,p1,p2);}

#endif // !defined(AFX_STDAFX_H__A9DB83DB_A9FD_11D0_BFD1_444553540000__INCLUDED_)
