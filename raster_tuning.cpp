#include "stdafx.h"
#include <sys/types.h>
#include <sys/stat.h>
#include "datainput.h"
#include "poisson_input.h"
#include "myopengl.h"

#ifdef _WINDOWS_USE_
unsigned int GetSeedFromTime()
{
  SYSTEMTIME sys_time;
  GetLocalTime(&sys_time);
  //QueryPerformanceCounter
  return (unsigned int)sys_time.wMilliseconds+(((sys_time.wDay*24+sys_time.wHour)*60+sys_time.wMinute)*60+sys_time.wSecond) * 1000;
}
#else
#include <sys/time.h>
unsigned int GetSeedFromTime()
{
  struct timeval tv;
  if (gettimeofday(&tv, NULL) != 0) {
    printf("Error: gettimeofday: Fail to get system time! (for random seed)");
    perror("");
    tv.tv_sec = time(NULL);
    if (tv.tv_sec == ((time_t) -1)) {
      perror("Tried to use time(), but failed.");
    }
  }
  return (unsigned int)(tv.tv_sec*1000000 + tv.tv_usec);
}
#endif

int RUN_DONE = 0;

GLfloat xdepth = 0.0;
GLfloat ydepth = 0.0;
GLfloat zdepth = -3.0;
int GLUT_WINDOW_ID = 0;
int GLOBAL_WINDOW_WIDTH = 640;
int GLOBAL_WINDOW_HEIGHT = 480;
int FIDDLE_PARAMETER = -3;
int STEPS_PER_DRAW = 128;        ///HERE
int STD_VIEW = 1;
int DRAW_FLAG = 0;
int FULL_SCREEN = 0;
int WINDOW_BIN_LENGTH = 512;
double SLIGHT_BIN = 0.5;          ///HERE SLIGHT_BIN = update_timestep of GLOBAL_STRA
int GRAYSCALE = 0;
double GLOBAL_TI = 0;
struct strobe **GLOBAL_STRA = NULL;

double Tstep = 1.0/32.0;
double g_comp_time = COMP_TIME;
double Rate_input = 0.0;
double Strength_Exinput = 0.0;
double Strength_Ininput = 0.0;
double Strength_CorEE = 0.0;
double Strength_CorIE = 0.0;
double Strength_CorII = 0.0;
double Strength_CorEI = 0.0;
double Current_0 = 0.00;
double Current_1 = 0.00;
double *phase = NULL;

struct raster spike_list;
struct raster RAS;

#if SMOOTH_CONDUCTANCE_USE
long initial_pertub_Ex_H;
long initial_pertub_In_H;
#endif
#if CORTICAL_STRENGTH_NONHOMO
long CORTICAL_SEED;
double** cortical_matrix = NULL;
#endif

struct neuron *neu = NULL;
struct vector *poisson_input = NULL;
long *initialseed_neuron = NULL;
double *last_input = NULL;
long initial_pertub_Vot = 0;
long initial_pertub_Ex = 0;
long initial_pertub_In = 0;
unsigned int initial_seed = 0;
long* ran_iy = NULL;
long** ran_iv = NULL;
int iset = 0;
double gset = 0.0;
double time_evolution = 0.0;
double last_time = 0.0;            // record time of last run
double g_simu_dt = 0;
double g_save_intervel = 0;
double g_RC_filter_ci = 0.5;          // filter coefficient
double g_RC_filter_co = 0.5;          // filter coefficient

bool g_no_graphic         = false;
bool g_b_verbose          = false;
bool g_b_quiet            = false;
bool g_b_save_while_cal   = true;
bool g_b_save_use_binary  = false;
bool g_b_RC_filter        = false;
bool g_b_auto_seed        = true;
bool g_b_RC_filter_coef_auto = true;

int g_num_neu_ex = Number_Exneuron;
int g_num_neu_in = Number_Inneuron;
int g_num_neu = Number_Exneuron + Number_Inneuron;
unsigned int initial_seed2 = 0;

double no_use_number = - 3.1416/2.71828*1.2345+1.4142/1.736;  // mark no change
double g_strength_corEE = no_use_number;
double g_strength_corEI = no_use_number;
double g_strength_corIE = no_use_number;
double g_strength_corII = no_use_number;
#if POISSON_INPUT_USE
double g_poisson_rate = no_use_number;
double g_poisson_strength = no_use_number;
double g_poisson_strength_in = no_use_number;
double *g_arr_poisson_rate = NULL;
double *g_arr_poisson_strength_E = NULL;
double *g_arr_poisson_strength_I = NULL;
double *arr_pr_tmp = NULL;
double *arr_ps_tmp = NULL;
double *arr_psi_tmp = NULL;
#endif

const int file_path_size = 1024;
char conf_file[file_path_size] = "";
char cor_mat_file[file_path_size] = "";
char g_staffsave_path[file_path_size] = "";
char g_ras_path[file_path_size] = "";
char g_conductance_path[file_path_size] = "";
char g_spike_interval_path[file_path_size] = "";
#if POISSON_INPUT_USE
char ps_file[file_path_size] = "";
char pr_file[file_path_size] = "";
char conf_file_default[] = "test2.txt";
#else
char conf_file_default[] = "test.txt";
#endif
char cor_mat_file_default[] = "cortical_matrix.txt";
const char g_staffsave_path_default[] = "data/staffsave.txt";

FILE *g_fout = NULL;
FILE *g_cond_out = NULL;

int RASTER_SIZE = (Number_Exneuron + Number_Inneuron)*5;
//int VECTOR_SIZE = (Number_Exneuron + Number_Inneuron)*5;

void ShowCLIHelp()
{
  printf("Usage: raster_tuning [OPTION1] [OPTION2] ...\n");
  printf("Simulate the Integrate-and-Fire model of neurons.\n");
  printf("\n");
  printf("  -ng           no GUI(graphical user interface)\n");
  printf("  -t TIME       set simulation time length, in millisecond\n");
  printf("                  Default: %f\n", g_comp_time);
  printf("  -n N [N2]     set number of neurons, N for Ex., N2 for In.\n");
  printf("                  Default: %d Ex. + %d In.\n", g_num_neu_ex, g_num_neu_in);
  printf("  -inf FILE     load configuration from FILE\n");
  printf("                  Default: \"%s\"\n", conf_file_default);
  printf("  -mat FILE     load cortical strength matrix from FILE\n");
  printf("                  Default: \"%s\"\n", cor_mat_file_default);
  printf("  -o FILE       save voltage data to FILE\n");
  printf("                  Default: \"%s\"\n", g_staffsave_path_default);
  printf("  --save-conductance FILE\n");
  printf("                save conductance to FILE\n");
  printf("  --save-spike FILE\n");
  printf("                save all spike events to FILE\n");
  printf("  --save-spike-interval FILE\n");
  printf("                save average spike interval of each neuron to FILE\n");
#if POISSON_INPUT_USE
  printf("  --read-pr FILE\n");
  printf("                read relative poisson input rate of each neuron from FILE\n");
  printf("  --read-ps FILE\n");
  printf("                read relative poisson input strength of each neuron from FILE\n");
#endif
  printf("  -s, --save-while-cal\n");
  printf("                save volt data while computing %s\n",
         g_b_save_while_cal?"(default)":"(set this if you want full data!)");
  printf("  --save-last   save latest (few) data at the end of computing %s\n",
         g_b_save_while_cal?"":"(default)");
  printf("  --bin-save    use raw double binary format to save the voltage data\n");
  printf("  -scee VALUE   set cortical strength (Ex. to Ex.)\n");
  printf("  -scei VALUE   set cortical strength (Ex. to In.)\n");
  printf("  -scie VALUE   set cortical strength (In. to Ex.)\n");
  printf("  -scii VALUE   set cortical strength (In. to In.)\n");
#if POISSON_INPUT_USE
  printf("  -pr VALUE     set poisson input rate\n");
  printf("  -ps VALUE     set poisson input strength for Ex. neuron\n");
  printf("  -psi VALUE    set poisson input strength for In. neuron\n");
  printf("  --pr-mul [VALUE ...] [VALUE@POSITION ...]\n");
  printf("                set relative poisson input rate for neurons\n");
  printf("  --ps-mul [VALUE ...] [VALUE@POSITION ...]\n");
  printf("                set relative poisson input strength for Ex. neurons\n");
  printf("  --psi-mul [VALUE ...] [VALUE@POSITION ...]\n");
  printf("                set relative poisson input strength for In. neurons\n");
#endif
  printf("  -dt VALUE     set time step. Default: %f (ms)\n", Tstep);
  printf("  --save-interval VALUE  or  -stv VALUE\n");
  printf("                set output time interval of data. Default: %f (ms)\n", SLIGHT_BIN);
  printf("  --seed VALUE  set seed for random number generator\n");
  printf("                Range: integer from 0 to 2^32-1 or float between 0 and 1\n");
  printf("  --seed-auto-on, --seed-auto-off      Default: %s\n",g_b_auto_seed?"on":"off");
  printf("                turn on or off the function of auto set random seed\n");
  printf("  --RC-filter   use RC low-pass filter befor sampling (instead of averaging)\n");
  printf("  --RC-filter ci co\n");
  printf("                same as --RC-filter, but let you set filter coefficients:\n");
  printf("                y[t] = co * y[t-dt] + ci * x[t] (x is input, y is output)\n");
  printf("                final output data is y[stv], y[2*stv], ... , y[k*stv]\n");
  printf("                stv means save time interval, see --save-interval\n");
  printf("  -v, --verbose show more information while executing\n");
  printf("  -q, --quiet   show only errors\n");
  printf("  -h, --help    show this help and exit\n");
  printf("  --version     output version information and exit\n");
  printf("\n");
  printf("Using options above will overwrite the corresponding parameters in configuration file. (which default is \"%s\")\n", conf_file_default);
  printf("\n");
  printf("Report raster_tuning bugs to xyy82148@sjtu.edu.cn or zdz@cims.nyu.edu\n");
  printf("\n");
}

void ShowCLIVersion()
{
  printf("raster_tuning 2.0.10 (alpha)\n");
  printf("Copyright: ZDZ etc.\n");
  printf("XYY branch version, based on \"clean IF code\" (19 Nov 2010)\n");
  printf("(Is this a free software?!)\n");
#ifdef _MSC_VER
  printf("Compile by Microsoft Visual Studio(C++): Version %d\n", _MSC_VER);
#elif __PATHCC__
  printf("Compile by PathScale Compiler: Version %s\n", __VERSION__);
#elif __GNUG__
  printf("Compile by GCC: Version %s\n", __VERSION__);
#elif __GNUC__
  printf("Compile by GNUC compatible compiler\n");
# ifdef __VERSION__
    printf("  Version: %s\n", __VERSION__);
# endif
#else
  printf("Compile by neither MS nor GNUC compiler\n");
# ifdef __VERSION__
    printf("__VERSION__ %s\n", __VERSION__);
# endif
#endif
  printf("%s, %s\n", __DATE__, __TIME__);
  printf("\n");
}

int ReadOneLongCmdPara(int argc, char *argv[], int pp, double *vec)
{
  char tmp_str[1024] = " ";
  while (++pp < argc) {
    if (argv[pp][0]=='-') {
      pp--;
      break;
    }
    if (argv[pp][0]=='@' || tmp_str[strlen(tmp_str)-1] == '@') {
      strcat(tmp_str, argv[pp]);
    } else {
      strcat(tmp_str, " ");
      strcat(tmp_str, argv[pp]);
    }
  }
  Str2Arr(tmp_str, vec);
  return pp;
}

int CheckDirAndCreate(const char *filepath)
{
  char path[file_path_size];
  strcpy(path, filepath);
  int l = (int)strlen(path);
  if (path[0]=='.' && path[1]=='/') {              // extract the path name
    l -= 2;
    for (int i=0; i<l; i++) { path[i] = path[i+2]; }
  }
  while (--l>=0) {
    if (path[l] != '\\' && path[l] != '/')
      path[l] = '\0';
    else
      break;
  }
  if (l < 0) return 0;
  if (l == 0 && path[0] == '.') return 0;
  path[l] = '\0';
  struct stat sb;
  if (stat(path, &sb) == -1) {
    printf("Seems no the \"%s\" directory. Trying to create one for you.\n", path);
    char cl[file_path_size];
#ifdef _WINDOWS_USE_
    sprintf(cl, "md %s", path);
#else
    sprintf(cl, "mkdir -p %s", path);
#endif
    int rt = system(cl);
    printf("  (Return value of mkdir: %d)\n", rt);
    return rt;
  } else {
    if ((sb.st_mode & S_IFMT) != S_IFDIR) {
      printf("Seems \"%s\" is not a directory!\n", path);
      return 1;
    }
  }
  return 0;
}

int main(int argc, char *argv[])
{
  int verbose_gl = 0;

  // Command line parameter translation
  g_comp_time = COMP_TIME;
  int pp = 0;                                   // pointer to parameter
  while (++pp < argc) {
    if (strcmp(argv[pp], "-ng")==0 || strcmp(argv[pp], "--no-graphic")==0) {
      g_no_graphic = true;                      // no openGL
      continue;
    }
    if (strcmp(argv[pp], "-t") == 0) {          // set the final time
      if (++pp >= argc) { pp--;  break; }
      double t_end = atof(argv[pp]);
      if (t_end <= 0 || t_end > 1e99) {
        printf("The time must between 0 and 1e99\n");
        return 1;
      }
      g_comp_time = t_end;
      continue;
    }
    if (strcmp(argv[pp], "-n")==0) {           // number of neurons
      if (++pp >= argc) break;
      int n_n = atoi(argv[pp]);
      if (n_n < 0 || n_n > NUM_NEU_MAX) {
        printf("Number of neurons must between 0 and %d\n", NUM_NEU_MAX);
        return 1;
      }
      g_num_neu_ex = n_n;
      g_num_neu_in = 0;
      if (pp+1 < argc) {
        if (argv[pp+1][0]!='-') {
          n_n = atoi(argv[++pp]);
          if (n_n < 0 || n_n > NUM_NEU_MAX) {
            printf("Number of neurons must between 0 and %d\n", NUM_NEU_MAX);
            return 1;
          }
          g_num_neu_in = n_n;
        }
      }
      if (g_num_neu_ex+g_num_neu_in <= 0) {
        printf("Nothing to calculate! (Number of neurons: %d Ex. + %d In.)\n",
               g_num_neu_ex, g_num_neu_in);
        return 0;
      }
      continue;
    }
    if (strcmp(argv[pp], "-inf")==0) {         // configuration file path
      if (++pp >= argc) break;
      strcpy(conf_file, argv[pp]);
      continue;
    }
    if (strcmp(argv[pp], "-mat")==0) {         // cortex strength matrix
      if (++pp >= argc) break;
      strcpy(cor_mat_file, argv[pp]);
      continue;
    }
    if (strcmp(argv[pp], "-o")==0) {           // set output file path
      if (++pp >= argc) break;
      strcpy(g_staffsave_path, argv[pp]);
      continue;
    }
    if (strcmp(argv[pp], "--save-conductance")==0) {  // save conductance
      if (++pp >= argc) break;
      strcpy(g_conductance_path, argv[pp]);
      continue;
    }
    if (strcmp(argv[pp], "--save-spike")==0) { // set spike-time-file path
      if (++pp >= argc) break;
      strcpy(g_ras_path, argv[pp]);
      continue;
    }
    if (strcmp(argv[pp], "--save-spike-interval")==0) {  // save average spike interval
      if (++pp >= argc) break;
      strcpy(g_spike_interval_path, argv[pp]);
      continue;
    }
#if POISSON_INPUT_USE
    if (strcmp(argv[pp], "--read-pr")==0) {
      if (++pp >= argc) break;
      strcpy(pr_file, argv[pp]);
      continue;
    }
    if (strcmp(argv[pp], "--read-ps")==0) {
      if (++pp >= argc) break;
      strcpy(ps_file, argv[pp]);
      continue;
    }
#endif
    if (strcmp(argv[pp], "-s")==0 || strcmp(argv[pp], "--save-while-cal")==0) {  // save data while calculating
      g_b_save_while_cal = true;
      continue;
    }
    if (strcmp(argv[pp], "--save-last")==0) {  // save only last few data
      g_b_save_while_cal = false;
      continue;
    }
    if (strcmp(argv[pp], "--bin-save")==0) {   // save data using binary format
      g_b_save_use_binary = true;
      continue;
    }
    if (strcmp(argv[pp], "-scee")==0) {        // set Strength_CorEE
      if (++pp >= argc) break;
      g_strength_corEE = atof(argv[pp]);
      continue;
    }
    if (strcmp(argv[pp], "-scei")==0) {        // set Strength_CorIE
      if (++pp >= argc) break;
      g_strength_corEI = atof(argv[pp]);
      continue;
    }
    if (strcmp(argv[pp], "-scie")==0) {        // set Strength_CorEI
      if (++pp >= argc) break;
      g_strength_corIE = atof(argv[pp]);
      continue;
    }
    if (strcmp(argv[pp], "-scii")==0) {        // set Strength_CorII
      if (++pp >= argc) break;
      g_strength_corII = atof(argv[pp]);
      continue;
    }
#if POISSON_INPUT_USE
    if (strcmp(argv[pp], "-pr")==0) {          // set poisson rate
      if (++pp >= argc) break;
      g_poisson_rate = atof(argv[pp]);
      continue;
    }
    if (strcmp(argv[pp], "-ps")==0) {          // set poisson strength for Ex.
      if (++pp >= argc) break;
      g_poisson_strength = atof(argv[pp]);
      continue;
    }
    if (strcmp(argv[pp], "-psi")==0) {         // set poisson strength for In.
      if (++pp >= argc) break;
      g_poisson_strength_in = atof(argv[pp]);
      continue;
    }
    if (strcmp(argv[pp], "--pr-mul")==0) {          // set poisson rate multiplier
      arr_pr_tmp = (double*)malloc(g_num_neu*sizeof(double));
      P_NULL_ERR(arr_pr_tmp, "Error: main: Allocation failed, command line interpretation.");
      for (int j=0; j<g_num_neu; j++) arr_pr_tmp[j] = 1;
      pp = ReadOneLongCmdPara(argc, argv, pp, arr_pr_tmp);
      continue;
    }
    if (strcmp(argv[pp], "--ps-mul")==0) {          // set poisson strength multiplier for Ex.
      arr_ps_tmp = (double*)malloc(g_num_neu*sizeof(double));
      P_NULL_ERR(arr_ps_tmp, "Error: main: Allocation failed, command line interpretation.");
      for (int j=0; j<g_num_neu; j++) arr_ps_tmp[j] = 1;
      pp = ReadOneLongCmdPara(argc, argv, pp, arr_ps_tmp);
      continue;
    }
    if (strcmp(argv[pp], "--psi-mul")==0) {          // set poisson strength multiplier for In.
      arr_psi_tmp = (double*)malloc(g_num_neu*sizeof(double));
      P_NULL_ERR(arr_psi_tmp, "Error: main: Allocation failed, command line interpretation.");
      for (int j=0; j<g_num_neu; j++) arr_psi_tmp[j] = 1;
      pp = ReadOneLongCmdPara(argc, argv, pp, arr_psi_tmp);
      continue;
    }
#endif
    if (strcmp(argv[pp], "-dt")==0) {          // set time step
      if (++pp >= argc) break;
      g_simu_dt = atof(argv[pp]);
      if (g_simu_dt<=0) {
        printf("Illegel time step!\nProgram terminated\n");
        return 2;
      }
      continue;
    }
    if (strcmp(argv[pp], "--save-interval")==0 || strcmp(argv[pp], "-stv")==0) {// set time interval of save
      if (++pp >= argc) break;
      g_save_intervel = atof(argv[pp]);
      if (g_save_intervel<g_simu_dt) {
        printf("Illegel time interval!\nProgram terminated\n");
        return 2;
      }
      continue;
    }
    if (strcmp(argv[pp], "--seed-auto-on")==0) {
      g_b_auto_seed = true;
      continue;
    }
    if (strcmp(argv[pp], "--seed-auto-off")==0) {
      g_b_auto_seed = false;
      continue;
    }
    if (strcmp(argv[pp], "-seed")==0) { // random number seed
      if (++pp >= argc) break;
      double s = atof(argv[pp]);
      if (0<s && s<1) s *= UINT_MAX;
      initial_seed2 = (unsigned int)s;
      g_b_auto_seed = false;
      continue;
    }
    if (strcmp(argv[pp], "--RC-filter")==0) {
      g_b_RC_filter = true;
      if (pp+2<argc && argv[pp+1][0]!='-' && argv[pp+2][0]!='-') {  // no error checking here
        g_b_RC_filter_coef_auto = false;
        g_RC_filter_ci = atof(argv[pp+1]);
        g_RC_filter_co = atof(argv[pp+2]);
        pp += 2;
        if (fabs(g_RC_filter_co)>=1) {
          printf("Unstable RC filter: y[t] = %e * y[t-1] + %e * x[t]\n", g_RC_filter_co, g_RC_filter_ci);
          printf("Absolute value of coefficient of y[t-1] must smaller than ONE.\n");
          printf("Program terminated\n");
          return 0;
        }
      } else {
        g_b_RC_filter_coef_auto = true;
      }
      continue;
    }
    if (strcmp(argv[pp], "-v")==0 || strcmp(argv[pp], "--verbose")==0) {
      g_b_verbose = true;
      continue;
    }
    if (strcmp(argv[pp], "-q")==0 || strcmp(argv[pp], "--quiet")==0) {
      g_b_quiet = true;
      continue;
    }
    if (strcmp(argv[pp], "-h")==0 || strcmp(argv[pp], "--help")==0) {
      ShowCLIHelp();
      return 0;
    }
    if (strcmp(argv[pp], "--version")==0) {
      ShowCLIVersion();
      return 0;
    }
    if (strcmp(argv[pp], "")==0) {
      continue;
    }
    break;
  }
  if (pp < argc) {
    printf("Invalid option: \"%s\"\n", argv[pp]);
    printf("Try `raster_tuning --help' for more information.\n");
    return 2;
  }
  if (g_b_quiet) g_b_verbose = false;
  verbose_gl = g_b_verbose;
  if (verbose_gl) { printf(" Initializing data\n"); fflush(stdout); }

  setglobals();               // assign space for variables

  if (g_ras_path[0] == '-')
    strcpy(g_ras_path, "data/ras.txt");
  if (conf_file[0] == '\0')
    strcpy(conf_file, conf_file_default);

  readinput(conf_file);       // read basic parameters from file

  if (g_simu_dt)        Tstep = g_simu_dt;
  if (g_save_intervel)  SLIGHT_BIN = g_save_intervel;
  // truncate SLIGHT_BIN to multiple of Tstep, due to strobeupdate().
  SLIGHT_BIN = floor(SLIGHT_BIN/Tstep)*Tstep;
  if (g_strength_corEE != no_use_number) Strength_CorEE = g_strength_corEE;
  if (g_strength_corEI != no_use_number) Strength_CorIE = g_strength_corEI;  // note: the name is inversed
  if (g_strength_corIE != no_use_number) Strength_CorEI = g_strength_corIE;
  if (g_strength_corII != no_use_number) Strength_CorII = g_strength_corII;
#if POISSON_INPUT_USE
  if (g_poisson_rate   != no_use_number) Rate_input = g_poisson_rate;
  if (g_poisson_strength != no_use_number)    Strength_Exinput = g_poisson_strength;
  if (g_poisson_strength_in != no_use_number) Strength_Ininput = g_poisson_strength_in;
#endif
  if (g_b_RC_filter) {
    if (g_b_RC_filter_coef_auto) {
      double u = M_PI*Tstep/SLIGHT_BIN;
      g_RC_filter_ci = u/(1+u);
      g_RC_filter_co = 1/(1+u);
    }
  }
  if (g_b_auto_seed) {        // use second & microsecond as random seed
    initial_seed = GetSeedFromTime();
  } else {
    if (initial_seed2) initial_seed = initial_seed2;
  }
#if CORTICAL_STRENGTH_NONHOMO
  if (cor_mat_file[0]=='\0')
    strcpy(cor_mat_file, cor_mat_file_default);
  read_cortical_matrix(cor_mat_file, cortical_matrix);
#endif

#if POISSON_INPUT_USE
  if (ReadPR(pr_file, g_arr_poisson_rate) < 0) {
    printf("Fail to read pr from file: %s\n",ps_file);
    printf("Program terminated.\n");
    return 1;
  }
  if (ReadPS(ps_file, g_arr_poisson_strength_E, g_arr_poisson_strength_I) < 0) {
    printf("Fail to read ps from file: %s\n",ps_file);
    printf("Program terminated.\n");
    return 1;
  }
  if (arr_pr_tmp) {
    for (int j=0; j<g_num_neu; j++) g_arr_poisson_rate[j] = arr_pr_tmp[j];
  }
  if (arr_ps_tmp) {
    for (int j=0; j<g_num_neu; j++) g_arr_poisson_strength_E[j] = arr_ps_tmp[j];
  }
  if (arr_psi_tmp) {
    for (int j=0; j<g_num_neu; j++) g_arr_poisson_strength_I[j] = arr_psi_tmp[j];
  }
  free(arr_pr_tmp);   arr_pr_tmp = NULL;
  free(arr_ps_tmp);   arr_ps_tmp = NULL;
  free(arr_psi_tmp);  arr_psi_tmp = NULL;
  for (int j=0; j<g_num_neu; j++) {
    g_arr_poisson_rate[j] *= Rate_input;
#if SMOOTH_CONDUCTANCE_USE
    g_arr_poisson_strength_E[j] *= Strength_Exinput;
    g_arr_poisson_strength_I[j] *= Strength_Ininput;
#else
    g_arr_poisson_strength_E[j] *= Strength_Exinput*2.0/(Time_ExCon);
    g_arr_poisson_strength_I[j] *= Strength_Ininput*5.0/(Time_InCon);
#endif
  }
  Rate_input = g_arr_poisson_rate[0];
#if SMOOTH_CONDUCTANCE_USE
  Strength_Exinput = g_arr_poisson_strength_E[0];
  Strength_Ininput = g_arr_poisson_strength_I[0];
#else
  Strength_Exinput = g_arr_poisson_strength_E[0]/(2.0/(Time_ExCon));
  Strength_Ininput = g_arr_poisson_strength_I[0]/(2.0/(Time_ExCon));
#endif
#endif

  input_initialization();

  if (g_staffsave_path[0] == '\0')
    strcpy(g_staffsave_path, g_staffsave_path_default);

  CheckDirAndCreate(g_staffsave_path);      // check the existence of the dir

  if (g_b_save_while_cal) {
    if (g_b_save_use_binary)
      g_fout = fopen(g_staffsave_path, "wb");
    else
      g_fout = fopen(g_staffsave_path, "w");
    if (g_fout == NULL) {
      printf("Fail to open \"%s\" for output!\n", g_staffsave_path);
      return 1;
    }
  }
  if (g_conductance_path[0]) {
    g_cond_out = fopen(g_conductance_path, "w");
    if (g_cond_out == NULL) {
      printf("Fail to open \"%s\" for conductance output\n", g_conductance_path);
      return 1;
    }
  }

  if (!g_b_quiet) {
    printf("Number of neurons: %d Ex. + %d In.\n", g_num_neu_ex, g_num_neu_in);
    printf("Input type: %s\n", POISSON_INPUT_USE?"poisson":"current");
    if (g_b_verbose && POISSON_INPUT_USE)
      printf("  Random seed: %u\n", initial_seed);
    printf("Cortical strength network type: %shomogeneous\n",
           CORTICAL_STRENGTH_NONHOMO?"non-":"");
    printf("Simulation time: %g ms,  dt = %f ms\n", g_comp_time, Tstep);
    printf("Data record interval: %f ms\n", SLIGHT_BIN);
    if (SMOOTH_CONDUCTANCE_USE) printf("Use smoothed conductance.\n");
  }
  if (g_b_verbose) {
    printf("Cortical strength of In. to In. neurons = %g\n", Strength_CorII);
    printf("Cortical strength of In. to Ex. neurons = %g\n", Strength_CorEI);
    printf("Cortical strength of Ex. to In. neurons = %g\n", Strength_CorIE);
    printf("Cortical strength of Ex. to Ex. neurons = %g\n", Strength_CorEE);
#if POISSON_INPUT_USE
    printf("Poisson input strength to In. neurons = %g\n", Strength_Ininput);
    printf("Poisson input strength to Ex. neurons = %g\n", Strength_Exinput);
    printf("Poisson input rate = %g\n", Rate_input);
    {
      int nn = g_num_neu;
      if (nn>20) {
        nn = 20;
        printf("  Only show first 20 neurons...\n");
      }
      printf("pr: ");
      for (int j=0; j<nn; j++) printf("%g, ", g_arr_poisson_rate[j]);
      printf("\nps: ");
      for (int j=0; j<nn; j++) printf("%g, ", g_arr_poisson_strength_E[j]);
      if (g_num_neu_in) {
        printf("\npsi: ");
        for (int j=0; j<nn; j++) printf("%g, ", g_arr_poisson_strength_I[j]);
      }
      printf("\n");
    }
#else
    printf("Current average:   %g\n", Current_0);
    printf("Current amplitude: %g\n", Current_1);
#endif
    printf("\n");
    printf("Read configuration from \"%s\"\n", conf_file);
    printf("Read cortical matrix from \"%s\"\n", cor_mat_file);
    printf("Save voltage data to \"%s\" %s %s\n", g_staffsave_path,
           g_b_save_while_cal?"":"(Save lastest data)",
           g_b_save_use_binary?"(binary mode)":"");
    if (g_ras_path[0])
      printf("Save firing data to \"%s\"\n", g_ras_path);
    if (g_conductance_path[0])
      printf("Save conductance data to \"%s\"\n", g_conductance_path);
    if (g_spike_interval_path[0])
      printf("Save average spike interval to \"%s\"\n", g_spike_interval_path);
    printf("\n");
    printf("Start computing...\n");
  }

  if (g_no_graphic) {       // non-GUI mode (non-interactive)
    double t1 = clock();
    int percent = 0;
    int len = 3;
    if (g_b_verbose) {
      len = printf("0%%");  fflush(stdout);
    }
    while (!RUN_DONE) {
      compute_perstep();
      if (g_b_verbose && time_evolution >= g_comp_time*0.01*(percent+1)) {
        percent = (int)(time_evolution/g_comp_time * 100);
        for (int j=0; j<len; j++) printf("\b"); // for matlab compatible (or just \r is OK)
        len = printf("%d%% time elapsed: %3.3f s     ",
                     percent, (clock()-t1)/CLOCKS_PER_SEC);
        fflush(stdout);
      }
    }
    if (g_b_verbose) {
      for (int j=0; j<len; j++) printf("\b");
      printf("\r100%%                                              ");
    }
    if (!g_b_quiet)
      printf("\nTime cost(in computation): %3.3f s\n", (clock()-t1)/CLOCKS_PER_SEC);
    LastRun();
    data_dump();
    if (!g_b_quiet)
      printf("done.\n");
    return 0;
  }

  // select display mode, 2xbuffer,rgba,alpha,depth_buffer
  if (verbose_gl) { printf(" calling glutInitDisplayMode\n"); fflush(stdout); }
//  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  // starts at upper left corner of screen
  if (verbose_gl) { printf(" ... WindowPosition\n"); fflush(stdout); }
  glutInitWindowPosition(0,0);
  // small window
  if (verbose_gl) { printf(" calling glutInitWindowSize\n"); fflush(stdout); }
  glutInitWindowSize(GLOBAL_WINDOW_WIDTH, GLOBAL_WINDOW_HEIGHT);
  // initialize glut
  if (verbose_gl) { printf(" calling glutInit\n"); fflush(stdout); }
  glutInit(&argc, argv);
  // open window
  if (verbose_gl) { printf(" ... CreateWindow\n"); fflush(stdout); }
  GLUT_WINDOW_ID = glutCreateWindow("Raster Tuning");
  // resizing function
  if (verbose_gl) { printf(" ReshapeFunc\n"); fflush(stdout); }
  glutReshapeFunc(ReSizeGLScene);
  // this does all the drawing for us
  if (verbose_gl) { printf(" DrawGlScene\n"); fflush(stdout); }
  glutDisplayFunc(DrawGLScene);
  // continuously redraw
  if (verbose_gl) { printf(" IdleFunc\n"); fflush(stdout); }
  glutIdleFunc(DrawGLScene);
  // input functions
  if (verbose_gl) { printf(" KeyboardFunc\n"); fflush(stdout); }
  glutKeyboardFunc(keyPressed);
  if (verbose_gl) { printf(" SpecialFunc\n"); fflush(stdout); }
  glutSpecialFunc(specialKeyPressed);
  // init window
  if (verbose_gl) { printf(" Finally Calling InitGL\n"); fflush(stdout); }
  InitGL(GLOBAL_WINDOW_WIDTH, GLOBAL_WINDOW_HEIGHT);

  // start processing
  if (verbose_gl) { printf(" Passing to main loop\n"); fflush(stdout); }
  glutMainLoop();

  return 0;
}
