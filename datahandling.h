///mark: modified on 2013/1/20
#ifndef _DATAHANDLING_H_
#define _DATAHANDLING_H_

/* This holds a strobe structure, this is used in open GL to show the
evolution of data on the screen. */
struct strobe {
  int length;        // the time length of saving data, the dimension is "ms"
  int tab;           // the next position in the structure to save new data
  double last_time;  // record the last updating time
  // the updating time step (not necessary equal to computing time step)
  double update_timestep;
  double *data;      // saving data in each cycle
  int cycle_bother;  // this determines whether to save data averaged over cycles
  int cyclenum;      // the number of cycles
  double *cycledata; // data averaged across cycles
};

// The following functions are related to assigning dynamic space
// same as the malloc(), calloc(), realloc() and free() functions in standard library
void * tmalloc(size_t size);
void * tcalloc(size_t nmemb,size_t size);
void * trealloc(void *p,size_t size);
void tfree(void *p);

/* Here are the strobe and strobetrace functions */
// initial setup for strobe structure, this determines the data length,
// updating time step and also whether doing cycle averaging or not
struct strobe * strobemake(int length,double update_timestep,int cycle_bother);
// update data in the strobe structure
void strobeupdateRCFilter(struct strobe *st,double t,double DT,double val);
void strobeupdate(struct strobe *st,double t,double DT,double val);
void strobe_veri(struct strobe *st,double t,double val);
void strobeupdate_old(struct strobe *st, double t, double val);
// release data structure
void strobetfree(struct strobe *st);

// this determines the value of color corresponding to the data value on the screen,
// the "coloringtype" determines the type of picture, 0 means greyscale, 7 means colorscale.
// "val" is the data value, "valmax" is the maximum value, "valmin" is the minimum value
// "rcolor","gcolor" and "bcolor" are the corresponding color values associated with
// the given data "val"
void colorscale(int coloringtype,double val,double valmax,double valmin,
                double *rcolor,double *gcolor,double *bcolor);

// this determines the statistical characteristics of the given data arrays.
// "s" is the data type, could be int or double, "ra" is the corresponding data array
// "length" is the data length, "max" and "min" are the maximum and minimum values of
// the data array "ra", "mean" and "stdev" are the mean value and standard deviations
void stats(const char *s,void *ra,int length,double *max,double *min,double *mean,double *stdev);

// change the color mode from HSV to RGB
void hsv2rgb(double h,double s,double v,double *r,double *g,double *b);

#endif
