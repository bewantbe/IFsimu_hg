#include "stdafx.h"
#include "datahandling.h"

void * tmalloc(size_t size)
{
  void *p;
  p = malloc(size);
  if (p==NULL) fprintf(stderr, "Error: Fail to allocate memory. (tmalloc)");
  return p;
}

void * tcalloc(size_t nmemb,size_t size)
{
  void *p;
  p = calloc(nmemb,size);
  if (p==NULL) fprintf(stderr, "Error: Fail to allocate memory. (tcalloc)");
  return p;
}

void * trealloc(void *p,size_t size)
{
  void *p2;
  p2 = realloc(p,size);
  if (p2==NULL) fprintf(stderr, "Error: Fail to allocate memory. (trealloc)");
  return p2;
}

void tfree(void *p)
{
  free(p);
}

struct strobe * strobemake(int length,double update_timestep,int cycle_bother)
{
  struct strobe *st=NULL;
  st = (struct strobe *) tmalloc(sizeof(struct strobe));
  st->length = length;
  st->tab = 0;
  st->last_time = GLOBAL_TI;
  st->update_timestep = update_timestep;
  st->data = (double *) tcalloc(st->length,sizeof(double));
  st->cycle_bother = cycle_bother;
  st->cyclenum = 0;
  st->cycledata = NULL;
  if (st->cycle_bother) {
    st->cycledata = (double *) tcalloc(st->length,sizeof(double));
  }
  if (st==NULL || st->data==NULL || (st->cycle_bother && st->cycledata==NULL)) {
    fprintf(stderr, "Error: strobemake: Allocation failed.");
    st=NULL;    // force to terminate.
  }
  return st;
}

// there is actually a small bug here: the first strobe data has only 4 data
// if you set st->update_timestep = 5*DT
void strobeupdate(struct strobe *st, double t, double DT, double val)
{
  /* this averages, rather than strobes, which is better for spike statistics
     in addition, time-steps are resolved accurately */
  double oldstep=0, newstep=0,
         nexttime = st->last_time + st->update_timestep;
  int newtab=0;
  if (t+DT < nexttime) {
    st->data[st->tab] += val * DT / st->update_timestep;
  } else { /* if (t+DT >= nexttime) */
    oldstep = minimum(st->update_timestep, maximum(0, nexttime-t));
    newstep = minimum(st->update_timestep, maximum(0, t+DT-nexttime));
    st->data[st->tab] += val * oldstep / st->update_timestep;
    newtab = st->tab+1;
    if (newtab == st->length) {
      newtab = 0;
      if (st->cycle_bother) {
        int i=0;
        for (i=0; i<st->length; i++) {
          st->cycledata[i] += st->data[i];
        }
        st->cyclenum++;
      }
    }
    st->data[newtab] = val * newstep / st->update_timestep;
    st->tab = newtab;
//    st->last_time = maximum(t+DT, st->last_time + st->update_timestep);
    st->last_time = nexttime;
  }
}

void strobe_veri(struct strobe *st,double t,double val)
{
    int t_int = t;
    if(t - t_int == 0){
        st->data[st->tab] = val;
        st->tab++;
    if (st->tab == st->length)
        st->tab = 0;
    }
}


// Only suitable for truncated update_timestep and dt
void strobeupdateRCFilter(struct strobe *st, double t, double DT, double val)
{
  double vo = g_RC_filter_co * st->data[st->tab] + g_RC_filter_ci * val;
  if (t+DT/2 >= st->last_time + st->update_timestep) {
    // due to floating point error, you can't set it to st->last_time + st->update_timestep
    st->last_time = t;
    if (++st->tab == st->length)
      st->tab = 0;
  }
  st->data[st->tab] = vo;
}

// now you can use cmd parameter "--RC-filter 0 1" to get the same results
//void strobeupdate_old(struct strobe *st, double t, double val)
//{
//  /* this strobes, as opposed to averaging, which is not good for discontinuous data
//	 like spike statistics */
//  if (t >= st->last_time + st->update_timestep){
//    st->data[st->tab] = val; st->tab++;
//    if (st->tab >= st->length){
//      st->tab=0;
//      if (st->cycle_bother) {
//        for (int i=0; i<st->length; i++)
//          st->cycledata[i] += st->data[i];
//        st->cyclenum++;
//      }
//    }
//    st->last_time = maximum(t,st->last_time + st->update_timestep);
//  }
//}

void strobetfree(struct strobe *st)
{
  tfree(st->data);
  if (st->cycle_bother) {
    tfree(st->cycledata);
  }
  tfree(st);
  st=NULL;
}


void colorscale(int coloringtype,double val,double valmax,double valmin,
                double *rcolor,double *gcolor,double *bcolor)
{
  double gamma=(coloringtype==5?1.0:0.5);
  int norm_flag=1;
  double norm=1;
  double v=0;
  double h=0;
  switch (coloringtype) {
    case 0: /* greyscale coloring */
      if (valmax<=valmin) {
        valmax=valmin+1;
      }
      v = (val-valmin)/(valmax-valmin);
      v = maximum(0,minimum(1,v));
      *rcolor = v;
      *gcolor = v;
      *bcolor = v;
      break;
    case 2: /* cyclical coloring */
      if (valmax<=valmin) {
        valmax=valmin+1;
      }
      v = 3*((val-valmin)/(valmax-valmin)-0.5);
      if (v < -0.5 || v > 0.5) { *rcolor = sqrt(fabs(v)-0.5);} else { *rcolor=0;}
      if (v > -0.5) { *gcolor = sqrt(1-fabs(v-0.5));} else { *gcolor=0;}
      if (v < 0.5) { *bcolor = sqrt(1-fabs(v+0.5));} else { *bcolor=0;}
      break;
    case 3: /* IC IS EC ES ALL coloring */
      switch ((int)val) {
        case 0: /* IC */
          *rcolor=0.25;
          *gcolor=0.75;
          *bcolor=1.00;
          break;
        case 1: /* IS */
          *rcolor=0.00;
          *gcolor=0.00;
          *bcolor=1.00;
          break;
        case 2: /* EC */
          *rcolor=1.00;
          *gcolor=0.75;
          *bcolor=0.25;
          break;
        case 3: /* ES */
          *rcolor=1.00;
          *gcolor=0.00;
          *bcolor=0.00;
          break;
        case 4: /* ALL */
          *rcolor=0.75;
          *gcolor=0.00;
          *bcolor=0.75;
          break;
        default:
          *rcolor=1;
          *gcolor=1;
          *bcolor=1;
          break;
      }
      break;
    case 4: /* coloring with gamma */
      if (valmax==valmin) { valmax++;}
      v = (val-valmin)/(valmax-valmin);
      if (v<0) { v=0;}
      if (v>1) { v=1;}
      v = pow(v,gamma);
      if (!GRAYSCALE) {
        if (v<0.5) {
          v = asin(2*v);
          *rcolor = 0;
          *gcolor = sin(v);
          *bcolor = cos(v);
        } else if (v>=0.5) {
          v = asin(2*(v-0.5));
          *rcolor = sin(v);
          *gcolor = cos(v);
          *bcolor = 0;
        }
      } else if (GRAYSCALE) {
        *rcolor = v;
        *gcolor = v;
        *bcolor = v;
      }
      break;
    case 5: /* super coloring with blue-cyan-green-yellow-red and gamma value */
      if (valmax==valmin) { valmax++;}
      v = (val-valmin)/(valmax-valmin);
      if (v<0) { v=0;}
      if (v>1) { v=1;}
      v = pow(v,gamma);
      if (!GRAYSCALE) {
        if (v>=0.0 && v<0.25) {
          *rcolor = 0;
          *gcolor = (v-0.0)/0.25;
          *bcolor = 1;
        }
        if (v>=0.25 && v<0.5) {
          *rcolor = 0;
          *gcolor = 1;
          *bcolor = 1-(v-0.25)/0.25;
        }
        if (v>=0.5 && v<0.75) {
          *rcolor = (v-0.5)/0.25;
          *gcolor = 1;
          *bcolor = 0;
        }
        if (v>=0.75 && v<=1) {
          *rcolor = 1;
          *gcolor = 1-(v-0.75)/0.25;
          *bcolor = 0;
        }
        if (norm_flag) {
          norm = sqrt(pow(*rcolor,2)+pow(*gcolor,2)+pow(*bcolor,2));
          *rcolor /= norm;
          *gcolor /= norm;
          *bcolor /= norm;
        }
      } else if (GRAYSCALE) {
        *rcolor = v;
        *gcolor = v;
        *bcolor = v;
      }
      break;
    case 6: /* use hsv2rgb with h in [0,360] */
      if (valmax<=valmin) { valmax=valmin+1;}
      v = crop(val,valmin,valmax);
      h = 359.99*(v-valmin)/(valmax-valmin+1);   // h = 359.99*(v-valmin)/(valmax-valmin)
      hsv2rgb(h,1,1,rcolor,gcolor,bcolor);
      break;
    case 7: /* use hsv2rgb with h in [240,-60] */
      if (valmax<=valmin) { valmax=valmin+1;}
      v = crop(val,valmin,valmax);
      h = 300*(valmax-v)/(valmax-valmin)-60;
      hsv2rgb(h,1,1,rcolor,gcolor,bcolor);
      break;
    default:
      if (valmax==valmin) { valmax++;}
      v = (val-valmin)/(valmax-valmin);
      if (v<0) { v=0;}
      if (v>1) { v=1;}
      if (!GRAYSCALE) {
        if (v<0.5) {
          v = asin(2*v);
          *rcolor = 0;
          *gcolor = sin(v);
          *bcolor = cos(v);
        } else if (v>=0.5) {
          v = asin(2*(v-0.5));
          *rcolor = sin(v);
          *gcolor = cos(v);
          *bcolor = 0;
        }
      } else if (GRAYSCALE) {
        *rcolor = v;
        *gcolor = v;
        *bcolor = v;
      }
      break;
  }
}

void stats(const char *s,void *ra,int length,double *max,double *min,
           double *mean,double *stdev)
{
  /* finds the stats of an array */
  int i=0;
  double max2=0,min2=0,mean2=0,stdev2=0;
  if (ra!=NULL) {
    if (strcmp(s,"double")==0) { // data is "double" type
      for (i=0; i<length; i++) {
        mean2 += ((double *)ra)[i];
      }
      mean2 /= (double)length;
      max2=mean2;
      min2=mean2;
      max2=((double *)ra)[0];
      min2=((double *)ra)[0];
      for (i=0; i<length; i++) {
        if (((double *)ra)[i]>max2) { max2=((double *)ra)[i];}
        if (((double *)ra)[i]<min2) { min2=((double *)ra)[i];}
        stdev2 += (((double *)ra)[i]-mean2)*(((double *)ra)[i]-mean2);
      }
      stdev2 /= (double)length;
      stdev2 = sqrt(stdev2);
    } else if (strcmp(s,"long")==0) { // data is "long int" type
      for (i=0; i<length; i++) {
        mean2 += ((long *)ra)[i];
      }
      mean2 /= (double)length;
      max2=mean2;
      min2=mean2;
      max2=((long *)ra)[0];
      min2=((long *)ra)[0];
      for (i=0; i<length; i++) {
        if (((long *)ra)[i]>max2) { max2=((long *)ra)[i];}
        if (((long *)ra)[i]<min2) { min2=((long *)ra)[i];}
        stdev2 += (((long *)ra)[i]-mean2)*(((long *)ra)[i]-mean2);
      }
      stdev2 /= (double)length;
      stdev2 = sqrt(stdev2);
    } else if (strcmp(s,"int")==0) { // data is "int" type
      for (i=0; i<length; i++) {
        mean2 += ((int *)ra)[i];
      }
      mean2 /= (double)length;
      max2=mean2;
      min2=mean2;
      max2=((int *)ra)[0];
      min2=((int *)ra)[0];
      for (i=0; i<length; i++) {
        if (((int *)ra)[i]>max2) { max2=((int *)ra)[i];}
        if (((int *)ra)[i]<min2) { min2=((int *)ra)[i];}
        stdev2 += (((int *)ra)[i]-mean2)*(((int *)ra)[i]-mean2);
      }
      stdev2 /= (double)length;
      stdev2 = sqrt(stdev2);
    }
    if (max!=NULL) {*max=max2;}
    if (min!=NULL) {*min=min2;}
    if (mean!=NULL) {*mean=mean2;}
    if (stdev!=NULL) {*stdev=stdev2;}
  }
}

void hsv2rgb(double h,double s,double v,double *r,double *g,double *b)
{
  int i;
  double f=0,p=0,q=0,t=0;
  if (s==0) { /* greyscale */
    *r=v;
    *g=v;
    *b=v;
    return;
  }
  h /= 60;
  i=floor(h); /* six sectors */
  f = h-i;
  p = v*(1-s);
  q = v*(1-s*f);
  t = v*(1-s*(1-f));
  switch (i) {
    case 0:
      *r=v;
      *g=t;
      *b=p;
      break;
    case 1:
      *r=q;
      *g=v;
      *b=p;
      break;
    case 2:
      *r=p;
      *g=v;
      *b=t;
      break;
    case 3:
      *r=p;
      *g=q;
      *b=v;
      break;
    case 4:
      *r=t;
      *g=p;
      *b=v;
      break;
    default: /* case 5 */
      *r=v;
      *g=p;
      *b=q;
      break;
  }
}
