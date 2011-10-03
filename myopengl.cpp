#include "stdafx.h"
#include "myopengl.h"
#include "datainput.h"
#include "datahandling.h"
#include "poisson_input.h"
#include <iomanip>
#include <iostream>
#include <sstream>

void ftexto(float x, float y, float z, const char *text)
{
  /* thanks to alex */
  const char *p = text;
  glRasterPos3f(x,y,z);
  for (; *p; p++) {
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *p);
  }
}

GLvoid InitGL(GLsizei Width, GLsizei Height)
{
  /* A general OpenGL initialization function.  Sets all of the initial parameters. */
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f); // This Will Clear The Background Color To Black
  glClearDepth(1.0);        // Enables Clearing Of The Depth Buffer
  glDepthFunc(GL_LESS);     // The Type Of Depth Test To Do
  glEnable(GL_DEPTH_TEST);      // Enables Depth Testing
  glShadeModel(GL_SMOOTH);      // Enables Smooth Color Shading
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();       // Reset The Projection Matrix
  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f); // Calculate The Aspect Ratio Of The Window
  glMatrixMode(GL_MODELVIEW);
}

GLvoid ReSizeGLScene(GLsizei Width, GLsizei Height)
{
  /* The function called when our window is resized */
  if (Height==0) { Height=1;} /* don't divide by zero */
  glViewport(0, 0, Width, Height); /* reset viewport & perspective */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);
  glMatrixMode(GL_MODELVIEW);
}

GLvoid glbox(double x,double y,double r,double rcolor,double gcolor,double bcolor)
{
  glColor3f(rcolor,gcolor,bcolor);
  glBegin(GL_QUADS);
  glVertex3f(x+r/2,y+r/2,0);
  glVertex3f(x+r/2,y-r/2,0);
  glVertex3f(x-r/2,y-r/2,0);
  glVertex3f(x-r/2,y+r/2,0);
  glEnd();
}

GLvoid Drawmenu(double side,double xoffset,double yoffset)
{
  char text[120];
  int menupos=0;
  double tcolor=0;
  const double b = 0.5;    // color level range (0, 1.0)

  menupos=menu_top-9;
  glColor3f(1,1,1);
  sprintf(text,"%d Ex. neurons, %d In. neurons",g_num_neu_ex,
          g_num_neu_in);
  ftexto(0*side+xoffset, menupos*side+yoffset, 0, text);

  menupos=menu_top-8;
  tcolor = b+(1-b)*(FIDDLE_PARAMETER==menupos);
  glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"T = %f, DT = 1/%d",time_evolution,(int)(1.0/Tstep));
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);

#if POISSON_INPUT_USE
  // use Poisson input **********************************
  menupos=menu_top-7;
  tcolor = b+(1-b)*(FIDDLE_PARAMETER==menupos);
  glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"Poisson input rate = %f",Rate_input);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);

  menupos=menu_top-6;
  tcolor = b+(1-b)*(FIDDLE_PARAMETER==menupos);
  glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"Poisson input strength to Ex. neurons = %f",Strength_Exinput);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);

  menupos=menu_top-5;
  tcolor = b+(1-b)*(FIDDLE_PARAMETER==menupos);
  glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"Poisson input strength to In. neurons = %f",Strength_Ininput);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
#else
  // use sinusoidal current input ***************************
  menupos=menu_top-7;
  tcolor = b+(1-b)*(FIDDLE_PARAMETER==menupos);
  glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"External cosinusoidal current rate = %f",Rate_input);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);

  menupos=menu_top-6;
  tcolor = b+(1-b)*(FIDDLE_PARAMETER==menupos);
  glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"External average current strength = %f",Current_0);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);

  menupos=menu_top-5;
  tcolor = b+(1-b)*(FIDDLE_PARAMETER==menupos);
  glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"External cosinusoidal current amplitude = %f",Current_1);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
#endif

  menupos=menu_top-4;
  tcolor = b+(1-b)*(FIDDLE_PARAMETER==menupos);
  glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"Cortical strength of Ex. to Ex. neurons = %f",Strength_CorEE);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);

  menupos=menu_top-3;
  tcolor = b+(1-b)*(FIDDLE_PARAMETER==menupos);
  glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"Cortical strength of Ex. to In. neurons = %f",Strength_CorIE);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);

  menupos=menu_top-2;
  tcolor = b+(1-b)*(FIDDLE_PARAMETER==menupos);
  glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"Cortical strength of In. to In. neurons = %f",Strength_CorII);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);

  menupos=menu_top-1;
  tcolor = b+(1-b)*(FIDDLE_PARAMETER==menupos);
  glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"Cortical strength of In. to Ex. neurons = %f",Strength_CorEI);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);

  menupos=menu_top;
  tcolor = b+(1-b)*(FIDDLE_PARAMETER==menupos);
  glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"STEPS_PER_DRAW = %d",STEPS_PER_DRAW);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);

  menupos=menu_top+1;
  tcolor = b+(1-b)*(FIDDLE_PARAMETER==menupos);
  glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"DRAW_FLAG = %d",DRAW_FLAG);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);

}

void keyPressed(unsigned char key, int x, int y)
{
  /* The function called whenever a normal key is pressed. */
  //Sleep(100); /* still not sure why this is called */
  switch (key) {
    case ESCAPE_KEY: /* stop everything */
      LastRun();       // pour data to file
      data_dump();
      glutDestroyWindow(GLUT_WINDOW_ID); /* shut down our window */
      if (!g_b_quiet)
        printf(RUN_DONE?"\nCalculation finished.\n":"\nComputation interrupted\n");
      fflush(stdout);
      exit(0);       /* exit the program...normal termination. */
      break;
    case SPACE_KEY: /* stop/restart temporarily */
      STEPS_PER_DRAW = !STEPS_PER_DRAW;
      if (!STEPS_PER_DRAW) {
        printf("PROCESS HALTED\n");
      }
      break;
    case ENTER_KEY:
      if (RUN_DONE) {              // continue calculating
        last_time += g_comp_time;
        RUN_DONE = 0;
      }
      break;
    case Q_KEY:
      if (FULL_SCREEN) {
        GLOBAL_WINDOW_WIDTH = 640;
        GLOBAL_WINDOW_HEIGHT = 480;
        FULL_SCREEN = 0;
      } else {
        GLOBAL_WINDOW_WIDTH = 1024;
        GLOBAL_WINDOW_HEIGHT = 768;
        FULL_SCREEN = 1;
      }
      break;
    case 'w':
      FIDDLE_PARAMETER++;
      if (FIDDLE_PARAMETER > menu_top+1) {
        // go back to the bottom line
        FIDDLE_PARAMETER = menu_top-8;
      }
      break;
    case 's':
      FIDDLE_PARAMETER--;
      if (FIDDLE_PARAMETER < menu_top-8) {
        // go back to the top line
        FIDDLE_PARAMETER = menu_top+1;
      }
      break;
    case 'd':
      switch (FIDDLE_PARAMETER) {
        case menu_top-8:
          Tstep *=2;
          break;
#if POISSON_INPUT_USE
        case menu_top-7:
          Rate_input *= 1.02;
          break;
        case menu_top-6:
          Strength_Exinput *= 1.02;
          break;
        case menu_top-5:
          Strength_Ininput *= 1.02;
          break;
#else
        case menu_top-7:
          Rate_input *= 1.02;
          break;
        case menu_top-6:
          Current_0 *= 1.02;
          break;
        case menu_top-5:
          Current_1 *= 1.02;
          break;
#endif
        case menu_top-4:
          Strength_CorEE *= 1.02;
          break;
        case menu_top-3:
          Strength_CorIE *= 1.02;
          break;
        case menu_top-2:
          Strength_CorII *= 1.02;
          break;
        case menu_top-1:
          Strength_CorEI *= 1.02;
          break;
        case menu_top:
          STEPS_PER_DRAW *= 2;
          if (STEPS_PER_DRAW == 0)
            STEPS_PER_DRAW = 1;
          if (STEPS_PER_DRAW > 65536)
            STEPS_PER_DRAW = 65536;
          break;
        case menu_top+1:
          DRAW_FLAG+=1;
          break;
        default:
          break;
      }
      break;
    case 'D':
      switch (FIDDLE_PARAMETER) {
        case menu_top-8:
          Tstep *=2;
          break;
#if POISSON_INPUT_USE
        case menu_top-7:
          Rate_input *= 1.02;
          break;
        case menu_top-6:
          Strength_Exinput *= 1.02;
          break;
        case menu_top-5:
          Strength_Ininput *= 1.02;
          break;
#else
        case menu_top-7:
          Rate_input *= 1.02;
          break;
        case menu_top-6:
          Current_0 *= 1.02;
          break;
        case menu_top-5:
          Current_1 *= 1.02;
          break;
#endif
        case menu_top-4:
          Strength_CorEE *= 1.02;
          break;
        case menu_top-3:
          Strength_CorIE *= 1.02;
          break;
        case menu_top-2:
          Strength_CorII *= 1.02;
          break;
        case menu_top-1:
          Strength_CorEI *= 1.02;
          break;
        case menu_top:
          STEPS_PER_DRAW *= 2;
          if (STEPS_PER_DRAW == 0)
            STEPS_PER_DRAW = 1;
          if (STEPS_PER_DRAW > 65536)
            STEPS_PER_DRAW = 65536;
          break;
        default:
          break;
      }
      break;
    case 'a':
      switch (FIDDLE_PARAMETER) {
        case menu_top-8:
          Tstep /=2;
          break;
#if POISSON_INPUT_USE
        case menu_top-7:
          Rate_input /= 1.02;
          break;
        case menu_top-6:
          Strength_Exinput /= 1.05;
          break;
        case menu_top-5:
          Strength_Ininput /= 1.05;
          break;
#else
        case menu_top-7:
          Rate_input /= 1.02;
          break;
        case menu_top-6:
          Current_0 /= 1.02;
          break;
        case menu_top-5:
          Current_1 /= 1.02;
          break;
#endif
        case menu_top-4:
          Strength_CorEE /= 1.02;
          break;
        case menu_top-3:
          Strength_CorIE /= 1.02;
          break;
        case menu_top-2:
          Strength_CorII /= 1.02;
          break;
        case menu_top-1:
          Strength_CorEI /= 1.02;
          break;
        case menu_top:
          STEPS_PER_DRAW /= 2;
          if (STEPS_PER_DRAW == 0) {
            STEPS_PER_DRAW = 1;
          }
          break;
        case menu_top+1:
          DRAW_FLAG-=1;
          break;
        default:
          break;
      }
      break;
    case 'A':
      switch (FIDDLE_PARAMETER) {
        case menu_top-8:
          Tstep /=2;
          break;
#if POISSON_INPUT_USE
        case menu_top-7:
          Rate_input /= 1.02;
          break;
        case menu_top-6:
          Strength_Exinput /= 1.02;
          break;
        case menu_top-5:
          Strength_Ininput /= 1.02;
          break;
#else
        case menu_top-7:
          Rate_input /= 1.02;
          break;
        case menu_top-6:
          Current_0 /= 1.02;
          break;
        case menu_top-5:
          Current_1 /= 1.02;
          break;
#endif
        case menu_top-4:
          Strength_CorEE /= 1.02;
          break;
        case menu_top-3:
          Strength_CorIE /= 1.02;
          break;
        case menu_top-2:
          Strength_CorII /= 1.02;
          break;
        case menu_top-1:
          Strength_CorEI /= 1.02;
          break;
        case menu_top:
          STEPS_PER_DRAW /= 2;
          if (STEPS_PER_DRAW == 0) {
            STEPS_PER_DRAW = 1;
          }
          break;
        default:
          break;
      }
      break;
  }
}

void specialKeyPressed(int key, int x, int y)
{
  /* The function called whenever a special key is pressed. */
  int mod=0;
  Sleep(100);
  switch (key) {
    case GLUT_KEY_PAGE_UP:
      zdepth /= 1.05;
      break;
    case GLUT_KEY_PAGE_DOWN:
      zdepth *= 1.05;
      break;
    case GLUT_KEY_UP:
      ydepth+= -0.05*zdepth;
      break;
    case GLUT_KEY_DOWN:
      ydepth-= -0.05*zdepth;
      break;
    case GLUT_KEY_LEFT:
      xdepth-= -0.05*zdepth;
      break;
    case GLUT_KEY_RIGHT:
      xdepth+= -0.05*zdepth;
      break;
    case GLUT_KEY_END:
      mod = glutGetModifiers();
      if (mod==GLUT_ACTIVE_SHIFT) {
        GLOBAL_WINDOW_WIDTH = 640;
        GLOBAL_WINDOW_HEIGHT = 480;
      } else {
        if (!STEPS_PER_DRAW) {
          compute_perstep();
        }
      }
      break;
    default:
      break;
  }
}

GLvoid Drawra(void *ra,const char *type,int rows,int cols,int coloringtype,int stduse,int gsuse,
              const char *text,double xside,double yside,double xoffset,double yoffset)
{
  /* draws a standard *ra */
  char text2[64];
  int i=0,j=0;
  double varmax=0,varmin=0,varmean=0,varstd=0,vMAX=0,vMIN=0;
  double *dra=NULL,*dra2=NULL;
  int *ira=NULL,*ira2=NULL;
  double side=maximum(xside,yside), bigside=maximum(rows,cols);
  double xord=0,yord=0;
  double rcolor=0,gcolor=0,bcolor=0;
  int clear_flag=0;

  if (strcmp(type,"double")==0) {
    dra = (double *)ra;
    if (gsuse && 0 /* put that "0" there so that only the "else" will be done */) {
      /* do nothing, because spacesmear is not implemented. */
    } else {
      dra2 = dra;
      clear_flag=0;
    }
    stats("double",dra2,rows*cols,&varmax,&varmin,&varmean,&varstd);
    if (stduse) {
      vMAX=varmean+STD_VIEW*varstd;
      vMIN=varmean-STD_VIEW*varstd;
    } else {
      vMAX=varmax;
      vMIN=varmin;
    }
    for (i=0; i<rows; i++) {
      for (j=0; j<cols; j++) {
        xord = (j+0.5)/bigside - (0.5*cols)/bigside;
        yord = -(i+0.5)/bigside + (0.5*rows)/bigside;;
        colorscale(coloringtype,dra2[i+j*rows],vMAX,vMIN,&rcolor,&gcolor,&bcolor);
        glbox(xord*xside+xoffset,
              yord*yside+yoffset,
              side/bigside,
              rcolor, gcolor, bcolor);
      }
    }
    if (clear_flag) {
      tfree(dra2);
      dra2=NULL;
    }
  } else if (strcmp(type,"int")==0) {
    ira = (int *)ra;
    ira2 = ira;
    clear_flag=0;
    stats("int",ira2,rows*cols,&varmax,&varmin,&varmean,&varstd);
    if (stduse) {
      vMAX=varmean+STD_VIEW*varstd;
      vMIN=varmean-STD_VIEW*varstd;
    } else {
      vMAX=varmax;
      vMIN=varmin;
    }
    for (i=0; i<rows; i++) {
      for (j=0; j<cols; j++) {
        xord = j/bigside;
        yord = -i/bigside;
        colorscale(coloringtype,ira2[i+j*rows],vMAX,vMIN,&rcolor,&gcolor,&bcolor);
        glbox(xord*xside+xoffset,yord*yside+yoffset,side/bigside,rcolor,gcolor,bcolor);
      }
    }
    if (clear_flag) {
      tfree(ira2);
      ira2=NULL;
    }
  }
  glColor3f(1,1,1);
  if (text!=NULL) {
    if (maximum(fabs(varmin),fabs(varmax))<0.005) {
      sprintf(text2,"%s-[%0.3e..%0.3e]",text,varmin,varmax);
    } else {
      sprintf(text2,"%s-[%0.2f..%0.2f]",text,varmin,varmax);
    }
  } else {
    if (maximum(fabs(varmin),fabs(varmax))<0.005) {
      sprintf(text2,"[%0.3e..%0.3e]",varmin,varmax);
    } else {
      sprintf(text2,"[%0.2f..%0.2f]",varmin,varmax);
    }
  }
  ftexto(xoffset, yoffset+0.13, 0, text2);
}


GLvoid Plotra(void *ra,char *type,int length,double max,double min,
              double xside,double yside,double xoffset,double yoffset,
              double rcolor,double gcolor,double bcolor)
{
  int nr=0;
  double max2=0,min2=0,mean=0,stdev=0;
  double xord=0,yord=0;
  double *dra=NULL;
  int *ira=NULL;
  if (strcmp(type,"double")==0) {
    dra = (double *)ra;
    if (max>min) {
      max2=max;
      min2=min;
    } else if (max<min) {
      stats(type,dra,length,&max2,&min2,&mean,&stdev);
    } else if (max==min) {
      stats(type,dra,length,&max2,&min2,&mean,&stdev);
      min2=0;
    }
    glColor3f(rcolor,gcolor,bcolor);
    glBegin(GL_LINE_STRIP); // lines which are not closed.
    for (nr=0; nr<length; nr++) {
      yord = (dra[nr]-min2)/(max2-min2);
      xord = ((double)nr+0.5)/(double)length;
      glVertex3f(xord*xside+xoffset,yord*yside+yoffset,0);
    }
    glEnd();
  } else if (strcmp(type,"int")==0) {
    ira = (int *)ra;
    if (max>min) {
      max2=max;
      min2=min;
    } else if (max<min) {
      stats(type,ira,length,&max2,&min2,&mean,&stdev);
    } else if (max==min) {
      stats(type,ira,length,&max2,&min2,&mean,&stdev);
      min2=0;
    }
    glColor3f(rcolor,gcolor,bcolor);
    glBegin(GL_LINE_STRIP);
    for (nr=0; nr<length; nr++) {
      yord = (ira[nr]-min2)/(max2-min2);
      xord = ((double)nr+0.5)/(double)length;
      glVertex3f(xord*xside+xoffset,yord*yside+yoffset,0);
    }
    glEnd();
  }
}

GLvoid Drawstra(struct strobe **stra,int ralength,double side,double xoffset,double yoffset,
                double max,double min,double rcolor,double gcolor,double bcolor,
                int drawvsplot)
{
  int sttab=stra[0]->tab, stlength=stra[0]->length,
      stcycle_bother=stra[0]->cycle_bother, stcyclenum=stra[0]->cyclenum;
  int na=0, nt=0, nt2=0;
  double xord=0, yord=0;
  double max2=max, min2=min, stdscale=0;
  double *meanra=NULL,*stdevra=NULL,mean=0,stdev=0,*ra=NULL;
  int coloringtype = rcolor >= 0 ? 0 : (int)gcolor;
  double rcolor2=0,gcolor2=0,bcolor2=0;
  int automatic_color=0;

  if (drawvsplot==0) {
    if (max<=min) {
      stdscale=17;
      meanra = (double*)tcalloc(ralength,sizeof(double));
      stdevra= (double*)tcalloc(ralength,sizeof(double));
      P_NULL_ERR(meanra, "Error: Drawstra: Allocation failed.");
      P_NULL_ERR(stdevra, "Error: Drawstra: Allocation failed.");
      //for (na=0; na<ralength; na++) {
      //  stats("double", stra[na]->data, stlength,
      //        NULL, NULL, meanra+na, stdevra+na);
      //}
      //stats("double", meanra, ralength,
      //      NULL,NULL, &mean, NULL);      // mean = average mean of each neuron
      //stats("double", stdevra, ralength,
      //      NULL, NULL, &stdev, NULL);    // stdev = average stdev of each neuron
      for (na=0; na<ralength; na++) {
        stats("double", stra[na]->data, stlength,
              meanra+na, stdevra+na, NULL, NULL);
      }
      stats("double", meanra, ralength,
            &max2, NULL, NULL, NULL);      // mean = average mean of each neuron
      stats("double", stdevra, ralength,
            NULL, &min2, NULL, NULL);     // stdev = average stdev of each neuron
      tfree(meanra);
      tfree(stdevra);
      //max2 = mean+stdscale*STD_VIEW*stdev;
      //min2 = mean-stdscale*STD_VIEW*stdev;
    }
    for (na=0; na<ralength; na++) {
      switch (coloringtype) {
        case -1: /* automatic coloring */
          automatic_color=1;
          break;
        case 0: /* specified coloring */
          rcolor2=rcolor;
          gcolor2=gcolor;
          bcolor2=bcolor;
          break;
        case 1:
          colorscale(7, na, ralength-1, 0,
                     &rcolor2,&gcolor2,&bcolor2);
          break;
        case 2:
          colorscale(6, na, ralength-1, 0,
                     &rcolor2,&gcolor2,&bcolor2);
          break;
        case 3:
          colorscale(3, na, 3, 1,
                     &rcolor2,&gcolor2,&bcolor2);
          break;
        default:
          rcolor2=rcolor;
          gcolor2=gcolor;
          bcolor2=bcolor;
          break;
      }
      for (nt=sttab; nt<sttab+stlength; nt++) {
        nt2 = periodize(nt,0,stlength);  // shift value in [o..stlength) because strobe is a "pseudo-loop"
        xord = (nt-sttab+0.5)/stlength-0.5;
        if (automatic_color) {
          colorscale(7, stra[na]->data[nt2], max2, min2,
                     &rcolor2,&gcolor2,&bcolor2);
          yord = (-na+0.5)/ralength;
          glbox(xord*side+xoffset, yord*side+yoffset,
                side/(double)minimum(ralength,stlength),
                rcolor2, gcolor2, bcolor2);
        } else {
          yord = (stra[na]->data[nt2] - min2)/(max2 - min2)-0.5;
          glbox(xord*side+xoffset, yord*side/5+yoffset,
                side/(double)minimum(100,stlength),
                rcolor2, gcolor2, bcolor2);
        }
      }
    }
    if (stcycle_bother && stcyclenum>0) {
      if (max<=min) {
        stdscale=17;
        meanra = (double*)tcalloc(ralength,sizeof(double));
        stdevra= (double*)tcalloc(ralength,sizeof(double));
        P_NULL_ERR(meanra, "Error: Drawstra: Allocation failed.");
        P_NULL_ERR(stdevra, "Error: Drawstra: Allocation failed.");
        for (na=0; na<ralength; na++) {
          stats("double",stra[na]->cycledata, stlength,
                NULL, NULL, meanra+na, stdevra+na);
        }
        stats("double",meanra,ralength,NULL,NULL,&mean,NULL);
        stats("double",stdevra,ralength,NULL,NULL,&stdev,NULL);
        tfree(meanra);
        tfree(stdevra);
        max2 = mean+stdscale*STD_VIEW*stdev;
        min2 = mean-stdscale*STD_VIEW*stdev;
      } else {
        max2 = max*stcyclenum;
        min2 = min*stcyclenum;
      }
      for (na=0; na<ralength; na++) {
        switch (coloringtype) {
          case -1:
            automatic_color=1;
            break;
          case 0:
            rcolor2=rcolor;
            gcolor2=gcolor;
            bcolor2=bcolor;
            break;
          case 1:
            colorscale(7,na,ralength-1,0,&rcolor2,&gcolor2,&bcolor2);
            break;
          case 2:
            colorscale(6,na,ralength-1,0,&rcolor2,&gcolor2,&bcolor2);
            break;
          case 3:
            colorscale(3,na,3,1,&rcolor2,&gcolor2,&bcolor2);
            break;
          default:
            rcolor2=rcolor;
            gcolor2=gcolor;
            bcolor2=bcolor;
            break;
        }
        for (nt=0; nt<stlength; nt++) {
          xord = (double)(nt+0.5)/(double)stlength;
          if (automatic_color) {
            colorscale(7,stra[na]->cycledata[nt]/stcyclenum,max2,min2,&rcolor2,&gcolor2,&bcolor2);
            yord = (double)(-na+0.5)/(double)ralength;
            glbox(xord*side+xoffset,
                  yord*side+yoffset,side/(double)minimum(ralength,stlength),
                  rcolor2, gcolor2, bcolor2);
          } else {
            yord = (stra[na]->cycledata[nt]-min2)/(max2-min2);
            glbox(xord*side+xoffset, yord*side+yoffset+side,
                  side/(double)minimum(100,stlength),
                  rcolor2, gcolor2, bcolor2);
          }
        }
      }
    }
  } else if (drawvsplot==1) {
    ra = (double *) tcalloc(ralength*stlength,sizeof(double));
    P_NULL_ERR(ra, "Error: Drawstra: Allocation failed.");
    for (na=0; na<ralength; na++) {
      for (nt=sttab; nt<sttab+stlength; nt++) {
        nt2 = periodize(nt,0,stlength);
        ra[na + (nt-sttab)*ralength] = stra[na]->data[nt2];
      }
    }
    Drawra(ra, "double", ralength, stlength,
           7,      // coloring type
           1,      // whether use standard var to determine color scale, 1 for use, 0 for not use
           0,      // gsuse
           NULL,   // text (char*)
           side,   // x side length
           side*100/g_num_neu,  // y side length
           xoffset, yoffset
          );
    if (stcycle_bother && stcyclenum>0) {
      for (na=0; na<ralength; na++) {
        for (nt=0; nt<stlength; nt++) {
          ra[na + nt*ralength] = stra[na]->cycledata[nt]/stcyclenum;
        }
      }
      Drawra(ra,"double",ralength,stlength,7,1,0,NULL,side,side*100/g_num_neu,xoffset,yoffset+side);
    }
    tfree(ra);
  }
}

/// Show information about firing time on OpenGL window
void DrawSpikeInfo()
{
  using std::setw;
  using std::setprecision;
  using std::ostringstream;
  double tlen   = 200;
  double tlast  = time_evolution;
  double tfirst = tlast - tlen;
  const double *rafrtm = RAS.array_firingtime;
  const int    *rafrid = RAS.array_index;
  int k = RAS.ras_index;
  double *frtmlast  = (double*)calloc(g_num_neu, sizeof(double));
  double *frtmfirst = (double*)calloc(g_num_neu, sizeof(double));
  size_t *frcnt     = (size_t*)calloc(g_num_neu, sizeof(size_t));
  if (frtmlast==NULL || frtmfirst==NULL || frcnt==NULL) {
    fprintf(stderr, "Error: DrawSpikeInfo: Allocation failed.'n");
    return ;
  }

  while (--k >= 0) {
    if (rafrtm[k] <= tfirst) break;
    if (frcnt[rafrid[k]]++ == 0)
      frtmlast[rafrid[k]] = rafrtm[k];
    else
      frtmfirst[rafrid[k]] = rafrtm[k];
  }
  int num_neu = g_num_neu;
  if (num_neu > 5) num_neu = 5;                  // only show first few of them
//  ostringstream sout1, sout2;
  ostringstream sout3, sout4;
//  sout1<<"fr time last  "<<std::fixed<<setprecision(1);
//  sout2<<"fr time frist "<<std::fixed<<setprecision(1);
  sout3<<"fr count "<<tlen<<"ms";
  sout4<<"fr interval   "<<std::fixed<<setprecision(1);
  for (int i=0; i<num_neu; i++) {
    if (frcnt[i] == 1) frtmfirst[i] = frtmlast[i];  // only one spike
//    sout1<<setw(10)<<frtmlast[i];
//    sout2<<setw(10)<<frtmfirst[i];
    sout3<<setw(10)<<frcnt[i];
    sout4<<setw(10)<<(frtmlast[i]-frtmfirst[i])/frcnt[i];
  }
//  ftexto(0*0.1-1.6, -17*0.1+0.8, 0, sout1.str().c_str());
//  ftexto(0*0.1-1.6, -18*0.1+0.8, 0, sout2.str().c_str());
  ftexto(0*0.1-1.6, -19*0.1+0.8, 0, sout3.str().c_str());
  ftexto(0*0.1-1.6, -20*0.1+0.8, 0, sout4.str().c_str());

  free(frtmfirst);
  free(frtmlast);
  free(frcnt);
}

void DrawGLScene()
{
  int index=0;
  char text[512];
  int LOCAL_DRAW_FLAG = size_neuronvar;

  for (index=0; index<STEPS_PER_DRAW; index++) {
    if (!RUN_DONE) {
      compute_perstep();
    }
  }

  /* draw stuff */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  glTranslatef(xdepth,ydepth,zdepth);

  if (RUN_DONE) {
    glColor3f(1,1,1);
    ftexto(0*0.1-1.6, -16*0.1+0.8, 0,
           "Program paused. Press Enter to continue. Press Esc to exit.");
    Sleep(200);                  // reduce the CPU usage a bit, seems useless
  }

  DrawSpikeInfo();
  Drawmenu(0.1, -1.6, 0.8);

  if (g_num_neu>0) {
    switch(abs(DRAW_FLAG%LOCAL_DRAW_FLAG)) {
      case 0:
        sprintf(text,"Voltage");
        break;
#if SMOOTH_CONDUCTANCE_USE
      case 1:
        sprintf(text,"gEX--GE");
        break;
      case 2:
        sprintf(text,"gEX--HE");
        break;
      case 3:
        sprintf(text,"gIn--GI");
        break;
      case 4:
        sprintf(text,"gIn--HI");
        break;
      case 5:
        sprintf(text,"refractory time");
        break;
#else
      case 1:
        sprintf(text,"gEX--GE");
        break;
      case 2:
        sprintf(text,"gIn--GI");
        break;
      case 3:
        sprintf(text,"refractory time");
        break;
#endif
      default:
        break;
    }

    // plot drawing style
    Drawstra(&(GLOBAL_STRA[
               0+abs(DRAW_FLAG%size_neuronvar)*g_num_neu
               ]),                   // struct strobe **stra
             minimum(10, g_num_neu),   // number of neurons (showed)
             1,                      // side length(unit)
             0.0, -0.5,              // x offset, y offset  0.2, 0.35
             0,    0,                // max, min, if (max<=min) auto determined
             -1,   1,   0,           // r, g, b
                                     // if r<0, coloringtype = g (auto set color)
             0                       // graph style: line plot
            );
    // drawra style
    Drawstra(&(GLOBAL_STRA[
               0+abs(DRAW_FLAG%size_neuronvar)*g_num_neu
               ]),                   // struct strobe **stra
             g_num_neu,                // number of neurons
             1,                      // side length(unit)
             0.0, -0.15,             // x offset, y offset  -0.8, -0.15,
             0, 0,                   // not use
             -1, 2, 0,               // not use
             1                       // graph style: color raster
            );
    ftexto(-0.8, 0.05, 0, text);
  }
  glutSwapBuffers();       // since double buffered
}
