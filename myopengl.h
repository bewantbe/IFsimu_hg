#ifndef _OPENGL_H_
#define _OPENGL_H_

/* The following function are used in open GL, some of them are not used in this
code, but may be used in the future. */

// this is used to print out some string "text" at the position (x,y,z) on the screen
void ftexto(float x, float y, float z, const char *text);

// open GL initial set up
GLvoid InitGL(GLsizei Width, GLsizei Height);

// reset the size of window shown on the screen
GLvoid ReSizeGLScene(GLsizei Width, GLsizei Height);

/* The following functions are used to draw different things, some of them are
used in the code, others are not used, but may be used in the future. */
// this is used to draw a box with center located at (x,y,z) and with color (r,g,b)
GLvoid glbox(double x,double y,double r,double rcolor,double gcolor,double bcolor);

// this is used to draw menu on the screen to test the range of parameters in network
// the left bottom point is located at (xoffset, yoffset), the size of each line
// in the menu is given by "side"
GLvoid Drawmenu(double side,double xoffset,double yoffset);

// these are used for corresponding changes of the picture shown on the screen when
// some keys on the keyboard are pressed
void keyPressed(unsigned char key, int x, int y);
void specialKeyPressed(int key, int x, int y);

// this is the most important function to draw several things on the screen
void DrawGLScene();

/* This is the most useful drawing functions for computational data. it draws a graph
with 2d-structure, for example, if we want to see the evolution of voltage for each
neuron in the network, we can use x-axis to indicate the time, y-axis indicates the
index of neuron, the color indicates the value of voltage. Therefore, "ra" corresponds
to the saved data (it does not matter how you save your data since "ra" is just a
pointer with any type), "type" means the data type (either "double" or "int" or
something else), "rows" is the length of x-axis (here is time length of window bin),
"cols" is the length of y-axis (here is the neuron index), "coloringtype" means the
type of the picture, 0 is greyscale, 7 is colorscale, other values between them are
seldomly used. The "stduse=1" means using "mean-STD_VIEW*stdev" as minimum and
"mean+STD_VIEW*stdev" as maximum to determine (r,g,b) for each data value of "ra",
otherwise (stduse=0) just use the absolute minimum and maximum to determine.
"guse=1" means smearing the picture, "guse=2" means not. Here, the "guse" are not
used.  "text" is the string you want to write in the picture such as the name of
picture, x coordinate and y-coordinate. "xside" and "yside" means the width and the
height of the picture on the screen, where x,y coordinates are rescaled to be the
coordinates of the picture. "xoffset" and "yoffset" means the starting position of
the openGL screen. */
GLvoid Drawra(void *ra,const char *type,int rows,int cols,int coloringtype,int stduse,int gsuse,
              const char *text,double xside,double yside,double xoffset,double yoffset);

/* This is for drawing normal picture with (y) vs (x) from data "ra". The
x-coordinate is the index of the sequence "ra" and y-coordinate is the data value of
"ra". "max" and "min" are either given or can be computed using "stats", they are
used to rescale the y-coordinate. "xside" and "yside" means the width and the height of
the picture on the screen, where x,y coordinates are rescaled to be the coordinates of
the picture. "xoffset" and "yoffset" are the starting position on the screen.
"rcolor", "gcolor" and "bcolor" is the assigned color. */
GLvoid Plotra(void *ra,char *type,int length,double max,double min,
              double xside,double yside,double xoffset,double yoffset,
              double rcolor,double gcolor,double bcolor);

/* This function draws the graph of "stra" which is a particular structure "strobe"
recording some data variable in a fixed time length (details can be seen from the
definition of strobe). Here, we save voltage, conductance, refractory residence time
for each neuron, therefore, total number of strobes is size_neuronvar*N (N is the toal
number of neurons, size_neuronvar is the total number of variables), therefore, the
first pointer chooses from 0 to ralength-1 which means the index of strobes, The
second pointer means each strobe structure with the evolution of variables stored.
"side" actually corresponds to the size of the picture on the screen, where the x,y
coordinates are rescaled to be the coordinates in the picture. "xoffset" and "yoffset"
means the starting position of the openGL screen. "max" and "min" are either given
or can be computed using "stats", "rcolor", "gcolor" and "bcolor" is the assigned
color. "drawvsplot" corresponds to the different graph style on the screen. */
GLvoid Drawstra(struct strobe **stra,int ralength,double side,double xoffset,double yoffset,
                double max,double min,double rcolor,double gcolor,double bcolor,
                int drawvsplot);

#endif
