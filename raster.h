///mark: modified on 2013/1/20
#ifndef _RASTER_H_
#define _RASTER_H_

// this is data structure to save raster (the X-axis is the firing time, and
// the Y-axis is the neuron index in the picture)
struct raster {
  int ras_size; // the length of both array_index and array_firingtime
  int ras_index; // the next position in the raster structure to save data

  int *array_index; // saves the index of firing neuron
  double *array_firingtime; // saves the firing time of the neuron
};

// initial setup for raster structure
void raster_initialize(raster &a);

// copy data from another raster structure
void raster_copy(raster &des, const raster &src);

// destroy the structure
void raster_destroy(raster &a);

// allocate space to save data (spike event)
void raster_allocate(raster &a, int dim);

// re-allocate space to save raster data, normally, this corresponds to the case
// that the raster data is full
void raster_reallocate(raster &a);

// insert one data (a spike event) in the structure
void raster_insert(raster &a, int index, double time);

// remove one data (a spike event) from the structure
void raster_erase(raster &a, int index);

// print out the raster information
void raster_output(raster &a);

// use fast sorting method to obtain a spike list with
// increasing order of firing time
void raster_quicksort(raster &spike_list);

// one part of functions for fast sorting method
void raster_q_sort(raster &spike_list, int left, int right);

#endif
