#ifndef _LOOP_H_
#define _LOOP_H_

// this is just vector structure
struct vector {
  int vect_size;  // the length of vector
  int vect_index; // the next position to save data in the vector
  double *vect_value;  // data value
};

// this is also a vector structure but with loop structure, meaning that if the
// structure is full, the next new data will be saved to the first postion in the
// vector and overwrite the original data
struct loop {
  int loop_index; // the next position to write new data
  vector loop_vect;
};

// initial setup for vector
void vector_initialize(vector &a);

// set data value in the vector
void vector_set_value(vector &a, int dim, double v);

// allocate space for the vector
void vector_allocate(vector &a, int dim);

// re-allocate space for the vector, this normally corresponds to the case
// that the structure is full
void vector_reallocate(vector &a);

// insert one data into the vector
void vector_insert(vector &a, double v);

// destroy the vector
void vector_destroy(vector &a);

// copy data value from another vector
void vector_copy(vector &des, const vector &src);

// compute the sum of data in the vector
double vector_sum(vector &a);

// print out the data in the vector
void vector_output(vector &a);

// initial setup for loop
void loop_initialize(loop &a);

// allocate space for the vector
void loop_allocate(loop &a, int dim);

// insert one data into the loop
void loop_insert(loop &a, double data);

// compute the sum of data in the loop
double loop_sum(loop &a);

// destroy the loop
void loop_destroy(loop &a);

// print out the data in the loop
void loop_output(loop &a);

#endif
