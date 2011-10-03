#include "stdafx.h"
#include "loop.h"

void vector_initialize(vector &a) // for a vector without being evaluated
{
  a.vect_size  = -1;
  a.vect_index = 0;
  a.vect_value = NULL;
}

void vector_set_value(vector &a, int dim, double v)
{
  free(a.vect_value);
  a.vect_size  = dim;
  a.vect_index = 0;
  a.vect_value = (double*)malloc(sizeof(double)*dim);
  P_NULL_ERR(a.vect_value, "Error: vector_set_value: allocation failed");
  for (int i = 0; i < dim; i++) {
    a.vect_value[i] = v;
  }
}

void vector_allocate(vector &a, int dim)
{
  free(a.vect_value);
  a.vect_index = 0;
  a.vect_size  = dim;
  a.vect_value = (double*)malloc(sizeof(double)*dim);
  P_NULL_ERR(a.vect_value, "Error: vector_set_value: Allocation failed!!");
}

void vector_reallocate(vector &a)
{
  if (a.vect_size<200)
    a.vect_size += VECTOR_ADD_SIZE;
  else                                 // changed by xyy
    a.vect_size *= 1.5;
  double *new_vect_value = (double*)malloc(a.vect_size*sizeof(double));
  if (new_vect_value == NULL) {
    fprintf(stderr, "Error: vector_reallocate: Allocation failed!!\n");
    return ;
  }
  for (int i=0; i<a.vect_index; i++) {
    new_vect_value[i] = a.vect_value[i];
  }
  free(a.vect_value);
  a.vect_value = new_vect_value;
}

void vector_insert(vector &a, double v)
{
  if (a.vect_index >= a.vect_size) {
    vector_reallocate(a);
  }

  a.vect_value[a.vect_index] = v;
  a.vect_index++;
}

void vector_destroy(vector &a)
{
  a.vect_size = -1;
  a.vect_index = 0;
  free(a.vect_value);
  a.vect_value = NULL;
}

void vector_copy(vector &des, const vector &src)
{
  free(des.vect_value);
  des.vect_size = src.vect_size;
  if (des.vect_size < 0) {
    P_ERR("Warning: vector_copy: The source vector is null!\n");
    return;
  }
  des.vect_value = (double*)malloc(sizeof(double)*src.vect_size);
  P_NULL_ERR(des.vect_value, "Error: vector_copy: Allocation fail.");
  memcpy(des.vect_value, src.vect_value, sizeof(double)*src.vect_size);
}

double vector_sum(vector &a)
{
  int i;
  double sum = 0.0;
  for (i=0; i<a.vect_size; i++) {
    sum += a.vect_value[i];
  }
  return sum;
}

void vector_output(vector &a)
{
  printf("vector output begins:\n");
  // a.vect_size implies whether the vector is evaluated or not
  for (int i=0; i<a.vect_size; i++) {
    printf("%g  ", a.vect_value[i]);
    if ((i+1)%6==0) {
      printf("\n");
    }
  }
  printf("vector output ends!\n");
}

void loop_initialize(loop &a)
{
  a.loop_index = -1;
  vector_initialize(a.loop_vect);
}

void loop_allocate(loop &a, int dim)
{
  a.loop_index = 0;
  vector_set_value(a.loop_vect, dim, 0);
}

void loop_insert(loop &a, double data)
{
//  a.loop_index = (a.loop_index % (a.loop_vect).vect_size);
  (a.loop_vect).vect_value[a.loop_index % (a.loop_vect).vect_size] = data;
  a.loop_index++;
}

double loop_sum(loop &a)
{
  return vector_sum(a.loop_vect);
}

void loop_destroy(loop &a)
{
  a.loop_index = -1;
  vector_destroy(a.loop_vect);
}

void loop_output(loop &a)
{
  printf("loop vector output begins:\n");

  int length = (a.loop_vect).vect_size;
  for (int i=0; i<length; i++) {
    printf("%g  ",(a.loop_vect).vect_value[(i+a.loop_index)%length]);
    if ((i+1)%6==0) {
      printf("\n");
    }
  }
  printf("loop vector output ends!\n");
}
