#include "stdafx.h"
#include "raster.h"

#define OLD_RASTER_USE 1

void raster_initialize(raster &a)
{
  a.ras_size = -1;
  a.ras_index = -1;

  a.array_index = NULL;
  a.array_firingtime = NULL;
}

void raster_copy(raster &des, const raster &src)
{
  int i;
  if (des.ras_size > 0) {
    raster_destroy(des);
  }

  raster_allocate(des, src.ras_size);
  des.ras_index = src.ras_index;
  for (i=0; i<des.ras_index; i++) {
    des.array_index[i] = src.array_index[i];
    des.array_firingtime[i] = src.array_firingtime[i];
  }
}

//void raster_copy(raster &des, const raster &src)
//{
//  if (des.ras_size < src.ras_index) {
//    raster_destroy(des);
//    raster_allocate(des, src.ras_size);
//  }
//  des.ras_index = src.ras_index;
//  for (int i=0; i<des.ras_index; i++) {
//    des.array_index[i] = src.array_index[i];
//    des.array_firingtime[i] = src.array_firingtime[i];
//  }
//}

void raster_destroy(raster &a)
{
  a.ras_index = -1;
  a.ras_size = -1;
  free(a.array_index);
  a.array_index = NULL;
  free(a.array_firingtime);
  a.array_firingtime = NULL;
}

void raster_allocate(raster &a, int dim)
{
  raster_destroy(a);

  a.ras_index = 0;
  a.ras_size = dim;

  a.array_index = (int *)malloc(dim*sizeof(int));
  a.array_firingtime = (double *)malloc(dim*sizeof(double));

  for (int i=0; i<dim; i++) {
    a.array_index[i] = -1;
    a.array_firingtime[i] = -1.0;
  }
}

void raster_clear(raster &a)
{
  a.ras_index = 0;
  for (int i=0; i<a.ras_size; i++) {
    a.array_index[i] = -1;
    a.array_firingtime[i] = -1.0;
  }
}

#if OLD_RASTER_USE
void raster_reallocate(raster &a)
{
  if (a.ras_size<200)
    a.ras_size += RASTER_ADD_SIZE;
  else
    a.ras_size *= 1.5;

  int *new_array_index = (int *)malloc(sizeof(int)*a.ras_size);
  double *new_array_firingtime = (double *)malloc(sizeof(double)*a.ras_size);
  if (new_array_index==NULL || new_array_firingtime==NULL) {
    printf("Out of memory! return by void raster_reallocate(raster &a)");
    return ;
  }

  for (int i=0; i<a.ras_index; i++) {
    new_array_index[i] = a.array_index[i];
    new_array_firingtime[i] = a.array_firingtime[i];
  }

  free(a.array_index);
  a.array_index = new_array_index;

  free(a.array_firingtime);
  a.array_firingtime = new_array_firingtime;
}

void raster_insert(raster &a, int index, double time)
{
  if ((a.ras_size-a.ras_index) < RASTER_ENDURE_SIZE) {
    raster_reallocate(a);
  }
  a.array_index[a.ras_index] = index;
  a.array_firingtime[a.ras_index] = time;
  a.ras_index ++;
}
#else // not saving all raster event, but only some last spike events
void raster_insert(raster &a, int index, double time)
{
  a.array_index[a.ras_index % a.ras_size] = index;
  a.array_firingtime[a.ras_index % a.ras_size] = time;

  a.ras_index++;
};
#endif // whether save all raster events

void raster_erase(raster &a, int index)
{
  if (index >= a.ras_index) {
    P_ERR("the erasing element is beyond raster's range!");
    return;
  }
  a.ras_index = a.ras_index - 1;
  for (int i=index; i<a.ras_index; i++) {
    a.array_index[i] = a.array_index[i+1];
    a.array_firingtime[i] = a.array_firingtime[i+1];
  }
  a.array_index[a.ras_index] = -1;
  a.array_firingtime[a.ras_index] = -1.0;
}

void raster_output(raster &a)
{
  P_ERR("output the firing events begin!\n");
  for (int i=0; i<a.ras_index; i++) {
    fprintf(stderr, "index: %d firing time: %g\n",
            a.array_index[i], a.array_firingtime[i]);
  }
  P_ERR("output the firing events end!\n");
}

// quick sort the spike_list in every time interval
void raster_quicksort(raster &spike_list)
{
  if (spike_list.ras_index > 0) {
    raster_q_sort(spike_list, 0, spike_list.ras_index - 1);
  }
}

void raster_q_sort(raster &spike_list, int left, int right)
{
  int index_pivot, l_hold, r_hold;
  double pivot;

  l_hold = left;
  r_hold = right;
  pivot = spike_list.array_firingtime[left];
  index_pivot = spike_list.array_index[left];

  while (left < right) {
    while ((spike_list.array_firingtime[right] >= pivot) && (left < right))
      right--;
    if (left != right) {
      spike_list.array_firingtime[left] = spike_list.array_firingtime[right];
      spike_list.array_index[left] = spike_list.array_index[right];
      left++;
    }
    while ((spike_list.array_firingtime[left] <= pivot) && (left < right))
      left++;
    if (left != right) {
      spike_list.array_firingtime[right] = spike_list.array_firingtime[left];
      spike_list.array_index[right] = spike_list.array_index[left];
      right--;
    }
  }
  spike_list.array_firingtime[left] = pivot;
  spike_list.array_index[left] = index_pivot;

  index_pivot = left;
  left = l_hold;
  right = r_hold;

  if (left < index_pivot)
    raster_q_sort(spike_list, left, index_pivot-1);
  if (right > index_pivot)
    raster_q_sort(spike_list, index_pivot+1, right);
}
