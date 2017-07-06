#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <math.h>
#include "flatKDTree.h"
#include "bitmap.h"

REAL *base_arr;

bool almost_eq(REAL A, REAL B) {
  double abs_diff = fabs(A - B);
  if (abs_diff <= EPSILON)
    return true;

  Real_t uA, uB;
  uA.f = A;
  uB.f = B;

  // check the signs
  int shift = 8*sizeof(REAL)-1;
  if ((uA.i >> shift ) != ( uB.i >> shift))
    return false;

  int ulpsDiff = labs(uA.i - uB.i);
  if (ulpsDiff <= MAX_ULP) 
    return true;

  return false;
}

int coord_comptor(const void * a, const void * b ) {
    size_t *ptr_a = *(size_t**)a;
    size_t *ptr_b = *(size_t**)b;

    REAL fa = base_arr[*ptr_a];
    REAL fb = base_arr[*ptr_b];

    if (almost_eq(fa,fb))
        return (size_t)(ptr_a-ptr_b);

    return (fa > fb)-(fa < fb);
}

// stable sort kd-array structure-of-arrays
void sort_rows(const REAL **data, const unsigned short k, const size_t n, size_t **idx) {
    
  size_t **idx_ptr = (size_t**)malloc(n * sizeof(size_t*));
  size_t  *idx_tmp =  (size_t*)malloc(n * sizeof(size_t));
  
  unsigned short dim;
  for (dim = 0; dim < k; ++dim) {

    size_t i;
    for (i = 0; i < n; ++i)
      idx_ptr[i] = &(*idx)[i];

    base_arr = (REAL*)data[dim];
    qsort(idx_ptr, n, sizeof(size_t*), coord_comptor);

    for (i = 0; i < n; ++i)
      idx_tmp[i] = *(idx_ptr[i]);

    for (i = 0; i < n; ++i)
      (*idx)[i] = idx_tmp[i];
  }
  free(idx_ptr);
  free(idx_tmp);
}


size_t unique_rows(REAL ***data, const unsigned short k, size_t n) {
  size_t n_unique;
  size_t i;
  size_t *coord_idx = (size_t*)malloc(n*sizeof(size_t));

  for (i=0; i<n; ++i)
    coord_idx[i] = i;
  sort_rows((const REAL**)(*data),k,n,&coord_idx);

  word_t *delete_map = create_bitmap(n);
  for (i=1;i<n;++i) {
    bool dup = true;
    unsigned short j;
    for(j=0;j<k;++j)
      dup = dup && almost_eq((*data)[j][coord_idx[i]],(*data)[j][coord_idx[i-1]]);

    if (dup)
      set_bit(delete_map, i);
  }
  n_unique = n;
  for (i=0;i<n;++i) {
    if (get_bit(delete_map,i)) {
      while (get_bit(delete_map,n_unique-1))
        --n_unique;
      size_t j;
      for(j=0;j<k;++j)
        (*data)[j][i] = (*data)[j][n_unique-1];
      --n_unique;
    }
    if (i>=n_unique)
      break;
  }
  free(coord_idx);
  free(delete_map);
  return n_unique;
}

int comptor( const void * a,const void * b ) {
  return (*(size_t*)a-*(size_t*)b);
}

void find_partition_boundaries(size_t n, unsigned int depth, size_t **partition_boundaries ) {
  unsigned int i, num_partitions;
  size_t j, mid;

  (*partition_boundaries)[0] = 0;
  (*partition_boundaries)[1] = n-1;
    
  for (i=0;i<depth;++i) {
    num_partitions = (unsigned int)pow(2.0,i);

    for (j=0;j<num_partitions;++j) {
      mid = (*partition_boundaries)[j*2] + 
        (((*partition_boundaries)[j*2+1] - (*partition_boundaries)[j*2] ) >> 1);

      (*partition_boundaries)[(2*num_partitions)+(2*j)+0] = mid;
      (*partition_boundaries)[(2*num_partitions)+(2*j)+1] = mid+1;
    }
    qsort(*partition_boundaries,2*(2*num_partitions),sizeof(int),comptor);
  }
}

int partition(REAL **kd_array, size_t **tree_index, unsigned int order_dim, size_t left, size_t right, size_t pivot_index) {
  size_t i;
  REAL pivot_value;
  REAL *data;

  data = kd_array[order_dim];

  pivot_value = data[(*tree_index)[pivot_index]];
  int store_index = left;

  SWAP((*tree_index)[pivot_index], (*tree_index)[right], size_t);

  for (i=left;i<right;++i) {
    if (data[(*tree_index)[i]]<pivot_value) {
      SWAP((*tree_index)[store_index], (*tree_index)[i], size_t);
      ++store_index;
    }
  }

  SWAP((*tree_index)[right], (*tree_index)[store_index], size_t);

  return store_index;
}

void nth_element(REAL **kd_array, size_t **tree_index, unsigned int order_dim, size_t left, size_t right, size_t n) {
  size_t pivot_index;
  size_t rnd_pivot_index;
  if (left==right)
    return;

  while (1) {
    assert(right >= left);
    rnd_pivot_index = left+(rand()%(right-left+1));
    pivot_index = partition(kd_array, tree_index, order_dim, left, right, rnd_pivot_index );

    if (n == pivot_index)
      return;
    else if (n < pivot_index)
      right = pivot_index-1;
    else
      left = pivot_index+1;
  }
}

unsigned int longest_dim(REAL **kd_array, size_t *tree_index, size_t left, size_t right) {
  REAL min[DIM], max[DIM], dist[DIM];
  size_t i;

  unsigned int dim;
  for (dim=0;dim<DIM;++dim)
    max[dim] = min[dim] = kd_array[dim][tree_index[left]];

  for (i=left+1;i<=right;++i) {
    for (dim=0;dim<DIM;++dim) {
      REAL val = kd_array[dim][tree_index[i]];
      max[dim] = max[dim] < val ? val : max[dim];
      min[dim] = min[dim] > val ? val : min[dim];
    }
  }

  for (dim=0;dim<DIM;++dim)
    dist[dim] = max[dim] - min[dim];

  REAL max_dist = dist[0];
  unsigned int max_dim = 0;
  for (dim=1;dim<DIM;++dim) {
    if(max_dist<dist[dim]) {
      max_dist = dist[dim];
      max_dim  = dim;
    }
  }
  return max_dim;
}

void tree_order(REAL **kd_array, size_t **tree_index, unsigned int order_dim, size_t left, size_t right ) {
  size_t mid;
  unsigned int long_dim;
  size_t len = right - left;

  // base case
  if (len == 1) 
    return;

  // divide
  mid = left + (len >> 1);

  long_dim = longest_dim(kd_array, *tree_index, left, right);
  nth_element(kd_array, tree_index, long_dim, left, right, mid);

  long_dim = longest_dim(kd_array, *tree_index, left, mid);
  tree_order(kd_array, tree_index, long_dim, left, mid);

  long_dim = longest_dim(kd_array, *tree_index, mid, right);
  tree_order(kd_array, tree_index, long_dim, mid, right);
    
}

/*
void shuffle_partitions( int **tree_data, int **partition_boundaries, int num_partitions, int** size_partition_sample ) {

    int i, j, r, left, right;

    for ( i = 0; i < num_partitions / 2; ++i ) {
        left =  (*partition_boundaries)[ i * 2 ];
        right = (*partition_boundaries)[ i * 2 + 1];

        for ( j = left; j < left + (*size_partition_sample)[ i ]; ++j ) {
            r = ( rand() % ( right + 1 - j ) ) + j;
            SWAP( (*tree_data)[ j ], (*tree_data)[ r ], int );
        }
    }
}

double fRand( double fMin, double fMax )
{
    double f = ( double )rand() / RAND_MAX;
    return fMin + f * ( fMax - fMin );        
}

void sub_sample_data( particle_data *particles, int *tree_data_index, int *partition_boundaries, int num_partitions, 
    int *size_partition_sample, double **sample_data ) {

    int pid = 0;

    int i, j, r, index;
    for ( i = 0; i < num_partitions / 2; ++i ) {
        for ( j = 0; j < size_partition_sample[ i ]; ++j ) {

            index = tree_data_index[ partition_boundaries[ i * 2 ] + j ];

            (*sample_data)[ pid * 4 + 0 ] = ( double )particles->x[ index ] + fRand( -1e-7, 1e-7 );
            (*sample_data)[ pid * 4 + 1 ] = ( double )particles->y[ index ] + fRand( -1e-7, 1e-7 );
            (*sample_data)[ pid * 4 + 2 ] = ( double )particles->z[ index ] + fRand( -1e-7, 1e-7 );
            ++pid;
        }
    }
}


// perform a 3-d rotation of the points
void rotate3d( double **sample_data, int num_points, double theta_x, double theta_y, double theta_z ) {

    int i, j, m;
    double x, y, z;
    double Qx[ 3 ][ 3 ];
    double Qy[ 3 ][ 3 ];
    double Qz[ 3 ][ 3 ];
    double Qt1[ 3 ][ 3 ];
    double Qt2[ 3 ][ 3 ];

    Qx[ 0 ][ 0 ] = 1.0f;            Qx[ 0 ][ 1 ] = 0.0f;           Qx[ 0 ][ 2 ] = 0.0f;
    Qx[ 1 ][ 0 ] = 0.0f;            Qx[ 1 ][ 1 ] = cos(theta_x);   Qx[ 1 ][ 2 ] = sin(theta_x);
    Qx[ 2 ][ 0 ] = 0.0f;            Qx[ 2 ][ 1 ] = -sin(theta_x);  Qx[ 2 ][ 2 ] = cos(theta_x);

    Qy[ 0 ][ 0 ] = cos(theta_y);    Qy[ 0 ][ 1 ] = 0.0f;            Qy[ 0 ][ 2 ] = -sin(theta_y); 
    Qy[ 1 ][ 0 ] = 0.0f;            Qy[ 1 ][ 1 ] = 1.0f;            Qy[ 1 ][ 2 ] = 0.0f; 
    Qy[ 2 ][ 0 ] = sin(theta_y);    Qy[ 2 ][ 1 ] = 0.0f;            Qy[ 2 ][ 2 ] = cos(theta_y);

    Qz[ 0 ][ 0 ] = cos(theta_z);    Qz[ 0 ][ 1 ] = sin(theta_z);    Qz[ 0 ][ 2 ] = 0.0f; 
    Qz[ 1 ][ 0 ] = -sin(theta_z);   Qz[ 1 ][ 1 ] = cos(theta_z);    Qz[ 1 ][ 2 ] = 0.0f; 
    Qz[ 2 ][ 0 ] = 0.0f;            Qz[ 2 ][ 1 ] = 0.0f;            Qz[ 2 ][ 2 ] = 1.0f;

    for ( i = 0; i < 3; ++i ) {
        for( j = 0; j < 3; ++j ) {
            Qt1[ i ][ j ] = 0;
            for( m = 0; m < 3; ++m ) 
                Qt1[ i ][ j ] += Qx[ i ][ m ] * Qy[ m ][ j ];
        }
    }

    for ( i = 0; i < 3; ++i ) {
        for( j = 0; j < 3; ++j ) {
            Qt2[ i ][ j ] = 0;
            for( m = 0; m < 3; ++m ) 
                Qt2[ i ][ j ] += Qt1[ i ][ m ] * Qz[ m ][ j ];
        }
    }

    for ( i = 0; i < num_points; ++i ) {
        x = (*sample_data)[ i * 4 + 0 ];
        y = (*sample_data)[ i * 4 + 1 ];
        z = (*sample_data)[ i * 4 + 2 ];
        (*sample_data)[ i * 4 + 0 ] = x * Qt2[ 0 ][ 0 ] + y * Qt2[ 0 ][ 1 ] + z * Qt2[ 0 ][ 2 ];
        (*sample_data)[ i * 4 + 1 ] = x * Qt2[ 1 ][ 0 ] + y * Qt2[ 1 ][ 1 ] + z * Qt2[ 1 ][ 2 ];
        (*sample_data)[ i * 4 + 2 ] = x * Qt2[ 2 ][ 0 ] + y * Qt2[ 2 ][ 1 ] + z * Qt2[ 2 ][ 2 ];
    }
}
*/
