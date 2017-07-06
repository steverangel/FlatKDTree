#ifdef __cplusplus
extern "C" {
#endif

#ifndef _FLATKDTREE_h
#define _FLATKDTREE_h

#define SWAP( a, b, type ) { type tmp; tmp = ( a ); ( a ) = ( b ); ( b ) = tmp; }

#ifndef DIM
#define DIM 3
#endif

#ifndef MAX_ULP
#define MAX_ULP 0
#endif

#ifdef SINGLE
#define REAL float
#define INT  int
#define EPSILON FLT_EPSILON
#endif

#ifdef DOUBLE
#define REAL double
#define INT  long long
#define EPSILON DBL_EPSILON
#endif

typedef enum { false, true } bool;    

typedef union Real_t{
  REAL f;
  INT  i;
} Real_t;

bool          almost_eq   (REAL   , REAL);
size_t        unique_rows (REAL***, const unsigned short, size_t); 
int           partition   (REAL** , size_t**            , unsigned int, size_t, size_t, size_t); 
void          nth_element (REAL** , size_t**            , unsigned int, size_t, size_t, size_t); 
unsigned int  longest_dim (REAL** , size_t*             , size_t      , size_t); 
void          tree_order  (REAL** , size_t**            , unsigned int, size_t, size_t);

#endif

#ifdef __cplusplus
}
#endif
