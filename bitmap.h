#ifndef _BITMAP_h
#define _BITMAP_h

#include <stdint.h>   // for uint32_t

typedef uint32_t word_t;    

// bitmap function prototypes
void set_bit(word_t *, size_t); 

void clear_bit(word_t *, size_t);

size_t get_bit(word_t *, size_t);

word_t* create_bitmap(size_t);

void delete_bitmap(word_t *);

#endif
