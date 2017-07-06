CC      = gcc
OPTS    = -O2
OPTS    = -g
INC     = -I './'
CFLAGS  = $(OPTS) $(INC) -DSINGLE
LDFLAGS = $(OPTS)

.c.o:
	$(CC) $(CFLAGS) -c $<

SRCS = $(wildcard *.c)
OBJS = $(SRCS:%.c=%.o)

all: driver
driver: $(OBJS)
	$(CC) $(OBJS) -o driver -lm

clean:
	rm -rf *.o driver 
