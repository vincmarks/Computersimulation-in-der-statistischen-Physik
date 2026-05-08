TARGET  = P1000
SUBTREE =

SRCs= melt.c \
      r250.c \
      sysin.c \
      sysout.c \
      wach.c \
      element.c \
      mcmove.c \
      cbgc.c \
      local.c \
      rep.c \
      pivot.c \
      pressure.c \
      knot.c \
      taylor.c \
      alexander_polynomial.c \
      numerical.c \
      weight.c

## OBJs:=$(subst .c,.o,$(SRCs))

OBJs= melt.o \
      r250.o \
      sysin.o \
      sysout.o \
      wach.o \
      element.o \
      mcmove.o \
      cbgc.o \
      local.o \
      rep.o \
      pivot.o \
      pressure.o \
      knot.o \
      taylor.o \
      alexander_polynomial.o \
      numerical.o \
      weight.o
 
HDRs= element.h \
      r250.h 

#CC        = env TARGET=cray-t3e cc
#CFLAGS    = $(STDCFLAGS) -O2 -X16 -I /usr/include/mpp
#CC         = icc
#CFLAGS     = $(STDCFLAGS) -O3 -xW -tpp7 -ipo -ipo_obj 
CC         = gcc
STDCFLAGS  = 
CFLAGS     = $(STDCFLAGS) -O3
#CFLAGS     = -q32 -O4
LDFLAGS    = -lm
DEBUG      = $(STDCFLAGS) -O0 -g
VERBOSE    =

.PHONY	: clean debug verbose backup


say	:
	echo $(SRCs)
	echo $(OBJs)

## Here is how to create the objects:

%.o	: %.c $(HDRs)
	$(CC) $(CFLAGS) -c $<

mc	: $(OBJs)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJs) $(LDFLAGS) 

generate : r250.o generate.o
	$(CC) $(CFLAGS) -o G r250.o generate.o $(LDFLAGS)

debug   : $(OBJs)
	$(CC) -o $(TARGET) $(OBJs) $(LDFLAGS) $(DEBUG)

verbose : $(SRCs)
	$(CC) -DVERBOSE -g -o bfl bdibfl.c $(STDCFLAGS) $(LDFLAGS)

backup  : $(SRCs)
	zoo auh melt.zoo $(SRCs)

distclean:
	rm -ef $(OBJs)

clean	:
	rm -f $(OBJs) $(TARGET)

new     :
	touch *.c
