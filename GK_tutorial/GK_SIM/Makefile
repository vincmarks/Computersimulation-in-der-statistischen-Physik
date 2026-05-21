TARGET  = LJ
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
      pressure.c \
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
      pressure.o \
      weight.o
 
HDRs= element.h \
      r250.h 

#CC        = env TARGET=cray-t3e cc
#CFLAGS    = $(STDCFLAGS) -O2 -X16 -I /usr/include/mpp
CC         = cc
STDCFLAGS  = 
CFLAGS     = $(STDCFLAGS) -O2
LDFLAGS    = -lm
DEBUG      = $(STDCFLAGS) -g
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
	rm -ef $(OBJs) $(TARGET)

new     :
	touch *.c
