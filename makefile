SHELL = /bin/sh
CC = g++
LIBS = -lgsl -lgslcblas -lfftw3 -lfftw3_omp -lm  -fopenmp
CPPFLAGS = -g -O2 -Wall

VPATH=%.h ./include
VPATH=%.o ./obj

OBJDIR = ./obj
INCLUDE_DIR = ./include
INCLUDES  := $(addprefix -I,$(INCLUDE_DIR))

objects = $(addprefix $(OBJDIR)/, test.o)

MY_APPS = zfgenerator

$(MY_APPS) : $(objects)
	${CC} -o ${MY_APPS} ${objects} ${CPPFLAGS} ${LIBS}

$(OBJDIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) ${LIBS} ${INCLUDES} $< -o $@


$(OBJDIR):
	mkdir $(OBJDIR)

.PHONY : clean
clean:
	rm -f ${MY_APPS}
	rm -f ${objects}
