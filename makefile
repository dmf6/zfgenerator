export CC = g++
export RM= rm -f
export MAKE= make
CPPFLAGS = -Wall -g -O
LIBS = -lgsl -lgslcblas -lfftw3  -L/usr/lib64/mysql -lmysqlcppconn 
INCLUDES = -I/usr/include/mysql
OBJECTS = test.o
MYAPP = zfgenerator

all: $(MYAPP)

${MYAPP} : $(OBJECTS)
	${CC} -o $(MYAPP) ${OBJECTS} ${CPPFLAGS} ${LIBS} 

%.o : %.cpp
	${CC} ${CPPFLAGS} ${LIBS} ${INCLUDES} -c $^ -o $@
clean:
	$(RM) ${OBJECTS} $(MYAPP)