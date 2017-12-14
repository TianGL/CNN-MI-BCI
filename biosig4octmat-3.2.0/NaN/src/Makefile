####################################################
# Copyright 2010, 2011,2012 Alois Schloegl 
# This is part of the NaN-toolbox - a statistics and machine learning toolbox for data with and without missing values. 
# http://pub.ist.ac.at/~schloegl/matlab/NaN/
####################################################

### modify directories according to your needs 

# Define non-default octave-version 
# Octave - global install  (e.g. from debian package)
# OCTAVE_VERSION= 
# Better alternative: define an OCTAVE_VERSION bash variable (or in .bashrc or .profile)
# OCTAVE_VERSION=-3.6.3

# Matlab configuration
#MATLABDIR = /usr/local/MATLAB/R2010b
# comment the following line if you use MATLAB on 32-bit operating system
MEX_OPTION += -largeArrayDims

# Mingw crosscompiler: available at http://www.nongnu.org/mingw-cross-env/
CROSS   = $(HOME)/src/mxe.github.schloegl/usr/bin/i686-pc-mingw32-
CROSS64 = $(HOME)/src/mxe.github.schloegl/usr/bin/x86_64-w64-mingw32-
# include directory for Win32-Matlab include
W32MAT_INC = $(HOME)/bin/win32/Matlab/R2010b/extern/include/
W64MAT_INC = $(HOME)/bin/win64/Matlab/R2010b/extern/include/
# path to GNUMEX libraries, available from here http://sourceforge.net/projects/gnumex/
GNUMEX   = $(HOME)/bin/win32/gnumex
GNUMEX64 = $(HOME)/bin/win64/gnumex
# building gnumex64 was difficult, these hints were quite useful: 
# http://sourceforge.net/mailarchive/forum.php?thread_name=AANLkTinZvxgC9ezp2P3UCX_a7TAUYuVsp2U40MQUV6qr%40mail.gmail.com&forum_name=gnumex-users
# Instead of building "mex shortpath.c" and "mex uigetpath.c", I used empty m-functions within argout=argin; 
####################################################

CC  	= gcc
CXX 	= g++
CFLAGS 	= -fopenmp -Wall -Wextra -Wconversion -O2 -fPIC 
OCTMEX	= mkoctfile$(OCTAVE_VERSION) --mex
RM      = rm

MEX_OPTION  += -lgomp CC\#$(CXX) CXX\#$(CXX) CFLAGS\#"$(CFLAGS) " CXXFLAGS\#"$(CFLAGS) "
MATMEX      = $(MATLABDIR)/bin/mex $(MEX_OPTION)

PROGS = histo_mex.mex covm_mex.mex kth_element.mex sumskipnan_mex.mex str2array.mex train.mex svmtrain_mex.mex svmpredict_mex.mex xptopen.mex 

### per default only the mex-files for octave are built
mex4o octave: $(PROGS)


### Matlab configuration - search for a matlab directory if not defined above
ifndef MATLABDIR
  ifneq ('',$(shell ls -1 /usr/local/MATLAB))
    # use oldest, typically mex-files a compatible with newer Matlab versions 
    MATLABDIR=$(shell ls -dt1 /usr/local/MATLAB/R* |tail -1)
    # alternatively, use latest matlab version
    #MATLABDIR=$(shell ls -dt1 /usr/local/MATLAB/R* |head -1)
  endif
endif

### if the MATLABDIR has been found or defined 
ifneq ('',MATLABDIR)
  ifneq ('',$(lshell ls -1 $(MATLABDIR)))
    MEX_EXT = $(shell $(MATLABDIR)/bin/mexext)
    mex4m matlab: $(patsubst %.mex, %.$(MEX_EXT), $(PROGS))
  endif
endif


mexw32 win32:  $(patsubst %.mex, %.mexw32, $(PROGS))
mexw64 win64:  $(patsubst %.mex, %.mexw64, $(PROGS))
all:    octave win32 win64 mex4m

clean:
	-$(RM) *.o *.obj *.o64 core octave-core *.oct *~ *.mex* 

#$(PROGS): Makefile

#########################################################
#	Octave, MATLAB on Linux 
#########################################################
%.oct: %.cc
	mkoctfile$(OCTAVE_VERSION) $<

%.mex: %.cpp
	$(OCTMEX) $<  -lgomp         ## Octave
%.$(MEX_EXT): %.cpp
	$(MATMEX) $<  -lgomp      ## Matlab

svmtrain_mex.mex svmpredict_mex.mex:     svmtrain_mex.cpp svm.h svm.cpp svm_model_matlab.c svmpredict_mex.cpp
	$(CXX) $(CFLAGS) -c svm.cpp
	#$(CC)  $(CFLAGS) $(INCOCTAVE) -c svm_model_matlab.c
	env CC=$(CXX) $(OCTMEX) svmtrain_mex.cpp svm.o svm_model_matlab.c
	env CC=$(CXX) $(OCTMEX) svmpredict_mex.cpp svm.o svm_model_matlab.c

svmtrain_mex.$(MEX_EXT) svmpredict_mex.$(MEX_EXT):     svmtrain_mex.cpp svm.h svm.cpp svm_model_matlab.c svmpredict_mex.cpp
	$(CXX) $(CFLAGS) -I $(MATLABDIR)/extern/include -c svm.cpp
	$(CC) $(CFLAGS) -I $(MATLABDIR)/extern/include -c svm_model_matlab.c
	$(MATMEX) svmtrain_mex.cpp svm.o svm_model_matlab.o
	$(MATMEX) svmpredict_mex.cpp svm.o svm_model_matlab.o

train.$(MEX_EXT) predict.$(MEX_EXT): train.c linear.h tron.o linear.o linear_model_matlab.c
	$(CC) $(CFLAGS) -I $(MATLABDIR)/extern/include -c linear_model_matlab.c
	$(MATMEX)  -lblas train.c tron.o linear.o linear_model_matlab.o 
	#$(MATMEX) -lblas predict.c tron.o linear.o linear_model_matlab.o

train.mex predict.mex: train.c linear.h tron.o linear.o linear_model_matlab.c
	env CC=$(CXX) $(OCTMEX) -lblas train.c tron.o linear.o linear_model_matlab.c
	#env CC=$(CXX) $(OCTMEX) -lblas predict.c tron.o linear.o linear_model_matlab.o

linear.o: linear.cpp linear.h
	$(CXX) $(CFLAGS) -c linear.cpp

tron.o: tron.cpp tron.h
	$(CXX) $(CFLAGS) -c tron.cpp


#########################################################
#	MATLAB/WIN32
#########################################################
%.obj: %.cpp
	$(CROSS)$(CXX) -fopenmp -c -DMATLAB_MEX_FILE -x c++ -o $@ -I$(W32MAT_INC) -O2 -DMX_COMPAT_32 $<
%.obj: %.c
	$(CROSS)$(CXX) -fopenmp -c -DMATLAB_MEX_FILE -x c++ -o $@ -I$(W32MAT_INC) -O2 -DMX_COMPAT_32 $< 

train.mexw32 predict.mexw32:  train.obj linear.obj linear_model_matlab.obj tron.obj
	$(CROSS)$(CXX) -shared $(GNUMEX)/mex.def -o $@ -L$(GNUMEX) -s $< linear_model_matlab.obj linear.obj tron.obj -llibmx -llibmex -llibmat -lcholmod -lblas

svmpredict_mex.mexw32 :  svmpredict_mex.obj svm.obj svm_model_matlab.obj 
	$(CROSS)$(CXX) -shared $(GNUMEX)/mex.def -o $@ -L$(GNUMEX) -s $< svm_model_matlab.obj svm.obj -llibmx -llibmex -llibmat -lcholmod
svmtrain_mex.mexw32 : svmtrain_mex.obj svm.obj svm_model_matlab.obj
	$(CROSS)$(CXX) -shared $(GNUMEX)/mex.def -o $@ -L$(GNUMEX) -s $< svm_model_matlab.obj svm.obj -llibmx -llibmex -llibmat -lcholmod

%.mexw32: %.obj
	$(CROSS)$(CXX) -shared $(GNUMEX)/mex.def -o $@ -L$(GNUMEX) -s $< -llibmx -llibmex -llibmat -lcholmod -lgomp -lpthread



#########################################################
#	MATLAB/WIN64
#########################################################

## ToDO: fix OpenMP support: currently -fopenmp causes Matlab to crash
%.o64: %.cpp
	$(CROSS64)$(CXX) -c -DMATLAB_MEX_FILE -x c++ -o $@ -I$(W64MAT_INC) -O2 $<
%.o64: %.c
	$(CROSS64)$(CXX) -c -DMATLAB_MEX_FILE -x c++ -o $@ -I$(W64MAT_INC) -O2 $< 

train.mexw64 predict.mexw64:  train.o64 linear.o64 linear_model_matlab.o64 tron.o64
	$(CROSS64)$(CXX) -shared $(GNUMEX64)/mex.def -o $@ -L$(GNUMEX64) -s $< linear_model_matlab.o64 linear.o64 tron.o64 -llibmx -llibmex -llibmat -lcholmod -lblas

svmpredict_mex.mexw64 :  svmpredict_mex.o64 svm.o64 svm_model_matlab.o64
	$(CROSS64)$(CXX) -shared $(GNUMEX64)/mex.def -o $@ -L$(GNUMEX64) -s $< svm_model_matlab.o64 svm.o64 -llibmx -llibmex -llibmat -lcholmod
svmtrain_mex.mexw64 : svmtrain_mex.o64 svm.o64 svm_model_matlab.o64
	$(CROSS64)$(CXX) -shared $(GNUMEX64)/mex.def -o $@ -L$(GNUMEX64) -s $< svm_model_matlab.o64 svm.o64 -llibmx -llibmex -llibmat -lcholmod

%.mexw64: %.o64
	$(CROSS64)$(CXX) -shared $(GNUMEX64)/mex.def -o $@ -L$(GNUMEX64) -s $< -llibmx -llibmex -llibmat -lcholmod -lgomp -lpthread

