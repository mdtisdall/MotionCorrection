#VERSION = DEBUG
VERSION = RELEASE

#CXXFLAGS += -I./
CXXFLAGS += -I/Users/dylan/Documents/Research/Students/Diana/MotionCorrection
CXXFLAGS += -I/usr/local/include
#CXXFLAGS += -I/Users/zyzdiana/GitHub/MotionCorrection/

LDFLAGS += -L/usr/local/lib
LDLIBS += -lfftw3f
LDLIBS += -lfftw3

ifeq ($(VERSION), DEBUG)                                                        
CXXFLAGS += -g
LDFLAGS += -g
endif                                                                           
ifeq ($(VERSION), RELEASE)                                                      
CXXFLAGS += -g
CXXFLAGS += -O3
LDFLAGS += -O3
endif  

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)                                                         
CXX = clang++
LDFLAGS += -framework Accelerate
CXXFLAGS += -DOSX
endif

ifeq ($(UNAME), Linux)
CXX = g++
CXXFLAGS += -DLINUX
endif


all: test test_read_file test_interp test_time_profiler test_rotate_coords

test_read_file:

test_interp:

test_time_profiler:

test_inverse:

test_rotate_coords:

test_tricubic:

test_gn:

fit_nav_reps: FFTWBuffer.o FFTOp.o

#test: BinaryFile_tests.o interp3D_tests.o
TESTOBJECTS += FFTWBuffer.o
TESTOBJECTS += FFTOp.o
TESTOBJECTS += BinaryFile_tests.o
TESTOBJECTS += Volume_tests.o
TESTOBJECTS += FFTWBuffer_tests.o
TESTOBJECTS += CentralDifferenceDifferentiator_tests.o
TESTOBJECTS += TrilinearInterpolator_tests.o
TESTOBJECTS += TricubicInterpolator_tests.o
TESTOBJECTS += CubicBSplineInterpolator_tests.o
TESTOBJECTS += FFTOp_tests.o
TESTOBJECTS += Gauss_Newton_Ref_Grad_tests.o
TESTOBJECTS += Gauss_Newton_New_Grad_tests.o
TESTOBJECTS += CircularMaskOp_tests.o
TESTOBJECTS += Weighted_Gauss_Newton_Ref_Grad_tests.o

test: $(TESTOBJECTS)

clean:
	rm -f *.o
	rm -f test
	rm -f fit_nav_reps
	rm -f test_read_file
	rm -f test_interp
	rm -f test_time_profiler
	rm -f test_tricubic
	rm -f test_gn

PHONY: .clean
