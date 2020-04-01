###########################################################
#
#		Fast Power Grid Simulator 
# 		(Netlist Parser and Nodal Voltage Solver)
#
###########################################################
#
#     $RCSfile: Makefile,v $
#     $Date: 2020/03/17 $
#     Authors: Duo Li, Sheldon Tan, and Han Zhou
#
#     Functions: Makefile
#
###########################################################

MAIN = pgsim_cmd

ITPP = /home/locker/EE/mscad/num_lib/itpp-4.3.1
CXSPARSE = /home/locker/EE/mscad/num_lib/CSparse
CXSPARSE = /home/locker/EE/mscad/num_lib/CXSparse
UMFPACK = /home/locker/EE/mscad/num_lib/umfpack/UMFPACK
UFCONFIG = /home/locker/EE/mscad/num_lib/umfpack/UFconfig
AMD = /home/locker/EE/mscad/num_lib/umfpack/AMD

INCFLAGS = -I$(ITPP) -I$(UMFPACK)/Include -I$(UMFPACK)/Source -I$(UFCONFIG) \
	-I$(CXSPARSE)/Include  -I$(AMD)/Include 

LIBFLAGS = -L$(ITPP)/lib -L$(UMFPACK)/Lib -L$(AMD)/Lib -Wl,-rpath,./ 
#LIBFLAGS = -L$(ITPP)/lib -L$(UMFPACK)/Lib -L$(AMD)/Lib -Wl,-rpath,$(ITPP)/lib

LIBDIR = /home/locker/EE/mscad/num_lib
LIBS = $(LIBDIR)/libumfpack.a $(LIBDIR)/libamd.a $(LIBDIR)/libcxsparse.a 

FLAGS_DEBUG = `/home/locker/EE/mscad/num_lib/itpp-4.3.1/bin/itpp-config --debug --cflags`
FLAGS_OPT   = `/home/locker/EE/mscad/num_lib/itpp-4.3.1/bin/itpp-config --cflags`

LIBS_DEBUG  = `/home/locker/EE/mscad/num_lib/itpp-4.3.1/bin/itpp-config --debug --static --libs`
LIBS_OPT    = `/home/locker/EE/mscad/num_lib/itpp-4.3.1/bin/itpp-config --static --libs`

SRCS = element.cpp mna.cpp parser.cpp function_ibm.cpp

MAINSRCS =em_cmd_ibm.cpp 

OBJS = $(addsuffix .o, $(basename $(SRCS)))
MAINOBJS = $(addsuffix .o, $(basename $(MAINSRCS)))

CPP = g++ 
#CPP = g++  -static 

THREADFLAGS = -pthread

DEBUGFLAGS = -g #-DNDEBUG -O3

$(MAIN): $(OBJS) $(MAINOBJS) $(TMAINOBJS) $(LIBS) $(CU_OBJS)  
	@echo "Link  em_cmd ...."
	@echo $(CPP) $(INCFLAGS) $(LIBFLAGS) $(FLAGS_OPT) $(DEBUGFLAGS) $(THREADFLAGS) -o em_cmd $(OBJS) $(MAINOBJS) $(LIBS) $(LIBS_OPT) $(CU_OBJS) $(cuLIB)
	@$(CPP) $(INCFLAGS) $(LIBFLAGS) $(FLAGS_OPT) $(DEBUGFLAGS) $(THREADFLAGS) -o em_cmd $(OBJS) $(MAINOBJS) $(LIBS) $(LIBS_OPT) $(CU_OBJS) $(cuLIB)

.cpp.o:
	$(CPP) $(INCFLAGS) $(LIBFLAGS) $(DEBUGFLAGS) -c $<

depend:
	makedepend $(INCFLAGS) -- $(SRCS) # DO NOT DELETE

clean:
	@rm -rf *.o *~

