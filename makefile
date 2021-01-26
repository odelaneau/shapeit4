#COMPILER MODE C++11
CXX=g++ -std=c++11

#HTSLIB LIBRARY [SPECIFY YOUR OWN PATHS]
HTSLIB_INC=$(EBROOTHTSLIB)/include/htslib
HTSLIB_LIB=$(EBROOTHTSLIB)/lib/libhts.a

#BOOST IOSTREAM & PROGRAM_OPTION LIBRARIES [SPECIFY YOUR OWN PATHS]
BOOST_INC=$(EBROOTBOOST)/include/boost
BOOST_LIB_IO=$(EBROOTBOOST)/lib/libboost_iostreams.a
BOOST_LIB_PO=$(EBROOTBOOST)/lib/libboost_program_options.a

#HTSLIB LIBRARY [SPECIFY YOUR OWN PATHS]
#HTSLIB_INC=/software/UHTS/Analysis/samtools/1.4/include
#HTSLIB_LIB=/software/UHTS/Analysis/samtools/1.4/lib64/libhts.a

#BOOST IOSTREAM & PROGRAM_OPTION LIBRARIES [SPECIFY YOUR OWN PATHS]
#BOOST_INC=/software/include
#BOOST_LIB_IO=/software/lib64/libboost_iostreams.a
#BOOST_LIB_PO=/software/lib64/libboost_program_options.a

#COMPILER & LINKER FLAGS

#Best performance is achieved with this. Use it if running on the same plateform you're compiling.
#CXXFLAG=-O3 -march=native

#Good performance and portable on most CPUs
CXXFLAG=-O3 -mavx2 -mfma 

#Portable version without avx2 (much slower)
#CXXFLAG=-O3

LDFLAG=-O3


#DYNAMIC LIBRARIES
DYN_LIBS=-lz -lbz2 -lm -lpthread -llzma -lcurl -lssl -lcrypto

#SHAPEIT SOURCES & BINARY
BFILE=bin/shapeit4.2
HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

#COMPILATION RULES
all: $(BFILE)

$(BFILE): $(OFILE)
	$(CXX) $(LDFLAG) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -o $@ $(DYN_LIBS)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAG) -c $< -o $@ -Isrc -I$(HTSLIB_INC) -I$(BOOST_INC)

clean: 
	rm -f obj/*.o $(BFILE)
