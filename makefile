#COMPILER MODE C++11
CXX=g++ -std=c++0x

#COMPILER FLAGS
CXXFLAG_REL=-O3
LIB_FLAGS=

#FILE LISTS
BFILE=bin/shapeit4
HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

#DEFAULT
all: default	

#UNIL & UNIGE DESKTOP
default: LIB_FLAGS=-lz -lbz2 -lm -lpthread -llzma 
default: HTSLD_INC=$(HOME)/Tools/htslib-1.9
default: HTSLD_LIB=$(HOME)/Tools/htslib-1.9
default: BOOST_INC=/usr/include
default: BOOST_LIB=/usr/lib/x86_64-linux-gnu
default: CXXFLAG=$(CXXFLAG_REL)
default: IFLAG=-Isrc -I$(HTSLD_INC) -I$(BOOST_INC)
default: LIB_FILES=$(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
default: LDFLAG=$(CXXFLAG_REL)
default: $(BFILE)

#COMPILATION RULES
$(BFILE): $(OFILE)
	$(CXX) $^ $(LIB_FILES) -o $@ $(LIB_FLAGS) $(LDFLAG)

obj/%.o: %.cpp $(HFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)

clean: 
	rm -f obj/*.o $(BFILE)
