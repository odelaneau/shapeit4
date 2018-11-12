##########################################
# SET CORRECTLY THESE 4 PATHS TO COMPILE #
##########################################
BOOST_LIB=
HTSLD_INC=
HTSLD_LIB=

#COMPILER MODE C++11
CXX=g++ -std=c++0x

#COMPILER FLAGS
CXXFLAG_REL=-O3
CXXFLAG_DBG=-g
LIB_FLAGS=

#FILE LISTS
BFILE=bin/shapeit4
HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

#DEFAULT VERSION (I.E. UNIGE DESKTOP RELEASE VERSION)
all: desktop	

#UNIGE DESKTOP RELEASE VERSION
desktop: LIB_FLAGS=-lz -lbz2 -lm -lpthread 
desktop: HTSLD_INC=$(HOME)/Tools/htslib-1.3
desktop: HTSLD_LIB=$(HOME)/Tools/htslib-1.3
desktop: BOOST_INC=/usr/include
desktop: BOOST_LIB=/usr/lib/x86_64-linux-gnu
desktop: CXXFLAG=$(CXXFLAG_REL)
desktop: IFLAG=-Isrc -I$(HTSLD_INC) -I$(BOOST_INC)
desktop: LIB_FILES=$(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
desktop: LDFLAG=$(CXXFLAG_REL)
desktop: $(BFILE)

#UNIGE DESKTOP RELEASE VERSION
cluster: LIB_FLAGS=-lz -lbz2 -lm -lpthread -lcurl -lcrypto -lssl -llzma
cluster: HTSLD_INC=/software/UHTS/Analysis/htslib/1.6/include
cluster: HTSLD_LIB=/software/UHTS/Analysis/htslib/1.6/lib
cluster: BOOST_INC=/software/include
cluster: BOOST_LIB=/software/lib64
cluster: CXXFLAG=$(CXXFLAG_REL)
cluster: IFLAG=-Isrc -I$(HTSLD_INC) -I$(BOOST_INC)
cluster: LIB_FILES=$(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
cluster: LDFLAG=$(CXXFLAG_REL)
cluster: $(BFILE)

#UNIL SERVER RELEASE VERSION
dcb2: LIB_FLAGS=-lz -lbz2 -lm -lpthread -lcurl -lcrypto -lssl -llzma 
dcb2: HTSLD_INC=/home/oldelan/lib/htslib-1.7
dcb2: HTSLD_LIB=/home/oldelan/lib/htslib-1.7
dcb2: BOOST_INC=/software/include
dcb2: BOOST_LIB=/software/lib64
dcb2: CXXFLAG=$(CXXFLAG_REL)
dcb2: IFLAG=-Isrc -I$(HTSLD_INC) -I$(BOOST_INC)
dcb2: LIB_FILES=$(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
dcb2: LDFLAG=$(CXXFLAG_REL)
dcb2: $(BFILE)

#UNIGE DESKTOP DEBUG VERSION
debug: LIB_FLAGS=-lz -lbz2 -lm -lpthread 
debug: HTSLD_INC=$(HOME)/Tools/htslib-1.3
debug: HTSLD_LIB=$(HOME)/Tools/htslib-1.3
debug: BOOST_INC=/usr/include
debug: BOOST_LIB=/usr/lib/x86_64-linux-gnu
debug: CXXFLAG=$(CXXFLAG_DBG)
debug: IFLAG=-Isrc -I$(HTSLD_INC) -I$(BOOST_INC)
debug: LIB_FILES=$(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
debug: LDFLAG=$(CXXFLAG_DBG)
debug: $(BFILE)

#COMPILATION RULES
$(BFILE): $(OFILE)
	$(CXX) $^ $(LIB_FILES) -o $@ $(LIB_FLAGS) $(LDFLAG)

obj/%.o: %.cpp $(HFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)

clean: 
	rm -f obj/*.o $(BFILE)
