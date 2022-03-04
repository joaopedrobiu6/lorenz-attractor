#### definitions

CXX := g++ -std=c++17
CXXFLAGS := $(shell root-config --cflags) -fPIC 

#### rules (target: dependencies -> actions)

SRCS := $(wildcard src/*.cpp)
OBJS := $(patsubst %.cpp, bin/%.o, $(notdir $(SRCS)) )
EXES := $(patsubst %.cpp, bin/%.exe, $(notdir $(wildcard main/*.cpp)) )
ROOT_INC := -I $(shell root-config --incdir)
ROOT_LIB := $(shell root-config --libs)
EIGEN_INC := -I /usr/include/eigen3

####
.PRECIOUS: $(OBJS)

VPATH = main:src

all: $(EXES)

lib: lib/libBIU.a

lib/%.a: $(OBJS)
	@echo making lib...[$^]
	ar ruv $@ $^ 
	ranlib $@

bin/%.exe: bin/%.o lib/libBIU.a
	@echo executable... $< [$@]
	$(CXX) $(CXXFLAGS) -o $@ $< -I src $(EIGEN_INC)  -L $(ROOT_LIB) -L lib -l BIU $(ROOT_LIB) 

bin/%.o: %.cpp
	@echo compiling... $< [$@]
	$(CXX) $(CXXFLAGS) -c -o $@ $< -I src $(EIGEN_INC) 

run:
	./bin/main.exe

py:
	python3.9 main/ODEdmd.py


clean:
	@echo cleaning...
	rm -f $(wildcard bin/*) $(wildcard lib/*)

clean_images:
	@echo cleaning...
	rm -f $(wildcard images/*.png)

dump:
	@echo SRCS...[$(SRCS)] [$(OBJS)] [$(EXES)]

