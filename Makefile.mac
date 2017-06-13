SHELL		=/bin/bash
CC		=c++ 
MAKEFILE	=Makefile
CFLAGS		=-std=c++11 -O3 


INCLUDE_MPI	=-I/usr/local/include
LIB_MPI		=-L/usr/local/lib
INCLUDE_PY	=-I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7
LIB_PY		=-L/System/Library/Frameworks/Python.framework/Versions/2.7
INCLUDE_NP	=-I/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/numpy/core/include
LIB_NP		=-L/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/numpy/core/lib


INCLUDES	=$(INCLUDE_MPI) $(INCLUDE_PY) $(INCLUDE_NP)
LIBS		=$(LIB_MPI) $(LIB_PY) $(LIB_NP)
LDFLAGS		=-lpython2.7 -lnpymath -lpthread -ldl -lutil -lmpi 
        
SRC		=src/
OBJ		=obj/
CPP_FILES	=$(wildcard $(SRC)*.cpp)
OBJ_FILES	=$(addprefix $(OBJ),$(notdir $(CPP_FILES:.cpp=.o))) 
OBJ_PY_FILES	=$(addprefix $(OBJ),$(notdir $(CPP_FILES:.cpp=.o.py))) 
SITE_PACKS	=/Library/Python/2.7/site-packages


mapp: prep $(OBJ_FILES) $(MAKEFILE) 
	$(CC) $(CFLAGS) $(OBJ_FILES) -o $@ $(LIBS) $(LDFLAGS)  

py: prep $(OBJ_PY_FILES) $(MAKEFILE)
	$(CC)  -bundle -undefined dynamic_lookup $(CFLAGS) $(OBJ_PY_FILES) $(LIBS) $(LDFLAGS) -o mapp.so

install:
	@sudo mv mapp.so $(SITE_PACKS)/mapp.so;

$(OBJ)%.o: $(SRC)%.cpp $(MAKEFILE)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@	


$(OBJ)%.o.py: $(SRC)%.cpp $(MAKEFILE)
	$(CC) -fPIC $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:  
	rm -rf $(OBJ)


prep:
	@mkdir -p $(OBJ);
