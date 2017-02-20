SHELL 	    = /bin/bash
SRC         = src/
PROGRAMS    = main
CC          = /opt/openmpi-intel/bin/mpic++ 
OBJ         = obj/
MAKEFILE    = Makefile
CFLAGS      = -std=c++11 -O3 
LIBS        = -L/home/sina/apps/python/lib/python2.7/site-packages/numpy/core/lib -L/home/sina/apps/python/lib  
INCLUDES    = -I/home/sina/apps/python/include/python2.7/ -I/home/sina/apps/python/lib/python2.7/site-packages/numpy/core/include
LDFLAGS	    = -lpython2.7 -lnpymath -lpthread -ldl -lutil 
        
CPP_FILES   = $(wildcard $(SRC)*.cpp)
H_FILES     = $(wildcard $(SRC)*.h)
OBJ_FILES   = $(addprefix $(OBJ),$(notdir $(CPP_FILES:.cpp=.o))) 


$(OBJ)%.o: $(SRC)%.cpp $(MAKEFILE)
	$(CC) -c $(CFLAGS) -o $@ $(INCLUDES) $<	

mapp:	prep $(OBJ_FILES) $(MAKEFILE) 
	$(CC) $(CFLAGS) $(OBJ_FILES) -o $@ $(LIBS) $(LDFLAGS)  

clean:  
	rm -rf $(OBJ)
prep:
	@mkdir -p $(OBJ); 


