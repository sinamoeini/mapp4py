#include "atoms.h"
#include "timer.h"
#include "memory.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor of timer
 --------------------------------------------*/
Timer::Timer()
{
    nmodes=0;
    curr_level=-1;
    level2mode_capacity=0;
    
    add_mode("comm ");
    add_mode("neigh");
    add_mode("force");
    add_mode("write");
    add_mode("other");
}
/*--------------------------------------------
 destructor of timer
 --------------------------------------------*/
Timer::~Timer()
{
    for(int imod=0;imod<nmodes;imod++)
        delete [] mode_names[imod];
    if(nmodes)
    {
        delete [] time;
        delete [] mode_names;
    }
    
    if(level2mode_capacity)
        delete [] level2mode;
}
/*--------------------------------------------
 destructor of timer
 --------------------------------------------*/
int Timer::add_mode(const char* mode_name)
{
    int imod=0;
    while(imod<nmodes && strcmp(mode_name,mode_names[imod]))
        imod++;

    if(imod!=nmodes)
        return imod;
    
    int lnght=static_cast<int>(strlen(mode_name))+1;
    GROW(time,nmodes,nmodes+1);
    GROW(mode_names,nmodes,nmodes+1);
    time[nmodes]=0.0;
    CREATE1D(mode_names[nmodes],lnght);
    memcpy(mode_names[nmodes],mode_name,lnght*sizeof(char));
    nmodes++;
    return nmodes-1;
}
/*--------------------------------------------
 start communication time;
 --------------------------------------------*/
void Timer::start(int mode)
{
    type0 t=MPI_Wtime();
    if(curr_level!=-1)
        time[level2mode[curr_level]]+=t;
    
    curr_level++;
    if(curr_level==level2mode_capacity)
    {
        GROW(level2mode,level2mode_capacity,level2mode_capacity+1);
        level2mode_capacity++;
    }
    level2mode[curr_level]=mode;
    
    time[mode]-=t;
    
}
/*--------------------------------------------
 start communication time;
 --------------------------------------------*/
void Timer::stop(int mode)
{
    /*-------temp_remove-------
    if(mode!=level2mode[curr_level])
        Error::abort("timer error");
    */
        
    type0 t=MPI_Wtime();
    time[mode]+=t;
    curr_level--;
    if(curr_level!=-1)
        time[level2mode[curr_level]]-=t;
}
/*--------------------------------------------
 start communication time;
 --------------------------------------------*/
type0 Timer::mode_time(int mode)
{
    return time[mode];
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Timer::init()
{
    type0 t=MPI_Wtime();
    tot_time=-t;
    for(int i=0;i<nmodes;i++)
        time[i]=0.0;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Timer::fin()
{
    type0 t=MPI_Wtime();
    tot_time+=t;

    t=tot_time;
    for(int i=0;i<nmodes;i++)
        t-=time[i];
    
    time[OTHER_TIME_mode]+=t;
}
/*--------------------------------------------
 print
 --------------------------------------------*/
void Timer::print_stats()
{
    /*-------temp_remove-------
    if(atoms->my_p==0)
    {
        fprintf(output,"time stats:\n");
        fprintf(output,"total time: %lf secs\n"
        ,tot_time);
        
        for(int imod=0;imod<nmodes;imod++)
        {
            fprintf(output,"%s time: %lf secs (%05.2lf%%)\n"
            ,mode_names[imod],time[imod],time[imod]*100.0/tot_time);
        }
    }
     */
}
/*--------------------------------------------
 print
 --------------------------------------------*/
type0 Timer::tst_start()
{
    tst_time=-MPI_Wtime();
    return 0.0;
}
/*--------------------------------------------
 print
 --------------------------------------------*/
type0 Timer::tst_stop()
{
    tst_time+=MPI_Wtime();
    return tst_time;
}
