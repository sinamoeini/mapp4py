#ifndef __MAPP__read__
#define __MAPP__read__
#include <mpi.h>
namespace MAPP_NS
{
    class Read
    {
    private:
    protected:
        MPI_Comm& world;
    public:
        Read(MPI_Comm&);
        virtual ~Read();

    };
}
#include <fstream>
namespace MAPP_NS
{
    class FileReader
    {
    private:
        std::ifstream fst;
        std::string str_line;
        MPI_Comm world;
        const int rank;
        size_t line_no;
    protected:
    public:
        char* file;
        bool finished;
        
        FileReader(class Communication*&,const char*);
        FileReader(class Communication*&,const std::string);
        FileReader(MPI_Comm,const char*);
        FileReader(MPI_Comm,const std::string);
        FileReader();
        ~FileReader();
        bool operator()(char*&,size_t&);
        bool skip(const int);
        size_t get_line_no();
        
        static size_t parse_line(char*,char**&,size_t&);
        static size_t hash_remover(char*&);
        static size_t concatenate(size_t,char**,char*&);
        
    };
}

#endif

