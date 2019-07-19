#ifndef __MAPP__logics__
#define __MAPP__logics__
#include "var_op.h"
enum {VERB,AUX_VERB};

/*------------------------------------------
 _       _____   _____   _   _____   _____  
| |     /  _  \ /  ___| | | /  ___| /  ___/ 
| |     | | | | | |     | | | |     | |___  
| |     | | | | | |  _  | | | |     \___  \ 
| |___  | |_| | | |_| | | | | |___   ___| | 
|_____| \_____/ \_____/ |_| \_____| /_____/
 ------------------------------------------*/
namespace MAPP_NS
{
    class Logics
    {
    private:
    protected:
        class Log;

        template<class VarClass> class LogUnary;
        template<class VarClass> class LogIS;
        
        template<class VarClass> class LogBinary;
        template<class VarClass> class LogGT;
        template<class VarClass> class LogGE;
        template<class VarClass> class LogEQ;
        template<class VarClass> class LogLE;
        template<class VarClass> class LogLT;
        class LogBinaryOp;
        class LogAND;
        class LogOR;
        class LogCond;
        class LogIF;
        class LogIFF;

        template<class VarClass>
        void creat_op(VarClass*,bool,VarClass*,bool,const char*);
        template<class VarClass>
        void creat_op(VarClass*,bool,const char*);
        
    public:
        Logics();
        template<class T0,class T1>
        Logics(T0&,const char*,T1&);
        template<class T0,class T1>
        Logics(T0&,const char*,T1&&);
        template<class T0>
        Logics(T0&,const char*);
        Logics(Logics&);
        Logics(Logics&&);
        ~Logics();
        
        
        
        
        Logics& operator = (Logics&);
        Logics& operator = (Logics&&);
        Logics&& operator + (Logics&&);
        Logics&& operator * (Logics&&);
        Logics&& operator - (Logics&&);
        Logics&& operator / (Logics&&);
        Logics& operator += (Logics&&);
        Logics& operator *= (Logics&&);
        Logics& operator -= (Logics&&);
        Logics& operator /= (Logics&&);
        
        Log* op;
        Logics* next;
        
        Logics& link();
        
        
        std::string print(void*);
        std::string operator()(void*);
        std::string operator()();
        static void var_log(Var*,Logics*);
        static bool var_log_nothrow(Var*,Logics*);
    };
}
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T0,class T1>
Logics::Logics(T0& var_l,const char* s,T1& var_r)
{
    bool l_alloc=false,r_alloc=false;
    Var* l=Var::find_var(var_l);
    Var* r=Var::find_var(var_r);

    if(!l)
    {
        l=new var<T0>(std::move(var_l));
        l_alloc=true;
    }
    if(!r)
    {
        r=new var<T1>(std::move(var_r));
        r_alloc=true;
    }
    creat_op(l,l_alloc,r,r_alloc,s);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T0,class T1>
Logics::Logics(T0& var_l,const char* s,T1&& var_r)
{
    bool l_alloc=false,r_alloc=false;
    Var* l=Var::find_var(var_l);
    Var* r=Var::find_var(var_r);
    
    if(!l)
    {
        l=new var<T0>(std::move(var_l));
        l_alloc=true;
    }
    if(!r)
    {
        r=new var<T1>(std::move(var_r));
        r_alloc=true;
    }
    creat_op(l,l_alloc,r,r_alloc,s);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T0>
Logics::Logics(T0& var_l,const char* s)
{
    bool l_alloc=false;
    Var* l=Var::find_var(var_l);
    /*
    if(!l)
    {
        l=new var<T0>(var_l);
        l_alloc=true;
    }*/
    creat_op(l,l_alloc,s);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
void Logics::creat_op(VarClass* l,bool l_alloc,const char* op_name)
{
    if(!strcmp(op_name,"set"))
        op=new LogIS<VarClass>(l,l_alloc,true);
    else if(!strcmp(op_name,"nset") || !strcmp(op_name,"!set"))
        op=new LogIS<VarClass>(l,l_alloc,false);
    
    next=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
void Logics::creat_op(VarClass* l,bool l_alloc,VarClass* r,bool r_alloc,const char* op_name)
{
    if(!strcmp(op_name,"gt"))
        op=new LogGT<VarClass>(l,l_alloc,r,r_alloc,true);
    else if(!strcmp(op_name,"ngt") || !strcmp(op_name,"!gt"))
        op=new LogGT<VarClass>(l,l_alloc,r,r_alloc,false);
    else if(!strcmp(op_name,"ge"))
        op=new LogGE<VarClass>(l,l_alloc,r,r_alloc,true);
    else if(!strcmp(op_name,"nge") || !strcmp(op_name,"!ge"))
        op=new LogGE<VarClass>(l,l_alloc,r,r_alloc,false);
    else if(!strcmp(op_name,"eq"))
        op=new LogEQ<VarClass>(l,l_alloc,r,r_alloc,true);
    else if(!strcmp(op_name,"neq") || !strcmp(op_name,"!eq"))
        op=new LogEQ<VarClass>(l,l_alloc,r,r_alloc,false);
    else if(!strcmp(op_name,"lt"))
        op=new LogLT<VarClass>(l,l_alloc,r,r_alloc,true);
    else if(!strcmp(op_name,"nlt") || !strcmp(op_name,"!lt"))
        op=new LogLT<VarClass>(l,l_alloc,r,r_alloc,false);
    else if(!strcmp(op_name,"le"))
        op=new LogLE<VarClass>(l,l_alloc,r,r_alloc,true);
    else if(!strcmp(op_name,"nle") || !strcmp(op_name,"!le"))
        op=new LogLE<VarClass>(l,l_alloc,r,r_alloc,false);
    next=NULL;
}
/*----------------------------------------------------
 _     _   _       _____   _____   _   _____   _____  
| |   / / | |     /  _  \ /  ___| | | /  ___| /  ___/ 
| |  / /  | |     | | | | | |     | | | |     | |___  
| | / /   | |     | | | | | |  _  | | | |     \___  \ 
| |/ /    | |___  | |_| | | |_| | | | | |___   ___| | 
|___/     |_____| \_____/ \_____/ |_| \_____| /_____/
 ----------------------------------------------------*/
namespace MAPP_NS
{
    class VLogics: public Logics
    {
    private:
    protected:
    public:
        template<typename T>
        VLogics(const char*,T&&);
        template<typename T>
        VLogics(const char*,T&);
        VLogics(const char* s);
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
VLogics::VLogics(const char* s,T& var_r):
Logics()
{
    bool r_alloc=false;
    Var* r=Var::find_var(var_r);
    if(!r)
    {
        r=new var<T>(std::move(var_r));
        r_alloc=true;
    }
    Var* l=NULL;
    creat_op(l,false,r,r_alloc,s);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class T>
VLogics::VLogics(const char* s,T&& var_r):
Logics()
{
    bool r_alloc=false;
    Var* r=Var::find_var(var_r);
    if(!r)
    {
        r=new var<T>(std::move(var_r));
        r_alloc=true;
    }
    Var* l=NULL;
    creat_op(l,false,r,r_alloc,s);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline VLogics::VLogics(const char* s):
Logics()
{
    Var* l=NULL;
    creat_op(l,false,s);
}
/*----------------------
 _       _____   _____ 
| |     /  _  \ /  ___|
| |     | | | | | |    
| |     | | | | | |  _ 
| |___  | |_| | | |_| |
|_____| \_____/ \_____/
 ----------------------*/
namespace MAPP_NS
{
    class Logics::Log
    {
    private:
    protected:
    public:
        static const std::string verb_aux;
        static const std::string negate;
        static const std::string verb;
        static const std::string s_verb;
        static const std::string p_verb;
        static const std::string s_n_verb;
        static const std::string p_n_verb;
        std::string op_name;
        const bool _is_;
        Log(std::string&&,const bool);
        virtual ~Log();
        void finish_sntnc(std::string&,std::string&,std::string&,std::string&,bool&,std::string&,int,bool&);
        void start_sntnc(std::string&,std::string&,std::string&,int,const std::string&,const std::string&);
        void process_sntnc(std::string&,std::string&,std::string&,std::string&,bool&,std::string&,int,bool&,const std::string&,const std::string&);
        virtual bool operator() ()=0;
        virtual bool operator() (void*)=0;
        virtual void print(std::string&,std::string&,std::string&,std::string&,bool&,std::string&,int,bool&)=0;
        virtual void print(std::string&,std::string&,std::string&,std::string&,bool&,std::string&,int,bool&,void*)=0;
        int eq(const char*,const char*);
        virtual Log* clone()=0;
    };
}
/*-------------------------------------------
 _   _   __   _       ___   _____   __    __
| | | | |  \ | |     /   | |  _  \  \ \  / /
| | | | |   \| |    / /| | | |_| |   \ \/ / 
| | | | | |\   |   / / | | |  _  /    \  /  
| |_| | | | \  |  / /  | | | | \ \    / /   
\_____/ |_|  \_| /_/   |_| |_|  \_\  /_/    
 -------------------------------------------*/
namespace MAPP_NS
{
    template<class VarClass>
    class Logics::LogUnary:public Log
    {
    public:
        VarClass* left;
        const bool left_alloc;
        LogUnary(VarClass*,const bool,std::string&&,const bool);
        virtual ~LogUnary();
        void print(std::string&,std::string&,std::string&,std::string&,bool&,std::string&,int,bool&);
        void print(std::string&,std::string&,std::string&,std::string&,bool&,std::string&,int,bool&,void*);
        virtual Log* clone()=0;
        virtual bool operator() ()=0;
        virtual bool operator() (void*);
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
Logics::LogUnary<VarClass>::LogUnary(VarClass* l,const bool l_alloc,std::string&& op_name,const bool _is_):
Log(std::move(op_name),_is_),
left(l),left_alloc(l_alloc)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
Logics::LogUnary<VarClass>::~LogUnary()
{
    if(left_alloc) delete left;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
void Logics::LogUnary<VarClass>::print(std::string& err_msg,std::string& L,std::string& op,std::string& R,bool& __is__,std::string& del,int v_flag,bool& p_flag)
{
    process_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag,left->name,std::string());
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
void Logics::LogUnary<VarClass>::print(std::string& err_msg,std::string& L,std::string& op,std::string& R,bool& __is__,std::string& del,int v_flag,bool& p_flag,void* v)
{
    int i=0;
    if(left==NULL)
    {
        
        left=static_cast<VarClass*>(v);
        i++;
    }
    process_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag,left->name,std::string());
    
    if(i) left=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
bool Logics::LogUnary<VarClass>::operator()(void* v)
{
    int i=0;
    if(left==NULL)
    {
        left=static_cast<VarClass*>(v);
        i++;
    }
    
    bool ans =this->operator()();

    
    if(i)
        left=NULL;
    return ans;
}
/*----------
 _   _____
| | /  ___/
| | | |___ 
| | \___  \
| |  ___| |
|_| /_____/
 ----------*/
namespace MAPP_NS
{
    
    template<class VarClass>
    class Logics::LogIS:public LogUnary<VarClass>
    {
    public:
        LogIS(VarClass*,const bool,const bool);
        bool operator() ();
        Log* clone();
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
Logics::LogIS<VarClass>::LogIS(VarClass* l,const bool l_alloc,const bool _is_):LogUnary<VarClass>(l,l_alloc,std::string("set"),_is_)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
bool Logics::LogIS<VarClass>::operator() ()
{
    if(Log::_is_)
        return (LogUnary<VarClass>::left->is_set());
    else
        return !(LogUnary<VarClass>::left->is_set());
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
Logics::Log* Logics::LogIS<VarClass>::clone()
{
    LogIS<VarClass>* log=new LogIS(LogUnary<VarClass>::left,LogUnary<VarClass>::left_alloc,Log::_is_);
    if(LogUnary<VarClass>::left_alloc)
        log->LogUnary<VarClass>::left=LogUnary<VarClass>::left->clone();
    return log;
}
/*-----------------------------------------------
 _____   _   __   _       ___   _____   __    __
|  _  \ | | |  \ | |     /   | |  _  \  \ \  / /
| |_| | | | |   \| |    / /| | | |_| |   \ \/ / 
|  _  { | | | |\   |   / / | | |  _  /    \  /  
| |_| | | | | | \  |  / /  | | | | \ \    / /   
|_____/ |_| |_|  \_| /_/   |_| |_|  \_\  /_/    
 -----------------------------------------------*/
namespace MAPP_NS
{
    template<class VarClass>
    class Logics::LogBinary:public Log
    {
    public:
        VarClass* left;
        bool left_alloc;
        VarClass* right;
        bool right_alloc;
        LogBinary(VarClass*,const bool,VarClass*,const bool,const char*,const bool);
        ~LogBinary();
        void print(std::string&,std::string&,std::string&,std::string&,bool&,std::string&,int,bool&);
        void print(std::string&,std::string&,std::string&,std::string&,bool&,std::string&,int,bool&,void*);
        virtual Log* clone()=0;
        virtual bool operator() ()=0;
        virtual bool operator() (void*);
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
Logics::LogBinary<VarClass>::LogBinary(VarClass* l,const bool l_alloc,VarClass* r,const bool r_alloc,const char* op_name,const bool _is_):
Log(op_name,_is_),left(l),left_alloc(l_alloc),right(r),right_alloc(r_alloc)
{}
/*--------------------------------------------
   
 --------------------------------------------*/
template<class VarClass>
Logics::LogBinary<VarClass>::~LogBinary()
{
    if(left_alloc) delete left;
    if(right_alloc) delete right;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
void Logics::LogBinary<VarClass>::print(std::string& err_msg,std::string& L,std::string& op,std::string& R,bool& __is__,std::string& del,int v_flag,bool& p_flag)
{
    process_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag,left->name,right->name);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
void Logics::LogBinary<VarClass>::print(std::string& err_msg,std::string& L,std::string& op,std::string& R,bool& __is__,std::string& del,int v_flag,bool& p_flag,void* v)
{
    int i=0;
    if(left==NULL)
    {
        left=static_cast<VarClass*>(v);
        i++;
    }
    if(right==NULL)
    {
        right=static_cast<VarClass*>(v);
        i+=2;
    }
    process_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag,left->name,right->name);
    
    if(i>1)
    {
        right=NULL;
        i-=2;
    }
    
    if(i)
        left=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
bool Logics::LogBinary<VarClass>::operator()(void* v)
{
    int i=0;
    if(left==NULL)
    {
        left=static_cast<VarClass*>(v);
        i++;
    }
    if(right==NULL)
    {
        right=static_cast<VarClass*>(v);
        i+=2;
    }
    bool ans =this->operator()();
    
    if(i>1)
    {
        right=NULL;
        i-=2;
    }
    
    if(i)
        left=NULL;
    return ans;
}
/*--------------
 _____   _____ 
/  ___| |_   _|
| |       | |  
| |  _    | |  
| |_| |   | |  
\_____/   |_|  
 --------------*/
namespace MAPP_NS
{
    template<class VarClass>
    class Logics::LogGT:public LogBinary<VarClass>
    {
    public:
        LogGT(VarClass*,const bool,VarClass*,const bool,const bool);
        bool operator() ();
        Log* clone();
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
Logics::LogGT<VarClass>::LogGT(VarClass* l,const bool l_alloc,VarClass* r,const bool r_alloc,const bool _is_):LogBinary<VarClass>(l,l_alloc,r,r_alloc,"greater than",_is_)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
bool Logics::LogGT<VarClass>::operator() ()
{
    if(Log::_is_)
        return  (*LogBinary<VarClass>::left > *LogBinary<VarClass>::right);
    else
        return !(*LogBinary<VarClass>::left > *LogBinary<VarClass>::right);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
Logics::Log* Logics::LogGT<VarClass>::clone()
{
    LogGT<VarClass>* log=new LogGT(LogBinary<VarClass>::left,LogBinary<VarClass>::left_alloc,LogBinary<VarClass>::right,LogBinary<VarClass>::right_alloc,Log::_is_);
    if(LogBinary<VarClass>::left_alloc)
        log->LogBinary<VarClass>::left=LogBinary<VarClass>::left->clone();
    if(LogBinary<VarClass>::right_alloc)
        log->LogBinary<VarClass>::right=LogBinary<VarClass>::right->clone();
    return log;
}
/*--------------
 _____   _____ 
/  ___| | ____|
| |     | |__  
| |  _  |  __| 
| |_| | | |___ 
\_____/ |_____|
 --------------*/
namespace MAPP_NS
{
    template<class VarClass>
    class Logics::LogGE:public LogBinary<VarClass>
    {
    public:
        LogGE(VarClass*,const bool,VarClass*,const bool,const bool);
        bool operator() ();
        Log* clone();
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
Logics::LogGE<VarClass>::LogGE(VarClass* l,const bool l_alloc,VarClass* r,const bool r_alloc,const bool _is_):LogBinary<VarClass>(l,l_alloc,r,r_alloc,"greater than or equal to",_is_)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
bool Logics::LogGE<VarClass>::operator() ()
{
    if(Log::_is_)
        return  (*LogBinary<VarClass>::left >= *LogBinary<VarClass>::right);
    else
        return !(*LogBinary<VarClass>::left >= *LogBinary<VarClass>::right);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
Logics::Log* Logics::LogGE<VarClass>::clone()
{
    Logics::LogGE<VarClass>* log=new LogGE(LogBinary<VarClass>::left,LogBinary<VarClass>::left_alloc,LogBinary<VarClass>::right,LogBinary<VarClass>::right_alloc,Log::_is_);
    if(LogBinary<VarClass>::left_alloc)
        log->LogBinary<VarClass>::left=LogBinary<VarClass>::left->clone();
    if(LogBinary<VarClass>::right_alloc)
        log->LogBinary<VarClass>::right=LogBinary<VarClass>::right->clone();
    return log;
}
/*----------------
 _____   _____   
| ____| /  _  \  
| |__   | | | |  
|  __|  | | | |  
| |___  | |_| |_ 
|_____| \_______|
 ----------------*/
namespace MAPP_NS
{
    template<class VarClass>
    class Logics::LogEQ:public LogBinary<VarClass>
    {
    public:
        LogEQ(VarClass*,const bool,VarClass*,const bool,const bool);
        bool operator() ();
        Log* clone();
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
Logics::LogEQ<VarClass>::LogEQ(VarClass* l,const bool l_alloc,VarClass* r,const bool r_alloc,const bool _is_):LogBinary<VarClass>(l,l_alloc,r,r_alloc,"equal to",_is_)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
bool Logics::LogEQ<VarClass>::operator() ()
{
    if(Log::_is_)
        return  (*LogBinary<VarClass>::left == *LogBinary<VarClass>::right);
    else
        return !(*LogBinary<VarClass>::left == *LogBinary<VarClass>::right);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
Logics::Log* Logics::LogEQ<VarClass>::clone()
{
    LogEQ<VarClass>* log=new LogEQ(LogBinary<VarClass>::left,LogBinary<VarClass>::left_alloc,LogBinary<VarClass>::right,LogBinary<VarClass>::right_alloc,Log::_is_);
    if(LogBinary<VarClass>::left_alloc)
        log->LogBinary<VarClass>::left=LogBinary<VarClass>::left->clone();
    if(LogBinary<VarClass>::right_alloc)
        log->LogBinary<VarClass>::right=LogBinary<VarClass>::right->clone();

    return log;
}
/*--------------
 _       _____ 
| |     | ____|
| |     | |__  
| |     |  __| 
| |___  | |___ 
|_____| |_____|
 --------------*/
namespace MAPP_NS
{
    template<class VarClass>
    class Logics::LogLE:public LogBinary<VarClass>
    {
    public:
        LogLE(VarClass*,const bool,VarClass*,const bool,const bool);
        bool operator() ();
        Log* clone();
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
Logics::LogLE<VarClass>::LogLE(VarClass* l,const bool l_alloc,VarClass* r,const bool r_alloc,const bool _is_):LogBinary<VarClass>(l,l_alloc,r,r_alloc,"less than or equal to",_is_)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
bool Logics::LogLE<VarClass>::operator() ()
{
    if(Log::_is_)
        return  (*LogBinary<VarClass>::left <= *LogBinary<VarClass>::right);
    else
        return !(*LogBinary<VarClass>::left <= *LogBinary<VarClass>::right);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
Logics::Log* Logics::LogLE<VarClass>::clone()
{
    LogLE<VarClass>* log=new LogLE(LogBinary<VarClass>::left,LogBinary<VarClass>::left_alloc,LogBinary<VarClass>::right,LogBinary<VarClass>::right_alloc,Log::_is_);
    if(LogBinary<VarClass>::left_alloc)
        log->LogBinary<VarClass>::left=LogBinary<VarClass>::left->clone();
    if(LogBinary<VarClass>::right_alloc)
        log->LogBinary<VarClass>::right=LogBinary<VarClass>::right->clone();
    return log;
}
/*--------------
 _       _____ 
| |     |_   _|
| |       | |  
| |       | |  
| |___    | |  
|_____|   |_|  
 --------------*/
namespace MAPP_NS
{
    template<class VarClass>
    class Logics::LogLT:public LogBinary<VarClass>
    {
    public:
        LogLT(VarClass*,const bool,VarClass*,const bool,const bool);
        bool operator() ();
        Log* clone();
    };
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
Logics::LogLT<VarClass>::LogLT(VarClass* l,const bool l_alloc,VarClass* r,const bool r_alloc,const bool _is_):LogBinary<VarClass>(l,l_alloc,r,r_alloc,"less than",_is_)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
bool Logics::LogLT<VarClass>::operator() ()
{
    if(Log::_is_)
        return  (*LogBinary<VarClass>::left < *LogBinary<VarClass>::right);
    else
        return !(*LogBinary<VarClass>::left < *LogBinary<VarClass>::right);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class VarClass>
Logics::Log* Logics::LogLT<VarClass>::clone()
{
    LogLT<VarClass>* log=new LogLT(LogBinary<VarClass>::left,LogBinary<VarClass>::left_alloc,LogBinary<VarClass>::right,LogBinary<VarClass>::right_alloc,Log::_is_);
    if(LogBinary<VarClass>::left_alloc)
        log->LogBinary<VarClass>::left=LogBinary<VarClass>::left->clone();
    if(LogBinary<VarClass>::right_alloc)
        log->LogBinary<VarClass>::right=LogBinary<VarClass>::right->clone();
    return log;
}
/*---------------------------------------------------------------
 _____   _   __   _       ___   _____   __    __  _____   _____
|  _  \ | | |  \ | |     /   | |  _  \  \ \  / / /  _  \ |  _  \
| |_| | | | |   \| |    / /| | | |_| |   \ \/ /  | | | | | |_| |
|  _  { | | | |\   |   / / | | |  _  /    \  /   | | | | |  ___/
| |_| | | | | | \  |  / /  | | | | \ \    / /    | |_| | | |    
|_____/ |_| |_|  \_| /_/   |_| |_|  \_\  /_/     \_____/ |_|
 ---------------------------------------------------------------*/
namespace MAPP_NS
{
    class Logics::LogBinaryOp:public Log
    {
    public:
        Log* left;
        Log* right;
        LogBinaryOp(Log*,Log*,const char*,const bool);
        ~LogBinaryOp();
        
        void print(std::string&,std::string&,std::string&,std::string&,bool&,std::string&,int,bool&);
        void print(std::string&,std::string&,std::string&,std::string&,bool&,std::string&,int,bool&,void*);
        virtual Log* clone()=0;
        virtual bool operator() ()=0;
        virtual bool operator() (void*)=0;
    };
}
/*-------------------------
     ___   __   _   _____
    /   | |  \ | | |  _  \
   / /| | |   \| | | | | |
  / / | | | |\   | | | | |
 / /  | | | | \  | | |_| |
/_/   |_| |_|  \_| |_____/
 -------------------------*/
namespace MAPP_NS
{
    class Logics::LogAND:public LogBinaryOp
    {
    public:
        LogAND(Log*,Log*);
        bool operator() ();
        bool operator() (void*);
        Log* clone();
    };
}
/*---------------
 _____   _____
/  _  \ |  _  \ 
| | | | | |_| | 
| | | | |  _  / 
| |_| | | | \ \ 
\_____/ |_|  \_\
 ---------------*/
namespace MAPP_NS
{
    class Logics::LogOR:public LogBinaryOp
    {
    public:
        LogOR(Log*,Log*);
        bool operator() ();
        bool operator() (void*);
        Log* clone();
    };
}
/*-------------------------------
 _____   _____   __   _   _____ 
/  ___| /  _  \ |  \ | | |  _  \
| |     | | | | |   \| | | | | |
| |     | | | | | |\   | | | | |
| |___  | |_| | | | \  | | |_| |
\_____| \_____/ |_|  \_| |_____/
 -------------------------------*/
namespace MAPP_NS
{
    class Logics::LogCond:public Log
    {
    public:
        Log* left;
        Log* right;
        LogCond(Log*,Log*,const char*,const bool);
        ~LogCond();
        
        virtual void print(std::string&,std::string&,std::string&,std::string&,bool&,std::string&,int,bool&)=0;
        virtual void print(std::string&,std::string&,std::string&,std::string&,bool&,std::string&,int,bool&,void*)=0;
        virtual Log* clone()=0;
        virtual bool operator() ()=0;
        virtual bool operator() (void*)=0;
    };
}
/*----------
 _   _____
| | |  ___|
| | | |__
| | |  __| 
| | | |    
|_| |_|    
 ----------*/
namespace MAPP_NS
{
    class Logics::LogIF:public LogCond
    {
    public:
        LogIF(Log*,Log*);
        bool operator() ();
        bool operator() (void*);
        void print(std::string&,std::string&,std::string&,std::string&,bool&,std::string&,int,bool&);
        void print(std::string&,std::string&,std::string&,std::string&,bool&,std::string&,int,bool&,void*);
        Log* clone();
    };
}
/*------------------
 _   _____   _____
| | |  ___| |  ___|
| | | |__   | |__  
| | |  __|  |  __| 
| | | |     | |    
|_| |_|     |_|
 ------------------*/
namespace MAPP_NS
{
    class Logics::LogIFF:public LogCond
    {
    public:
        LogIFF(Log*,Log*);
        bool operator() ();
        bool operator() (void*);
        void print(std::string&,std::string&,std::string&,std::string&,bool&,std::string&,int,bool&);
        void print(std::string&,std::string&,std::string&,std::string&,bool&,std::string&,int,bool&,void*);
        Log* clone();
    };
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/









#endif /* logics_hpp */
