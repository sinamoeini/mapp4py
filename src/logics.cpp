#include "logics.h"

using namespace MAPP_NS;
const std::string Logics::Log::verb_aux =std::string("should");
const std::string Logics::Log::negate=std::string("not");
const std::string Logics::Log::verb=std::string("be");
const std::string Logics::Log::s_verb=std::string("is");
const std::string Logics::Log::p_verb=std::string("are");
const std::string Logics::Log::s_n_verb=std::string("are not");
const std::string Logics::Log::p_n_verb=std::string("should not be");
/*------------------------------------------
 _       _____   _____   _   _____   _____  
| |     /  _  \ /  ___| | | /  ___| /  ___/ 
| |     | | | | | |     | | | |     | |___  
| |     | | | | | |  _  | | | |     \___  \ 
| |___  | |_| | | |_| | | | | |___   ___| | 
|_____| \_____/ \_____/ |_| \_____| /_____/
 ------------------------------------------*/
Logics::Logics():
op(NULL),
next(NULL)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::~Logics()
{
    delete op;
    delete next;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::Logics(Logics& other)
{
    if(other.op!=NULL)
        this->op=other.op->clone();
    else
        this->op=NULL;
    
    if(other.next!=NULL)
    {
        this->next=new Logics(*(other.next));
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::Logics(Logics&& other)
{
    this->op=other.op;
    other.op=NULL;
    
    this->next=other.next;
    other.next=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Logics::operator=(Logics& other)
{
    delete this->op;
    if(other.op!=NULL)
        this->op=other.op->clone();
    else
        this->op=NULL;
    
    if(other.next!=NULL)
    {
        delete this->next;
        this->next=new Logics(*(other.next));
    }
    
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Logics::operator=(Logics&& other)
{
    
    delete this->op;
    this->op=other.op;
    other.op=NULL;
    
    
    delete this->next;
    this->next=other.next;
    other.next=NULL;
    
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics&& Logics::operator+(Logics&& other)
{
    if(this->op==NULL)
    {
        delete this->op;
        this->op=other.op;
        other.op=NULL;
        return std::move(*this);
    }
    Log* op_=new LogOR(this->op,other.op);
    op=op_;
    other.op=NULL;
    return std::move(*this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics&& Logics::operator*(Logics&& other)
{
    if(this->op==NULL)
    {
        delete this->op;
        op=other.op;
        other.op=NULL;
        return std::move(*this);
    }
    
    Log* op_=new LogAND(this->op,other.op);
    op=op_;
    other.op=NULL;
    return std::move(*this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics&& Logics::operator/(Logics&& other)
{
    
    if(this->op==NULL)
    {
        delete this->op;
        this->op=other.op;
        other.op=NULL;
        return std::move(*this);
    }
    
    
    Log* op_=new LogIF(this->op,other.op);
    op=op_;
    other.op=NULL;
    return std::move(*this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics&& Logics::operator-(Logics&& other)
{
    if(this->op==NULL)
    {
        delete this->op;
        op=other.op;
        other.op=NULL;
        return std::move(*this);
    }
    
    Log* op_=new LogIFF(this->op,other.op);
    op=op_;
    other.op=NULL;
    return std::move(*this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Logics::operator+=(Logics&& other)
{
    if(this->op==NULL)
    {
        delete this->op;
        op=other.op;
        other.op=NULL;
        return *this;
    }
    
    Log* op_=new LogOR(this->op,other.op);
    op=op_;
    other.op=NULL;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Logics::operator*=(Logics&& other)
{
    if(this->op==NULL)
    {
        op=other.op;
        other.op=NULL;
        return *this;
    }
    
    Log* op_=new LogAND(this->op,other.op);
    op=op_;
    other.op=NULL;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Logics::operator-=(Logics&& other)
{
    if(this->op==NULL)
    {
        delete this->op;
        op=other.op;
        other.op=NULL;
        return *this;
    }
    
    Log* op_=new LogIFF(this->op,other.op);
    op=op_;
    other.op=NULL;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Logics::operator/=(Logics&& other)
{
    if(this->op==NULL)
    {
        op=other.op;
        other.op=NULL;
        return *this;
    }
    
    Log* op_=new LogIF(this->op,other.op);
    op=op_;
    other.op=NULL;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::var_log(Var* v,Logics* logs)
{
    
    std::string err_msg=logs->operator()(v);
    if(!err_msg.empty()) throw err_msg;
    
    if(v->rank)
    {
        for(size_t i=0;i<v->size;i++)
        {
            try
            {
                var_log(&(v->operator[](i)),logs-1);
            }
            catch(std::string& err_msg)
            {
                throw err_msg;
            }
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Logics::var_log_nothrow(Var* v,Logics* logs)
{
    
    if(logs->op && !logs->op->operator()(v))
        return false;
    
    if(v->rank)
        for(size_t i=0;i<v->size;i++)
            if(!var_log_nothrow(&(v->operator[](i)),logs-1))
                return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
std::string Logics::print(void* v)
{
    int v_flag=AUX_VERB;
    if(op==NULL)
        return NULL;
    
    std::string err_msg;
    std::string del;
    std::string l;
    std::string r;
    std::string oper;
    std::string _verb_;
    bool p_flag;
    bool __is__;
    
    op->print(err_msg,l,oper,r,__is__,del,v_flag,p_flag,v);
    op->finish_sntnc(err_msg,l,oper,r,__is__,del,v_flag,p_flag);
    err_msg+=".";
    
    if(next)
    {
        std::string next_msg=next->print(v);
        if(!next_msg.empty())
            err_msg+=" "+next_msg;
    }
    
    return err_msg;
}
/*--------------------------------------------
 
 --------------------------------------------*/
std::string Logics::operator()(void* v)
{
    
    if(op==NULL)
        return std::string();
        
    if(!op->operator()(v))
    {
        return print(v);
    }
    
    if(next!=NULL)
        return next->operator()(v);
    
    return std::string();
}
/*--------------------------------------------
 
 --------------------------------------------*/
std::string Logics::operator()()
{
    return operator()(NULL);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Logics::link()
{
    
    if(next==NULL)
    {
        next=new Logics();
    }
    
    return *next;
}
/*----------------------
 _       _____   _____ 
| |     /  _  \ /  ___|
| |     | | | | | |    
| |     | | | | | |  _ 
| |___  | |_| | | |_| |
|_____| \_____/ \_____/
 ----------------------*/
Logics::Log::Log(std::string&& op_name_,const bool _is):
op_name(op_name_),
_is_(_is)
{

}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::Log::~Log()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
int Logics::Log::eq(const char* st0,const char* st1)
{
    if(st0==NULL && st1==NULL)
        return 1;
    if(st0==NULL)
        return 0;
    if(st1==NULL)
        return 0;
    if(strcmp(st0,st1))
        return 0;
    
    return 1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::Log::finish_sntnc(std::string& err_msg,std::string& L,std::string& op,std::string& R,bool& __is__,std::string& del,int v_flag,bool& p_flag)
{
    if(L.empty())
        return;
    err_msg+=L+" ";
    
    if(v_flag==VERB)
    {
        if(__is__ && p_flag)
            err_msg+=p_verb+" ";
        else if(__is__ && !p_flag)
            err_msg+=s_verb+" ";
    
        else if(!__is__ && p_flag)
            err_msg+=p_n_verb+" ";
        else
            err_msg+=s_n_verb+" ";
    }
    else
    {
        if(__is__)
            err_msg+=verb_aux+" "+verb+" ";
        else
            err_msg+=verb_aux+" "+negate+" "+verb+" ";
    }
    err_msg+=op;
    
    if(!R.empty())
        err_msg+=" "+R;

    
    L.clear();
    op.clear();
    R.clear();
    p_flag=false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::Log::start_sntnc(std::string& L,std::string& op,std::string& R,int v_flag,const std::string& l_name,const std::string& r_name)
{
    L=l_name;
    op=op_name;
    R=r_name;
}
/*--------------------------------------------
 3: (!L ,  O ,  R ) -> ( L0+L1 , O        , R     )
 4: ( L , !O , !R ) -> ( L     , O0+R0+O1 , R1    )
 5: ( L , !O ,  R ) -> ( L     , O0+O1    , R     )
 6: ( L ,  O , !R ) -> ( L     , O        , R0+R1 )
 
 2: (!L ,  O , !R ) ->   x
 1: (!L , !O ,  R ) ->   x
 0: (!L , !O , !R ) ->   x
 
 an example:
 
 #define mul(A,B) A.print(NULL).c_str(),B.print(NULL).c_str(),(A*B).print(NULL).c_str()
 int l0,l1,r0,r1;
 var<int>L0(l0,"L0");
 var<int>L1(l1,"L1");
 var<int>R0(r0,"R0");
 var<int>R1(r1,"R1");
 const char* op1="eq";
 const char* op0="!gt";
 const char* form="case %s\n\t%s\n+\t%s\n=\t%s\n"
 "--------------------------------------------------------------------------\n";
 const char* c1="1: (!L , !O ,  R ) ->   x";
 const char* c2="2: (!L ,  O , !R ) ->   x";
 const char* c3="3: (!L ,  O ,  R ) -> ( L0+L1 , O        , R     )";
 const char* c4="4: ( L , !O , !R ) -> ( L     , O0+R0+O1 , R1    )";
 const char* c5="5: ( L , !O ,  R ) -> ( L     , O0+O1    , R     )";
 const char* c6="6: ( L ,  O , !R ) -> ( L     , O        , R0+R1 )";
 printf(form,c4,mul(Logics(l0,op0,r0),Logics(l0,op1,r1)));
 printf(form,c2,mul(Logics(l0,op0,r0),Logics(l1,op0,r1)));
 printf(form,c1,mul(Logics(l0,op0,r0),Logics(l1,op1,r0)));
 printf(form,c3,mul(Logics(l0,op0,r0),Logics(l1,op0,r0)));
 printf(form,c5,mul(Logics(l0,op0,r0),Logics(l0,op1,r0)));
 printf(form,c6,mul(Logics(l0,op0,r0),Logics(l0,op0,r1)));
 
 const char* op2="!eq";
 printf(form,c2,mul(Logics(l0,op0,r0),Logics(l1,op2,r1)));
 printf(form,c3,mul(Logics(l0,op0,r0),Logics(l1,op2,r0)));
 printf(form,c6,mul(Logics(l0,op0,r0),Logics(l0,op2,r1)));
 --------------------------------------------*/
void Logics::Log::process_sntnc(std::string& err_msg,std::string& L,std::string& op,std::string& R,bool& __is__,std::string& del,int v_flag,bool& p_flag,const std::string& l_name,const std::string& r_name)
{
    if(L.empty())
    {
        start_sntnc(L,op,R,v_flag,l_name,r_name);
        __is__=_is_;
        p_flag=false;
        return;
    }
    
    int case_no=4*(L==l_name)+2*(__is__==_is_ && op==op_name)+(R==r_name);
    
    if(case_no<3 || case_no==7)
    {
        finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
        err_msg+=", "+del+" ";
        print(err_msg,L,op,R,__is__,del,v_flag,p_flag);
        return;
    }
 
    /*--------------------------------------------------
     3: (!L ,  O ,  R ) -> ( L0+L1 , O        , R     )
     --------------------------------------------------*/
    if(case_no==3)
    {
        p_flag=true;
        L+=" "+del+" "+l_name;
    }
    /*--------------------------------------------------
     4: ( L , !O , !R ) -> ( L     , O0+R0+O1 , R1    )
     4: ( L , !O , !R ) -> ( L     , O1+R1+O0 , R0    )
     --------------------------------------------------*/
    else if(case_no==4)
    {
        if(_is_!=__is__)
        {
            if(!__is__)
            {
                if(R.empty())
                    op=op_name+" "+del+" "+negate+" "+op;
                else
                    op=op_name+" "+r_name+" "+del+" "+negate+" "+op;
                __is__=true;
            }
            else
            {
                if(R.empty())
                    op+=" "+del+" "+negate+" "+op_name;
                else
                    op+=" "+R+" "+del+" "+negate+" "+op_name;
                R=r_name;
            }
        }
        else
        {
            if(R.empty())
                op+=" "+del+" "+op_name;
            else
                op+=" "+R+" "+del+" "+op_name;
            R=r_name;
        }
    }
    /*--------------------------------------------------
     5: ( L , !O ,  R ) -> ( L     , O0+O1    , R     )
     5: ( L , !O ,  R ) -> ( L     , O1+O0    , R     )
     --------------------------------------------------*/
    else if(case_no==5)
    {
        
        if(_is_!=__is__)
        {
            if(!__is__)
            {
                op=op_name+" "+del+" "+negate+" "+op;
                __is__=true;
            }
            else
                op+=" "+del+" "+negate+" "+op_name;
        }
        else
            op+=" "+del+" "+op_name;
    }
    /*--------------------------------------------------
     6: ( L ,  O , !R ) -> ( L     , O        , R0+R1 )
     --------------------------------------------------*/
    else if(case_no==6)
    {
        if(!r_name.empty())
            R+=" "+del+" "+r_name;
    }
}
/*---------------------------------------------------------------
 _____   _   __   _       ___   _____   __    __  _____   _____
|  _  \ | | |  \ | |     /   | |  _  \  \ \  / / /  _  \ |  _  \
| |_| | | | |   \| |    / /| | | |_| |   \ \/ /  | | | | | |_| |
|  _  { | | | |\   |   / / | | |  _  /    \  /   | | | | |  ___/
| |_| | | | | | \  |  / /  | | | | \ \    / /    | |_| | | |    
|_____/ |_| |_|  \_| /_/   |_| |_|  \_\  /_/     \_____/ |_|
 ---------------------------------------------------------------*/
Logics::LogBinaryOp::LogBinaryOp(Log* l,Log* r,const char* op_name,const bool _is_):
Log(op_name,_is_),left(l),right(r)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::LogBinaryOp::~LogBinaryOp()
{
    delete left;
    delete right;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogBinaryOp::print(std::string& err_msg,std::string& L,std::string& op,std::string& R,bool& __is__,std::string& del,int v_flag,bool& p_flag)
{
    left->print(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    del=op_name;
    right->print(err_msg,L,op,R,__is__,del,v_flag,p_flag);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogBinaryOp::print(std::string& err_msg,std::string& L,std::string& op,std::string& R,bool& __is__,std::string& del,int v_flag,bool& p_flag,void* v)
{
    left->print(err_msg,L,op,R,__is__,del,v_flag,p_flag,v);
    del=op_name;
    right->print(err_msg,L,op,R,__is__,del,v_flag,p_flag,v);
}
/*-------------------------
     ___   __   _   _____
    /   | |  \ | | |  _  \
   / /| | |   \| | | | | |
  / / | | | |\   | | | | |
 / /  | | | | \  | | |_| |
/_/   |_| |_|  \_| |_____/
 -------------------------*/
Logics::LogAND::LogAND(Log* l,Log* r):LogBinaryOp(l,r,"and",true)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Logics::LogAND::operator()()
{return ((*left)() && (*right)());}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Logics::LogAND::operator()(void* v)
{return ((*left)(v) && (*right)(v));}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::Log* Logics::LogAND::clone()
{
    Log* l=left->clone();
    Log* r=right->clone();
    return new LogAND(l,r);
}
/*---------------
 _____   _____
/  _  \ |  _  \ 
| | | | | |_| | 
| | | | |  _  / 
| |_| | | | \ \ 
\_____/ |_|  \_\
 ---------------*/
Logics::LogOR::LogOR(Log* l,Log* r):LogBinaryOp(l,r,"or",true)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Logics::LogOR::operator()()
{return ((*left)() || (*right)());}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Logics::LogOR::operator()(void* v)
{return ((*left)(v) || (*right)(v));}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::Log* Logics::LogOR::clone()
{
    Log* l=left->clone();
    Log* r=right->clone();
    return new LogOR(l,r);
}
/*-------------------------------
 _____   _____   __   _   _____ 
/  ___| /  _  \ |  \ | | |  _  \
| |     | | | | |   \| | | | | |
| |     | | | | | |\   | | | | |
| |___  | |_| | | | \  | | |_| |
\_____| \_____/ |_|  \_| |_____/
 -------------------------------*/
Logics::LogCond::LogCond(Log* l,Log* r,const char* op_name,const bool _is_):
Log(op_name,_is_),left(l),right(r)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::LogCond::~LogCond()
{
    delete left;
    delete right;
}
/*----------
 _   _____
| | |  ___|
| | | |__
| | |  __| 
| | | |    
|_| |_|    
 ----------*/
Logics::LogIF::LogIF(Log* l,Log* r):LogCond(l,r,"if",true)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Logics::LogIF::operator()()
{
    if(!(*left)())
        return true;
    if((*right)())
        return true;
    
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Logics::LogIF::operator()(void* v)
{
    
    if(!(*left)(v))
        return true;
    
    if((*right)(v))
        return true;
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::Log* Logics::LogIF::clone()
{
    Log* l=left->clone();
    Log* r=right->clone();
    return new LogIF(l,r);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogIF::print(std::string& err_msg,std::string& L,std::string& op,std::string& R,bool& __is__,std::string& del,int v_flag,bool& p_flag)
{
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    if(!err_msg.empty())
        err_msg+=". ";
    
    v_flag=AUX_VERB;
    right->print(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);

    err_msg+=", when ";
    
    v_flag=VERB;
    left->print(err_msg,L,op,R,__is__,del,v_flag,p_flag);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogIF::print(std::string& err_msg,std::string& L,std::string& op,std::string& R,bool& __is__,std::string& del,int v_flag,bool& p_flag,void* v)
{
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    if(!err_msg.empty())
        err_msg+=". ";
    
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    if(!err_msg.empty())
        err_msg+=". ";
    
    v_flag=AUX_VERB;
    right->print(err_msg,L,op,R,__is__,del,v_flag,p_flag,v);
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    
    err_msg+=", when ";
    
    v_flag=VERB;
    left->print(err_msg,L,op,R,__is__,del,v_flag,p_flag,v);
}
/*------------------
 _   _____   _____
| | |  ___| |  ___|
| | | |__   | |__  
| | |  __|  |  __| 
| | | |     | |    
|_| |_|     |_|
 ------------------*/
Logics::LogIFF::LogIFF(Log* l,Log* r):LogCond(l,r,"iff",true)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Logics::LogIFF::operator()()
{
    if((*left)())
    {
        if((*right)())
            return true;
    }
    else
        if(!(*right)())
            return true;

    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Logics::LogIFF::operator()(void* v)
{
    if((*left)(v))
    {
        if((*right)(v))
            return true;
    }
    else
        if(!(*right)(v))
            return true;
    
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::Log* Logics::LogIFF::clone()
{
    Log* l=left->clone();
    Log* r=right->clone();
    return new LogIFF(l,r);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogIFF::print(std::string& err_msg,std::string& L,std::string& op,std::string& R,bool& __is__,std::string& del,int v_flag,bool& p_flag)
{
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    if(!err_msg.empty())
        err_msg+=". ";

    v_flag=AUX_VERB;
    right->print(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    err_msg+=", if and only if ";

    v_flag=VERB;
    left->print(err_msg,L,op,R,__is__,del,v_flag,p_flag);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogIFF::print(std::string& err_msg,std::string& L,std::string& op,std::string& R,bool& __is__,std::string& del,int v_flag,bool& p_flag,void* v)
{
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    if(!err_msg.empty())
        err_msg+=". ";
    
    v_flag=AUX_VERB;
    right->print(err_msg,L,op,R,__is__,del,v_flag,p_flag,v);
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    err_msg+=", if and only if ";
    
    v_flag=VERB;
    left->print(err_msg,L,op,R,__is__,del,v_flag,p_flag,v);
}
