// newran.cpp -----------------------------------------------------------

// NEWRAN02B - 22 July 2002

#define WANT_MATH
#include "../Utility/include.h"

#include "newran.h"

#ifdef use_namespace
namespace NEWRAN {
#endif



//**************************** utilities ******************************
inline Real square(Real x) { return x*x; }
inline ExtReal square(const ExtReal& x) { return x*x; }

static void ErrorNoSpace() { throw myException("Newran: out of space@ ErrorNoSpace()"); }

//************************* end of definitions ************************


Real Random::Raw()                           // get new uniform random number
{
   // m = 2147483647 = 2^31 - 1; a = 16807;
   // 127773 = m div a; 2836 = m mod a
   long iseed = (long)seed;
   long hi = iseed / 127773L;                 // integer division
   long lo = iseed - hi * 127773L;            // modulo
   iseed = 16807 * lo - 2836 * hi;
   if (iseed <= 0) iseed += 2147483647L;
   seed = (double)iseed; return seed*4.656612875e-10;
}

Real Random::Density(Real) const
{ throw myException("density function not defined@Random::Density"); return 0.0; }

#ifdef _MSC_VER
static void DoNothing(int) {}
#endif

Real Random::Next()                          // get new mixed random number
{
   if (!seed)
      throw myException("Random number generator not initialised@Random::Next");
   int i = (int)(Raw()*128);               // 0 <= i < 128
#ifdef _MSC_VER
   DoNothing(i); DoNothing(i);
#endif
   Real f = Buffer[i]; Buffer[i] = Raw();  // Microsoft release gets this wrong
   return f;

   // return Mother(&iseed);
}

double Random::Get()                  // get random number seed
{ return seed/2147483648UL; }

Random::Random(double s)            // set random number seed
                                      // s must be between 0 and 1
{
   if (s>=1.0 || s<=0.0)
      throw myException("Newran: seed out of range@Random::Random");
   //iseed = 2147483648L * s;         // for Mother
   seed = (long)(s*2147483648UL);
   for (int i = 0; i<128; i++) Buffer[i] = Raw();
}


PosGen::PosGen(double s):Random(s)                             // Constructor
{
   #ifdef MONITOR
      cout << "constructing PosGen\n";
   #endif
   NotReady=true;
}

PosGen::~PosGen()
{
   if (!NotReady)
   {
      #ifdef MONITOR
	 cout << "freeing PosGen arrays\n";
      #endif
      delete [] sx; delete [] sfx;
	  sx=sfx=0;
   }
   #ifdef MONITOR
      cout << "destructing PosGen\n";
   #endif
}

void PosGen::Build(bool sym)                 // set up arrays
{
   #ifdef MONITOR
      cout << "building PosGen arrays\n";
   #endif
   int i;
   NotReady=false;
   sx=new Real[60]; sfx=new Real[60];
   if (!sx || !sfx) ErrorNoSpace();
   Real sxi=0.0; Real inc = sym ? 0.01 : 0.02;
   for (i=0; i<60; i++)
   {
      sx[i]=sxi; Real f1=Density(sxi); sfx[i]=f1;
      if (f1<=0.0) goto L20;
      sxi+=inc/f1;
   }
   throw myException("Newran: area too large@PosGen::Build");
L20:
   if (i<50) throw myException("Newran: area too small@PosGen::Build");
   xi = sym ? 2*i : i;
   return;
}

Real PosGen::Next()
{
   Real ak,y; int ir;
   if (NotReady) Build(false);
   do
   {
      Real r1=Random::Next();
      ir = (int)(r1*xi); Real sxi=sx[ir];
      ak=sxi+(sx[ir+1]-sxi)*Random::Next();
      y=sfx[ir]*Random::Next();
   }
   while ( y>=sfx[ir+1] && y>=Density(ak) );
   return ak;
}

Real SymGen::Next()
{
   Real s,ak,y; int ir;
   if (NotReady) Build(true);
   do
   {
      s=1.0;
      Real r1=Random::Next();
      if (r1>0.5) { s=-1.0; r1=1.0-r1; }
      ir = (int)(r1*xi); Real sxi=sx[ir];
      ak=sxi+(sx[ir+1]-sxi)*Random::Next();
      y=sfx[ir]*Random::Next();
   }
   while ( y>=sfx[ir+1] && y>=Density(ak) );
   return s*ak;
}

AsymGen::AsymGen(Real modex,double s):Random(s)                 // Constructor
{
   #ifdef MONITOR
      cout << "constructing AsymGen\n";
   #endif
   mode=modex; NotReady=true;
}

void AsymGen::Build()                        // set up arrays
{
   #ifdef MONITOR
      cout << "building AsymGen arrays\n";
   #endif
   int i;
   NotReady=false;
   sx=new Real[121]; sfx=new Real[121];
   if (!sx || !sfx)  ErrorNoSpace();
   Real sxi=mode;
   for (i=0; i<120; i++)
   {
      sx[i]=sxi; Real f1=Density(sxi); sfx[i]=f1;
      if (f1<=0.0) goto L20;
      sxi+=0.01/f1;
   }
   throw myException("Newran: area too large (a)@AsymGen::Build");
L20:
   ic=i-1; sx[120]=sxi; sfx[120]=0.0;
   sxi=mode;
   for (; i<120; i++)
   {
      sx[i]=sxi; Real f1=Density(sxi); sfx[i]=f1;
      if (f1<=0.0) goto L30;
      sxi-=0.01/f1;
   }
   throw myException("Newran: area too large (b)@AsymGen::Build");
L30:
   if (i<100)  throw myException("Newran: area too small@AsymGen::Build");
   xi=i;
   return;
}

Real AsymGen::Next()
{
   Real ak,y; int ir1;
   if (NotReady) Build();
   do
   {
      Real r1=Random::Next();
      int ir=(int)(r1*xi); Real sxi=sx[ir];
      ir1 = (ir==ic) ? 120 : ir+1;
      ak=sxi+(sx[ir1]-sxi)*Random::Next();
      y=sfx[ir]*Random::Next();
   }
   while ( y>=sfx[ir1] && y>=Density(ak) );
   return ak;
}

AsymGen::~AsymGen()
{
   if (!NotReady)
   {
      #ifdef MONITOR
	 cout << "freeing AsymGen arrays\n";
      #endif
      delete [] sx; delete [] sfx;
   }
   #ifdef MONITOR
      cout << "destructing AsymGen\n";
   #endif
}

Normal::Normal(double s):SymGen(s),count(0)
{
   if (count) { NotReady=false; xi=Nxi; sx=Nsx; sfx=Nsfx; }
   else { Build(true); Nxi=xi; Nsx=sx; Nsfx=sfx; }
   count++;
}

Normal::~Normal()
{
   count--;
   if (count) NotReady=true;                     // disable freeing arrays
}

Real Normal::Density(Real x) const               // normal density
{ return (fabs(x)>8.0) ? 0 : 0.398942280 * exp(-x*x / 2); }

Real  Normal::NextNonStand(Real rmean,Real rvariance){

	Real X = Next();
	Real stddev = sqrt( rvariance );

	return rmean + stddev * X;
}

Real Cauchy::Density(Real x) const               // Cauchy density function
{ return (fabs(x)>1.0e15) ? 0 : 0.31830988618 / (1.0+x*x); }

Real  Cauchy::NextNonStand(Real rmean,Real rvariance){

	Real X = Next();
	Real stddev = sqrt( rvariance );

	return rmean + stddev * X;
}

Real Exponential::Density(Real x) const          // Negative exponential
{ return  (x > 40.0 || x < 0.0) ? 0.0 : exp(-x); }


// Identification routines for each class - may not work on all compilers?

char* Random::Name()            { return (char*)"Random";           }
char* Uniform::Name()           { return (char*)"Uniform";          }
char* Constant::Name()          { return (char*)"Constant";         }
char* PosGen::Name()            { return (char*)"PosGen";           }
char* SymGen::Name()            { return (char*)"SymGen";           }
char* AsymGen::Name()           { return (char*)"AsymGen";          }
char* Normal::Name()            { return (char*)"Normal";           }
char* Cauchy::Name()            { return (char*)"Cauchy";           }
char* Exponential::Name()       { return (char*)"Exponential";      }
char* Levy::Name()		{ return (char*)"Levy";      }


Levy::Levy(Real c,double s):AsymGen(c/3,s){
	sc=c;
}
Real Levy::Density(Real x) const{
	if(x<=0.0) return 0.0;
	Real y;
	y=sqrt(0.5*sc/3.14159265358979323846)*exp(-0.5*sc/x)/pow(x,1.5);
	return (fabs(x)>1.0e15) ? 0 :y;
}



#ifdef use_namespace
}
#endif


