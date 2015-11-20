// newran.h ------------------------------------------------------------

// NEWRAN02B - 22 July 2002

#ifndef NEWRAN_LIB
#define NEWRAN_LIB 0

//******************* utilities and definitions *************************

#include "../Utility/include.h"
#include "../Utility/myexcept.h"
#include "extreal.h"

#ifdef use_namespace
namespace NEWRAN { using namespace RBD_COMMON; }
namespace RBD_LIBRARIES { using namespace NEWRAN; }
namespace NEWRAN {
#endif

typedef Real (*PDF)(Real);                // probability density function

Real ln_gamma(Real);                      // log gamma function

//**************** uniform random number generator **********************



class Random                              // uniform random number generator
{
   double seed;                    // seed
   //unsigned long iseed;          // for Mother
   Real Buffer[128];               // for mixing random numbers
   Real Raw();                     // unmixed random numbers
   void operator=(const Random&) {}       // private so can't access
public:
   Random(double s);             // set seed (0 < seed < 1)
   double Get();                   // get seed
   virtual Real Next();                   // get new value
   virtual char* Name();                  // identification
   virtual Real Density(Real) const;      // used by PosGen & Asymgen
   Random() {}                            // do nothing
   virtual ~Random() {}                   // make destructors virtual
   virtual ExtReal Mean() const { return 0.5; }
                                          // mean of distribution
   virtual ExtReal Variance() const { return 1.0/12.0; }
					  // variance of distribution
   virtual void tDelete() {}              // delete components of sum
   virtual int nelems() const { return 1; }


};


//****************** uniform random number generator *********************

class Uniform : public Random
{
   void operator=(const Uniform&) {}      // private so can't access

public:

   char* Name();                          // identification
   Uniform(double s):Random(s){}                           // set value
   Real Next() {  
		return Random::Next(); 
   }
   ExtReal Mean() const { return 0.5; }
   ExtReal Variance() const { return 1.0/12.0; }
   Real Density(Real x) const { return (x < 0.0 || x > 1.0) ? 0 : 1.0; }
};


//************************* return constant ******************************

class Constant : public Random
{
   void operator=(const Constant&) {}     // private so can't access
   Real value;                            // value to be returned

public:
   char* Name();                          // identification
   Constant(Real v) { value=v; }          // set value
//   Real Next();
   Real Next() {  return value; }
   ExtReal Mean() const { return value; }
   ExtReal Variance() const { return 0.0; }
};

//**************** positive random number generator **********************

class PosGen : public Random              // generate positive rv
{
   void operator=(const PosGen&) {}       // private so can't access

protected:
   Real xi, *sx, *sfx;
   bool NotReady;
   void Build(bool);                      // called on first call to Next

public:
   char* Name();                          // identification
   PosGen(double s);           // constructor
   ~PosGen();                             // destructor
   Real Next();                           // to get a single new value
   ExtReal Mean() const { return (ExtReal)Missing; }
   ExtReal Variance() const { return (ExtReal)Missing; }
};

//**************** symmetric random number generator **********************

class SymGen : public PosGen              // generate symmetric rv
{
   void operator=(const SymGen&) {}       // private so can't access

public:
	SymGen(double s):PosGen(s){};
   char* Name();                          // identification
   Real Next();                           // to get a single new value
};

//**************** normal random number generator **********************

class Normal : public SymGen              // generate standard normal rv
{
   void operator=(const Normal&) {}       // private so can't access

   Real Nxi, *Nsx, *Nsfx;          // so we need initialise only once
   long count;                     // assume initialised to 0

public:
   char* Name();                          // identification
   Normal(double s);
   ~Normal();
   Real Density(Real) const;              // normal density function
   ExtReal Mean() const { return 0.0; }
   ExtReal Variance() const { return 1.0; }
   Real  NextNonStand(Real rmean,Real rvariance);
};

//**************** Cauchy random number generator **********************

class Cauchy : public SymGen              // generate standard cauchy rv
{
   void operator=(const Cauchy&) {}       // private so can't access

public:
	Cauchy(double s):SymGen(s){};
   char* Name();                          // identification
   Real Density(Real) const;              // Cauchy density function
   ExtReal Mean() const { return Indefinite; }
   ExtReal Variance() const { return PlusInfinity; }
   Real NextNonStand(Real rmean,Real rvariance);
};

//**************** Exponential random number generator **********************

class Exponential : public PosGen         // generate standard exponential rv
{
   void operator=(const Exponential&) {}  // private so can't access

public:
	Exponential(double s):PosGen(s){};
   char* Name();   // identification
   Real Density(Real) const;              // Exponential density function
   ExtReal Mean() const { return 1.0; }
   ExtReal Variance() const { return 1.0; }
};

//**************** asymmetric random number generator **********************

class AsymGen : public Random             // generate asymmetric rv
{
   void operator=(const AsymGen&) {}      // private so can't access
   Real xi, *sx, *sfx; int ic;
   bool NotReady;
   void Build();                          // called on first call to Next

protected:
   Real mode;

public:
   char* Name();                          // identification
   AsymGen(Real,double s);                         // constructor (Real=mode)
   ~AsymGen();                            // destructor
   Real Next();                           // to get a single new value
   ExtReal Mean() const { return (ExtReal)Missing; }
   ExtReal Variance() const { return (ExtReal)Missing; }
};

class  Levy : public AsymGen{
	void operator=(const Levy&) {}  // private so can't access

public:
   Real sc;
   Levy(Real,double s);
   char* Name();                    // identification
   Real Density(Real) const;              // Levy density function
   ExtReal Mean() const { return Indefinite; }
   ExtReal Variance() const { return Indefinite; }
};


#ifdef use_namespace
}
#endif

#endif

// body file: newran.cpp



