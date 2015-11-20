#ifndef FSCHWEFEL_2_6_H
#define FSCHWEFEL_2_6_H

#include "BenchmarkFunction.h"
class FSchwefel_2_6 : public BenchmarkFunction
{
    public:
		FSchwefel_2_6(ParamMap &v);
        FSchwefel_2_6(const int rId, const int rDim, string& rName);
        virtual ~FSchwefel_2_6();
    protected:
        void initialize();
        void loadData();
        void evaluate__(double const *x,vector<double>& obj);
    private:
          int ** mpp_a;
          double * mp_b;

};

#endif // FSCHWEFEL_2_6_H
