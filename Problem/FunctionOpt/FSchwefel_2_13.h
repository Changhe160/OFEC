#ifndef FSCHWEFEL_2_13_H
#define FSCHWEFEL_2_13_H

#include "BenchmarkFunction.h"

class FSchwefel_2_13 : public BenchmarkFunction
{
    public:
		FSchwefel_2_13(ParamMap &v);
        FSchwefel_2_13(const int rId, const int rDimNumber, string& rName);
        virtual ~FSchwefel_2_13();
    protected:
        void initialize();
        void loadData();
        void evaluate__(double const *x,vector<double>& obj);
    private:
          int ** mpp_a;
          int ** mpp_b;
          double * mp_alpha;
};

#endif // FSCHWEFEL_2_13_H
