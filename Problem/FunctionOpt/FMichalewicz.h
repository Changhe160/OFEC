#ifndef FMICHALEWICZ_H
#define FMICHALEWICZ_H

#include "BenchmarkFunction.h"


class FMichalewicz : public BenchmarkFunction
{
    public:
		FMichalewicz(ParamMap &v);
        FMichalewicz(const int rId, const int rDimNumber, string& rName);
        virtual ~FMichalewicz();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
		int m_m;
    private:
};

#endif // FMICHALEWICZ_H