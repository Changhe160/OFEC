#ifndef FBRANINRCOS_H
#define FBRANINRCOS_H

#include "BenchmarkFunction.h"


class FBraninRCOS : public BenchmarkFunction
{
    public:
		FBraninRCOS(ParamMap &v);
        FBraninRCOS(const int rId, const int rDimNumber, string& rName);
        virtual ~FBraninRCOS();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};

#endif // FBRANINRCOS_H