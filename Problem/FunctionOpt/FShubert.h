#ifndef FSHUBERT_H
#define FSHUBERT_H

#include "BenchmarkFunction.h"


class FShubert : public BenchmarkFunction
{
    public:
		FShubert(ParamMap &v);
        FShubert(const int rId, const int rDimNumber, string& rName);
        virtual ~FShubert();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};

#endif // FSHUBERT_H