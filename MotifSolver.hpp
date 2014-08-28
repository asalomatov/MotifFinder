#if !defined(_MOTIFSOLVER_H)
#define _MOTIFSOLVER_H

#include "DESolver.h"



class MotifSolver : public DESolver
{
public:
    MotifSolver(int dim,int pop);
    ~MotifSolver(void);
    void MapPars( const double *pars);
    double EnergyFunction(double *trial, bool &bAtSolution);
    int Optimize();
    void SetMaxGenerations(const int n){MAX_GENERATIONS=n;}
private:
    int MAX_GENERATIONS;
    int count;
};

#endif
