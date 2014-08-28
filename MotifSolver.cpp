#include "MotifSolver.hpp"

MotifSolver::MotifSolver(int dim,int pop) : DESolver(dim,pop), count(0)
{	
}
MotifSolver::~MotifSolver(void)
{
}
double MotifSolver::EnergyFunction(double *trial,bool &bAtSolution)
{
    double result=0;
    MapPars(trial);
    if(!ValidCandidate())
        return 10.0;
    result = -1.0 * EvaluateCandidate();
    if(GetCount() >= GetCountCutoff())
        AppendMotif();
    if (count++ % nPop == 0)
        printf("%d %E\n",count / nPop + 1,Energy());
    if((unsigned) result==GetNSeqRead()) {bAtSolution=true;}
    return(result);
}
void MotifSolver::MapPars( const double *pars)
{
    BSeq bs;
    //bs.a = ""; bs.b = "";
    double sA_star = (double)GetSeqAMissMatch()/GetSeqALength(); 
    double sA_O = sA_star+(1-sA_star)*.5; 
    double sA_T = sA_O+(1-sA_star)*.25; 
    double sB_star = (double)GetSeqBMissMatch()/GetSeqBLength(); 
    double sB_O = sA_star+(1-sB_star)*.5; 
    double sB_T = sA_O+(1-sB_star)*.25; 
    for(int i=0;i<Dimension();++i)
    {
        if(i<(int)GetSeqALength())
        {
            if(Double1IsLessThanDouble2(pars[i], sA_star, EPS))
                bs.a.append("*");
            else
                if(Double1IsGreaterThanDouble2(pars[i], sA_star, EPS) && 
                    Double1IsLessThanDouble2(pars[i], sA_O, EPS))
                    bs.a.append("O");
                else
                    if(Double1IsGreaterThanDouble2(pars[i], sA_O, EPS) && 
                        Double1IsLessThanDouble2(pars[i], sA_T, EPS))
                        bs.a.append("T");
                    else
                        if(Double1IsGreaterThanDouble2(pars[i], sA_T, EPS)) 
                            bs.a.append("C");
        }
        else
        {
            if(Double1IsLessThanDouble2(pars[i], sB_star, EPS))
                bs.b.append("*");
            else
                if(Double1IsGreaterThanDouble2(pars[i], sB_star, EPS) && 
                    Double1IsLessThanDouble2(pars[i], sB_O, EPS))
                    bs.b.append("O");
                else
                    if(Double1IsGreaterThanDouble2(pars[i], sB_O, EPS) && 
                        Double1IsLessThanDouble2(pars[i], sB_T, EPS))
                        bs.b.append("T");
                    else
                        if(Double1IsGreaterThanDouble2(pars[i], sB_T, EPS)) 
                            bs.b.append("C");
        }
    }
//    std::cout<<"Candidate is : "<<bs.a+"."+bs.b<<std::endl;
    SetCandidate(bs);
}
int MotifSolver::Optimize()
{
    std::fstream fout;
    double* min=new double[Dimension()];
    double* max=new double[Dimension()];
    int ii;

    for(ii=0;ii<Dimension();ii++)
    {
            max[ii] =  1.0;
            min[ii] =  0.0;
    }
    Setup(min,max,stRand1Bin,0.8,0.3);
//	pick from the list below
//stBest1Exp			0
//stRand1Exp			1
//stRandToBest1Exp		2
//stBest2Exp			3
//stRand2Exp			4
//stBest1Bin			5
//stRand1Bin			6
//stRandToBest1Bin		7
//stBest2Bin			8
//stRand2Bin	        	9
    printf("Calculating...\n\n");
    Solve(MAX_GENERATIONS);

    MapPars(Solution());
    std::cout << "Best Energy inner opt: "<< Energy() << std::endl<< std::endl;
    std::cout<< "\n\nBest Coefficients:\n";
    BSeq bs = GetCandidate();
    std::cout<<bs.a+"."+bs.b<<std::endl;;

    if (max) delete max;
    if (min) delete min;
    
    return 0;
}
