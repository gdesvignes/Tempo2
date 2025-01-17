#include <tempo2.h>

void t2Fit(pulsar *psr,unsigned int npsr, const char *covarFuncFile);

unsigned int t2Fit_getFitData(pulsar *psr, double* x, double* y,
        double* e, int* ip);

void t2Fit_fillGlobalFitInfo(pulsar* psr, unsigned int npsr,FitInfo &OUT);

void t2Fit_fillFitInfo(pulsar* psr, FitInfo &OUT);

void t2Fit_buildDesignMatrix(pulsar* psr,int ipsr,double x, int ipos, double* afunc);
void t2Fit_buildConstraintsMatrix(pulsar* psr,int ipsr, int iconstraint, double* afunc);

void t2Fit_updateParameters(pulsar *psr,int ipsr,double *val,double *error);
