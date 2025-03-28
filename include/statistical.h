#ifndef _STATISTICAL_H
#define _STATISTICAL_H

#include "mp.h"

Real Boot_Mean(Real Boot[], int Nboot);
Real Boot_Mean(PrecVec Boot, int Nboot);
Real Boot_Sigma(Real Boot[], int Nboot);
Real Boot_Sigma(PrecVec Boot, int Nboot);

int indice(int ia, int ib, int ic, int id, int ie, int na, int nb, int  nc, int nd, int ne);
int indice(int ia, int ib, int ic, int na, int nb, int  nc);
int indice(int ia, int ib, int na, int nb);

#endif