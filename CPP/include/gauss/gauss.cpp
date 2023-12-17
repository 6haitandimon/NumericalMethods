#include "gauss.h"

// void Zeros(std::vector<double> &result, int &size){
//     result.clear();
//     for(int index = 0; index < size; index++)
//         result.push_back(0);
// }

// void MatrixMultiplication(std::vector<std::vector<double> > A, std::vector<double> x, std::vector<double> &matrixMultiplication){
//     int matrixSize = A.size();
//     Zeros(matrixMultiplication, matrixSize);
//     for(int index = 0; index < matrixSize; index++){
//         for(int jIndex = 0; jIndex < matrixSize; jIndex++)
//             matrixMultiplication[index] += A[index][jIndex] * x[jIndex];
//     }
//     return;
// }


bool GaussElimination(double *An, double **X, int x)
{
	for (int k = 0; k < x; k++)
	{
		double max = fabs(X[k][k]);
		int remeber = k;		
		for (int i = k + 1; i < x; i++)
		{
			if (max < fabs(X[i][k]))
			{
				max = fabs(X[i][k]);
				remeber = i;
			}
		}

		if (fabs(max - 0) < 1e-6)
		{
			return 0;
		}

		if (k != remeber)				
		{
			double *temp = X[k];
			X[k] = X[remeber];
			X[remeber] = temp;
		}


		double lead = X[k][k];			
		for (int r = k; r < x + 1; r++)
		{
			X[k][r] /= lead;
		}
		
		for (int i = k + 1; i < x; i++)
		{
			double temp = X[i][k];
			for (int j = k; j < x + 1; j++)
			{
				X[i][j] -= X[k][j] * temp;
			}
		}
		
	}

	An[x - 1] = X[x - 1][x + 1 - 1];				
	for (int i = x - 2; i >= 0; i--)
	{
		An[i] = X[i][x + 1 - 1];
		for (int j = i + 1; j < x + 1 - 1; j++)
		{
			An[i] -= X[i][j] * An[j];
		}
	}
	return 1;
}