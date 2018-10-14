#include <cmath>
#include <float.h>
#include <iostream>
#include "MultiplicativePRNG.h"


using namespace std;


double calcSingleDefiniteIntegral(double lowerBound, double upperBound, int nodeNum, double (*calcFunction)(double),
                                  PRNG* prng);
double calcSingleIntegralFunction(double x);
double* calcResultsSingleIntegral(double lowerBound, double upperBound, double (*calcFunction)(double), PRNG* prng,
                                  bool isNeedResetPRNG, int calcsNumber, int nodeNumStep, int* nodeNums);

double calcDoubleDefiniteIntegral(double lowerBoundX, double upperBoundX, double (*calcLowerBoundY)(double),
                                  double (*calcUpperBoundY)(double), double areaSquare,
                                  double (*calcFunction)(double, double), int num, PRNG* prng);
double calcDoubleIntegralFunction(double x, double y);
double calcDoubleIntegralLowerBoundY(double x);
double calcDoubleIntegralUpperBoundY(double x);
double* calcResultsDoubleIntegral(double lowerBoundX, double upperBoundX, double (*calcLowerBoundY)(double),
                                  double (*calcUpperBoundY)(double), double areaSquare,
                                  double (*calcFunction)(double, double), PRNG* prng,
                                  bool isNeedResetPRNG, int calcsNumber, int nodeNumStep, int* nodeNums);

double* calcIntegralResiduals(double exactSolution, double *results, int num);



void calcSolutions()
{
    const double singleIntegralLowerBound = -1.5;
    const double singleIntegralUpperBound = +1.5;

    const double doubleIntegralLowerBoundX = -1.0;
    const double doubleIntegralUpperBoundX = 1.0;

    const int calcsNum = 20;
    const int nodeNumStep = 10;

    const double exactSolutionSingleIntegral = 1.85544;
    const double exactSolutionDoubleIntegral = 6.28319;

    const long long module = 2LL << 30;
    const long long seed = 79507;
    const int multiplier = 79507;

    PRNG* prng = new MultiplicativePRNG(module, seed, multiplier);

    auto nodeNumsSingleIntegral = new int[calcsNum];
    auto resultsSingleIntegral = calcResultsSingleIntegral(singleIntegralLowerBound, singleIntegralUpperBound,
            calcSingleIntegralFunction, prng, true, calcsNum, nodeNumStep, nodeNumsSingleIntegral);
    auto residualsSingleIntegral = calcIntegralResiduals(exactSolutionSingleIntegral, resultsSingleIntegral, calcsNum);

    auto nodeNumsDoubleIntegral = new int[calcsNum];
    auto resultsDoubleIntegral = calcResultsDoubleIntegral(doubleIntegralLowerBoundX, doubleIntegralUpperBoundX,
                                                           calcDoubleIntegralLowerBoundY, calcDoubleIntegralUpperBoundY,
                                                           M_PI, calcDoubleIntegralFunction, prng, true, calcsNum,
                                                           nodeNumStep, nodeNumsDoubleIntegral);
    auto residualsDoubleIntegral = calcIntegralResiduals(exactSolutionDoubleIntegral, resultsDoubleIntegral, calcsNum);

    FILE* outputFile;
    fopen_s(&outputFile, "output.txt", "w");

    fprintf(outputFile, "Single integral:\n%s;%s\n", "nodeNum", "residual");

    for (int i = 0; i < calcsNum; i++)
    {
//        printf("Nodes num: %d\n", nodeNumsSingleIntegral[i]);
//        printf("Result: %lf\n", resultsSingleIntegral[i]);
//        printf("Residual: %lf\n", residualsSingleIntegral[i]);
        fprintf(outputFile, "%d;%lf\n", nodeNumsSingleIntegral[i], residualsSingleIntegral[i]);
    }

    fprintf(outputFile, "\nDouble integral:\n%s;%s\n", "nodeNum", "residual");

    for (int i = 0; i < calcsNum; i++)
    {
//        printf("Nodes num: %d\n", nodeNumsDoubleIntegral[i]);
//        printf("Result: %lf\n", resultsDoubleIntegral[i]);
//        printf("Residual: %lf\n", residualsDoubleIntegral[i]);
        fprintf(outputFile, "%d;%lf\n", nodeNumsDoubleIntegral[i], residualsDoubleIntegral[i]);
    }

    fclose(outputFile);

    delete[] nodeNumsSingleIntegral;
    delete[] resultsSingleIntegral;
    delete[] residualsSingleIntegral;
    delete[] nodeNumsDoubleIntegral;
    delete[] resultsDoubleIntegral;
    delete[] residualsDoubleIntegral;
    delete prng;
}

double calcSingleDefiniteIntegral(double lowerBound, double upperBound, int nodeNum, double (*calcFunction)(double),
                                  PRNG* prng)
{
    double sum = 0;
    double len = (upperBound - lowerBound);

    for (int i = 0; i < nodeNum; i++)
    {
        sum += calcFunction(lowerBound + len * prng->next());
    }

    return len * sum / nodeNum;
}

double calcSingleIntegralFunction(double x)
{
    return exp(-pow(x, 6.0));
}

double* calcResultsSingleIntegral(double lowerBound, double upperBound, double (*calcFunction)(double), PRNG* prng,
                                  bool isNeedResetPRNG, int calcsNumber, int nodeNumStep, int* nodeNums)
{
    int nodeNum = 0;
    auto results = new double[calcsNumber];

    for (int i = 0; i < calcsNumber; i++)
    {
        nodeNum += nodeNumStep;
        nodeNums[i] = nodeNum;
        results[i] = calcSingleDefiniteIntegral(lowerBound, upperBound, nodeNum, calcFunction, prng);
        if (isNeedResetPRNG)
        {
            prng->reset();
        }
    }

    return results;
}

double* calcIntegralResiduals(double exactSolution, double *results, int num)
{
    auto residuals = new double[num];

    for (int i = 0; i < num; i++)
    {
        residuals[i] = abs(exactSolution - results[i]);
    }

    return residuals;
}

double calcDoubleDefiniteIntegral(double lowerBoundX, double upperBoundX, double (*calcLowerBoundY)(double),
                                  double (*calcUpperBoundY)(double), double areaSquare,
                                  double (*calcFunction)(double, double), int num, PRNG* prng)
{
    double result = 0;
    double prngX = 0;
    double prngY = 0;

    for (int i = 0; i < num; i++)
    {
        prngX = lowerBoundX + (upperBoundX - lowerBoundX) * prng->next();
        prngY = calcLowerBoundY(prngX) + (calcUpperBoundY(prngX) - calcLowerBoundY(prngX)) * prng->next();
        result += calcFunction(prngX, prngY);
    }

    return areaSquare / num * result;
}

double calcDoubleIntegralFunction(double x, double y)
{
    return 1 / sqrt(x * x + y * y);
}

double calcDoubleIntegralLowerBoundY(double x)
{
    return -sqrt(1 - x * x);
}

double calcDoubleIntegralUpperBoundY(double x)
{
    return sqrt(1 - x * x);
}

double* calcResultsDoubleIntegral(double lowerBoundX, double upperBoundX, double (*calcLowerBoundY)(double),
                                  double (*calcUpperBoundY)(double), double areaSquare,
                                  double (*calcFunction)(double, double), PRNG* prng,
                                  bool isNeedResetPRNG, int calcsNumber, int nodeNumStep, int* nodeNums)
{
    int nodeNum = 0;
    auto results = new double[calcsNumber];

    for (int i = 0; i < calcsNumber; i++)
    {
        nodeNum += nodeNumStep;
        nodeNums[i] = nodeNum;
        results[i] = calcDoubleDefiniteIntegral(lowerBoundX, upperBoundX, calcLowerBoundY, calcUpperBoundY, areaSquare,
                                                calcFunction, nodeNum, prng);
        if (isNeedResetPRNG)
        {
            prng->reset();
        }
    }

    return results;
}
