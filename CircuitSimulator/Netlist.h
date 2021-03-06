#ifndef NETLIST_H
#define NETLIST_H

#include <vector>
using std::vector;

#include <string>
using std::string;

#include "Component.h"
#include "Utilities.h"

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <Eigen/Sparse>
#include <Eigen/src/Core/IO.h>

using Eigen::SparseMatrix;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::Triplet;
using Eigen::ColMajor;

enum solutionType {FREQ, FREQMOR, DC, TIMEEULER, TIMETRAPEZOIDAL};

class Netlist {
public:
	Netlist();
	void readNetlist(string inputFilename);
	void populateMatrices();
	void simulate(string outputFilename);
	vector<int> nodes;
	void prettyPrintNetlist();
	void prettyPrintMatrices();

private:
	vector<Component*> circuitElements;
	vector<solutionType> solutions;

	int numVoltages = 0;
	int numCurrents = 0;
	int nodeToTrack = 0;
	int newIndexCounter;
	int highestNode = 0;

	SparseMatrix<double, ColMajor> G;
	Matrix<double, Dynamic, 1> X;
	SparseMatrix<double, ColMajor> C;
	Matrix<double, Dynamic, 1> B;

	vector<T> GVals;
	vector<T> XVals;
	vector<T> CVals;
	vector<T> BVals;

	void incrVoltages(int amount);
	void incrCurrents(int amount);
	void resizeMatrices();
	void calculateNewIndiciesNonDistributed();
	void calculateNewIndiciesDistributed();
	void correctDistributedNodes();
	void solveFrequency(string& filename);
	void solveFrequencyMOR(string& filename);
	void solveTimeBackwardEuler(string& filename);
	void solveTimeTrapezoidalRule(string& filename);
	double calcMagnitude(int row, Matrix<std::complex<double>, Dynamic, 1>& matrix);
	double calcPhase(int row, Matrix<std::complex<double>, Dynamic, 1>& matrix);
};
#endif