#define _USE_MATH_DEFINES
#define _SCL_SECURE_NO_WARNINGS
#include <ctime>
#include "Netlist.h"

#include "Resistor.h"
#include "Capacitor.h"
#include "Inductor.h"
#include "IndVoltageSource.h"
#include "IndCurrentSource.h"
#include "VCVS.h"
#include "VCCS.h"
//TODO: implement these components
//#include "MutualInductance.h"
//#include "OpAmp.h"
#include "VAC.h"
#include "VPulse.h"
#include "VStep.h"
#include "CNT.h"

#include <Eigen/SparseLU>
#include <Eigen/QR>
#include <Eigen/Eigen>
using Eigen::SparseLU;
using Eigen::HouseholderQR;
using Eigen::ColMajor;
using Eigen::COLAMDOrdering;

Netlist::Netlist() : GVals(0), XVals(0), CVals(0), BVals(0)
{
	nodes.push_back(0);
}

void Netlist::incrVoltages(int amount)
{
	numVoltages += amount;
}

void Netlist::incrCurrents(int amount)
{
	numCurrents += amount;
}

void Netlist::readNetlist(string inputFilename)
{
	// attempt to open file
	ifstream fileReader(inputFilename);
	if (fileReader.fail())
		Utilities::Error("problem opening the input file");

	vector<string> tokens;
	string nextLine;
	std::getline(fileReader, nextLine);
	while (!(nextLine == ".end")) {
		tokens = Utilities::parseLine(nextLine);

		if (tokens.at(0).at(0) == '.') {
			if (tokens.at(0).substr(1) == "freq") {
				solutions.push_back(solutionType::FREQ);
				std::istringstream(tokens.at(1)) >> nodeToTrack;
				nodeToTrack -= 1;
			}
			else if (tokens.at(0).substr(1) == "freqmor") {
				solutions.push_back(solutionType::FREQMOR);
				std::istringstream(tokens.at(1)) >> nodeToTrack;
				nodeToTrack -= 1;
			}
			else if (tokens.at(0).substr(1) == "timeeuler") {
				solutions.push_back(solutionType::TIMEEULER);
				std::istringstream(tokens.at(1)) >> nodeToTrack;
				nodeToTrack -= 1;
			}
			else if (tokens.at(0).substr(1) == "timetrapezoidal") {
				solutions.push_back(solutionType::TIMETRAPEZOIDAL);
				std::istringstream(tokens.at(1)) >> nodeToTrack;
				nodeToTrack -= 1;
			}
			else if (tokens.at(0).substr(1) == "dc") {
				solutions.push_back(solutionType::DC);
			}
			else {
				Utilities::Error("Unrecognized solution type.");
				exit(-1);
			}
		}
		else if (tokens.at(0).at(0) == '#') {
			//do nothing, it's a comment
		}
		else {
			int nodeOne = stoi(tokens.at(1));
			int nodeTwo = stoi(tokens.at(2));

			if (nodeOne > highestNode)
				highestNode = nodeOne;
			if (nodeTwo > highestNode)
				highestNode = nodeTwo;

			Component* newComponent = Utilities::parseComponent(tokens);
			circuitElements.push_back(newComponent);
		}

		std::getline(fileReader, nextLine);
	}
	fileReader.close();
	cout << "finished reading netlist" << endl;
	//correctDistributedNodes();
}

void Netlist::resizeMatrices()
{
	int dimension = numVoltages + numCurrents;
	cout << "Voltages: " << numVoltages << endl;
	cout << "Currents: " << numCurrents << endl;
	cout << "Dimension: " << dimension << endl;
	G.resize(dimension, dimension);
	G.setZero();
	X.resize(dimension, 1);
	X.setZero();
	C.resize(dimension, dimension);
	C.setZero();
	B.resize(dimension, 1);
	B.setZero();
}

void Netlist::calculateNewIndiciesNonDistributed()
{
	// this calculates the new indices for components that need to augment a matrix
	// calculation of the new indices for components that need to use a 
	// current equation had to wait until the number of voltage equations was known
	// newIndex is really numVoltages + 1 to start at new rows/columns augmented on matrices
	// then we subtract 1 to offset the fact that matrix indices start at 0
	// this is why we post increment the value newIndex

	for (Component* c : circuitElements) {
		if (Inductor* p = dynamic_cast<Inductor*> (c)) {
			p->setNewIndex(newIndexCounter++);
		}
		else if (IndVoltageSource* p = dynamic_cast<IndVoltageSource*> (c)) {
			p->setNewIndex(newIndexCounter++);
		}
		else if (VCVS* p = dynamic_cast<VCVS*> (c)) {
			p->setNewIndex(newIndexCounter++);
		}
		else if (VPulse* p = dynamic_cast<VPulse*> (c)) {
			p->setNewIndex(newIndexCounter++);
		}
		else if (VStep* p = dynamic_cast<VStep*> (c)) {
			p->setNewIndex(newIndexCounter++);
		}
		else if (VAC* p = dynamic_cast<VAC*> (c)) {
			p->setNewIndex(newIndexCounter++);
		}
	}
}

void Netlist::calculateNewIndiciesDistributed() {

	cout << "correcting indices of inductors from CNT netlist" << endl;

	for (Component* c : circuitElements) {
		if (CNT* p = dynamic_cast<CNT*> (c)) {
			for (Component* c1 : p->subComponents) {
				if (Inductor* p1 = dynamic_cast<Inductor*> (c1)) {
					p1->setNewIndex(newIndexCounter++);
				}
			}
		}
	}
}

void Netlist::correctDistributedNodes() {

	cout << "correcting distributed node values" << endl;

	for (Component* c : circuitElements) {
		if (CNT* p = dynamic_cast<CNT*> (c)) {
			p->setNodeOffset(highestNode - 1);
			highestNode = p->modifyNodeValues();
		}
	}
}

void Netlist::populateMatrices()
{
	for (Component* c : circuitElements) {
		incrVoltages(c->numVoltagesToAdd(nodes));
		incrCurrents(c->numCurrentsToAdd());
	}
	newIndexCounter = numVoltages;

	calculateNewIndiciesNonDistributed();
	//calculateNewIndiciesDistributed();
	resizeMatrices();
	for (Component* c : circuitElements) {
		c->insertStamp(GVals, XVals, CVals, BVals);
	}
	G.setFromTriplets(GVals.begin(), GVals.end());
	C.setFromTriplets(CVals.begin(), CVals.end());

	for (T t : XVals) {
		int row = t.row();
		int col = t.col();
		double val = t.value();
		X(row, col) = X(row, col) + val;
	}

	for (T t : BVals) {
		int row = t.row();
		int col = t.col();
		double val = t.value();
		B(row, col) = B(row, col) + val;
	}
}

void Netlist::simulate(string outputFilename)
{
	for (solutionType s : solutions) {
		switch (s) {
		case FREQ:
			solveFrequency(outputFilename);
			break;
		case FREQMOR:
			solveFrequencyMOR(outputFilename);
			break;
		case TIMEEULER:
			solveTimeBackwardEuler(outputFilename);
			break;
		case TIMETRAPEZOIDAL:
			solveTimeTrapezoidalRule(outputFilename);
			break;
		case DC:
			//TODO: implement DC analysis(transient at a single time point)
			break;
		}
	}
}

void Netlist::solveFrequency(string& outputFilename)
{
	int VACNewIndex = 0;
	double VACAmplitude = 0;
	double stepSize = 0;
	double currentFreq = 0;
	double numSteps = 0;
	double magnitude = 0;
	double phase = 0;

	//can only sweep one VAC Source at the moment
	for (Component* c : circuitElements) {
		if (VAC* p = dynamic_cast<VAC*> (c)) {
			VACNewIndex = p->newIndex;
			VACAmplitude = p->amplitude;
			stepSize = (p->maxFrequency - p->minFrequency) / p->numSteps;
			currentFreq = p->minFrequency;
			numSteps = p->numSteps;
			break;
		}
	}

	// attempt to open file
	ofstream fileWriter(outputFilename);
	if (fileWriter.fail())
		Utilities::Error("problem opening the output file");

	SparseMatrix<std::complex<double>, ColMajor> GPlusSC;
	double wSweep;
	std::complex<double> s;
	const int size = numVoltages + numCurrents;
	Matrix<std::complex<double>, Dynamic, 1> BNew(size, 1), XNew(size, 1);

	GPlusSC.setZero();
	BNew.setZero();
	XNew.setZero();

	SparseLU<SparseMatrix<std::complex<double>, ColMajor>, COLAMDOrdering<int>> solver;

	for (int i = 0; i < numSteps; i++, currentFreq += stepSize) {
		//remove dc components from B matrix
		BNew.setZero();
		//set AC component magnitude to amplitude
		BNew(VACNewIndex, 0) = VACAmplitude;

		//create static A matrix
		wSweep = 2 * M_PI * currentFreq;
		s.imag(wSweep);			
		GPlusSC = G.cast<std::complex<double>>() + (s * C.cast<std::complex<double>>());

		//set up solver
		GPlusSC.makeCompressed();
		solver.compute(GPlusSC);

		//solve system
		XNew = solver.solve(BNew);

		magnitude = calcMagnitude(nodeToTrack, XNew);
		phase = calcPhase(nodeToTrack, XNew);

		//print to file
		fileWriter << currentFreq << "\t" << magnitude << "\t" << phase << endl;
	}
	fileWriter.close();
}

void Netlist::solveFrequencyMOR(string& outputFilename)
{
	cout << "Using Model Order Reduction" << endl;

	const int size = numVoltages + numCurrents;
	int numberOfMoments = 100;
	int numberOfFreqHops = 1;
	int newSize = numberOfMoments * numberOfFreqHops;

	int VACNewIndex = 0;
	double VACAmplitude = 0;
	double stepSize = 0;
	double currentFreq = 0;
	double numSteps = 0;
	double magnitude = 0;
	double phase = 0;

	double maxFreq = 0;

	//can only sweep one VAC Source at the moment
	for (Component* c : circuitElements) {
		if (VAC* p = dynamic_cast<VAC*> (c)) {
			VACNewIndex = p->newIndex;
			VACAmplitude = p->amplitude;
			stepSize = (p->maxFrequency - p->minFrequency) / p->numSteps;
			currentFreq = p->minFrequency;
			numSteps = p->numSteps;
			maxFreq = p->maxFrequency;
			break;
		}
	}

	/*
	//setup B matrix to be normalized
	Matrix<double, Dynamic, 1> BNew(size, 1), initialMoment(size, 1);
	//remove dc components from B matrix
	BNew.setZero();
	initialMoment.setZero();

	//set AC component magnitude to amplitude
	BNew(VACNewIndex, 0) = VACAmplitude;

	// attempt to open file
	ofstream fileWriter(outputFilename);
	if (fileWriter.fail())
		Utilities::Error("problem opening the output file");

	cout << "Generating Orthogonal Q Vectors" << endl;

	SparseMatrix<double, ColMajor> negC = -C;
	Matrix<double, Dynamic, 1> curVec;
	
	//set up solver
	SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int>> solver;
	G.makeCompressed();
	solver.compute(G);

	//solve for initial moment
	initialMoment = solver.solve(BNew);

	SparseMatrix<double, ColMajor> qMatrix(size, numberOfMoments);
	Matrix<double, Dynamic, 1> curVector(size, 1), curVectorPrime(size, 1);
	//generate first Q vector
	qMatrix.col(0) = initialMoment.sparseView() / initialMoment.norm();
	//generate rest of Q vectors
	for (int vectorIndex = 1; vectorIndex < numberOfMoments; vectorIndex++) {
		//tempVector = A * prev Q vector
		curVector = solver.solve(negC * qMatrix.col(vectorIndex - 1));
		//calculate qVectorPrime
		curVectorPrime = curVector;
		for (int index = 0; index < vectorIndex; index++) {
			curVectorPrime -= ((curVector.dot(qMatrix.col(index).toDense())) * qMatrix.col(index));
		}
		//create qVector
		qMatrix.col(vectorIndex) = (curVectorPrime.sparseView() / curVectorPrime.norm());
	}

	cout << "Done Making Q Vectors" << endl;

	cout << "Generating Hat Matrices" << endl;
	SparseMatrix<double, ColMajor> GHat, CHat;
	Matrix<double, Dynamic, 1> BHat;
	GHat = (qMatrix.transpose() * G) * qMatrix;
	CHat = (qMatrix.transpose() * C) * qMatrix;
	BHat = qMatrix.transpose() * B.sparseView();

	cout << "Done making Hat Matrices" << endl;

	SparseLU<SparseMatrix<std::complex<double>, ColMajor>, COLAMDOrdering<int>> complexSolver;
	SparseMatrix<std::complex<double>, ColMajor> GPlusSC;
	double wSweep;
	std::complex<double> s;
	Matrix<std::complex<double>, Dynamic, 1> BModifiable(numberOfMoments, 1), Z(numberOfMoments, 1), extractSolution(size, 1);

	GPlusSC.setZero();
	BModifiable.setZero();
	Z.setZero();

	cout << "Solving" << endl;
	for (int i = 0; i < numSteps; i++, currentFreq += stepSize) {
		//reset BNew matrix
		BModifiable = BHat;

		//create static A matrix
		wSweep = 2 * M_PI * currentFreq;
		s.imag(wSweep);
		GPlusSC = GHat.cast<std::complex<double>>() + (s * CHat.cast<std::complex<double>>());

		//set up solver
		GPlusSC.makeCompressed();
		complexSolver.compute(GPlusSC);
		//solve system
		Z = complexSolver.solve(BModifiable);
		extractSolution = qMatrix * Z;
		magnitude = calcMagnitude(nodeToTrack, extractSolution);
		phase = calcPhase(nodeToTrack, extractSolution);

		//print to file
		fileWriter << currentFreq << "\t" << magnitude << "\t" << phase << endl;
	}
	cout << "Done Solving" << endl;
	fileWriter.close();
	*/

	// attempt to open file
	ofstream fileWriter(outputFilename);
	if (fileWriter.fail())
		Utilities::Error("problem opening the output file");

	cout << "Generating Orthogonal Q Vectors" << endl;

	SparseMatrix<double, ColMajor> negC = -C;
	Matrix<double, Dynamic, 1> curVec;

	//set up solver
	SparseLU<SparseMatrix<std::complex<double>, ColMajor>, COLAMDOrdering<int>> centerSolver;
	double  freqStep = maxFreq / numberOfFreqHops;
	int curCenter = 0;
	SparseMatrix<std::complex<double>, ColMajor> GPlusSCCenter(size, size);

	double wCenter;
	std::complex<double> sCenter;

	SparseMatrix < std::complex<double>, ColMajor > qMatrix(size, newSize);

	//setup B matrix to be normalized
	Matrix < std::complex<double>, Dynamic, 1 > BNew(size, 1), initialMoment(size, 1);

	for (int curStep = 0; curStep < numberOfFreqHops; curStep++) {
		//remove dc components from B matrix
		BNew.setZero();
		initialMoment.setZero();

		//set AC component magnitude to amplitude
		BNew(VACNewIndex, 0) = VACAmplitude;

		if (curStep == 0) {
			wCenter = 0;
		}
		else {
			wCenter = 2 * M_PI * curStep * freqStep;
		}

		sCenter.imag(-wCenter);
		GPlusSCCenter = G.cast<std::complex<double>>() + (sCenter * C.cast<std::complex<double>>());
		GPlusSCCenter.makeCompressed();
		centerSolver.compute(GPlusSCCenter);

		int indexZero = curStep * numberOfMoments;

		//solve for initial moment
		initialMoment = centerSolver.solve(BNew);

		Matrix<std::complex<double>, Dynamic, 1> curVector(size, 1), curVectorPrime(size, 1);
		//generate first Q vector
		qMatrix.col(indexZero) = initialMoment.sparseView() / initialMoment.norm();
		//generate rest of Q vectors
		for (int vectorIndex = indexZero + 1; vectorIndex < indexZero + numberOfMoments; vectorIndex++) {
			//tempVector = A * prev Q vector
			curVector = centerSolver.solve(negC.cast<std::complex<double>>() * qMatrix.col(vectorIndex - 1));
			//calculate qVectorPrime
			curVectorPrime = curVector;
			for (int index = indexZero; index < vectorIndex; index++) {
				curVectorPrime -= ((curVector.dot(qMatrix.col(index).toDense())) * qMatrix.col(index));
			}
			//create qVector
			qMatrix.col(vectorIndex) = (curVectorPrime.sparseView() / curVectorPrime.norm());
		}
	}

	cout << "Done Making Q Vectors" << endl;

	cout << "Generating Hat Matrices" << endl;
	SparseMatrix<std::complex<double>, ColMajor> GHat(newSize, newSize), CHat(newSize, newSize);
	Matrix<std::complex<double>, Dynamic, 1> BHat(newSize, 1);
	GHat = (qMatrix.transpose() * G.cast<std::complex<double>>()) * qMatrix;
	CHat = (qMatrix.transpose() * C.cast<std::complex<double>>()) * qMatrix;
	BHat = qMatrix.transpose() * B.cast<std::complex<double>>().sparseView();

	cout << "Done making Hat Matrices" << endl;

	SparseLU<SparseMatrix<std::complex<double>, ColMajor>, COLAMDOrdering<int>> complexSolver;
	SparseMatrix<std::complex<double>, ColMajor> GPlusSC;
	double wSweep;
	std::complex<double> s;
	Matrix<std::complex<double>, Dynamic, 1> BModifiable(newSize, 1), Z(newSize, 1), extractSolution(size, 1);

	GPlusSC.setZero();
	BModifiable.setZero();
	Z.setZero();

	cout << "Solving" << endl;
	for (int i = 0; i < numSteps; i++, currentFreq += stepSize) {
		//reset BNew matrix
		BModifiable = BHat;

		//create static A matrix
		wSweep = 2 * M_PI * currentFreq;
		s.imag(wSweep);
		GPlusSC = GHat.cast<std::complex<double>>() + (s * CHat.cast<std::complex<double>>());

		//set up solver
		GPlusSC.makeCompressed();
		complexSolver.compute(GPlusSC);
		//solve system
		Z = complexSolver.solve(BModifiable);
		extractSolution = qMatrix * Z;
		magnitude = calcMagnitude(nodeToTrack, extractSolution);
		phase = calcPhase(nodeToTrack, extractSolution);

		//print to file
		fileWriter << currentFreq << "\t" << magnitude << "\t" << phase << endl;
	}
	cout << "Done Solving" << endl;
	fileWriter.close();
}

void Netlist::solveTimeBackwardEuler(string& outputFilename)
{
	int vNewIndex = 0;
	double stepSize = 0;
	double currentTime = 0;
	double numSteps = 0;
	double amplitude = 0;
	double riseTime = 0;
	double pulseWidth = 0;
	double magnitude = 0;
	double vPulseValue = 0;

	//can only sweep for on VPulse/VStep source at the moment
	for (Component* c : circuitElements) {
		if (VPulse* p = dynamic_cast<VPulse*> (c)) {
			vNewIndex = p->newIndex;
			stepSize = (p->maxTime - p->minTime) / p->numSteps;
			currentTime = p->minTime;
			numSteps = p->numSteps;
			amplitude = p->amplitude;
			riseTime = p->riseTime;
			pulseWidth = p->pulseWidth;
		}
		else if (VStep* p = dynamic_cast<VStep*> (c)) {
			vNewIndex = p->newIndex;
			stepSize = (p->maxTime - p->minTime) / p->numSteps;
			currentTime = p->minTime;
			numSteps = p->numSteps;
			amplitude = p->amplitude;
			riseTime = p->riseTime;
			pulseWidth = p->maxTime;
		}
	}

	// attempt to open file
	ofstream fileWriter(outputFilename);
	if (fileWriter.fail())
		Utilities::Error("problem opening the output file");
	
	SparseMatrix<double, ColMajor> COverH(C), AStatic(G);

	COverH /= stepSize;
	AStatic += COverH;

	const int size = numVoltages + numCurrents;
	Matrix<double, Dynamic, 1> BCurrentTimePoint(size, 1), XPreviousTimePoint(size, 1), XCurrentTimePoint(size, 1);

	BCurrentTimePoint.setZero();
	XCurrentTimePoint.setZero();
	XPreviousTimePoint.setZero();

	SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int>> solver;

	AStatic.makeCompressed();
	solver.compute(AStatic);

	for (int i = 0; i < numSteps; i++, currentTime += stepSize) {
		BCurrentTimePoint.setZero();

		//calculate VPulseValue
		if (0 <= currentTime && currentTime <= riseTime) {
			vPulseValue = ((amplitude / riseTime) * currentTime);
		}
		else if (riseTime <= currentTime && currentTime <= riseTime + pulseWidth) {
			vPulseValue = amplitude;
		}
		else if (riseTime + pulseWidth <= currentTime && currentTime <= 2 * riseTime + pulseWidth) {
			vPulseValue = (amplitude - ((amplitude/riseTime) * (currentTime - (riseTime + pulseWidth))));
		}
		else {
			vPulseValue = 0;
		}
		
		BCurrentTimePoint(vNewIndex, 0) = vPulseValue;
		
		BCurrentTimePoint = COverH * XPreviousTimePoint + BCurrentTimePoint;

		XCurrentTimePoint = solver.solve(BCurrentTimePoint);

		magnitude = XCurrentTimePoint(nodeToTrack, 0);

		fileWriter << currentTime << "\t" << magnitude << endl;

		XPreviousTimePoint = XCurrentTimePoint;
	}

	fileWriter.close();
}

void Netlist::solveTimeTrapezoidalRule(string& outputFilename)
{
	int vNewIndex = 0;
	double stepSize = 0;
	double currentTime = 0;
	double numSteps = 0;
	double amplitude = 0;
	double riseTime = 0;
	double pulseWidth = 0;
	double magnitude = 0;
	double vPulseValue = 0;

	//can only sweep for on VPulse/VStep source at the moment
	for (Component* c : circuitElements) {
		if (VPulse* p = dynamic_cast<VPulse*> (c)) {
			vNewIndex = p->newIndex;
			stepSize = (p->maxTime - p->minTime) / p->numSteps;
			currentTime = p->minTime;
			numSteps = p->numSteps;
			amplitude = p->amplitude;
			riseTime = p->riseTime;
			pulseWidth = p->pulseWidth;
		}
		else if (VStep* p = dynamic_cast<VStep*> (c)) {
			vNewIndex = p->newIndex;
			stepSize = (p->maxTime - p->minTime) / p->numSteps;
			currentTime = p->minTime;
			numSteps = p->numSteps;
			amplitude = p->amplitude;
			riseTime = p->riseTime;
			pulseWidth = p->maxTime;
		}
	}

	// attempt to open file
	ofstream fileWriter(outputFilename);
	if (fileWriter.fail())
		Utilities::Error("problem opening the output file");

	SparseMatrix<double, ColMajor> COverH(C), GOverTwo(G), AStatic, BStatic;

	COverH /= stepSize;
	GOverTwo /= 2;

	AStatic = COverH + GOverTwo;
	BStatic = COverH - GOverTwo;

	const int size = numVoltages + numCurrents;
	Matrix<double, Dynamic, 1> BPreviousTimePoint(size, 1), BCurrentTimePoint(size, 1), XPreviousTimePoint(size, 1), XCurrentTimePoint(size, 1), BDynamic(size, 1);

	BPreviousTimePoint.setZero();
	BCurrentTimePoint.setZero();
	XPreviousTimePoint.setZero();
	XCurrentTimePoint.setZero();

	SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int>> solver;

	AStatic.makeCompressed();
	solver.compute(AStatic);

	for (int i = 0; i < numSteps; i++, currentTime += stepSize) {
		BCurrentTimePoint.setZero();

		//calculate VPulseValue
		if (0 <= currentTime && currentTime <= riseTime) {
			vPulseValue = ((amplitude / riseTime) * currentTime);
		}
		else if (riseTime <= currentTime && currentTime <= riseTime + pulseWidth) {
			vPulseValue = amplitude;
		}
		else if (riseTime + pulseWidth <= currentTime && currentTime <= 2 * riseTime + pulseWidth) {
			vPulseValue = (amplitude - ((amplitude / riseTime) * (currentTime - (riseTime + pulseWidth))));
		}
		else {
			vPulseValue = 0;
		}
		//set BCurrentTimePoint with the correct voltage value for the pulse input
		BCurrentTimePoint(vNewIndex, 0) = vPulseValue;
		//BDynamic = (B(tk)+B(tk-1))/2 + ((C/h)-(G/2))*X(tk-1)
		BDynamic = ((BCurrentTimePoint + BPreviousTimePoint) / 2) + (BStatic * XPreviousTimePoint);

		XCurrentTimePoint = solver.solve(BDynamic);

		magnitude = XCurrentTimePoint(nodeToTrack, 0);

		fileWriter << currentTime << "\t" << magnitude << endl;

		XPreviousTimePoint = XCurrentTimePoint;
		BPreviousTimePoint = BCurrentTimePoint;
	}

	fileWriter.close();
}

double Netlist::calcMagnitude(int row, Matrix<std::complex<double>, Dynamic, 1>& matrix)
{
	double real = matrix(row, 0).real();
	double imag = matrix(row, 0).imag();
	return sqrt((real * real) + (imag * imag));
}

double Netlist::calcPhase(int row, Matrix<std::complex<double>, Dynamic, 1>& matrix)
{
	double real = matrix(row, 0).real();
	double imag = matrix(row, 0).imag();
	return atan(imag/real);
}

void Netlist::prettyPrintNetlist()
{
	//TODO: add toString functions to each Component derived class
}

void Netlist::prettyPrintMatrices()
{
	//TODO: print the matrices to a excel spreadsheet or something
	cout << "\n\n\nG:\n" << G << endl;
	cout << "\n\n\nX:\n" << X << endl;
	cout << "\n\n\nC:\n" << C << endl;
	cout << "\n\n\nB:\n" << B << endl;
}