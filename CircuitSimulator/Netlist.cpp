#include "Netlist.h"
#include <math.h>

#include "Resistor.h"
#include "Capacitor.h"
#include "Inductor.h"
#include "IndVoltageSource.h"
#include "IndCurrentSource.h"
#include "VCVS.h"
#include "VCCS.h"
//#include "MutualInductance.h"
//#include "OpAmp.h"
//#include "VAC.h"
#include "VPulse.h"
#include "VStep.h"

#include <Eigen/SparseLU>
#include <Eigen/Eigen>
using Eigen::SparseLU;
using Eigen::PartialPivLU;
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

void Netlist::readNetlist(string inputFilename)
{
	// attempt to open file
	ifstream fileReader(inputFilename);
	if (fileReader.fail())
		Error("problem opening the input file");

	string nextLine;
	std::getline(fileReader, nextLine);
	while (!(nextLine == ".end")) {
		parseLine(nextLine);
		std::getline(fileReader, nextLine);
	}
	fileReader.close();
	populateMatrices();
}

void Netlist::parseLine(string& nextLine)
{
	std::transform(nextLine.begin(), nextLine.end(), nextLine.begin(), ::tolower);

	string token;
	vector<string> tokens;
	while (token != nextLine) {
		token = nextLine.substr(0, nextLine.find_first_of(" "));
		nextLine = nextLine.substr(nextLine.find_first_of(" ") + 1);
		tokens.push_back(token);
	}

	if (tokens.at(0).at(0) == '.') {
		if (tokens.at(0).substr(1) == "freq") {
			solutions.push_back(solutionType::FREQ);
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
			Error("Unrecognized solution type.");
			exit(-1);
		}
	}
	else if (tokens.at(0).at(0) == '#') {
		//do nothing, it's a comment
	}
	else {
		circuitElements.push_back(parseComponent(tokens));
		//parseComponent(tokens);
	}
}

/*
Letter codes for different circuit elements

Resistor R
Capacitor C
Inductor L
Independent voltage source & stimulus V
Independent current source & stimulus I
Voltage-controlled voltage source E
Voltage-controlled current source G
Mutual Inductance/Coupling K
*/
Component* Netlist::parseComponent(vector<string>& tokens)
{	
	Component* newComponent = nullptr;

	int nodeOne;
	int nodeTwo;
	nodeOne = std::stoi(tokens.at(1));
	nodeTwo = std::stoi(tokens.at(2));

	switch (tokens.at(0).at(0)) {
	case 'r':
		newComponent = new Resistor(tokens.at(0), nodeOne, nodeTwo, convert(tokens.at(3)));
		break;
	case 'c':
		newComponent = new Capacitor(tokens.at(0), nodeOne, nodeTwo, convert(tokens.at(3)));
		break;
	case 'l':
		newComponent = new Inductor(tokens.at(0), nodeOne, nodeTwo, convert(tokens.at(3)));
		break;
	case 'v':
		newComponent = new IndVoltageSource(tokens.at(0), nodeOne, nodeTwo, convert(tokens.at(3)));
		break;
	case 'i':
		newComponent = new IndCurrentSource(tokens.at(0), nodeOne, nodeTwo, convert(tokens.at(3)));
		break;
	case 'e':
		newComponent = new VCVS(tokens.at(0), nodeOne, nodeTwo, std::stoi(tokens.at(3)), std::stoi(tokens.at(4)), convert(tokens.at(5)));
		break;
	case 'g':
		newComponent = new VCCS(tokens.at(0), nodeOne, nodeTwo, std::stoi(tokens.at(3)), std::stoi(tokens.at(4)), convert(tokens.at(5)));
		break;
	/*
	case 'k':

		break;
	case 'o':

		break;
	case 'a':

		break;
	*/
	case 'p':
		newComponent = new VPulse(tokens.at(0), nodeOne, nodeTwo, convert(tokens.at(3)), convert(tokens.at(4)), 
			convert(tokens.at(5)), convert(tokens.at(6)), convert(tokens.at(7)), std::stoi(tokens.at(8)));
		break;
	case 's':
		newComponent = new VStep(tokens.at(0), nodeOne, nodeTwo, convert(tokens.at(3)), convert(tokens.at(4)),
			convert(tokens.at(5)), convert(tokens.at(6)), std::stoi(tokens.at(7)));
		break;
	default:
		Error("Unknown component.");
		exit(-1);
		break;
	}

	incrVoltages(newComponent->numVoltagesToAdd(nodes));
	incrCurrents(newComponent->numCurrentsToAdd());

	return newComponent;
}

double Netlist::convert(string& token)
{
	if (isANumber(token)) {
		return std::stod(token);
	}
	else {
		/*
		there should be a trailing modifier after the number
		F	E-15	femto
		P	E-12	pico
		N	E-9		nano
		U	E-6		micro
		M	E-3		milli
		K	E+3		kilo
		MEG E+6 	mega
		G 	E+9 	giga
		T 	E+12 	tera
		*/
		string* value = splitString(token);
		double baseNum = stod(value[0]);

		if (value[1] == "a") {
			return baseNum *= pow(10, -18);
		}
		else if (value[1] == "f") {
			return baseNum *= pow(10, -15);
		}
		else if (value[1] == "p") {
			return baseNum *= pow(10, -12);
		}
		else if (value[1] == "n") {
			return baseNum *= pow(10, -9);
		}
		else if (value[1] == "u") {
			return baseNum *= pow(10, -6);
		}
		else if (value[1] == "m") {
			return baseNum *= pow(10, -3);
		}
		else if (value[1] == "k") {
			return baseNum *= pow(10, 3);
		}
		else if (value[1] == "meg") {
			return baseNum *= pow(10, 6);
		}
		else if (value[1] == "g") {
			return baseNum *= pow(10, 9);
		}
		else if (value[1] == "t") {
			return baseNum *= pow(10, 12);
		}
		else {
			cerr << "Error: " << token << endl;
			Error("Unknown value modifier in netlist.");
			exit(-1);
		}
	}
	return 0.0;
}

void Netlist::calculateNewIndicies()
{
	// this calculates the new indices for components that need to augment a matrix
	// calculation of the new indices for components that need to use a 
	// current equation had to wait until the number of voltage equations was known
	// newIndex is really numVoltages + 1 to start at new rows/columns augmented on matrices
	// then we subtract 1 to offset the fact that matrix indices start at 0
	// this is why we post increment the value newIndex

	int newIndex = numVoltages;
	for (Component* c : circuitElements) {
		if (Inductor* p = dynamic_cast<Inductor*> (c)) {
			p->setNewIndex(newIndex++);
		}
		else if (IndVoltageSource* p = dynamic_cast<IndVoltageSource*> (c)) {
			p->setNewIndex(newIndex++);
		}
		else if (VCVS* p = dynamic_cast<VCVS*> (c)) {
			p->setNewIndex(newIndex++);
		}
		else if (VPulse* p = dynamic_cast<VPulse*> (c)) {
			p->setNewIndex(newIndex++);
		}
		else if (VStep* p = dynamic_cast<VStep*> (c)) {
			p->setNewIndex(newIndex++);
		}
	}
}

void Netlist::populateMatrices()
{
	resizeMatrices();
	calculateNewIndicies();
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
		case TIMEEULER:
			solveTimeBackwardEuler(outputFilename);
			break;
		case TIMETRAPEZOIDAL:
			solveTimeTrapezoidalRule(outputFilename);
			break;
		case DC:
			//TODO implement DC analysis(transient at a single time point)
			break;
		}
	}
}

void Netlist::solveFrequency(string& outputFilename)
{

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
		else {
			//free memory from unused Component objects
			c->~Component();
		}
	}

	// attempt to open file
	ofstream fileWriter(outputFilename);
	if (fileWriter.fail())
		Error("problem opening the output file");
	
	SparseMatrix<double, Eigen::ColMajor> COverH(C), AStatic(G);

	COverH /= stepSize;
	AStatic += COverH;


	const int size = numVoltages + numCurrents;
	Matrix<double, Dynamic, 1> BDynamic(size, 1), BCurrentTimePoint(size, 1), XPreviousTimePoint(size, 1), XCurrentTimePoint(size, 1);

	BDynamic.setZero();
	BCurrentTimePoint.setZero();
	XCurrentTimePoint.setZero();
	XPreviousTimePoint.setZero();

	SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int>> solver;

	AStatic.makeCompressed();
	solver.compute(AStatic);

	for (int i = 0; i < numSteps; i++) {
		BCurrentTimePoint = B;

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
		
		BDynamic = COverH * XPreviousTimePoint;

		BCurrentTimePoint = BCurrentTimePoint + BDynamic;

		XCurrentTimePoint = solver.solve(BCurrentTimePoint);

		magnitude = XCurrentTimePoint(nodeToTrack, 0);

		fileWriter << currentTime << "\t" << magnitude << endl;

		XPreviousTimePoint = XCurrentTimePoint;
		//XCurrentTimePoint = X;
		currentTime += stepSize;
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
		else {
			//free memory from unused Component objects
			c->~Component();
		}
	}

	// attempt to open file
	ofstream fileWriter(outputFilename);
	if (fileWriter.fail())
		Error("problem opening the output file");

	SparseMatrix<double, Eigen::ColMajor> COverH(C), GOverTwo(G), AStatic, BStatic;

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

	for (int i = 0; i < numSteps; i++) {
		BCurrentTimePoint = B;

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
		//XCurrentTimePoint = X;
		currentTime += stepSize;
	}

	fileWriter.close();
}

double Netlist::calcMagnitude(int row, int col)
{
	return 0.0;
}

double Netlist::calcPhase(int row, int col)
{
	return 0.0;
}

void Netlist::prettyPrintNetlist()
{
	//TODO add toString functions to each Component derived class
}

void Netlist::prettyPrintMatrices()
{
	cout << "\n\n\nG:\n" << G << endl;
	cout << "\n\n\nX:\n" << X << endl;
	cout << "\n\n\nC:\n" << C << endl;
	cout << "\n\n\nB:\n" << B << endl;
}

bool Netlist::isANumber(string& token)
{
	for (char c : token) {
		if (!isdigit(c) && c != '.' && c != 'e' && c != '-') {
			return false;
		}
	}
	return true;
}

string* Netlist::splitString(string& token)
{
	int index = 0;
	for (char c : token) {
		if (!isdigit(c) && c != '.') {
			break;
		}
		else
			index++;
	}

	string* result = new string[2];
	result[0] = token.substr(0, index);
	result[1] = token.substr(index);
	return result;
}
