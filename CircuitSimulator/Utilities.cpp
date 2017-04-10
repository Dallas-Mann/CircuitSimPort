#include "Utilities.h"

#include "Resistor.h"
#include "Capacitor.h"
#include "Inductor.h"
#include "IndVoltageSource.h"
#include "IndCurrentSource.h"
#include "VCVS.h"
#include "VCCS.h"
//TODO: implement these three components
//#include "MutualInductance.h"
//#include "OpAmp.h"
#include "VAC.h"
#include "VPulse.h"
#include "VStep.h"
#include "CNT.h"

// Print the correct usage in case of user syntax error.
void Utilities::Usage(char* arg, string message) {
	cerr << endl << "Message: " << message << endl;
	cerr << "Usage:   " << arg << " <netlist.txt> <soultion.txt>" << endl;
	exit(-1);
}

void Utilities::Error(string message) {
	cerr << endl << "Message: " << message << endl;
	exit(-1);
}

bool Utilities::areSame(double a, double b){
	return fabs(a - b) < 0.1;
}

vector<string> Utilities::parseLine(string& nextLine)
{
	std::transform(nextLine.begin(), nextLine.end(), nextLine.begin(), ::tolower);

	string token;
	vector<string> tokens;
	while (token != nextLine) {
		token = nextLine.substr(0, nextLine.find_first_of(" "));
		nextLine = nextLine.substr(nextLine.find_first_of(" ") + 1);
		tokens.push_back(token);
	}

	return tokens;
}

Component* Utilities::parseComponent(vector<string>& tokens)
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
		*/
	case 'a':
		newComponent = new VAC(tokens.at(0), nodeOne, nodeTwo, convert(tokens.at(3)), convert(tokens.at(4)),
			convert(tokens.at(5)), convert(tokens.at(6)), std::stoi(tokens.at(7)));
		break;
	case 'p':
		newComponent = new VPulse(tokens.at(0), nodeOne, nodeTwo, convert(tokens.at(3)), convert(tokens.at(4)),
			convert(tokens.at(5)), convert(tokens.at(6)), convert(tokens.at(7)), std::stoi(tokens.at(8)));
		break;
	case 's':
		newComponent = new VStep(tokens.at(0), nodeOne, nodeTwo, convert(tokens.at(3)), convert(tokens.at(4)),
			convert(tokens.at(5)), convert(tokens.at(6)), std::stoi(tokens.at(7)));
		break;
	case 't':
		newComponent = new CNT(tokens.at(0), nodeOne, nodeTwo, tokens.at(3));
		break;
	default:
		Utilities::Error("Unknown component.");
		exit(-1);
		break;
	}

	return newComponent;
}

double Utilities::convert(string& token)
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
			Utilities::Error("Unknown value modifier in netlist.");
			exit(-1);
		}
	}
}

bool Utilities::isANumber(string& token)
{
	//TODO: should check if number is a regular expression instead of this
	for (char c : token) {
		if (!isdigit(c) && c != '.' && c != 'e' && c != '-') {
			return false;
		}
	}
	return true;
}

string* Utilities::splitString(string& token)
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