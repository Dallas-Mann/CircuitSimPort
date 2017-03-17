#include "CNT.h"
#include "Netlist.h"
#include "Resistor.h"
#include "Capacitor.h"
#include "Inductor.h"

CNT::CNT(string id, int nodeOne, int nodeTwo, string sourceFile) :
	id{ id }, nodeOne{ nodeOne }, nodeTwo{ nodeTwo }, sourceFile{ sourceFile } {
	//populate list of subComponents to get the number of voltages/currents
	readNetlist(sourceFile);
}

void CNT::insertStamp(vector<T>& GVals, vector<T>& XVals, vector<T>& CVals, vector<T>& BVals) {
	for (Component* c : subComponents) {
		c->insertStamp(GVals, XVals, CVals, BVals);
	}
}

int CNT::numVoltagesToAdd(vector<int> &nodes) {
	int val = 0;
	for (Component* c : subComponents) {
		val += c->numVoltagesToAdd(nodes);
	}
	cout << "VALUE: " << val << endl;
	return val;
}

int CNT::numCurrentsToAdd() {
	int val = 0;
	for (Component* c : subComponents) {
		val += c->numCurrentsToAdd();
	}
	return val;
}

void CNT::setNodeOffset(int offset) {
	nodeOffset = offset;
}

void CNT::readNetlist(string inputFilename)
{
	// attempt to open file
	ifstream fileReader(inputFilename);
	if (fileReader.fail())
		Utilities::Error("problem opening the input file");

	int counter = 0;
	bool firstLine = true;

	vector<string> tokens;
	string nextLine;
	std::getline(fileReader, nextLine);
	while (!(nextLine == ".end")) {
		tokens = Utilities::parseLine(nextLine);
		std::getline(fileReader, nextLine);

		if (tokens.at(0).at(0) == '#') {
			//do nothing, it's a comment
		}
		else {
			Component* newComponent = Utilities::parseComponent(tokens);
			subComponents.push_back(newComponent);
		}
		firstLine = false;
	}
	fileReader.close();
}

int CNT::modifyNodeValues() {
	//clear list of subComponents to read them in again and modify the node values
	for (Component* c : subComponents) {
		c->~Component();
	}
	subComponents.clear();
	//modify node one and node two values to account for the nodeOffset, only if they are not equal to the CNT nodeone/nodetwo values

	// attempt to open file
	ifstream fileReader(sourceFile);
	if (fileReader.fail())
		Utilities::Error("problem opening the input file");

	int highestNode = 0;
	bool firstLine = true;

	vector<string> tokens;
	string nextLine;
	std::getline(fileReader, nextLine);
	
	cout << "**************CNT****************" << endl;

	while (!(nextLine == ".end")) {
		tokens = Utilities::parseLine(nextLine);
		std::getline(fileReader, nextLine);

		if (tokens.at(0).at(0) == '#') {
			//do nothing, it's a comment
		}
		else {
			int currentNodeOne = stoi(tokens.at(1));
			int currentNodeTwo = stoi(tokens.at(2));

			//default to the CNT node values
			int newNodeOne = currentNodeOne;
			int newNodeTwo = currentNodeTwo;

			//change node values based on node Offset, and keep incrementing with counter
			//this maps the first and last node in the netlist file to the actual nodes of the CNT in the general circuit
			if (firstLine && currentNodeOne != 0) {
				newNodeOne = nodeOne;
			}
			else if (!firstLine && currentNodeOne != 0) {
				newNodeOne = currentNodeOne + nodeOffset;
			}

			if (nextLine != ".end" && currentNodeTwo != 0) {
				newNodeTwo = currentNodeTwo + nodeOffset;
			}
			else if (nextLine == ".end" && currentNodeTwo != 0) {
				newNodeTwo = nodeTwo;
			}

			/*
			cout << "Old Node One: " << currentNodeOne << endl;
			cout << "New Node One: " << newNodeOne << endl;
			cout << "Old Node Two: " << currentNodeTwo << endl;
			cout << "New Node Two: " << newNodeTwo << endl;
			cout << "\n" << endl;
			*/

			if (newNodeOne > highestNode)
				highestNode = newNodeOne;
			if (newNodeTwo > highestNode)
				highestNode = newNodeTwo;

			//convert back to strings and create the component objects
			tokens.at(1) = std::to_string(newNodeOne);
			tokens.at(2) = std::to_string(newNodeTwo);

			Component* newComponent = Utilities::parseComponent(tokens);
			subComponents.push_back(newComponent);
		}
		firstLine = false;
	}
	fileReader.close();
	return highestNode;
}