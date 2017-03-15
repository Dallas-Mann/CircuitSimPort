#include "Resistor.h"

Resistor::Resistor(string id, int nodeOne, int nodeTwo, double resistance) : 
	id{ id }, nodeOne{ nodeOne }, nodeTwo{ nodeTwo }, resistance{ resistance } {
}

void Resistor::insertStamp(vector<T>& GVals, vector<T>& XVals, vector<T>& CVals, vector<T>& BVals) {
	double conductance = 1.0 / resistance;
	// 0th node is ground node, and thus not implemented in our matrices
	// because of this we need to offset all the matrix indices by -1
	int indexOne = nodeOne - 1;
	int indexTwo = nodeTwo - 1;
	if (nodeOne == 0) {
		GVals.push_back({indexTwo, indexTwo, conductance});
	}
	else if (nodeTwo == 0) {
		GVals.push_back({ indexOne, indexOne, conductance });
	}
	else {
		GVals.push_back({ indexOne, indexOne, conductance });
		GVals.push_back({ indexTwo, indexTwo, conductance });
		GVals.push_back({ indexOne, indexTwo, -conductance });
		GVals.push_back({ indexTwo, indexOne, -conductance });
	}
}

int Resistor::numVoltagesToAdd(vector<int> &nodes) {
	int val = 0;
	if (std::find(nodes.begin(), nodes.end(), nodeOne) == nodes.end()) {
		nodes.push_back(nodeOne);
		val++;
	}
	if (std::find(nodes.begin(), nodes.end(), nodeTwo) == nodes.end()) {
		nodes.push_back(nodeTwo);
		val++;
	}
	return val;
}

int Resistor::numCurrentsToAdd() {
	return 0;
}
