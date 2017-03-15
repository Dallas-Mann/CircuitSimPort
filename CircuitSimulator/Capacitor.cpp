#include "Capacitor.h"

Capacitor::Capacitor(string id, int nodeOne, int nodeTwo, double capacitance) : 
	id{ id }, nodeOne{ nodeOne }, nodeTwo{ nodeTwo }, capacitance{ capacitance } {
}

void Capacitor::insertStamp(vector<T>& GVals, vector<T>& XVals, vector<T>& CVals, vector<T>& BVals) {
	// 0th node is ground node, and thus not implemented in our matrices
	// because of this we need to offset all the matrix indices by -1
	int indexOne = nodeOne - 1;
	int indexTwo = nodeTwo - 1;
	if (nodeOne == 0) {
		CVals.push_back({ indexTwo, indexTwo, capacitance });
	}
	else if (nodeTwo == 0) {
		CVals.push_back({ indexOne, indexOne, capacitance });
	}
	else {
		CVals.push_back({ indexOne, indexOne, capacitance });
		CVals.push_back({ indexTwo, indexTwo, capacitance });
		CVals.push_back({ indexOne, indexTwo, -capacitance });
		CVals.push_back({ indexTwo, indexOne, -capacitance });
	}
}

int Capacitor::numVoltagesToAdd(vector<int> &nodes) {
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

int Capacitor::numCurrentsToAdd() {
	return 0;
}
