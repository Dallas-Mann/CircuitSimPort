#include "VCCS.h"

VCCS::VCCS(string id, int nodeOne, int nodeTwo, int nodeThree, int nodeFour, double gain) :
	id{ id }, nodeOne{ nodeOne }, nodeTwo{ nodeTwo }, nodeThree{ nodeThree }, nodeFour{ nodeFour }, gain{ gain } {
}

void VCCS::insertStamp(vector<T>& GVals, vector<T>& XVals, vector<T>& CVals, vector<T>& BVals) {
	// 0th node is ground node, and thus not implemented in our matrices
	// because of this we need to offset all the matrix indices by -1
	int indexOne = nodeOne - 1;
	int indexTwo = nodeTwo - 1;
	int indexThree = nodeThree - 1;
	int indexFour = nodeFour - 1;
	if (!(nodeOne == 0 || nodeThree == 0)) {
		GVals.push_back({ indexThree, indexOne, gain });
	}
	if (!(nodeTwo == 0 || nodeThree == 0)) {
		GVals.push_back({ indexThree, indexTwo, -gain });
	}
	if (!(nodeOne == 0 || nodeFour == 0)) {
		GVals.push_back({ indexFour, indexOne, -gain });
	}
	if (!(nodeTwo == 0 || nodeFour == 0)) {
		GVals.push_back({ indexFour, indexTwo, gain });
	}
}

int VCCS::numVoltagesToAdd(vector<int> &nodes) {
	int val = 0;
	if (std::find(nodes.begin(), nodes.end(), nodeOne) == nodes.end()) {
		nodes.push_back(nodeOne);
		val++;
	}
	if (std::find(nodes.begin(), nodes.end(), nodeTwo) == nodes.end()) {
		nodes.push_back(nodeTwo);
		val++;
	}
	if (std::find(nodes.begin(), nodes.end(), nodeOne) == nodes.end()) {
		nodes.push_back(nodeThree);
		val++;
	}
	if (std::find(nodes.begin(), nodes.end(), nodeTwo) == nodes.end()) {
		nodes.push_back(nodeFour);
		val++;
	}
	return val;
}

int VCCS::numCurrentsToAdd() {
	return 0;
}
