#include "VCVS.h"

VCVS::VCVS(string id, int nodeOne, int nodeTwo, int nodeThree, int nodeFour, double gain) :
	id{ id }, nodeOne{ nodeOne }, nodeTwo{ nodeTwo }, nodeThree{ nodeThree }, nodeFour{ nodeFour }, newIndex{ -1 }, gain{ gain } {
}

void VCVS::insertStamp(vector<T>& GVals, vector<T>& XVals, vector<T>& CVals, vector<T>& BVals) {
	// 0th node is ground node, and thus not implemented in our matrices
	// because of this we need to offset all the matrix indices by -1
	int indexOne = nodeOne - 1;
	int indexTwo = nodeTwo - 1;
	int indexThree = nodeThree - 1;
	int indexFour = nodeFour - 1;
	if (!nodeOne == 0) {
		GVals.push_back({ newIndex, indexOne, -gain });
	}
	if (!nodeTwo == 0) {
		GVals.push_back({ newIndex, indexTwo, gain });
	}
	if (!nodeThree == 0) {
		GVals.push_back({ newIndex, indexThree, 1 });
		GVals.push_back({ indexThree, newIndex, 1 });
	}
	if (!nodeFour == 0) {
		GVals.push_back({ newIndex, indexFour, -1 });
		GVals.push_back({ indexFour, newIndex, -1 });
	}
}

int VCVS::numVoltagesToAdd(vector<int> &nodes) {
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

int VCVS::numCurrentsToAdd() {
	return 1;
}

void VCVS::setNewIndex(int val)
{
	newIndex = val;
}
