// CircuitSimulator.cpp : Defines the entry point for the console application.
#include "Netlist.h"

int main(int argc, char* argv[]) {
	// check for the correct number of arguments
	if (argc != 3)
		Utilities::Usage(argv[0], "Wrong number of arguments.");

	Netlist *netlist = new Netlist();
	netlist->readNetlist(argv[1]);
	netlist->populateMatrices();
	//netlist->prettyPrintMatrices();
	netlist->simulate(argv[2]);
    return 0;
}