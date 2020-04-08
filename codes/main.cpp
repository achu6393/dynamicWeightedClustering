/*
 * MainMin.cpp
 *
 *  Created on: 17 May 2019
 *      Author: apandiya
 */


#include "GAClustering.hpp"

#include <conio.h>

#include <fstream>
#include <string>



using namespace std;

int main(){
	// Sensor Network Details
	int nSensor=2000;
	vector<SensorNode> sSet(nSensor);

	// Output File
	if(ofile.std::ofstream::is_open()) ofile.close();
	ofile.open("OutputFile.csv", std::ofstream::out | std::ofstream::app);
	ofile<<"Pecentage of nodes charged:"<<",";
	// Power Delivery Vehicle Details
	PDV pdv;
	bool flightStatus = true;

	//Check flight energy
	if (pdv.ePDV < 0.30*C_max){
		printf("Energy Low. Flight to BS");
		flightStatus = false;
		break;
	}

	do{
		pdv.readEnergy(sSet);
		if (pdv.eCheck(sSet)) {
			flightStatus = false;
			pdv.ePDV = C_max;
		}
	}while (flightStatus);

		return (0);

}

