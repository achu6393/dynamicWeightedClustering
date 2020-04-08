/*
 * Constant.h
 *
 *  Created on: 28 May 2019
 *      Author: apandiya
 */

#ifndef CONSTANT_H_
#define CONSTANT_H_



namespace constants{

//constants
const int inf = 99999999;
const signed negInf = -99999999;

const double pi = 3.1415926535897932384626;



//System Parameters
double alpha_mat = 3.21;
double nAcous= 0.00858;


//Circuit Efficiencies
float cEffInd = 0.50;			//IPT Efficiency[%] - worst case
float cEffPD = 0.90;			//Piezo Driver Efficiency[%]
float cEffAC = 0.98;			//Discharge
float cEffPDAC = 0.98;			//Acoustic to DC Efficiency[%]


//Power Delivery Vehicle
float C_max = 187;						//Battery Energy Capacity [Wh]
float I_PDV = 15.96;						//Current Intake [A]
float V_PDV = 22.8;							//Voltage Intake[V]
float pFlight = 363.888;						// Power Rating of the PDV [W]
float fSpeed = 2160000;						//Flight Speed[cm/h]

//Sensor System Parameter
float Ton = 1;					//Ton [ms]
float Toff = 0.5;				//Toff [ms]
double mu = 0.10;				//SC Leakage[%]
int allowableFails =5;			//Allowable Fails per Sensor


//Acoustic Limits
double acousE = 12;					//Energy Sending [J]
double dAcou = 70;					//Maximum Acoustic distance [cm]
double g_k_ = 0.01;					//Channel Loss
double N_o = 21.5;					//Far Field range
double fAcous = 2*pi*47500;			//angular frequency of Acoustic waves[Hz]

//Acoustic Conditioning
float V_send = 150;					//Sending Voltage
float V_d = 0.9;					//Voltage drop cross diode
float tPulse = 0.004;

//Algorithm Constants
float wInc = 0.1;
int NP = 5 + 1;
int CR = 50;
int pop = 50;
int gen = 50;

int wDistance = 1;
int wEnergy = 1000000;


}

#endif /* CONSTANT_H_ */
