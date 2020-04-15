/*
 * SensorNode.hpp
 *
 *  Created on: 21 May 2019
 */

#ifndef SENSORNODE_HPP_
#define SENSORNODE_HPP_

#include "pointCircle.hpp"
#include "constants.h"
#include <fstream>
std::ofstream ofile;


using namespace constants;

struct cluster;
bool Psensor_flag = false;

double get_E();
bool mySort (const point* a,const point* b);
int randIndGenerator(int size);


struct Sense{
public:
	mutable double SC_C,SC_R, SC_E ,  SC_V;			//Sleep consumption not included as included in
	double V_Sense, I_Sense, I_Idle, P_ont;
	int TCycle;											//information from the energy status.
	double SC_VMax, SC_VMin , SC_VCritical;
	double SCycle, ICycle,CCycle;
	mutable int T;
};

struct Harvest{
public:
	mutable double effE1, effE2, effAC, effD;
};

struct SensorNode : public Harvest, public Sense{
public:
	point position;
	mutable bool E_flag = false;
	mutable int f_Count;

	//Member Functions
	void Initialize_Param();
	void acousTransfer(double);
	void T_Sensor();
	void P_Sensor();
	float updateWeight(const double v);
	void updateV();
	void updateE() { this->SC_E = this->SC_C * 0.5 * this->SC_V * this->SC_V;};
	void readEnergy(double t);

	//Constructors
	SensorNode() : position(), f_Count(0) {Initialize_Param();};
	SensorNode(int a) : position(a), f_Count(0) {Initialize_Param();};
	SensorNode(double a, double b) : position(a,b)  , f_Count(0) {Initialize_Param();};
};

struct cluster final{
public:
	std::vector<SensorNode*> pContains;
	point* O;


	int maxPopulation;

	cluster() : O(), maxPopulation(NP) {};
};

struct PDV: public Sense{
public:
	point fPosition;
	double tTime, ePDV;
	int tCycle;
	std::vector<point> flightPath;
	std::vector<std::vector<cluster>> targetC;
	std::vector<double> targetM;
	std::vector<std::vector<cluster>> trailC;
	std::vector<double> trailM;
	std::vector<std::vector<cluster>> solutionC;
	std::vector<double> solutionM;

	friend struct cluster;
	//Constructors
	PDV() : fPosition(0,0,15), tTime(0.0), tCycle(0){
		ePDV = C_max;};

	//Memeber Functions
	bool eCheck(std::vector<SensorNode> &S);
	void updatePosition(const point &sPos);
	void updateEnergy(double t) {ePDV -= (pFlight*(t+0.0056));};	//Including take off and landing 5 + 15
	double checkEnergy(double t) {return (pFlight*(t+0.0056));};
	void readEnergy(std::vector<SensorNode> &S);
	void geneticInitialisation(std::vector<cluster> &initC, int k,std::vector<SensorNode> &S,std::vector<SensorNode*> P);
	void geneticCrossover(std::vector<std::vector<cluster>> targetC, std::vector<std::vector<cluster>> &trailC);
	void GeneticClustering(std::vector<SensorNode> &S,std::vector<SensorNode*> P);
	int bestSolution();
	double geneticFitnessFunction(std::vector<SensorNode>S, std::vector<cluster> C, SensorNode* S0);
	void OnlyInductive(std::vector<SensorNode> &S, std::vector<SensorNode*> P);
};

int randIndGenerator(int size){
	std::uniform_int_distribution<int> distribution_ind(0,size);
	return distribution_ind(generator_int);
}

SensorNode* get_ptr(point* b){
    static typename std::aligned_storage<sizeof(SensorNode),alignof(SensorNode)>::type buffer;

    SensorNode* p=static_cast<SensorNode*>(static_cast<void*>(&buffer));
    ptrdiff_t const offset=static_cast<char*>(static_cast<void*>(&p->position))-static_cast<char*>(static_cast<void*>(&buffer));
    return static_cast<SensorNode*>(static_cast<void*>(static_cast<char*>(static_cast<void*>(b))-offset));
}

void SensorNode::Initialize_Param(){
	//Randomly placing diffferent type of nodes
	if (distribution_int(generator_int) % 2 ==0)
		P_Sensor();
	else
		T_Sensor();
	//Harvest
	effE1 = cEffInd;
	effD = cEffPD;
	effAC = cEffAC;
	effE2 = cEffPDAC;
}

void SensorNode::T_Sensor(){
	//Sense
	V_Sense = 1.5;
	I_Sense = 5.4 * pow(10.0, -6.0);
	I_Idle = 4 * pow(10.0, -9.0);
	P_ont = 0.002;
	TCycle = 10; //[x 10sec]
	T = 10;
	SCycle = 0.002;
	ICycle = 99.498;
	CCycle = 0.500;

	//Supercap - HB Supercapacitors
	SC_C = 6;							//Capacitance [F]
	SC_R =0.10;						//ESR  [Ohm]
	SC_VMin = 1.75;						//Voltage [V]
	SC_VCritical = 1.50;				//[V]
	SC_VMax = 2.50;						//[V]
	SC_V = get_E()*0.6;
	updateE();
}

void SensorNode::P_Sensor(){
	//Sense
	V_Sense = 3.30;
	I_Sense = 1.2 * pow(10.0, -3.0);
	I_Idle = 4 * pow(10.0, -6.0);
	P_ont = 0.002;
	TCycle = 1;  //[x 10sec]
	SCycle = 0.002;
	ICycle = 9.498;
	CCycle = 0.500;
	T = 1;

	//Supercap
	SC_C = 3;
	SC_R = 0.2;
	SC_VMin = 3.50;
	SC_VCritical = 3.30;
	SC_VMax = 5.00;
	SC_V = get_E();
	updateE();

	//flag for PDV
	Psensor_flag = true;

}

void PDV::OnlyInductive(std::vector<SensorNode> &S, std::vector<SensorNode*> P){
	//Step 1: Set the Max. K
	double tempS = 0.5 * S[4].SC_C * S[4].SC_VMax *S[4].SC_VMax;	//T Sensor
	if (Psensor_flag) tempS = 37.5;									//P Sensor
	int needS = P.size(),unchargedS=0;

	int k =  floor(((ePDV-6)*3600)/(tempS + (1.8*3600)));
	if (k > static_cast<int>(P.size())) k = P.size(); //urgent or no cycle
	double dPDV =0;
	for (int j = 0; j < static_cast<int>(P.size()); j++){
			this->flightPath.push_back(P[j]->position);
		}
//	ofile<<flightPath.size() << "," ;

	do{
		std::vector<double> t = fPosition.calcDistance(flightPath);
		int fInd = std::distance(t.begin(),std::min_element(t.begin(),t.end()));

		//Checking if enough energy for PDV visit
		if(checkEnergy(flightPath[fInd].calcDistance(origin)/fSpeed)
				+ checkEnergy(fPosition.calcDistance(flightPath[fInd])/fSpeed) > this->ePDV)
			break;

		dPDV +=  fPosition.calcDistance(flightPath[fInd]);
		updatePosition(flightPath[fInd]);
		
		//Calculate eWSN
		//Update Center Node

			SensorNode* thisS = P[fInd];
			//Charging Center node to 100%
			double pckt = ((0.5*thisS->SC_C) * (pow(thisS->SC_VMax,2)- pow(thisS->SC_V,2)))/cEffInd;
			this->ePDV -= (pckt/3600)*6;			//Joules -> Wh;
			this->tTime += pckt/5; 		//Charging rate of 5J/s
			thisS->SC_V = thisS->SC_VMax;
			thisS->updateE();


			thisS->E_flag =false;
			thisS->f_Count = 0;

		flightPath.erase(flightPath.begin() + fInd);P.erase(P.begin() + fInd);
		k--;
	}while(k !=0 && this->ePDV > 0.3*C_max);
	dPDV += fPosition.calcDistance(origin); //Flight back to Base
	updatePosition(origin);

	for(int i = 0; i < static_cast<int>(S.size()); i++){
		if (S[i].SC_V < S[i].SC_VMin) unchargedS++;
	}

	float perC = ((needS - unchargedS)*100)/needS;
	ofile<<perC<< ",";

//	std::cout<<perC<< '\t';
	this->tTime = ((C_max - this->ePDV)/pFlight)*3600;
												// XX Flight time + tafe off, landing time[s]
	//Update all Sensor Energy
	for(int i=0; i<static_cast<int>(S.size()); i++){
		if(S[i].SC_V > S[i].SC_VCritical){
			S[i].readEnergy(this->tTime);
			S[i].updateV();
		}
	}
};

float SensorNode::updateWeight(const double v){
	if (this->SC_VMax - v <0.25) return 10;
	else if (this->SC_VMax - v > 0.25 && this->SC_VMax - v <1.0) return 9;
	else if (this->SC_VMax - v > 1.0 && this->SC_VMax - v <1.25) return 8;
	else if (this->SC_VMax - v > 1.25 && this->SC_VMax - v <1.5) return 7;
	else if (this->SC_VMax - v > 1.5 && this->SC_VMax - v <1.75) return 4;
	else
		return 3;
}

void SensorNode::updateV(){
	double temp = (2* this->SC_E)/this->SC_C;
	this->SC_V = sqrt(temp);
}

void SensorNode::acousTransfer(double d){
	d = d/100;						// cm to m
	double g_d =  exp(pow(fAcous,nAcous) * d * -1 *  alpha_mat);

	this->SC_E += this->effD*this->effD*this->effE2*g_d*acousE;
	this->updateV();
	this->SC_V = std::min(this->SC_V, SC_VMax);updateE();
}

bool PDV :: eCheck(std::vector<SensorNode> &S){
	std::vector<SensorNode*> P;		//Less Energy SN
//	std::vector<SensorNode*> pointerS;		//Charged SN


	for(unsigned int i = 0; i< S.size(); i++){
		if(S[i].SC_V < S[i].SC_VMin){
			if (S[i].E_flag && S[i].f_Count == 0)
				S[i].position.w = S[i].updateWeight(S[i].SC_V);
			else if (S[i].E_flag && S[i].f_Count >= allowableFails)
				S[i].position.w = 2;
			else {
				S[i].E_flag = true;
				S[i].position.w = 8;
			}
			P.push_back(&S[i]);
//			pointerS.push_back(&S[i]);
		}
//		else pointerS.push_back(&S[i]);
	};
	if(P.size() > 20) {
		/*Comment the method of WPT not used*/
		GeneticClustering(S, P);
//		OnlyInductive(S, P);
		return true;
	}
	return false;
}

void PDV::updatePosition(const point &sPos){
	updateEnergy(fPosition.calcDistance(sPos)/fSpeed);		//t = Distance/speed
	fPosition.x = sPos.x;
	fPosition.y = sPos.y;
}

void SensorNode::readEnergy(double t){
	int senseCycle = floor(t/(3*this->TCycle));

	SC_E -= ((senseCycle*SCycle*V_Sense*I_Sense)+(senseCycle*ICycle*I_Idle*V_Sense) + (senseCycle*CCycle*0.00245)); // J Energy for
																													// Communication and Micrprocessor
	this->updateV();
	int d=0;d++;
}

void PDV::readEnergy(std::vector<SensorNode> &S){

	for(unsigned int i = 0; i< S.size(); i++){
		S[i].T--;
		//Ping at every 1 seconds
		//If sense cycle
		if (S[i].T == 0){
			//Check if SN has enough energy
			if(S[i].SC_V > S[i].SC_VCritical) S[i].readEnergy(6);
			//Else increase failure count and move point to critical
			else {
				S[i].f_Count++;
				S[i].E_flag = true;   //y??
			}
			S[i].T = S[i].TCycle;
		}
		//Idle cycle
		else{
			S[i].SC_E -= S[i].SC_V *  S[i].ICycle;
			S[i].updateV();
		}

	}
}

point* constton(const point* p){
	return const_cast<point*> (p);
}

bool mySort (const point* a, const point* b){
	return (a->w < b->w);
}

// Optimisaion Algorithm

void PDV::geneticInitialisation(std::vector<cluster> &initC, int k, std::vector<SensorNode> &S, std::vector<SensorNode*> P){
	//Select cluster centers randomly
	for(int j =0; j <k; j++){
		initC.emplace_back();
		int thisInd = randIndGenerator(P.size()-1);
		initC[j].O = constton(&P[thisInd]->position);
		initC[j].pContains.clear();
		initC[j].pContains.push_back(P[thisInd]);
		P.erase(P.begin() + thisInd);
	}
	//Assign the end nodes
	for(int j = 0; j < static_cast<int>(S.size()); j++){
		int thisInd = 0;
		for(int l=1; l<k; l++){
			if(S[j].position.calcDistance(*initC[l-1].O) >  S[j].position.calcDistance(*initC[l].O))
				thisInd = l;
		}
		//If the nearest center is within Acoustic limits
		if (S[j].position.calcDistance(*initC[thisInd].O) < dAcou && S[j].position.calcDistance(*initC[thisInd].O) != 0)
			initC[thisInd].pContains.push_back(&S[j]);
	}
}

void deleteDuplicates(cluster &C){
	for(int i =0; i<static_cast<int>(C.pContains.size()); i++){
		for(int j= i+1; j < static_cast<int>(C.pContains.size()); j++){
			if(C.pContains[i] == C.pContains[j]){
				C.pContains.erase(C.pContains.begin() + j); //If matches, remove
				j--;
			};
		}
	}
}
void PDV::geneticCrossover(std::vector<std::vector<cluster>> targetC, std::vector<std::vector<cluster>> &trailC){

	//Two-point crossover
	for(int i=0; i<pop; i = i+2){
		int r1 = randIndGenerator(pop-1);
		int r2 = randIndGenerator(pop-1);
		trailC[i] = targetC[r1];
		trailC[i+1] = targetC[r2];
//
		for(int k = 0; (k < static_cast<int>(trailC[i].size())) && CR < randIndGenerator(100); k++){
			bool flag =false;
			for(int m=0; m<static_cast<int>(std::min(trailC[i][k].pContains.size(),trailC[i + 1][k].pContains.size())); m++){
				if (CR <randIndGenerator(100)) {flag = true; break;}
				trailC[i][k].pContains[m] = trailC[i+1][k].pContains[m];
				trailC[i+1][k].pContains[m] = trailC[i][k].pContains[m];
			}if (flag) break;
		}

	//Swap Mutation
	int randInd1 = randIndGenerator(trailC[i].size() - 1);
	int randInd2 = randIndGenerator(trailC[i].size() - 1);
	if (randInd1 != randInd2) //No point swapping in the same cluster
		std::swap(trailC[i][randInd1].pContains[randIndGenerator(trailC[i][randInd1].pContains.size()- 1)],
				trailC[i][randInd2].pContains[randIndGenerator(trailC[i][randInd2].pContains.size()- 1)]);


	randInd1 = randIndGenerator(trailC[i+1].size()- 1);
	randInd2 = randIndGenerator(trailC[i+1].size()- 1);
	if (randInd1 != randInd2)
		std::swap(trailC[i+1][randInd1].pContains[randIndGenerator(trailC[i+1][randInd1].pContains.size()- 1)],
				trailC[i+1][randInd2].pContains[randIndGenerator(trailC[i+1][randInd2].pContains.size()- 1)]);

	//Delete duplicates
	for(int j =0; j <  static_cast<int>(trailC[i].size()); j++){ //For every cluster
		deleteDuplicates(trailC[i][j]);
	}
	//again for i+1
	for(int j =0; j <  static_cast<int>(trailC[i+1].size()); j++){ //For every cluster
		deleteDuplicates(trailC[i+1][j]);
	}

	}
};


int PDV::bestSolution(){
	int bestInd= 0;
	for(int i =1; i < pop; i++){
		if(this->targetM[bestInd] < this->targetM[i]) bestInd = i;
	}
	return bestInd;
}

void penalityFunction(std::vector<cluster> &C, std::vector<SensorNode> &S){
	for (int k=0; k< static_cast<int>(C.size()); k++){
		if (static_cast<int>(C[k].pContains.size()) > NP){
			C.emplace_back();
			do{
				int randInd = randIndGenerator(C[k].pContains.size()-1);
				C.emplace_back();
				C[C.size()].pContains.push_back(C[k].pContains[randInd]);
				C[k].pContains.erase(C[k].pContains.begin() + randInd);
			}while(static_cast<int>(C[k].pContains.size()) != NP);
			C[C.size()].O = &C[C.size()].pContains[0]->position;
		}
		if (static_cast<int>(C[k].pContains.size()) < NP){
			for(int j = 0; j < static_cast<int>(S.size()); j++){
				if(C[k].pContains[0]->position.calcDistance(S[j].position) < dAcou && &(C[k].pContains[0]->position) != &S[j].position)
					C[k].pContains.push_back(&S[j]);
			}
		};
		deleteDuplicates(C[k]);
	}
};

void kOptimisation(std::vector<cluster> &C){
	for (int k=0; k< static_cast<int>(C.size())-1; k++){
		for(int j = k+1; j< static_cast<int>(C.size()); j++){
			if(C[k].pContains[0]->position.calcDistance(C[j].pContains[0]->position) < 0.5 * dAcou){
				C[k].pContains.push_back(C[j].pContains[0]);
				C.erase(C.begin() + j);
			}
		}
	}
}

void PDV::GeneticClustering(std::vector<SensorNode> &S, std::vector<SensorNode*> P){

	//Intialisation
	//Step 1: Set the Max. K
	double tempS = 0.5 * S[4].SC_C * S[4].SC_VMax *S[4].SC_VMax;	//T Sensor
	if (Psensor_flag) tempS = 37.5;									//P Sensor
	int needS = P.size(),unchargedS=0;


	int k =  floor(((ePDV-6)*3600)/(tempS + (1.8*3600)));
	if (k > static_cast<int>(P.size())) k = P.size(); //urgent or no cycle


	//Step 2: Intialise the cluster
	for (int i=0; i<pop; i++){
		//Initialise Target Vector
		targetC.emplace_back();
		geneticInitialisation(targetC[i],k,S,P);
		//Compute Target Metric
		targetM.emplace_back();
		this->targetM[i] = geneticFitnessFunction(S, targetC[i], &S[0]);
		this->ePDV = C_max;this->tTime = 0.0;

		//Initialise the First target vector
		trailC.emplace_back();
		this->geneticInitialisation(this->trailC[i], k, S, P);
		//Compute trail Metric
		trailM.emplace_back();
		trailM[i] = geneticFitnessFunction(S, trailC[i], &S[0]);
		this->ePDV = C_max;this->tTime = 0.0;
	}

	int thisGen = 0;
	do{
		geneticCrossover(targetC, trailC);

		//Filter trail by degree

		//check if other centres are in acoustic distance


		for (int i=0; i< pop; i++){
//			penalityFunction(targetC[i],S); //Comment when includling angular penality function
			targetM[i] = geneticFitnessFunction(S, targetC[i], &S[0]);
			this->ePDV = C_max; this->tTime = 0.0;
//			kOptimisation(trailC[i]);   //Comment when not optimising for clusters
//			penalityFunction(trailC[i],S);
			trailM[i] = geneticFitnessFunction(S, trailC[i], &S[0]);
			this->ePDV = C_max; this->tTime = 0.0;
		}
		//Selection
		for(int i =0;  i <pop; i++){
			if (trailM[i] > targetM[i]){
				targetC[i] = trailC[i];
				targetM[i] = trailM[i];
			}
		}
		//Print all metric
//		ofile<<"Generation:"<<thisGen<<std::endl;
//		for(int i =0;  i <pop; i++){
//			ofile<<targetM[i] << "," << trailM[i]<<std::endl;
//		}

		thisGen++;
	} while (thisGen < gen);

	//Best Solution
	int bestInd = bestSolution();
//	ofile<<"Best Solution"<<std::endl;
//	ofile<<targetM[bestInd]<< "," << targetC[bestInd].size()<<std::endl;


	//Do every energy cal for the solution

	//Optimised Path
	double dPDV = 0;
	for (int j = 0; j < static_cast<int>(targetC[bestInd].size()); j++){
		this->flightPath.push_back(targetC[bestInd][j].pContains[0]->position);
	}
//	ofile<<flightPath.size() << "," ;
	do{
		std::vector<double> t = this->fPosition.calcDistance(this->flightPath);
		int fInd = std::distance(t.begin(),std::min_element(t.begin(),t.end()));

		//Checking if enough energy for PDV visit
		if(checkEnergy(flightPath[fInd].calcDistance(origin)/fSpeed)
				+ checkEnergy(fPosition.calcDistance(flightPath[fInd])/fSpeed) > this->ePDV)
			break;

//		std::cout<<flightPath[fInd].x << "," << flightPath[fInd].y << "," << t[fInd]<< "," << this->ePDV<<std::endl;
//		ofile<<flightPath[fInd].x << "," << flightPath[fInd].y << "," << t[fInd]<< "," << this->ePDV<<std::endl;
		dPDV +=  fPosition.calcDistance(flightPath[fInd]);
		updatePosition(flightPath[fInd]);
		//Calculate eWSN
		//Update Center Node
		{

			SensorNode* thisS = targetC[bestInd][fInd].pContains[0];
			//Charging Center node to 100%
			double pckt = ((0.5*thisS->SC_C) * (pow(thisS->SC_VMax,2)- pow(thisS->SC_V,2)))/cEffInd;
			this->ePDV -= (pckt/3600);			//Joules -> Wh
			thisS->SC_V = thisS->SC_VMax;
			thisS->updateE();

			//Discharging node to E2 - 8J for Power transfer
			thisS->SC_E -= acousE;
			thisS->updateV();

			thisS->E_flag =false;
			thisS->f_Count = 0;
		}
		//Update End Nodes
		for(int k = 1; k < static_cast<int>(targetC[bestInd][fInd].pContains.size()); k++){
			auto thisS = targetC[bestInd][fInd].pContains[k] - &S[0];
//			ofile<< "END NODE CHARGING From" << "," << S[thisS].SC_V<< ",";
			S[thisS].acousTransfer(S[thisS].position.calcDistance(targetC[bestInd][fInd].pContains[0]->position));
//			ofile<< S[thisS].SC_V<<std::endl;

			S[thisS].E_flag =false;
			S[thisS].f_Count = 0;
		}
		flightPath.erase(flightPath.begin() + fInd);
		targetC[bestInd].erase(targetC[bestInd].begin() + fInd);

	}while(flightPath.size());
	dPDV += fPosition.calcDistance(origin); //Flight back to Base
	updatePosition(origin);

	for(int i = 0; i < static_cast<int>(S.size()); i++){
		if (S[i].SC_V < S[i].SC_VMin) unchargedS++;
	}

	float perC = ((needS - unchargedS)*100)/needS;
	ofile<<perC<< ",";
	this->tTime = ((C_max - this->ePDV)/pFlight)*3600;
												// XX Flight time + tafe off, landing time[s]
	//Update all Sensor Energy
	for(int i=0; i<static_cast<int>(S.size()); i++){
		if(S[i].SC_V > S[i].SC_VCritical){
			S[i].readEnergy(this->tTime);
			S[i].updateV();
		}
	}
};

double  PDV::geneticFitnessFunction(std::vector<SensorNode> S, std::vector<cluster> C, SensorNode* S0){

	std::vector<point> thisFlightPath;
	for (int j = 0; j < static_cast<int>(C.size()); j++){
		thisFlightPath.push_back(C[j].pContains[0]->position);
	}

	double dPDV =0, eWSN =0;
	// Calculate dPDV
	while (thisFlightPath.size()){
		std::vector<double> t = fPosition.calcDistance(thisFlightPath);
		int fInd = std::distance(t.begin(),std::min_element(t.begin(),t.end()));
		dPDV +=  fPosition.calcDistance(thisFlightPath[fInd]);
		updatePosition(thisFlightPath[fInd]);

		//Calculate eWSN
		//Update Center Node
		{
			auto thisS = C[fInd].pContains[0] - S0;
			//Charging Center node to 100%
			double pckt = ((0.5*S[thisS].SC_C) * (pow(S[thisS].SC_VMax,2)- pow(S[thisS].SC_V,2)))/cEffInd;
			this->ePDV -= (pckt/3600);			//Joules -> Wh
			S[thisS].SC_V = S[thisS].SC_VMax;
			S[thisS].updateE();

			//Discharging node to E2 - 8J for Power transfer
			S[thisS].SC_E -= acousE;
			S[thisS].updateV();
		}
		//Update End Nodes
		for(int k = 1; k < static_cast<int>(C[fInd].pContains.size()); k++){
			auto thisS = C[fInd].pContains[k] - S0;

			S[thisS].acousTransfer(S[thisS].position.calcDistance(C[fInd].pContains[0]->position));
		}
		thisFlightPath.erase(thisFlightPath.begin() + fInd);
		C.erase(C.begin() + fInd);
	}
	dPDV += fPosition.calcDistance(origin); //Flight back to Base
	dPDV = dPDV/100; // from cm to m
	updatePosition(origin);

	this->tTime = ((C_max - this->ePDV)/pFlight)*3600;
												// XX Flight time + tafe off, landing time[s]

	//Update all Sensor Energy
	for(int j = 0; j < static_cast<int>(S.size()); j++){
		if(S[j].SC_V > S[j].SC_VCritical) {
			S[j].readEnergy(this->tTime);
			S[j].updateV();
		}
	}

	for(int j = 0; j < static_cast<int>(S.size()); j++){
		eWSN += S[j].SC_E;
	}


	//check pdv energy
	if (this->ePDV < 0) return -1;
	else
		return (wEnergy*eWSN) + (1/(wDistance*dPDV));
}

#endif /* SENSORNODE_HPP_ */
