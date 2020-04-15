/*
 * MinCircle.hpp
 *
 * Smallest encircling cluster for single and multi cluster sensor network.
 *
 * 		Created on: 17 May 2019
 */


#ifndef MINCIRCLE_HPP_
#define MINCIRCLE_HPP_

#include <cmath>
#include <valarray>
#include <random>
#include <chrono>
#include <iostream>
#include <stdio.h>


// Setting seed for random number generator
unsigned seed = static_cast<int> (std::chrono::system_clock::now().time_since_epoch().count());

// Uniform-Real random number generator
#include "constants.h"
static std::default_random_engine generator(seed);
std::uniform_real_distribution<double> distribution(3.0,5.0); //min and max value of V in a network


// Uniform-integer random number generator
static std::default_random_engine generator_int(seed);
std::uniform_int_distribution<int> distribution_int(0,200000); //Distance of the network in cm

int getIndInt(int limit){
	std::uniform_int_distribution<int> distribution_indInt(0,limit);
	return distribution_indInt(generator_int);
};

struct point{
public:
	double x,y;
	int w;

	//Constructors
	point () : x(distribution_int(generator_int)),y(distribution_int(generator_int)), w(10) {};
//	point () : x(5),y(5), w(10) {};
	point (int a) :x(getIndInt(a)),y(getIndInt(a)), w(10) {};
	point (double a, double b) : x(a), y(b), w(10) {};
	point (double a, double b, double c) : x(a),y(b), w(c){};

	//Member Functions
//	float updateWeight(double);
	point subtract(const point &p) const;
	double calcDistance(const point &p) const;
	std::vector<double> calcDistance(const std::vector<point> &P) const;
	double cross(const point &p) const;
	point add(const point &p) const;
	bool isEqual(const point &p) const;
};


const point origin(0,0,15);

//float point::updateWeight(const double SC){
//	if (round(SC_max/SC) < 2.5) return 0;
//	return round(SC_max/SC);
//}

point point::subtract(const point &p) const{
	return point(x-p.x, y-p.y);
}

double point::calcDistance(const point &p) const{
	return std::hypot(x-p.x, y-p.y);
}

double point::cross(const point &p) const{
	return x* p.y - y*p.x;
}

point point::add(const point &p) const{
	return point(x+p.x, y+p.y);
}

std::vector<double> point::calcDistance(const std::vector<point> &P) const{
	std::vector<double> d;
	d.reserve(P.size());
	for(const point &thisP : P){
		d.push_back(calcDistance(thisP));
	}
	return d;
}

bool point::isEqual(const point &p) const{
	if(x == p.x && y == p.y)
		return true;
	return false;
}

struct circle {
public:
	point O;
	double r;

	circle() : O(), r(-1) {};
	circle(point p, double a) : O(p), r(a) {};

	//Member functions

	bool contains(point &p) const;
	bool contains(std::vector<point> &p);
};

bool circle::contains(point &p) const{
	return O.calcDistance(p) <= r;
};

bool circle::contains(std::vector<point> &p){
	for(point &P : p){
		if(!contains(P)) return false;
	}
	return true;
}

double get_E(){
	return distribution(generator);
};

circle intilizeC(point &p, double a) {
	return circle(p,a);
};


#endif /* MINCIRCLE_HPP_ */
