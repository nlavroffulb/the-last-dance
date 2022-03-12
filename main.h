#pragma once
#include <iostream>
#include<cmath>
#include<vector>
#include<string>
#include <time.h>
const double d{ 1 };


std::vector<double> last_dance(std::vector<double> initial, std::vector<double> N_position, double new_end_to_end_jump) {
	double L{ vector_modulus(vector_subtraction(N_position,initial)) };
	double l{ L - new_end_to_end_jump };

	double alpha{ cosine_rule_angle(d,L,l) };

	double c{ d * cos(alpha) }, r{ d * sin(alpha) };

	//get a vector on the plan containing the circle
	//using the fact that dot product to the normal is 0.
	//normal is the vector from initial to N. 
	std::vector<double> normal{ vector_subtraction(N_position,initial) };
	std::vector<double> unit_normal{ normalize(normal) };

	//random angle on circle
	double phi{ rand2(0,1) * atan(1) * 4 };
	double x{ r * cos(phi) }, y{ r * sin(phi) };
	double z{ -(1 / unit_normal[2]) * (unit_normal[0] * x + unit_normal[1] * y) };
	
	//vector going from centre of circle to a point on the circle
	std::vector<double> r_vector{ x,y,z };
	//vector going from centre of coordinate system to centre of circle
	std::vector<double> c_vector{ c * unit_normal[0],c * unit_normal[1],c * unit_normal[2] };

	//vector of the ith monomer
	std::vector<double> new_position{ vector_addition(c_vector,r_vector) };
	new_position = vector_addition(new_position, initial);

	//resultant end_to_end distance
	double r12{ vector_modulus(vector_subtraction(N_position,new_position)) };

	std::cout << "Generated r12: " << r12 << std::endl;
	std::cout << "Actual r12 " << l << std::endl;

	return new_position;
}
