//#pragma once
//
//#pragma once
//
//#include <iostream>
//#include<cmath>
//#include<vector>
//#include<string>
//#include <time.h>
//
////custom header files
////#include"rejection_sampling.h"
//#include "MathFunctions.h"
//#include "ran2.h"
//#include "RejectionSampling.h"
////#include "polymer_monomer.h"
//
////const double pi = 3.14159;
////const double a = 1.0;
//
//
////flexible chain, fixed length
//double yamakawa(double& distance, int& number_of_segments) {//could uninclude monomer separation if the chain length is given in terms of it
//
//	double sum_value{ 0 };
//	int n{ 0 };
//	for (int k{ 0 }; k <= (number_of_segments - (distance / a)) / 2; k++) {
//		sum_value = sum_value + pow(-1, k) * choose(number_of_segments, k) * pow(number_of_segments - 2.0 * k - (distance / a), number_of_segments - 2);
//	}
//
//	return sum_value / (pow(2, number_of_segments + 1) * factorial(number_of_segments - 2) * pi * pow(a, 2) * distance);
//
//}
//
////yamakawa max
//double yamakawa_max(double& initial_end_end_distance, int& segments) {
//
//	double r_min_after_jump{ initial_end_end_distance - a };
//
//	return yamakawa(r_min_after_jump, segments);
//
//}
//
//double rand2(double min_value, double max_value) {
//
//	return double(((max_value - min_value) * rand()) / RAND_MAX) + min_value;
//}
//
//std::vector<double> random_direction(double& radius, bool& cartesian) {
//	double theta_rand{ acos(1 - 2 * rand2(0, 1)) }, phi_rand{ 2 * pi * rand2(0,1) };
//	std::vector<double> spherical_coordinates{ radius,theta_rand,phi_rand };
//
//	if (cartesian) {
//		return spherical_to_cartesian(spherical_coordinates);
//	}
//	else {
//		return { radius, theta_rand, phi_rand };
//	}
//}
//
//std::vector<double> random_direction(bool& cartesian) {
//	double theta_rand{ acos(1 - 2 * rand2(0, 1)) }, phi_rand{ 2 * pi * rand2(0,1) };
//	std::vector<double> spherical_coordinates{ 1,theta_rand,phi_rand };
//
//	if (cartesian) {
//		return spherical_to_cartesian(spherical_coordinates);
//	}
//	else {
//		return { 1, theta_rand, phi_rand };
//	}
//}
//std::vector<double> rejection_sample(double initial_end_distance, std::vector<double> initial_position, std::vector<double>& monomer_N_position, int& number_of_segments, std::vector<int>& rejections_vector) {
//	bool cart{ true };
//	std::vector<double> jump{ random_direction(cart) };
//	double new_end_distance{ vector_modulus(vector_addition(initial_position,jump)) };
//
//	while (new_end_distance > initial_end_distance) {
//		std::cout << "rejection 1" << std::endl;
//		return rejection_sample(initial_end_distance, initial_position, monomer_N_position, number_of_segments, rejections_vector);
//	}
//
//	//YAMAKAWA MAX
//	double initial_r12{ vector_modulus(vector_subtraction(monomer_N_position,initial_position)) };
//
//
//	double pdf_max{ yamakawa_max(initial_r12,number_of_segments) };
//
//	double random_num_Y{ rand2(0,pdf_max) };
//
//	while (random_num_Y > yamakawa(new_end_distance, number_of_segments)) {
//		std::cout << "rejection2" << std::endl;
//		return rejection_sample(initial_end_distance, initial_position, monomer_N_position, number_of_segments, rejections_vector);
//	}
//
//	return vector_addition(initial_position, jump);
//}
