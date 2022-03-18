#pragma once

#include <iostream>
#include<cmath>
#include<vector>
#include<string>
#include <time.h>


//custom header files
//#include"rejection_sampling.h"
#include "MathFunctions.h"
#include "ran2.h"
//#include "polymer_monomer.h"

//const double pi = 3.14159;
//const double a = 1.0;


//flexible chain, fixed length
double yamakawa(double& distance, int& number_of_segments) {//could uninclude monomer separation if the chain length is given in terms of it

	double sum_value{ 0 };
	int n{ 0 };
	for (int k{ 0 }; k <= (number_of_segments - (distance / a)) / 2; k++) {
		sum_value = sum_value + pow(-1, k) * choose(number_of_segments, k) * pow(number_of_segments - 2.0 * k - (distance / a), number_of_segments - 2);
	}

	return sum_value / (pow(2, number_of_segments + 1) * factorial(number_of_segments - 2) * pi * pow(a, 2) * distance);

}

//yamakawa max
double yamakawa_max(double& initial_end_end_distance, int& segments) {

	double r_min_after_jump{ initial_end_end_distance - a };

	return yamakawa(r_min_after_jump, segments);

}

double rand2(double min_value, double max_value) {

	return double(((max_value - min_value) * rand()) / RAND_MAX) + min_value;
}

std::vector<double> random_direction(double& radius, bool& cartesian) {
	double theta_rand{ acos(1 - 2 * rand2(0, 1)) }, phi_rand{ 2 * pi * rand2(0,1) };
	std::vector<double> spherical_coordinates{ radius,theta_rand,phi_rand };

	if (cartesian) {
		return spherical_to_cartesian(spherical_coordinates);
	}
	else {
		return { radius, theta_rand, phi_rand };
	}
}

std::vector<double> sign_of_jump(std::vector<double> previous_position, std::vector<double> Nth_position, std::vector<double> new_position){
	std::vector<double> R{ (vector_subtraction(Nth_position,previous_position)) };
	double mod_R{ vector_modulus(vector_subtraction(Nth_position,previous_position)) };
	for (int i{ 0 }; i < 3; i++) {
		if (abs(Nth_position[i] - new_position[i]) > abs(R[i])) {
			new_position[i] = -new_position[i];
		}

	}

	return new_position;

}

std::vector<double> generate_new_position(double accepted_r12, std::vector<double> previous_position, std::vector<double> Nth_position) {
	double r12{ accepted_r12 };

	//vector from {i-1}th monomer to Nth monomer.
	double mod_R{ vector_modulus(vector_subtraction(Nth_position,previous_position)) };

	std::vector<double> R{ (vector_subtraction(Nth_position,previous_position)) };

	//calculate angle between new position and R. 
	double alpha{ acos(((a * a) - (r12 * r12) + (mod_R * mod_R)) / (2 * a * mod_R)) };

	//generate random theta component of new position in spherical polar coordinates.
	double theta{ rand2(0,2 * atan(1) * 4) };

	//parameters of our formula 
	double a{ R[1] * sin(theta) };
	double b{ R[0] * sin(theta) };
	double c{ mod_R * cos(alpha) - R[2] * cos(theta) };
	//need to check that the discriminant is positive

	while ((a * a) + (b * b) - (c * c) < 0) {
		theta = rand2(0,2 * atan(1) * 4);
		a = R[1] * sin(theta);
		b = R[0] * sin(theta);
		c = mod_R * cos(alpha) - R[2] * cos(theta);

	}
	double d{ sqrt((a * a) + (b * b) - (c * c)) };


	//the formula. gives 2 possible values we just calculate one.
	double phi{ atan((b * d + a * c) / (b * c - a * d)) };

	//new position vector
	std::vector<double> new_position{ cos(phi) * sin(theta),sin(phi) * sin(theta),cos(theta) };
	std::cout << "position" << std::endl;
	for (auto i : new_position) {
		std::cout << i << std::endl;
	}


	new_position = vector_addition(new_position, previous_position);
	std::cout << "position" << std::endl;
	for (auto i : new_position) {
		std::cout << i << std::endl;
	}
	new_position = sign_of_jump(previous_position, Nth_position, new_position);
	std::cout << "position" << std::endl;
	for (auto i : new_position) {
		std::cout << i << std::endl;
	}

	return new_position;
}

std::vector<double> random_walk_rejection_sample2(std::vector<double> old_position, std::vector<double>& monomer_N_position, int& number_of_segments, std::vector<int>& rejections_vector) {

	//YAMAKAWA MAX
	double initial_r12{ vector_modulus(vector_subtraction(monomer_N_position,old_position)) };
	double jump{ rand2(0,a) };

	double trial_r12{ initial_r12 - jump };

	double pdf_max{ yamakawa_max(initial_r12,number_of_segments) };

	double random_num_Y{ rand2(0,pdf_max) };

	while (random_num_Y > yamakawa(trial_r12,number_of_segments)) {
		//std::cout << trial_r12 << std::endl;
		return random_walk_rejection_sample2(old_position, monomer_N_position, number_of_segments, rejections_vector);
	}

	double rand_phi{ rand2(0,2 * atan(1) * 4) };
	std::cout << "initial r12 " << initial_r12 << std::endl;

	double theta{ acos(((a * a) - (trial_r12 * trial_r12) + (initial_r12 * initial_r12)) / (2 * a * initial_r12)) };
	std::cout << "theta " << theta<<std::endl;
	//std::vector<double> new_position{ {a * cos(rand_phi) * sin(theta),a * sin(rand_phi) * sin(theta),a * cos(theta)} };
	//std::vector<double> new_position{ {a * cos(rand_phi),a * sin(rand_phi),-1 * a * cos(theta)} };
	std::vector<double> new_position{ generate_new_position(trial_r12,old_position,monomer_N_position) };

	//double compare_r12{ vector_modulus(vector_subtraction(monomer_N_position,new_position)) };
	double compare_r12{ vector_modulus(vector_subtraction(monomer_N_position,new_position)) };



	std::cout << "Generated r12: " << trial_r12 << std::endl;
	std::cout << "Compare with generated position: " << compare_r12 << std::endl;

	return new_position;
	//random number between 0 and a. We add that to the R12
	//random number between R12 for monomer {i-1} and R12 corresponding to moving a distance of a towards the end monomer (ie the position of the i-1 monomer + distance a towards end monomer). the random number is what we'll rejection sample. =X
	//generate random number between 0 and Yamakawa max. =Y

	//

	//R12 is the distance between the endpoint and the ith monomer. 

	//if the random R12 we picked is accepted then we transform into a direction using polar coordinates. or maybe not. TBD

}



std::vector<double> random_walk_rejection_sample(std::vector<double> old_position, std::vector<double>& monomer_N_position, int monomer_id, int& total_number_of_monomers, std::vector<int>& rejections_vector) {

	int N{ total_number_of_monomers };
	int i{ monomer_id };
	double r{ a };
	bool cartesian_coords{ true };
	double pdf_max{ 0 };

	static int rejections{ 0 };
	rejections++;

	if (rejections > 3000) {
		std::cout << "stack overflow";
		exit(EXIT_FAILURE);
	}
	else {
		//std::vector<int> *rejections_per_iteration{ rejections_vector };

		//random jump direction
		std::vector<double> jump_direction{ random_direction(r,cartesian_coords) };
		//std::cout << "random direction {x,y,z}" << "\n";
		//for (auto i : jump_direction) {
		//	//std::cout << "random direction\n";
		//	std::cout << i << ' ' << "\n";

		//}

		//end-to-end distance before jump
		std::vector<double> end_to_end_vector{ vector_subtraction(old_position,monomer_N_position) };
		double mod_end_to_end_vector{ vector_modulus(end_to_end_vector) };
		double initial_end_to_end{ mod_end_to_end_vector };

		//maximum value of the yamakawa function
		pdf_max = yamakawa_max(initial_end_to_end, N);


		//new randomly generated position
		std::vector<double> new_position{ vector_addition(old_position,jump_direction) };

		//end-to-end distance after jump
		std::vector<double> new_end_to_end_vector{ vector_subtraction(new_position,monomer_N_position) };
		double r_12{ vector_modulus(new_end_to_end_vector) };
		std::cout << "r12: " << r_12 << "\n" << std::endl;

		//if the distance
		if (r_12 == 0.0) {
			std::cout << "rejected: exception\n" << std::endl;
			return random_walk_rejection_sample(old_position, monomer_N_position, i, N, rejections_vector);
		}
		//this condition will guarantee that the polymer will not be overstretched at any point. 
		else if (r_12 > (total_number_of_monomers * 1.0 - 1.0) * a) {
			std::cout << "rejected: would result in overstretching\n" << std::endl;
			return random_walk_rejection_sample(old_position, monomer_N_position, i, N, rejections_vector);

		}
		else {

			//std::cout << "yamakawa_max: " << pdf_max << "\n";



			//r_min and r_max
			double r_max{ distance_to_end(old_position, monomer_N_position) + a }, r_min{ distance_to_end(old_position, monomer_N_position) };


			//X is the end_to_end distance. Lies between mod(r)+a and mod(r)-a.
			double Y{ rand2(0,pdf_max) };

			std::cout << "old position" << "\n";
			for (auto i : old_position) {

				std::cout << i << "\n";
			}
			int number_of_segments{ N - i };
			if (Y <= yamakawa(r_12, number_of_segments)) {

				std::cout << "new position {x,y,z}" << "\n";

				for (auto i : new_position) {

					std::cout << i << "\n";
				}

				std::cout << "r12: " << r_12 << "\n" << std::endl;

				std::cout << "monomer number: " << monomer_id << "\n" << std::endl;
				std::cout << "total number of rejections" << rejections << "\n" << std::endl;
				//rejections_per_iteration.push_back(&rejections);
				rejections_vector.push_back(rejections);
				return new_position;
			}
			else {
				std::cout << "rejected: rejection sampling method\n";
				return random_walk_rejection_sample(old_position, monomer_N_position, i, N, rejections_vector);
			}

		}

	}
}



double test(std::vector<double> N, std::vector<double> z) {
	double jump{ 0.75 };


	std::vector<double> R{ vector_subtraction(N,z) };
	double mod_R{ vector_modulus(R) };

	double r{ mod_R - jump };

	double alpha{ acos((1 +(mod_R*mod_R)-(r*r)) / (2*mod_R)) };
	//std::cout << "Alpha " << alpha << std::endl;

	//random theta
	double theta{ 1 };

	//parameters of our formula 
	double l{ R[1] * sin(theta) };
	double b{ R[0] * sin(theta) };
	double c{ mod_R * cos(alpha) - R[2] * cos(theta) };
	double d{ sqrt((l*l) + (b * b) - (c * c)) };

	//formula
	double phi1{ atan((b * d + l * c) / (b * c - l * d)) };
	double phi2{ atan((b * d - l * c) / (b * c + l * d)) };

	std::vector<double> b1{ z };
	b1[0] = b1[0] + cos(phi1) * sin(theta);
	b1[1] = b1[1] + sin(phi1) * sin(theta);
	b1[2] = cos(theta);

	std::vector<double> b2{ z };
	b2[0] = b2[0] + cos(phi2) * sin(theta);
	b2[1] = b2[1] + sin(phi2) * sin(theta);
	b2[2] = cos(theta);


	double generated_r1{ vector_modulus(vector_subtraction(N,b1)) };
	double generated_r2{ vector_modulus(vector_subtraction(N,b2)) };

	std::cout << phi1 << std::endl;
	std::cout << phi2 << std::endl;
	std::cout << "Mod R " << mod_R << std::endl;
	std::cout << "Target r " << r << std::endl;
	std::cout << "1st Generated r"<<generated_r1<<std::endl;
	std::cout << "2nd Generated r"<<generated_r2<<std::endl;

	return 0;

}

double test2(std::vector<double> N, std::vector<double> z) {
	double jump{ 0.75 };


	std::vector<double> R{ vector_subtraction(N,z) };
	double mod_R{ vector_modulus(R) };

	double r{ mod_R - jump };

	double alpha{ acos((1 + (mod_R * mod_R) - (r * r)) / (2 * mod_R)) };
	//std::cout << "Alpha " << alpha << std::endl;

	//random theta
	double theta{ 1 };

	//parameters of our formula 
	double l{ R[1] * sin(theta) };
	double b{ R[0] * sin(theta) };
	double c{ mod_R * cos(alpha) - R[2] * cos(theta) };
	double d{ sqrt((l * l) + (b * b) - (c * c)) };

	//formula
	double phi1{ asin(c/sqrt((l*l)+(b*b))) };
	double phi2{ atan((b * d - l * c) / (b * c + l * d)) };

	std::vector<double> b1{ z };
	b1[0] = b1[0] + cos(phi1) * sin(theta);
	b1[1] = b1[1] + sin(phi1) * sin(theta);
	b1[2] = cos(theta);

	std::vector<double> b2{ z };
	b2[0] = b2[0] + cos(phi2) * sin(theta);
	b2[1] = b2[1] + sin(phi2) * sin(theta);
	b2[2] = cos(theta);


	double generated_r1{ vector_modulus(vector_subtraction(N,b1)) };
	double generated_r2{ vector_modulus(vector_subtraction(N,b2)) };

	std::cout << phi1 << std::endl;
	std::cout << phi2 << std::endl;
	std::cout << "Mod R " << mod_R << std::endl;
	std::cout << "Target r " << r << std::endl;
	std::cout << "1st Generated r" << generated_r1 << std::endl;
	std::cout << "2nd Generated r" << generated_r2 << std::endl;

	return 0;

}

