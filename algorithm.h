#pragma once

#include <iostream>
#include<cmath>
#include<vector>
#include<string>
#include <time.h>
#include "math_functions.h"
#include "main.h"


std::vector<std::vector<double>> polymer_generator(std::vector<double> start_position, std::vector<double> end_position, int number_of_segments) {
	int n{ number_of_segments };
	std::vector<std::vector<double>> polymer;
	std::vector<double> new_position;
	std::vector<double> old_position{ start_position };

	int number_of_monomers{ n + 1 };

	for (int i{ 0 }; i < n; i++) {

		//std::cout << "monomer " << i << "\n";
		new_position = rejection_sample(old_position, end_position, number_of_monomers);
		polymer.push_back(new_position);

		//for (auto i : new_position) {
		//	std::cout << "INSERTING INTO POLYMER " << i << "\n";
		//}

		old_position = new_position;
		if (i > 0) {
			double distance;
			std::vector<double> distance_vector{ vector_subtraction(polymer[i], polymer[i - 1]) };
			distance = abs(vector_modulus(distance_vector));

			std::cout << "difference between final position vectors; " << distance << "\n";

		}
	}

	for (int i = 0; i < polymer.size(); i++)
	{
		std::cout << "{";
		for (int j = 0; j < polymer[i].size(); j++)
		{
			std::cout << polymer[i][j] << " ";
		}
		std::cout << "}\n";

	}
	std::cout << "{";
	for (auto i : end_position) {
		std::cout << i << " ";
	}
	std::cout << "}" << std::endl;
	return polymer;
}
