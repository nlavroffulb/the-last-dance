#pragma once

#include <iostream>
#include<cmath>
#include<vector>
#include<string>
#include <time.h>

//vector addition
std::vector<double> vector_addition(std::vector<double> vector1, std::vector<double> vector2) {

    return { vector1[0] + vector2[0],vector1[1] + vector2[1],vector1[2] + vector2[2] };
}
std::vector<double> vector_subtraction(std::vector<double> vector1, std::vector<double> vector2) {
    return { vector1[0] - vector2[0],vector1[1] - vector2[1],vector1[2] - vector2[2] };

}

//MODULUS
double vector_modulus(std::vector<double> vector1) {
    return sqrt(pow(abs(vector1[0]), 2) + pow(abs(vector1[1]), 2) + pow(abs(vector1[2]), 2));
}

//DOT PRODUCT
double dot_product(std::vector<double> vector1, std::vector<double> vector2) {
    return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2];
}
double angle_between_vectors(std::vector<double> vector1, std::vector<double> vector2) {

    return acos(dot_product(vector1, vector2) / (vector_modulus(vector1) * vector_modulus(vector2)));
}


int quadrant_of_a_2d_vector(std::vector<double> two_d_vector) {
    if (two_d_vector[0] == 0.0 || two_d_vector[1] == 0) { return 0; }
    if (two_d_vector[0] > 0.0 && two_d_vector[1] > 0.0) { return 1; }
    if (two_d_vector[0] > 0.0 && two_d_vector[1] < 0.0) { return 4; }
    if (two_d_vector[0] < 0.0 && two_d_vector[1] > 0.0) { return 2; }
    if (two_d_vector[0] < 0.0 && two_d_vector[1] < 0.0) { return 3; }
}


std::vector<double> normalize(std::vector<double> unnormalized_vector) {
    std::vector<double> unit_vector;
    for (int i{ 0 }; i < unnormalized_vector.size(); i++) {
        unit_vector.push_back(unnormalized_vector[i] / vector_modulus(unnormalized_vector));
    }

    return unit_vector;
}


double rand2(double min_value, double max_value) {

    return double(((max_value - min_value) * rand()) / RAND_MAX) + min_value;
}

double cosine_rule_angle(double a, double b, double c) {
    return acos(((a * a) + (b * b) - (c * c)) / (2 * a * b));
}