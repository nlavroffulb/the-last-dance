#pragma once
#pragma once

#include <iostream>
#include<cmath>
#include<vector>
#include<string>
#include <time.h>
const double pi = 3.14159;

double const a{ 1.0 };

const std::vector<double> x_hat{ 1.0,0,0 };
const std::vector<double> y_hat{ 0,1.0,0 };
const std::vector<double> z_hat{ 0,0,1.0 };




//********************************************************************************************//
//********************************************************************************************//
//********************************************************************************************//
int factorial(int n) {
	if (n == 0) return 1;
	else if (n > 50) {
		std::cout << "Overload from too big a factorial";
		exit(0);
	}
	else {
		return n * factorial(n - 1);
	}

}
//********************************************************************************************//
//********************************************************************************************//
//********************************************************************************************//

int choose(int n, int  k) {
	if (k == 0) return 1;
	return (n * choose(n - 1, k - 1)) / k;

}

// Create a vector of evenly spaced numbers.
std::vector<double> evenly_spaced_range(double min, double max, size_t N) {
	std::vector<double> range;
	double delta = (max - min) / double(N - 1);
	for (int i = 0; i < N; i++) {
		range.push_back(min + i * delta);
	}
	return range;
}

double distance_to_end(std::vector<double> position_i, std::vector<double> position_N) {
	return sqrt(pow(position_i[0] - position_N[0], 2) + pow(position_i[1] - position_N[1], 2) + pow(position_i[2] - position_N[2], 2));
}

//********************************************************************************************//
//********************************************************************************************//
//********************************************************************************************//

std::vector<double>spherical_to_cartesian(std::vector<double>& spherical_coordinates) {
	std::vector<double> cartesian_coordinates;

	cartesian_coordinates.push_back(spherical_coordinates[0] * cos(spherical_coordinates[2]) * sin(spherical_coordinates[1]));
	cartesian_coordinates.push_back(spherical_coordinates[0] * sin(spherical_coordinates[2]) * sin(spherical_coordinates[1]));
	cartesian_coordinates.push_back(spherical_coordinates[0] * cos(spherical_coordinates[1]));

	return cartesian_coordinates;
}

//********************************************************************************************//
//********************************************************************************************//
//********************************************************************************************//

double rand2(double min_value, double max_value) {

    return double(((max_value - min_value) * rand()) / RAND_MAX) + min_value;
}


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

//MULTIPLICATION BY SCALAR
std::vector<double> multiplication_by_scalar(double constant, std::vector<double> vector) {
    return { constant * vector[0],constant * vector[1],constant * vector[2] };
}
//void multiplication_by_scalar(double constant, std::vector<double> vector) {
//    vector =  {constant * vector[0],constant * vector[1],constant * vector[2] };
//}

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



void set_small_values_to_zero(double& element) {
    if (abs(element) < 1e-6) {
        element = 0;
    }
}

std::vector<double> clockwise_rotation_about_z_axis(std::vector<double> point, double angle) {
    double theta{ angle };
    std::vector<double> P{ point };
    std::vector<std::vector<double>> z_axis_rotation_matrix{ {cos(theta),-1.0 * sin(theta),0}, {sin(theta),cos(theta),0}, {0,0,1} };

    std::vector<double> P_prime{ 0,0,0 };

    for (int i{ 0 }; i < z_axis_rotation_matrix.size(); i++) {
        for (int j{ 0 }; j < P.size(); j++) {
            P_prime[i] += z_axis_rotation_matrix[i][j] * P[j];
            set_small_values_to_zero(P_prime[i]);
        }
    }
    for (auto i : P_prime) {
        std::cout << i << "\n";
    }
    return P_prime;

}

std::vector<double> clockwise_rotation_about_y_axis(std::vector<double> point, double angle) {
    double theta{ angle };
    std::vector<double> P{ point };

    std::vector<std::vector<double>> y_axis_rotation_matrix{ {cos(theta),0,sin(theta)}, {0,1,0}, {-sin(theta),0,cos(theta)} };

    std::vector<double> P_prime{ 0,0,0 };

    for (int i{ 0 }; i < y_axis_rotation_matrix.size(); i++) {
        for (int j{ 0 }; j < P.size(); j++) {
            P_prime[i] += y_axis_rotation_matrix[i][j] * P[j];
            set_small_values_to_zero(P_prime[i]);

        }
    }
    //std::cout << "p prime: \n";
    //for (auto i : P_prime) {
    //    std::cout << i << "\n";
    //}

    return P_prime;
}

std::vector<double> clockwise_rotation_about_x_axis(std::vector<double> point, double angle) {
    double theta{ angle };
    std::vector<double> P{ point };

    std::vector<std::vector<double>> x_axis_rotation_matrix{ {1,0,0}, {0,cos(theta),-sin(theta),0},{0,sin(theta),cos(theta),0} };

    std::vector<double> P_prime{ 0,0,0 };

    for (int i{ 0 }; i < x_axis_rotation_matrix.size(); i++) {
        for (int j{ 0 }; j < P.size(); j++) {
            P_prime[i] += x_axis_rotation_matrix[i][j] * P[j];
            set_small_values_to_zero(P_prime[i]);

        }
    }
    for (auto i : P_prime) {
        std::cout << i << "\n";
    }

    return P_prime;
}


std::vector<double> point_translation(std::vector<double> point, std::vector<double> new_origin) {
    return vector_addition(point, new_origin);
}

//angles between axes

// angle with x_axis for example: dot product (x_hat,projection(transformed_point)) / mod(transformed_point)

void clockwise_rotation_angle(std::vector<double> point, double& angle) {
    if (quadrant_of_a_2d_vector(point) == 2 || quadrant_of_a_2d_vector(point) == 3) {
        angle = -angle;
    }
    else if (quadrant_of_a_2d_vector(point) == 0) {
        if (point[0] > 1) {
            angle = 0;
        }
        else if (point[0] < 1) {
            angle = -pi;
        }
        else if (point[1] > 1) {
            angle = 3 * pi / 2;
        }
        else if (point[1] < 1) {
            angle = pi / 2;
        }
    }
}


double angle_with_xy_plane(std::vector<double> point) {
    std::vector<double> P{ point };

    //new p is projection
    P[2] = 0;
    std::vector<double> projection{ P };

    //angle between projection and x_axis
    double theta;
    theta = angle_between_vectors({ 0,1,0 }, projection);
    clockwise_rotation_angle({ P[0],P[1] }, theta);

    return theta;
}

double angle_with_yz_plane(std::vector<double> point) {
    std::vector<double> P{ point };

    //new p is projection
    P[0] = 0;
    std::vector<double> projection{ P };

    //angle between projection and x_axis
    double theta;
    theta = angle_between_vectors({ 0,0,1 }, projection);
    clockwise_rotation_angle({ P[1],P[2] }, theta);

    return theta;

}

std::vector<double> forward_rotation_angles_for_specific_transformation(std::vector<double> new_z_axis) {
    std::vector<double> rotation_angles;
    rotation_angles.push_back(angle_with_xy_plane(new_z_axis));

    new_z_axis = clockwise_rotation_about_z_axis(new_z_axis, angle_with_xy_plane(new_z_axis));
    rotation_angles.push_back(angle_with_yz_plane(new_z_axis));


    return rotation_angles;
}

std::vector<double> backward_rotation_angles_for_specific_transformation(std::vector<double> new_z_axis) {
    std::vector<double> forward_angles{ forward_rotation_angles_for_specific_transformation(new_z_axis) };
    return { -forward_angles[1],-forward_angles[0] };
}

std::vector<double> backward_double_rotation_and_translation_transformation(std::vector<double> point, std::vector<double> new_origin, std::vector<double> rotation_angles) {
    std::vector<double> P{ point };


    P = clockwise_rotation_about_y_axis(P, rotation_angles[0]);
    P = clockwise_rotation_about_z_axis(P, rotation_angles[1]);

    P = point_translation(P, new_origin);

    return P;
}


std::vector<double> point_on_cone(std::vector<double> penultimate_monomer_position, std::vector<double> ultimate_monomer_position) {
    std::vector<double> P1{ penultimate_monomer_position }, P2{ ultimate_monomer_position };

    double r12{ vector_modulus(vector_subtraction(P1,P2)) };

    std::vector<double> coordinates;

    double cone_radius{ sqrt(1 - pow(r12 / 2,2)) };
    double r{ cone_radius };

    return { r * cos(2 * pi * rand2(0,1)),r * sin(2 * pi * rand2(0,1)),r12 / 2 };

}

std::vector<double> new_position(std::vector<double> point, std::vector<double> new_origin, std::vector<double> new_axis) {
    std::vector<double> rotation_angles{ backward_rotation_angles_for_specific_transformation(new_axis) };

    std::vector<double> final_result{ backward_double_rotation_and_translation_transformation(point,new_origin,rotation_angles) };

    return final_result;

}

std::vector<double> random_position_in_circle_around_point(std::vector<double> translation_vector, std::vector<double> z_axis) {

    double t{ 2 * pi * rand2(0,1) };
    std::vector<double> random_position{ 1.75 * cos(t), 1.75 * sin(t),0 };
    for (auto i : random_position) {
        std::cout << i << " ";
    }
    std::cout << std::endl;

    std::vector<double> P1{ z_axis };
    std::vector<double> P2{ random_position };
    double distance{ vector_modulus(vector_subtraction(P1,P2)) };
    std::cout << "distance pre transformation: " << distance << std::endl;

    double theta{ angle_with_xy_plane(P1) };

    P1 = clockwise_rotation_about_z_axis(P1, theta);
    P2 = clockwise_rotation_about_z_axis(P2, theta);

    double distance_post_transf{ vector_modulus(vector_subtraction(P1,P2)) };
    std::cout << "distance post transformation: " << distance_post_transf << std::endl;

    double beta{ angle_with_yz_plane(P1) };

    std::cout << "angle with yz: " << beta << std::endl;

    P1 = clockwise_rotation_about_x_axis(P1, beta);
    P2 = clockwise_rotation_about_x_axis(P2, beta);

    double d1{ vector_modulus(vector_subtraction(P1,P2)) };
    std::cout << "distance post 2nd trans: " << d1 << std::endl;

    P2 = point_translation(P2, translation_vector);
    double d2{ vector_modulus(vector_subtraction(P2,translation_vector)) };
    std::cout << "final distance " << d2 << std::endl;

    return P2;
}
