// Stage 1 ULB.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Polymer_Generator.h"
#include "ran2.h"
#include "rejection_sampling_new_try.h"
#include<algorithm>
int main()
{
    srand(time(NULL));
    
    //std::vector<int> rejections_per_iteration;

    //polymer_generator({ 0,0,0 }, { 6,0,0 }, 8, rejections_per_iteration);
    std::vector<double> N_pos{6.5,0,0 };
    int n{ 8 };
    //random_walk_rejection_sample2({0,0,0}, N_pos, n, rejections_per_iteration);
    

    //test2({ 2,5,6 }, { 1,1,1 });

    //rejection_sample(5, { 0,0,1 }, N_pos, n, rejections_per_iteration);

    ran2_state_t obj;
    ran2_set(&obj, rand());
    double rn{ 0 };

    //for (int i{ 0 }; i < 10; i++) {
    //    rn = ran2_get_double(&obj);
    //    std::cout << rn << std::endl;

    //}
    std::cout << generate_rn() << std::endl;

    //void* vstate{ 0 };
    //double y{ 1.0 };
    //double* x;
    //x = &y;

    //std::vector<double> random_numbers;


    //for (int i{ 0 }; i < 1000; i++) {
    //    double* z{ &y };
    //    random_numbers.push_back(ran2_get_double(z));
    //    //std::cout << ran2_get_double(z)<<std::endl;

    //}
    //double sum{ 0 };
    //for (auto i : random_numbers) {
    //    std::cout << i << std::endl;
    //    sum += i;
    //}
    //double average{ sum / random_numbers.size() };
    //std::cout << "The average is " << average << std::endl;
    //std::cout << ran2_get_double(x)<< std::endl;

    //delete x;


}

