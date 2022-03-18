// PolymerMonomerObject Tests.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include "polymer class.h"
//#include "maths_functions.h"
//#include "maths_functions.h"

//need to print out monomers in a structure (functiond definition is already there. 







int main()
{

    //polymer p({ 0,0,0 }, { 0,0,9 });

    //p.add_monomer({ 0,0,3 });
    //p.add_monomer({ 0,0,4 });
    //p.add_monomer({ 0,0,5 });

    //p.link_monomers(0, 3);//1 structure
    //p.number_of_structures();
    //std::cout<<p(0)->structure_length() << "\n";
    //p.link_monomers(1, 4);//2 structures
    //p.number_of_structures();
    //std::cout<<p(1)->structure_length()<<"\n";

    //p.link_monomers(2, 4);//2 structures, previous one is broken
    //p.number_of_structures();
    //std::cout << p(1)->structure_length() << "\n";

    //p.structures_monomer_is_founding_member_of(p[0]);
    ////std::cout<<p(2)->structure_length();


    polymer kissing_hairpin({ 0,0,0 }, { 9,0,0 });

    kissing_hairpin.add_monomer({ 1,0,0 },"G");
    kissing_hairpin.add_monomer({ 2,1,0 },"C");
    kissing_hairpin.add_monomer({ 1,2,0 },"A");
    kissing_hairpin.add_monomer({ 2,3,0 },"A");
    kissing_hairpin.add_monomer({ 3,3,0 },"G");
    kissing_hairpin.add_monomer({ 4,2,0 },"C");
    kissing_hairpin.add_monomer({ 3,1,0 },"A");
    kissing_hairpin.add_monomer({ 4,0,0 },"A");
    kissing_hairpin.add_monomer({ 5,0,0 },"U");
    kissing_hairpin.add_monomer({ 6,1,0 },"G");
    kissing_hairpin.add_monomer({ 5,2,0 },"C");
    kissing_hairpin.add_monomer({ 6,3,0 },"A");
    //kissing_hairpin.add_monomer({ 7,3,0 },"C");
    //kissing_hairpin.add_monomer({ 8,2,0 },"G");
    //kissing_hairpin.add_monomer({ 7,1,0 },"A");
    //kissing_hairpin.add_monomer({ 8,0,0 },"U");
    //kissing_hairpin.add_monomer({ 9,0,0 },"G");

    std::vector < std::vector<int>> c_matrix;
    c_matrix = kissing_hairpin.compatibility_matrix();

    for (int i{ 0 }; i < c_matrix.size(); i++) {
        for (int k{ 0 }; k < c_matrix[i].size(); k++) {
            std::cout << c_matrix[i][k] << " ";
        }
        std::cout << "\n";
    }

    //std::vector<int> v;
    /*v = kissing_hairpin.compatible_region(c_matrix);
    for (auto i : v) {
        std::cout << i << "\n";
    }*/

    std::vector<int> compatible_monomers;
    compatible_monomers = kissing_hairpin.compatible_region(c_matrix);
    for (auto i : compatible_monomers) {
        std::cout << i << " ";
    }

    kissing_hairpin.print_bases();
    //1st test of kissing hairpin
    //kissing_hairpin.link_monomers(2, 7);
    //kissing_hairpin.link_monomers(6, 11);
    //kissing_hairpin.link_monomers(10, 15);

    //2nd test of kissing hairpin
    //kissing_hairpin.link_monomers(1, 4);
    //kissing_hairpin.link_monomers(3, 6);
    //kissing_hairpin.link_monomers(5, 8);
    //kissing_hairpin[1]->linked_to()->get_position();

    kissing_hairpin.helix_vectors(compatible_monomers);
    //kissing_hairpin[6]->get_position();
    
    //kissing_hairpin.structure_type(kissing_hairpin(0));
    //kissing_hairpin(0)->print_structure_type();
    //kissing_hairpin.structure_type(kissing_hairpin(1));
    //kissing_hairpin(1)->print_structure_type();
    //kissing_hairpin.structure_type(kissing_hairpin(2));
    //kissing_hairpin(2)->print_structure_type();


    //kissing_hairpin.output_for_ovito(kissing_hairpin.generate_vector_of_monomer_positions());

    std::vector<double> v{ {0,0,3/ sqrt(6)} };
    std::vector<double> u{ {0.75/ sqrt(2),0,0} };
    std::vector<double> x{ 0,0,0 };
    std::vector<std::vector<double>> s_x;
    std::vector<std::vector<double>> s_y;

    kissing_hairpin.generate(3, v, u, x, s_x, s_y);
    kissing_hairpin.print_for_ovito(s_x, s_y);

    std::vector<double> vvv;
    vvv = kissing_hairpin.rotate_helix(u, v, 0.8934184156629714);

}

