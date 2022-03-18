#include "polymer class.h"
#include<iostream>
#include <string>
#include <time.h>
#include<fstream>
#include<algorithm>
#include "maths_functions.h"

//#define PI;
polymer::polymer(int num) {
    N = num;

    std::vector<monomer*> M;
    for (int i{ 0 }; i < N; i++) {
        M.push_back(new monomer());
    }
    chain = M;
}

polymer::polymer(std::vector<double> beginning_position, std::vector<double> end_position) {
    N = 2;
    chain.push_back(new monomer(beginning_position,0,"U"));
    chain.push_back(new monomer(end_position,1,"G"));
}


//general functions of the polymer and given monomer in the chain
int polymer::chain_length() {
    std::cout << chain.size() << "\n";
    return chain.size();
}

void polymer::get_monomer_position(int id) {
    chain[id]->get_position();
}

void polymer::add_monomer(std::vector<double> position, std::string base_type) {
    auto iterator = chain.end() - 1;
    chain[chain.size() - 1]->id++;
    chain.insert(iterator, new monomer(position,chain.size()-1,base_type));
    
    N++;
}

std::vector<std::vector<double>> polymer::generate_vector_of_monomer_positions()
{
    std::vector<std::vector<double>> temp;
    for (int i{ 0 }; i < chain.size(); i++) {
        temp.push_back(chain[i]->get_position());
    }

    return temp;
}


monomer* polymer::operator[](int i) {
    return chain[i];
}

structure* polymer::operator()(int i)
{
    return vector_of_structures[i];
}

//structure related functions
int polymer::number_of_structures() {
    std::cout << "There are " << vector_of_structures.size() << " structures in this polymer.\n";
    return vector_of_structures.size();
}
structures polymer::structures_monomer_is_member_of(monomer* monomer_x)
{

    structures structures_list;

    for (int i{ 0 }; i < vector_of_structures.size(); i++) {
        for (int j{ 0 }; j < vector_of_structures[i]->structure_length(); j++) {
            if (vector_of_structures[i]->return_monomer(j) == monomer_x) {
                structures_list.push_back(vector_of_structures[i]);
            }
        }
    }
    //std::cout << "part of " << structures_list.size() << " structures\n";

    if (structures_list.size() > 0) {

        return structures_list;
    }
    else exit;

}

structures polymer::structures_monomer_is_founding_member_of(monomer* monomer_x)
{
    structures founding_member_structure;
    structures structures_list{ structures_monomer_is_member_of(monomer_x) };
    for (int i{ 0 }; i < structures_list.size(); i++) {
        if (structures_list[i]->return_monomer(0) == monomer_x || structures_list[i]->return_monomer(structures_list[i]->structure_length() - 1) == monomer_x) {
            founding_member_structure.push_back(structures_list[i]);
        }
    }
    //std::cout << "Founding member of " << founding_member_structure.size() << " structures.\n";
    return founding_member_structure;
}

int polymer::index_of_structure_in_vector(structure* structure_x) {
    for (int i{ 0 }; i < vector_of_structures.size(); i++) {
        if (vector_of_structures[i] == structure_x) {
            return i;
        }
    }
    return -1;
}

//*********************************************************************

void polymer::define_structure(int start_monomer, int end_monomer) {
    vector_of_structures.push_back(new structure());

    for (int i{ start_monomer }; i <= end_monomer; i++) {
        vector_of_structures[vector_of_structures.size() - 1]->add_monomer(chain[i]);
    }

}
void polymer::print_monomers_in_structure(int structure_id) {
}

void polymer::print_structure(int structure_id)
{
}

//need to remove
void polymer::link_monomers(int monomer_i, int monomer_j) {

    int i{ monomer_i}, j{ monomer_j};
    
    //need to erase the structures that monomers i & j are 'founding members' of i.e they are the linked monomers that enclose the rest of the structure.
    for (int k{ 0 }; k < structures_monomer_is_founding_member_of(chain[i]).size(); k++) {
        try {
            int structure_index{ index_of_structure_in_vector(structures_monomer_is_founding_member_of(chain[i])[k]) };
            vector_of_structures.erase(vector_of_structures.begin() + structure_index - 1);
            std::cout << "Structure DESTROYED\n";

        }
        catch(...){std::cout << "Vector subscript out of bounds.\n"; }
    }
    for (int m{ 0 }; m < structures_monomer_is_founding_member_of(chain[j]).size(); m++) {
        try {
            int structure_index{ index_of_structure_in_vector(structures_monomer_is_founding_member_of(chain[j])[m]) };
            vector_of_structures.erase(vector_of_structures.begin() + structure_index - 1);
            std::cout << "Structure DESTROYED\n";

        }
        catch (...) {std::cout << "Vector subscript out of bounds.\n";}


    }



    //cannot link adjacent monomers
    if (abs(i - j) > 1) {

        //if either of the monomers were linked to another monomer we want that other monomer to be unlinked (point to nullptr). 
        if (chain[i]->linked_to() != nullptr) {
            chain[i]->linked_to()->link_to(nullptr);
        }
        else if (chain[j]->linked_to() != nullptr) {
            chain[j]->linked_to()->link_to(nullptr);
        }
        chain[i]->link_to(chain[j]);



        //also need to delete the old structure (if it existed)

        define_structure(monomer_i, monomer_j);

    }
    else { std::cout << "Cannot link those monomers. They are adjacent.\n" << std::endl; }
}



void polymer::structure_type(structure* structure_x)
{
    //by default, the structures we create are hairpins. will encode more as we go along but for now we don't need to code specifically for hairpins.

    //do need to identify kissing hairpins. these join two hairpins.
    
    //want to loop through all the monomers in our potential kissing hairpin structure. if there are two monomers within the structure (so not the beginning
    // or final monomer of the potential kissing hairpin) that are linked (NOT to each other) then we have a kissing hairpin
    bool condition_1{ false }, condition_2{ false };
    
    for (int i{ 0 }; i < structure_x->structure_length(); i++) {
        if (structure_x->return_monomer(i)->linked_to() != nullptr) {
            if (structure_x->return_monomer(i)->linked_to()->id < structure_x->return_monomer(i)->id) {
                condition_1 = true;
            }
            if (condition_1 == true && structure_x->return_monomer(i)->linked_to()->id > structure_x->return_monomer(i)->id) {
                condition_2 = true;
            }

        }
    }

    if (condition_1 == true && condition_2 == true) {
        structure_x->set_structure_type("Kissing hairpin");
    }
}

//*****************************************************************************///
//*****************************************************************************///
//*****************************************************************************///
void polymer::output_for_ovito(std::vector<std::vector<double>> polymer)
{
    std::string file_name{ "polymer " + std::to_string(rand()) + ".txt" };
    std::ofstream ovito_output(file_name);

    ovito_output << polymer.size() << "\n\n";
    for (int i{ 0 }; i < polymer.size(); i++) {
        std::string line{ "" };
        if (i % 2 == 0) {
            line += "C";
        }
        else if (i % 2 != 0) {
            line += "O";
        }
        for (int j{ 0 }; j < polymer[i].size(); j++) {
            line += " " + std::to_string(polymer[i][j]);
        }
        if (i == polymer.size() - 1) {
            //std::cout << "this line ran";
            ovito_output << line << std::endl;;
        }
        else {
            ovito_output << line << "\n";
        }
    }

    ovito_output.close();
}
void polymer::print_bases()
{
    std::cout << std::endl;
    for (int i{ 0 }; i < chain.size(); i++) {
        std::cout << chain[i]->get_base() << " - ";
    }
    std::cout << std::endl;
}
//*****************************************************************************///
//*****************************************************************************///
//*****************************************************************************///

bool polymer::compatible_bases(monomer* monomer_i,monomer* monomer_j) {
    if (monomer_i->get_base() == "U" && monomer_j->get_base() == "A" || monomer_i->get_base() == "A" && monomer_j->get_base() == "U") { return true; }
    else if (monomer_i->get_base() == "G" && monomer_j->get_base() == "C" || monomer_i->get_base() == "C" && monomer_j->get_base() == "G") { return true; }
    else { return false; }
}

std::vector<std::vector<int>> polymer::compatibility_matrix()
{
    std::vector<std::vector<int>> v;

    for (int i{ 0 }; i < chain.size(); i++) {
        std::vector<int> monomer_compatabilities;

        for (int j{ 0 }; j < chain.size(); j++) {

            if (compatible_bases(chain[i], chain[j]) == true) {
                monomer_compatabilities.push_back(1);
            }
            else { monomer_compatabilities.push_back(0); }

        }
        v.push_back(monomer_compatabilities);
    }
    return v;
}



std::vector<int> polymer::compatible_region(std::vector<std::vector<int>> matrix) {


    std::vector<std::vector<int>> compatible_regions;
    std::vector<std::vector<int>> vector_of_diagonals;

    //row diagonals
    for (int i{ 5 }; i < matrix.size(); i++) {
        std::vector<int> diagonal;

        for (int k{ 0 }; k <= i; k++) {
            diagonal.push_back(matrix[i - k][k]);

        }
        vector_of_diagonals.push_back(diagonal);

    }

    //column diagonals
    for (int j{ 1 }; j < matrix.size() - 5; j++) {
        std::vector<int> diagonal;

        for (int k{ 0 }; k < matrix.size() - j; k++) {
            diagonal.push_back(matrix[matrix.size() - 1 - k][j + k]);

        }
        vector_of_diagonals.push_back(diagonal);
    }

    //find compatible region. 3 or more consecutive 1s
    for (int i{ 0 }; i < vector_of_diagonals.size(); i++) {
        int counter{ 0 };
        std::vector<int> temp;

        for (int j{ 0 }; j < vector_of_diagonals[i].size(); j++) {
            if (j + 1 < vector_of_diagonals[i].size()) {
                if (vector_of_diagonals[i][j] == 1 && vector_of_diagonals[i][j + 1] == 1) {
                    temp.push_back(5 + i - j);
                    temp.push_back(j);

                    counter++;
                }

            }

        }
        if (counter >= 3) {

            //sort highest to lowest and delete duplicates
            std::sort(temp.begin(), temp.end());
            temp.erase(std::unique(temp.begin(), temp.end()), temp.end());
            compatible_regions.push_back(temp);

        }

    }
    std::cout << "\n";
    for (int i{ 0 }; i < compatible_regions.size(); i++) {
        if (compatible_regions[i].size() == 6) {
            return compatible_regions[i];
        }
    }
}



void polymer::create_structure_from_compatible_region(std::vector<int> compatible_region)
{

}

void polymer::helix_vectors(std::vector<int> monomers_in_region)
{
    //v, helix direction
    std::vector<double> v;
    //u, versor joining the centre of the helix to the first base in the first base pair
    std::vector<double> u;

    v = vector_subtraction(chain[monomers_in_region[0]]->get_position(), chain[monomers_in_region[2]]->get_position());

    std::vector<double> monomer_N_position{ random_position_in_circle_around_point(chain[monomers_in_region[0]]->get_position(),v) };

    u = multiplication_by_scalar(0.5, (vector_subtraction(chain[monomers_in_region[0]]->get_position(), monomer_N_position)));
    
}

//major groove
const double Mg = 0.47;
//minor groove
const double mg = 1.08;
//diameter
const double dm = 2.5;

double the = acos((2 * dm * dm / 4 - mg * mg) / (2 * dm * dm / 4.));
//print('the=', the, the * 180 / np.pi)


//rotathe u around v of theta
std::vector<double> polymer::rotate_helix(std::vector<double> u, std::vector<double> v, double the) {
    double cthe = cos(the);
    double sthe = sin(the);
    double v_dot_u = dot_product(u,v);
    //v_vec_u = np.zeros(3)
    std::vector<double> v_vec_u;

    v_vec_u.push_back(v[1] * u[2] - v[2] * u[1]);
    v_vec_u.push_back(v[2] * u[0] - v[0] * u[2]);
    v_vec_u.push_back(v[0] * u[1] - v[1] * u[0]);

    std::cout << std::endl;
    std::vector<double> wx;
    wx = multiplication_by_scalar(cthe, u);
    wx = vector_addition(wx, multiplication_by_scalar(sthe, v_vec_u));
    wx = vector_addition(wx, multiplication_by_scalar((1.0 - cthe) * v_dot_u, v));


    return wx;

}
//
//generate an A - RNA helix of n bp
//v: helix direction
//x : euclidean position of the starting base pair backbone
//u : versor joining the center of the helix to x; uand v are orthogonal
void polymer::generate(int n, std::vector<double> v, std::vector<double> u, std::vector<double> x, std::vector<std::vector<double>>& strand_1, std::vector<std::vector<double>>& strand_2) {
    std::vector<std::vector<double>> strand_x(n);
    std::vector<std::vector<double>> strand_y(n);
    double theta = 0.75 * pi;
    double phi = 33 * pi / 180;
    double dv{ 0.28 }, radius{ 2.5 / 2 };

    std::vector<double> running_centre{ vector_subtraction(x,multiplication_by_scalar(radius,u)) };

    std::vector<double> w{ rotate_helix(u,v,theta) };

    std::vector<double> y = vector_addition(running_centre, multiplication_by_scalar(radius, w));

    strand_x[0] = x;
    strand_y[n-1] = y;
    for (auto i : strand_y[n-1]) {
        std::cout << i << " ";
    }

    int running_n{ 1 };
    while (running_n < n) {
        running_centre = vector_addition(running_centre, multiplication_by_scalar(0.28, v));
        std::vector<double> up = rotate_helix(u, v, phi);

        std::vector<double> wp = rotate_helix(w, v, phi);

        u = up;

        w = wp;
        std::vector<double> running_x{ vector_addition(running_centre,multiplication_by_scalar(radius,u)) };

        std::vector<double> running_y{ vector_addition(running_centre,multiplication_by_scalar(radius,w)) };

        strand_x[running_n] = running_x;
        strand_y[n - running_n - 1] = running_y;
        running_n++;

    }
    for (int m{ 0 }; m < strand_y.size(); m++) {
        for (int k{ 0 }; k < strand_y[0].size(); k++) {
            std::cout << strand_y[m][k] << " ";
        }
        std::cout << std::endl;
    }
    strand_1 = strand_x;
    strand_2 = strand_y;
}

void polymer::print_for_ovito(std::vector<std::vector<double>> strand_1, std::vector<std::vector<double>> strand_2) {
    std::string file_name{ "polymer " + std::to_string(rand()) + ".txt" };
    std::ofstream ovito_output(file_name);


    int i = 0;
    while (i < strand_1.size()) {
        std::cout << "C ";
        for (auto k : strand_1[i]) {
            std::cout << k<< " ";
        }
        std::cout << std::endl;
        i++;


    }
    int j{ 0 };
    while (j < strand_2.size()){
        std::cout << "C ";
        for (auto k : strand_2[j]) {
            std::cout << k << " ";
        }
        std::cout << std::endl;

        //print('C', s_y[i, 0], s_y[i, 1], s_y[i, 2]);
        j++;

    }
}