#pragma once
#include "monomer class.h"
//#include "polymer class.h"
#include <vector>

const std::string structure_type[2] = { "Hairpin", "Kissing hairpin" };

class structure{
    //friend class polymer;
private:
    monomers some_structure;
    int id{ 0 };
    std::string type{ structure_type[0]};
public:
    structure() = default;

    //structure(monomer* monomer_in_structure, int structure_id) {
    //    id = structure_id;
    //    some_structure.push_back(monomer_in_structure);
    //}
    void add_monomer(monomer* some_monomer);

    int structure_length() {
        //std::cout << "There are " << some_structure.size() << " in this structure.\n";
        return some_structure.size();
    }
    
    void set_structure_type(std::string s_type) { type = s_type; }
    void print_structure_type() { std::cout << "This structure is a " << type << "\n"; }

    monomer* operator[](int i);

    monomer* return_monomer(int i) {
        return some_structure[i];
    }

    ~structure() { std::cout << "Structure destructor called.\n"; }
};

