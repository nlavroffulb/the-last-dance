#include"structure class.h"
#include <iostream>

void structure::add_monomer(monomer* some_monomer) {
    some_structure.push_back(some_monomer);
}

monomer* structure::operator[](int i) {
    return some_structure[i];
}

