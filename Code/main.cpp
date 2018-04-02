//
//  main.cpp
//  Adjoint-MC
//
//  Created by Patrick Jenny on 14/02/18.
//  Copyright © 2018 Patrick Jenny. All rights reserved.
//

#include <iostream>
#include "driver.h"

using namespace std;

int main(int argc, char** argv) {
    Driver d;
    d.solve_Burger();
//d.solve_NS    ();
    cout << "\ndone\n";
    return 0;
}