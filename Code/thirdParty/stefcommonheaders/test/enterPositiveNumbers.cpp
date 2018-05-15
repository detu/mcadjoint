#include <iostream>
#include <stefCommonHeaders/checkedNumber.hpp>


stefCommonHeaders::PositiveNumber<double> checkIfCopyConstructorWorking( stefCommonHeaders::PositiveNumber<double> num) {
    return num + 1;
}

int main() {
    stefCommonHeaders::PositiveNumber<double> num = 1;
    std::cout << "Enter positive numbers\n";

    for (;;) {
        std::cin >> num;
        if (!std::cin.good()) {
            break;
        }
        std::cout << "Entered " << num << "\n";
    }

    return 0;
}
