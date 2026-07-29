//
// Created by LiuRuichen on 2026/7/29.
//

# include <iostream>
using namespace std;

void FizzBuzz(int max_number) {
    int current_number = 1;
    while (current_number <= max_number) {
        if (current_number % 15 == 0) {cout << "FizzBuzz\n";}
        else if (current_number % 3 == 0) {cout << "Fizz\n";}
        else if (current_number % 5 == 0) {cout << "Buzz\n";}
        else {cout << current_number << "\n";}
        current_number += 1;
    }

}


int main() {
    int max_number {100};
    FizzBuzz(max_number);
}
