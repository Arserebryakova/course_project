#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "solver.hpp"
#include "problem.hpp"

int main(int argc, char** argv)
{
    Problem equation(argc, argv);
    EulerSolver solver;
    Result answer = solver.solve(equation);

    std::cout << std::fixed << std::setprecision(equation.accuracy());
    std::cout << "First moment = " << answer.N << std::endl;
    std::cout << "C(0) = " << answer.getC0() << std::endl;

    std::ofstream out(equation.path());
    out << std::fixed << std::setprecision(equation.accuracy());
    for (size_t i = 0; i < equation.nodes(); ++i)
        out << equation.origin() + i * equation.step() << " " << answer.C[i] << std::endl;
    return 0;
}
