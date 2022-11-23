// DemoForGeneticAlgorithm.cpp: 定义应用程序的入口点。
//

#include "DemoForGeneticAlgorithm.h"
#include "GeneicAlgorithm.hpp"

#include <fstream>
#include <cmath>

using namespace std;

int main()
{
	GeneicAlgorithm ga(10);

	auto eval_function = [](std::vector<double> v) {
		double result = 0; std::for_each(v.begin(), v.end(), [&](double x) {result += pow(std::ceil(x+0.5),2); }); return result; };

	auto bounds_of_variables = std::vector<struct GeneicAlgorithm::Bound>{};
	for (auto i = 0; i < 10; i++)
		bounds_of_variables.push_back(GeneicAlgorithm::Bound { -100, 100, 0.01 });
	auto max_evolution_time = 500;
	auto population_size = 50;//较大比较好，收敛很快
	auto probability_crossover = 0.4;//不可以过大！否则很难收敛
	auto probability_mutation = 0.0001;//不可以过大！否则很难收敛
	auto print_statistics = true;

	auto [statistics_result,value_result, adaptibility_result]=
		ga.execute(eval_function,bounds_of_variables,
		max_evolution_time,population_size,probability_crossover,
			probability_mutation,print_statistics);

	std::cout << "\n result is: point (";
	for (auto val : value_result)
	{
		std::cout << val << "  ";
	}
	std::cout << ")" << std::endl << " adaptibility is " << adaptibility_result << std::endl;

	if (not print_statistics)
		return 0;

	std::fstream output_file("C:/Users/Lhuna/Downloads/1.csv");
	if (not output_file.is_open())
	{
		std::cerr << "output file does not exist!" << std::endl;
		return -1;
	}
	for (auto item : statistics_result)
	{
		output_file << item << std::endl;
	}

	return 1;
}
