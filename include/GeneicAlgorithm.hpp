#include "DemoForGeneticAlgorithm.h"

#include <exception>
#include <vector>
#include <cmath>
#include <string>
#include <functional>
#include <algorithm>
#include <random>

class BaseAlgorithm
{
public:
	BaseAlgorithm();
	BaseAlgorithm(BaseAlgorithm&&) = delete;
	BaseAlgorithm(const BaseAlgorithm&) = delete;
};

class GeneicAlgorithm
{
public:
	GeneicAlgorithm(int i) {};
	bool tester()
	{
		auto bounds = std::vector<struct GeneicAlgorithm::Bound>{ {-5,5,0.0001}, { -5,5,0.0001 } };
		auto result = encoder(2, std::vector<double>{-1.6707,-4.5202},bounds);
		auto result2 = decoder(2, bounds, result);
		eval_function = [](std::vector<double> v) {
			double result = 0.1; std::for_each(v.begin(), v.end(), [&](double x) {result *= x; }); return result; };
		std::vector<GeneicAlgorithm::Chromosome> chromosome{
			{"11",2,std::vector<double>{1,1},std::vector<size_t>()},
			{"22",2,std::vector<double>{2,2},std::vector<size_t>()},
			{"33",2,std::vector<double>{3,3},std::vector<size_t>()},
		};
		decltype(chromosome) result3 = roulette_selection(chromosome);
		return true;
	}

public:
	//struct and complex types

	//interval args
	struct Bound
	{
		double LowerBound;
		double UpperBound;
		double Precision;
	};

private:
	//struct and complex types
	
	//Chromosome
	class Chromosome
	{
	public:
		std::string str;//chromosome in string form
		unsigned int dimension;//variable number
		std::vector<double> values;//variable values
		std::vector<size_t> sub_length;//length of each variable in chromosome
	};


private:
	//methods

	GeneicAlgorithm::Chromosome encoder(unsigned int, const std::vector<double>&, const std::vector<struct GeneicAlgorithm::Bound>&);
	std::vector<double> decoder(unsigned int, const std::vector<struct GeneicAlgorithm::Bound>&, const std::string&);
	std::vector<GeneicAlgorithm::Chromosome> roulette_selection(std::vector<GeneicAlgorithm::Chromosome>& x);
	

private:
	//variables and constants
	
	//encoding related
	std::vector<GeneicAlgorithm::Chromosome> chromosomes;//store all chromosome
	std::vector< GeneicAlgorithm::Bound> bounds_of_variables;//store all boundary of variables

	//evaluation related
	std::function<double(std::vector<double>)> eval_function;//evaluate function

protected:
	//utility methods

	std::vector<double> get_random_numbers(unsigned int N, struct GeneicAlgorithm::Bound={0,1,0.0000000001});

};


