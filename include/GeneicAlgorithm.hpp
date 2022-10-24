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
		//auto result = encoder(2, std::vector<double>{-1.6707,-4.5202},bounds);

		eval_function = [](std::vector<double> v) {
			double result = 0.01; std::for_each(v.begin(), v.end(), [&](double x) {result *= x; }); return result; };
		
		std::vector<GeneicAlgorithm::Chromosome> chromosome;
		for (auto i = 1.0; i < 5; i+=0.1)
		{
			chromosome.push_back(std::move(encoder(2, { (double)i,(double)i}, bounds)));
		}

		decltype(chromosome) result3 = roulette_selection(chromosome);

		auto ret_val=crossover(result3);
		//ret_val=mutation(result3);

		return true;
	}

	bool run();

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
		std::vector<struct Bound> bounds;//boundary of variable

	};


private:
	//methods

	GeneicAlgorithm::Chromosome encoder(unsigned int, const std::vector<double>&, const std::vector<struct GeneicAlgorithm::Bound>&);
	std::vector<double> decoder(unsigned int, const std::vector<struct GeneicAlgorithm::Bound>&, const std::string&);
	std::vector<GeneicAlgorithm::Chromosome> roulette_selection(std::vector<GeneicAlgorithm::Chromosome>& x);//roulette selection function : choose offspring chromosomes
	int crossover(std::vector<GeneicAlgorithm::Chromosome>&);//do crossover on chromosome set
	int mutation(std::vector<GeneicAlgorithm::Chromosome>&);//do mutation on chromosome set

private:
	//variables and constants
	
	//encoding related
	std::vector<GeneicAlgorithm::Chromosome> chromosomes;//store all chromosome
	std::vector< GeneicAlgorithm::Bound> bounds_of_variables;//store all boundary of variables

	//evaluation related
	std::function<double(std::vector<double>)> eval_function;//evaluate function

	//crossover and mutation related
	double probability_crossover = 0.8;//probability of crossover on each pair of chromosome 
	double probality_mutation=0.01;//probality of mutation on each point of gene

protected:
	//utility methods

	std::vector<double> get_random_numbers(unsigned int, struct GeneicAlgorithm::Bound={ 0, 1, 0.0000000001 });

};


