#include "GeneicAlgorithm.hpp"

GeneicAlgorithm::Chromosome GeneicAlgorithm::encoder(unsigned int Dimension,
	const std::vector<double>& Variables, const std::vector<struct GeneicAlgorithm::Bound>& Bounds)
{
	if (Dimension < 1)
		throw std::logic_error("Variable dimension cannot smaller than 1!");
	if (Variables.size() != Dimension or Bounds.size() != Dimension)
		throw std::runtime_error("Variable counts do not equal to dimension!");

	GeneicAlgorithm::Chromosome chromosome;
	chromosome.dimension = Dimension;
	chromosome.values = Variables;

	size_t length = 0;// length of chromosome
	chromosome.sub_length.resize(Dimension);// length of each variable in chromosome
	for (int i = 0; i < Dimension; i++)
	{
		auto length_interval = Bounds[i].UpperBound - Bounds[i].LowerBound;
		chromosome.sub_length[i] = std::ceil(std::log(length_interval / Bounds[i].Precision) / std::log(2));
		length += chromosome.sub_length[i];
	}

	chromosome.str.reserve(length);

	//lambda converting decimal to binary string
	auto to_binary_string = [](long decimal, size_t length){
		std::string result(length,'0');
		auto pt = result.rbegin();
		while (decimal)
		{
			*pt++ = decimal % 2 == 1 ? '1' : '0';
			decimal /= 2;

		}
		return result;
	};

	for (int i = 0; i < Dimension; i++)
	{
		long sub_chromosome = (Variables[i] - Bounds[i].LowerBound) / Bounds[i].Precision;
		auto binary_string = to_binary_string(sub_chromosome, chromosome.sub_length[i]);
		chromosome.str += binary_string;
	}

	chromosome.adaptibility = eval_function(chromosome.values);
	return chromosome;
}

std::vector<double> GeneicAlgorithm::decoder(unsigned int Dimension,
	const std::vector<struct GeneicAlgorithm::Bound>& Bounds, const std::string& Chromosome)
{
	if (Dimension < 1)
		throw std::logic_error("Variable dimension cannot smaller than 1!");
	if (Bounds.size() != Dimension)
		throw std::runtime_error("Variable bounds counts do not equal to dimension!");

	std::vector<size_t> sub_length(Dimension);// length of each variable in chromosome
	for (int i = 0; i < Dimension; i++)
	{
		auto length_interval = Bounds[i].UpperBound - Bounds[i].LowerBound;
		sub_length[i] = std::ceil(std::log(length_interval / Bounds[i].Precision) / std::log(2));
	}

	//lambda converting binary string to decimal
	auto to_decimal = [&](size_t begin,size_t end){
		double val = 0;
		for (int i = 0; i < end - begin; i++)
		{
			val += pow(2, i) * ('1' == Chromosome[end - i - 1] ? 1 : 0);
		}
		return val;
	};
	std::vector<double> variable_return(Dimension,0);
	size_t begin_pt = 0;//pointing at string in begin
	size_t end_pt = 0;//pointing at string in end
	for (int i = 0; i < Dimension; i++)
	{
		end_pt += sub_length[i];
		variable_return[i] = to_decimal(begin_pt, end_pt) * Bounds[i].Precision + Bounds[i].LowerBound;
		begin_pt += sub_length[i];
	}
	return variable_return;
}

std::vector<GeneicAlgorithm::Chromosome> GeneicAlgorithm::roulette_selection(std::vector<GeneicAlgorithm::Chromosome>& Chromosome_Sets)
{
	double sum_of_adaptability = 0;//i.e. F=sigma(eval(v_k))
	std::vector<double> accumulative_probility(Chromosome_Sets.size(),0);//i.e. q_k, accumulative probility of each individual v_k
	

	for (auto i = 0; i < Chromosome_Sets.size(); i++)
	{
		sum_of_adaptability += Chromosome_Sets[i].adaptibility;
		std::for_each(accumulative_probility.begin(), accumulative_probility.begin() + i + 1, [&](double &item) {
			item += Chromosome_Sets[i].adaptibility;
			}
		);
	}

	std::for_each(accumulative_probility.begin(), accumulative_probility.end(), [&](double &item) {
		item /= sum_of_adaptability;
		}
	);

	auto probilities = get_random_numbers(Chromosome_Sets.size());
	std::vector<GeneicAlgorithm::Chromosome> return_set(Chromosome_Sets.size());

	for (auto i = 0; i < probilities.size(); i++)
	{
		decltype(i) j;
		for (j = 1; j < accumulative_probility.size(); j++)
		{
			if (accumulative_probility[j] < probilities[i])
			{
				return_set[i] = Chromosome_Sets[j-1];
				break;
			}
		}
		if (accumulative_probility.size() == j)
			return_set[i] = *(Chromosome_Sets.end()-1);

	}

	return return_set;
}

int GeneicAlgorithm::crossover(std::vector<GeneicAlgorithm::Chromosome>& chromosomeSet)
{
	//Generate a random number and determine whether to cross according to the size of the random number and the probability
	//Generate a random number and use this random number as the crossover point
	//Use the decoder to decode after the crossover to check whether the value of the daughter chromosome is legal. If it is not legal, re-determine the crossover site
	//Crossover all chromosomes in pairs to obtain offspring

	int crossover_count = 0;

	auto try_cross = [&](auto i) ->bool{
		GeneicAlgorithm::Chromosome& chromosome_left = chromosomeSet[i * 2];
		GeneicAlgorithm::Chromosome& chromosome_right = chromosomeSet[i * 2 + 1];
		
		//select crossover point
		struct GeneicAlgorithm::Bound bound { 1, chromosome_left.str.length()-1, 1 };
		auto cross_point = get_random_numbers(1, bound);
		cross_point[0] = std::round(cross_point[0]);
		
		//do crossover
		std::swap_ranges(chromosome_left.str.begin() + cross_point[0], chromosome_left.str.end(), chromosome_right.str.begin() + cross_point[0]);
		
		//get value from chromosome string
		auto val_left = decoder(chromosome_left.dimension, bounds_of_variables, chromosome_left.str);
		auto val_right = decoder(chromosome_right.dimension, bounds_of_variables, chromosome_right.str);
		
		//begin value check
		for (int j = 0; j < val_left.size(); j++)
		{
			if (val_left[j]<bounds_of_variables[j].LowerBound or val_left[j]>bounds_of_variables[j].UpperBound)
			{
				//restore crossover
				std::swap_ranges(chromosome_left.str.begin() + cross_point[0], chromosome_left.str.end(), chromosome_right.str.begin() + cross_point[0]);
				return false;
			}
			if (val_right[j]<bounds_of_variables[j].LowerBound or val_right[j]>bounds_of_variables[j].UpperBound)
			{
				//restore crossover
				std::swap_ranges(chromosome_left.str.begin() + cross_point[0], chromosome_left.str.end(), chromosome_right.str.begin() + cross_point[0]);
				return false;
			}
		}

		//update values
		chromosome_left.values = val_left;
		chromosome_right.values = val_right;
		chromosome_left.adaptibility = eval_function(chromosome_left.values);
		chromosome_right.adaptibility = eval_function(chromosome_right.values);
		return true;
	};

	auto probility_of_crossover_or_not = get_random_numbers(chromosomeSet.size()/2);//generate N/2 random numbers to ensure each pair of chromosome does crossover or not
	for (auto i = 0; i < probility_of_crossover_or_not.size(); i++)
	{
		if (probility_of_crossover_or_not[i] < probability_crossover)
		{
			crossover_count++;

			//keep crossover until success
			while (not try_cross(i));

		}
	}

	return crossover_count;
}

int GeneicAlgorithm::mutation(std::vector<GeneicAlgorithm::Chromosome>& chromosomeSet)
{
	//Generate random numbers for each site to determine whether to mutate.
	//Judge the legitimacy of the offspring chromosomes.
	
	int mutation_count = 0;

	//generate M random number to ensure mutation or not
	std::vector<double> probility_of_mutation_or_not;

	auto try_mutate = [&](int i) {
		probility_of_mutation_or_not = get_random_numbers(chromosomeSet[0].str.length());
		std::string str_before = chromosomeSet[i].str;//store chromosome string for restore
		int mutation_count_temp = 0;
		for (auto j = 0; j < chromosomeSet[i].str.size(); j++)
		{
			if (probility_of_mutation_or_not[j] < probality_mutation)
				chromosomeSet[i].str[j] = chromosomeSet[i].str[j] == '1' ? '0' : '1', mutation_count_temp++;
		}

		auto val = decoder(chromosomeSet[i].dimension, bounds_of_variables, chromosomeSet[i].str);
		for (auto k = 0; k < val.size(); k++)
		{
			if (val[k]<bounds_of_variables[k].LowerBound or val[k]>bounds_of_variables[k].UpperBound)
			{
				chromosomeSet[i].str = str_before;
				return false;
			}
		}

		chromosomeSet[i].values = val;
		chromosomeSet[i].adaptibility = eval_function(chromosomeSet[i].values);
		mutation_count += mutation_count_temp;
		return true;
	};

	for (auto i = 0; i < chromosomeSet.size(); i++)
	{
		while (not try_mutate(i));
	}
	return mutation_count;
}

std::vector<double> GeneicAlgorithm::get_random_numbers(unsigned int N, struct GeneicAlgorithm::Bound bound)
{
	std::random_device rd;  // Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(bound.LowerBound,bound.UpperBound);

	std::vector<double> result(N);
	std::for_each(result.begin(), result.end(), [&](double& item) {
		item = dis(gen);
		}
	);

	return result;
}

bool GeneicAlgorithm::execute()
{

	//only for debug
	eval_function = [](std::vector<double> v) {
		double result = 0; std::for_each(v.begin(), v.end(), [&](double x) {result += x*x; }); return result; };
	bounds_of_variables = std::vector<struct GeneicAlgorithm::Bound>{ {-100,100,0.01}, { -100,100,0.01}
	, { -100,100,0.01} , { -100,100,0.01}, { -100,100,0.01}, { -100,100,0.01}, { -100,100,0.01}, { -100,100,0.01}
	, { -100,100,0.01} , { -100,100,0.01} };
	dimension = 10;
	max_evolution_time = 500;
	population_size = 50;
	probability_crossover = 0.4;
	probality_mutation = 0.01;

	std::vector<GeneicAlgorithm::Chromosome> chromosomeSet(population_size);

	//init chromosome randomly
	for (auto& item : chromosomeSet)
	{
		std::vector<double> variables(dimension);
		for (auto i = 0; i < dimension; i++)
		{
			variables[i] = get_random_numbers(1, bounds_of_variables[i])[0];
		}

		item = encoder(dimension, variables, bounds_of_variables);
	}

	//lambda seeking best item in set and show
	auto show_best_item = [&](int turn) {
		auto best_element = std::max_element(chromosomeSet.begin(), chromosomeSet.end(), 
			[](Chromosome item1, Chromosome item2) {
				return item1.adaptibility < item2.adaptibility;
			});
		std::cout << "round " << turn<<' ';
		(*best_element).show_values();
	};

	for (auto i = 0; i < max_evolution_time; i++)
	{
		chromosomeSet = roulette_selection(chromosomeSet);
		auto ret_val = crossover(chromosomeSet);
		ret_val = mutation(chromosomeSet);
		
		show_best_item(i);
	}

	show_best_item(max_evolution_time);

	return true;
}