#include "GeneicAlgorithm.hpp"

GeneicAlgorithm::Chromosome GeneicAlgorithm::encoder(unsigned int Dimension,
	const std::vector<double>& Variables, const std::vector<struct GeneicAlgorithm::Bound>& Bounds)
{
	if (Dimension < 1)
		throw std::logic_error("Variable dimension cannot smaller than 1!");
	if (Variables.size() != Dimension or Bounds.size() != Dimension)
		throw std::runtime_error("Variable counts do not equal to dimension!");

	size_t length = 0;// length of chromosome
	std::vector<size_t> sub_length(Dimension);// length of each variable in chromosome
	for (int i = 0; i < Dimension; i++)
	{
		auto length_interval = Bounds[i].UpperBound - Bounds[i].LowerBound;
		sub_length[i] = std::ceil(std::log(length_interval / Bounds[i].Precision) / std::log(2));
		length += sub_length[i];
	}
	GeneicAlgorithm::Chromosome chromosome;
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
		auto binary_string = to_binary_string(sub_chromosome, sub_length[i]);
		chromosome.str += binary_string;
	}

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
		begin_pt += end_pt;
	}
	return variable_return;
}

std::vector<GeneicAlgorithm::Chromosome> GeneicAlgorithm::roulette_selection(std::vector<GeneicAlgorithm::Chromosome>& Chromosome_Sets)
{
	double sum_of_adaptability = 0;//i.e. F=sigma(eval(v_k))
	std::vector<double> adaptabilities(Chromosome_Sets.size());//i.e. p_k, adaptability of each individual v_k
	decltype(adaptabilities) accumulative_probility(Chromosome_Sets.size(),0);//i.e. q_k, accumulative probility of each individual v_k
	

	for (auto i = 0; i < Chromosome_Sets.size(); i++)
	{
		adaptabilities[i] = eval_function(Chromosome_Sets[i].values);
		sum_of_adaptability += adaptabilities[i];
		std::for_each(accumulative_probility.begin(), accumulative_probility.begin() + i + 1, [&](double &item) {
			item += adaptabilities[i];
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

std::vector<double> GeneicAlgorithm::get_random_numbers(unsigned int N, struct GeneicAlgorithm::Bound bound= { 0,1,0.0000000001 })
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