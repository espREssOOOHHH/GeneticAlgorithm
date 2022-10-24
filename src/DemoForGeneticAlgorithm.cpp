// DemoForGeneticAlgorithm.cpp: 定义应用程序的入口点。
//

#include "DemoForGeneticAlgorithm.h"
#include "GeneicAlgorithm.hpp"

using namespace std;

int main()
{
	GeneicAlgorithm ga(1);

	for (auto i = 0; i < 10000; i++)
	{
		ga.tester();
	}

	return 0;
}
