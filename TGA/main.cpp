//main.cpp
#include <iostream>
#include "Pruning_Genetic.h"

extern Individual Population[GROUP_SCALE + 1];
#ifdef HighThree
extern Individual High[];
#endif // HighThree
extern ofstream out;

int main()
{
	int Xnration;
	int i;
	int seed = 142;

	showTime();
	initGroup(seed);
	evaluate();
	selectBest();

	for (Xnration = 0; Xnration < MAX_GENS; Xnration++)
	{
		selector(seed);
		crossover(seed);
		mutate(seed);
		report(Xnration);
		evaluate();
		elitist();
	}

	cout << "\n";
	cout << "  Best member after " << MAX_GENS << " Xnrations:\n";
	cout << "\n";

	for (i = 0; i < N_VARS; i++)
	{
		cout << "  X(" << i + 1 << ") = " << Population[GROUP_SCALE].Xn[i] << "\n";
	}
	cout << "\n";
	cout << "  Best Fitness = " << Population[GROUP_SCALE].Fitness << "\n";

#ifdef HighThree
	for (int j = TIME - 1; j > -1; j--) {
		cout << "第" << TIME - j << "好的方案为：适应度" << High[j].Fitness << endl;
	}
	for (int j = TIME - 1; j > -1; j--) {
		out << "第" << TIME - j << "好的方案为：适应度" << High[j].Fitness << endl;
		showCutting(High[j].Xn);
	}
#else
	showCutting(Population->Xn);
#endif //HighThree

	showTime();
	while (1);
	return 0;
}

