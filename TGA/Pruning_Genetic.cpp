//Pruning_Genetic.cpp
#define _CRT_SECURE_NO_WARNINGS

# include <cstdlib>
# include <iostream>
# include <iomanip>

# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
#include "Pruning_Genetic.h"
#include <queue>


//申请种群内存，其中多加1个是放置上一代中最优秀个体
struct Individual Population[GROUP_SCALE + 1];
#ifdef HighThree
Individual High[TIME];
#endif // HighThree



ofstream out = ofstream("out.txt");
#ifdef PROFIT
double Profit[P_Num] = {
	23.0,23.0,19.9,19.9,21,21,16,16
};

#endif // PROFIT

int P_order[P_Num/2] = { 1,3,4,2 };

struct GCode Genetics[P_Num] = {		
									{74973,P1_L,P1_W},{74973,P1_W,P1_L},//1
									{92974,P2_L,P2_W},{92974,P2_W,P2_L},//3
									{69975,P4_L,P4_W},{69975,P4_W,P4_L},//4
									{134514,P3_L,P3_W},{134514,P3_W,P3_L}//2
									
								
};

X_Range  XnRange[N_VARS] = {
							{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},
							{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},
							{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},
							{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},
							{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},
							{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},
							{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},
							{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME},{0,TRY_TIME}

};

struct HLine
{
	int Width;//可用的宽度在变化
	int Base;//左下角横坐标
	int Height;//高度
	
};

//越低的越靠上
bool operator<(HLine a, HLine b) {
	if (a.Height > b.Height) return true;
	return false;    
}

void swapReference(int* a, int i, int j) {
	int tmp = a[i];
	a[i] = a[j];
	a[j] = tmp;
}

void randperm(int N,int*a) {
	
	for (int i = 0; i < N; i++)
		a[i] = i;
	for (int i = 1; i < N; i++)
		swapReference(a, i, randT(0, i));
	
}





//打印排列结果
int showCutting(int* Xn) {
	vector<int> sum;
	for (int j = 0; j < P_Num; j++) {
		sum.push_back(0);
	}


	
	//std::priority_queue实现的大根堆
	priority_queue<HLine> products;
	products.push(HLine{ B_W,0,0 });

	int i = 0;
	do {

		HLine lowline = products.top();
		products.pop();
		if (B_H - lowline.Height < P_MIN_WIDTH || lowline.Width < P_MIN_WIDTH) {//边界剪枝
			continue;
		}
		
		// 1 宽度装 0 不能装
		if (CanBePut(lowline.Width, Genetics[(int)Xn[i]].Length, Genetics[(int)Xn[i]].Width, lowline.Height) == 1) {
			//可以放下那么
			//1.水平线宽度改变
			//2.加入新水平线
			int newwidth, newheight,oldbase;
			
				Genetics[(int)Xn[i]].Length;
				newheight = Genetics[(int)Xn[i]].Length;
				newwidth = Genetics[(int)Xn[i]].Width;
				lowline.Width -= newwidth;
				oldbase = lowline.Base;
				lowline.Base += newwidth;

			products.push(lowline);
			products.push(HLine{ newwidth,oldbase, newheight + lowline.Height });
			sum[Xn[i]]++;
			out << "第" << i << "个零件P"<< P_order[Xn[i] / 2] << " : " << "左下角坐标为(" << oldbase << "," << lowline.Height << "),放置角度是" << ((Xn[i] + 1)%2 == 0 ? "竖直." : "水平.") << endl;
			i++;

		}
		else {

			
			//不可以放下那么
			//1.对其他零件进行判断，取最接近的，可以就改变
			//2.都不可以就换水平线(继续循环)
			int best_width = lowline.Width;
			int best_Pid = -1;

			//int a[TRY_TIME];
			//memset(a, 0, TRY_TIME * sizeof(int));
			//randperm(TRY_TIME, a);

			for (int j = 0; j < TRY_TIME; j++) {//此处当零件过多时，改为随机算法，随机取N个
				int k = j;


				if (CanBePut(lowline.Width, Genetics[k].Length, Genetics[k].Width, lowline.Height) == 1) {
					int local_width = lowline.Width - Genetics[k].Width;
					if (local_width < best_width) {
						best_width = local_width;
						best_Pid = k;
					}
				}
			}

			if (best_Pid == -1) {//没有可以替换的
				continue;
			}
			else {//有可以替换的
				Xn[i] = best_Pid;

				int newwidth, newheight;

				Genetics[(int)Xn[i]].Length;
				newheight = Genetics[(int)Xn[i]].Length;
				newwidth = Genetics[(int)Xn[i]].Width;
				lowline.Width -= newwidth;
				lowline.Base += newwidth;

				products.push(lowline);
				products.push(HLine{ newwidth,lowline.Base - newwidth, newheight + lowline.Height });
				sum[Xn[i]]++;
				out << "第" << i << "个零件P" << P_order[Xn[i]/2 ] << " : " << "左下角坐标为(" << lowline.Base - newwidth << "," << lowline.Height << "),放置角度是" << ((Xn[i] + 1) % 2 == 0 ? "竖直." : "水平.") << endl;
				i++;
			}
		}
	} while (!products.empty());
	int Psum = 0;
	for (int j = 0; j < P_Num/2; j+=2) {
		out << "P" << P_order[j/2 ] << "零件共使用" << sum[j]+sum[j+1] << "个" << endl;
		Psum += sum[j]+sum[j + 1];
	}
	out << "总共使用零件" << Psum << "个。"<<endl;
	
	return i;
}



//改进的最低水平线算法计算对于Xn[]能放多少零件，返回能放的零件个数。
int Cutting(int* Xn) {
	//std::priority_queue实现的大根堆
	priority_queue<HLine> products;
	products.push(HLine{ B_W,0,0 });

	int i = 0;
	do {

		HLine lowline = products.top();
		products.pop();
		if (B_H - lowline.Height < P_MIN_WIDTH || lowline.Width < P_MIN_WIDTH) {//剪枝
			continue;
		}
		
		// 1 宽度装 0 不能装
		if (CanBePut(lowline.Width, Genetics[(int)Xn[i]].Length, Genetics[(int)Xn[i]].Width, lowline.Height) == 1) {
			//可以放下那么
			//1.水平线宽度改变
			//2.加入新水平线
			int newwidth, newheight;
		
			
				Genetics[(int)Xn[i]].Length;
				newheight = Genetics[(int)Xn[i]].Length;
				newwidth = Genetics[(int)Xn[i]].Width;
				lowline.Width -= newwidth;
				lowline.Base += newwidth;
			
			products.push(lowline);
			products.push(HLine{ newwidth,lowline.Base-newwidth, newheight + lowline.Height });
			i++;

		}
		else {
			
			//不可以放下那么
			//1.对其他零件进行判断，取最接近的，可以就改变
			//2.都不可以就换水平线(继续循环)
			int best_width = lowline.Width;
			int best_Pid = -1;


			//int a[TRY_TIME];
			//memset(a, 0, TRY_TIME * sizeof(int));
			//randperm(TRY_TIME, a);

			for (int j = 0; j < TRY_TIME; j++) {//此处当零件过多时，改为随机算法，随机取N个
				int k = j;
								
				
				if (CanBePut(lowline.Width, Genetics[k].Length, Genetics[k].Width, lowline.Height) == 1) {
					int local_width = lowline.Width - Genetics[k].Width;
					if (local_width < best_width) {
						best_width = local_width;
						best_Pid = k;
					}
				}
			}

			if (best_Pid == -1) {//没有可以替换的
				continue;
			}
			else {//有可以替换的
				Xn[i] = best_Pid;
				
				int newwidth, newheight;
				
					Genetics[(int)Xn[i]].Length;
					newheight = Genetics[(int)Xn[i]].Length;
					newwidth = Genetics[(int)Xn[i]].Width;
					lowline.Width -= newwidth;
					lowline.Base += newwidth;
				
				products.push(lowline);
				products.push(HLine{ newwidth,lowline.Base - newwidth, newheight + lowline.Height });
				i++;
			}


		}
	} while (!products.empty());
	return i;
}

/*
inline int CanBePut(int Lwidth,int Plength,int Pwidth,int Lheight) {
	if (Lwidth >= Plength && (Lheight + Pwidth) <= B_H) //能装下长度
		return 2;
	else //不能装下长度
		if (Lwidth >= Pwidth && (Lheight + Plength) <= B_H) //能装下宽度
			return 1;
		else
			return 0;
}
*/

inline int CanBePut(int Lwidth, int Plength, int Pwidth, int Lheight) {
		if (Lwidth >= Pwidth && (Lheight + Plength) <= B_H) //能装下宽度
			return 1;
		else
			return 0;
}


//有交配权的所有父代进行交叉
void crossover(int& seed)
{
	const double a = 0.0;
	const double b = 1.0;
	int mem;
	int one;
	int first = 0;
	double x;

	for (mem = 0; mem < GROUP_SCALE; ++mem)
	{
		x = randT(0.0, 1.0);
		//x = r8_uniform_ab(a, b, seed);//产生交配概率

		if (x < P_MATING)
		{
			++first;

			if (first % 2 == 0)//交配
			{
				Xover(one, mem, seed);
			}
			else
			{
				one = mem;
			}

		}
	}
	return;
}

//对最差的染色体和最优的染色体的处理，起到优化的目的
void elitist()
{
	int i;
	double best;
	int best_mem;
	double worst;
	int worst_mem;

	best = Population[0].Fitness;
	worst = Population[0].Fitness;

	for (i = 0; i < GROUP_SCALE - 1; ++i)
	{
#ifdef HighThree

		for (int j = TIME-1; j >-1; j--) {
			if (Population[i].Fitness >= High[j].Fitness) {
				High[j] = Population[i];
				break;
			}
		}
		
		

#endif // HighThree

		if (Population[i + 1].Fitness < Population[i].Fitness)
		{

			if (best <= Population[i].Fitness)
			{
				best = Population[i].Fitness;
				best_mem = i;
			}

			if (Population[i + 1].Fitness <= worst)
			{
				worst = Population[i + 1].Fitness;
				worst_mem = i + 1;
			}

		}
		else
		{

			if (Population[i].Fitness <= worst)
			{
				worst = Population[i].Fitness;
				worst_mem = i;
			}

			if (best <= Population[i + 1].Fitness)
			{
				best = Population[i + 1].Fitness;
				best_mem = i + 1;
			}

		}

	}

	//对于当前代的最优值的处理，如果当前的最优值小于上一代则将上一代的值最优个体取代当前的最弱个体
	//基因保留
	if (Population[GROUP_SCALE].Fitness <= best)
	{
		for (i = 0; i < N_VARS; i++)
		{
			Population[GROUP_SCALE].Xn[i] = Population[best_mem].Xn[i];
		}
		Population[GROUP_SCALE].Fitness = Population[best_mem].Fitness;
	}
	else
	{
		for (i = 0; i < N_VARS; i++)
		{
			Population[worst_mem].Xn[i] = Population[GROUP_SCALE].Xn[i];
		}
		Population[worst_mem].Fitness = Population[GROUP_SCALE].Fitness;
	}
	return;
}





//计算适应度值
void evaluate()
{
	int member;
	int i;

	
	
#ifdef PROFIT
	for (member = 0; member < GROUP_SCALE; member++)
	{
		int sum = Cutting(Population[member].Xn);
		double sum_pro = 0;
		for (int i = 0; i < sum; i++) {
			sum_pro += Profit[(size_t)Population[member].Xn[i]];
		}
		Population[member].Fitness = sum_pro;
	}

#else
	for (member = 0; member < GROUP_SCALE; member++)
	{
		int sum = Cutting(Population[member].Xn);
		double sum_area = 0;
		for (int i = 0; i < sum; i++) {
			sum_area += Genetics[(size_t)Population[member].Xn[i]].Area;
		}
		Population[member].Fitness = sum_area / B_AREA;
	}

#endif // PROFIT

		
	return;

	

}


//产生整形的随机数
int i4_uniform_ab(int a, int b, int& seed)
{
	int c;
	const int i4_huge = 2147483647;
	int k;
	float r;
	int value;

	if (seed == 0)
	{
		cerr << "\n";
		cerr << "I4_UNIFORM_AB - Fatal error!\n";
		cerr << "  Input value of SEED = 0.\n";
		exit(1);
	}
	//保证a小于b
	if (b < a)
	{
		c = a;
		a = b;
		b = c;
	}

	k = seed / 127773;
	seed = 16807 * (seed - k * 127773) - k * 2836;

	if (seed < 0)
	{
		seed = seed + i4_huge;
	}

	r = (float)(seed) * 4.656612875E-10;
	//
	//  Scale R to lie between A-0.5 and B+0.5.
	//
	r = (1.0 - r) * ((float)a - 0.5)
		+ r * ((float)b + 0.5);
	//
	//  Use rounding to convert R to an integer between A and B.
	//
	value = round(r);//四舍五入
	//保证取值不越界
	if (value < a)
	{
		value = a;
	}
	if (b < value)
	{
		value = b;
	}

	return value;
}

//初始化种群个体
void initGroup(int& seed)

{
	int i;
	int j;
	double lbound;
	double ubound;
	// 
	//  initGroup variables within the bounds 
	//
	for (i = 0; i < N_VARS; i++)
	{
		//input >> lbound >> ubound;

		for (j = 0; j < GROUP_SCALE; j++)
		{
			Population[j].Fitness = 0;
			Population[j].ReFitness = 0;
			Population[j].SumFitness = 0;
			Population[j].Xn[i] = randT(XnRange[i].Lower, XnRange[i].Upper);
			//Population[j].Xn[i] = r8_uniform_ab(XnRange[i].Lower, XnRange[i].Upper, seed);
		}
	}

	return;
}


//挑选出最大值，保存在种群数组的最后一个位置
void selectBest()
{
	int cur_best;
	int mem;
	int i;

	cur_best = 0;

	for (mem = 0; mem < GROUP_SCALE; mem++)
	{
		if (Population[GROUP_SCALE].Fitness < Population[mem].Fitness)
		{
			cur_best = mem;
			Population[GROUP_SCALE].Fitness = Population[mem].Fitness;
		}
	}

	for (i = 0; i < N_VARS; i++)
	{
		Population[GROUP_SCALE].Xn[i] = Population[cur_best].Xn[i];
	}

	return;
}

//个体变异
void mutate(int& seed)
{
	const double a = 0.0;
	const double b = 1.0;
	int i;
	int j;
	double lbound;
	double ubound;
	double x;

	for (i = 0; i < GROUP_SCALE; i++)
	{
		for (j = 0; j < N_VARS; j++)
		{
			//x = r8_uniform_ab(a, b, seed);
			x = randT(a, b);//突变概率
			if (x < P_MUTATION)
			{
				lbound = XnRange[j].Lower;
				ubound = XnRange[j].Upper;
				Population[i].Xn[j] = randT(lbound, ubound);
				//Population[i].Xn[j] = r8_uniform_ab(lbound, ubound, seed);
			}
		}
	}

	return;
}

//模板函数，用于生成各种区间上的数据类型
template<typename T>
T randT(T Lower, T Upper)
{
	return rand() / (double)RAND_MAX * (Upper - Lower) + Lower;
}


vector<int>& RandomSeq(){

	vector<int> A;
	A.push_back(1);
	A.push_back(0);
	if (rand() % 2 != 0) {
		A[0] = 0;
		A[1] = 1;
	}
	return A;
}



//产生小数随机数
double r8_uniform_ab(double a, double b, int& seed)

{
	int i4_huge = 2147483647;
	int k;
	double value;

	if (seed == 0)
	{
		cerr << "\n";
		cerr << "R8_UNIFORM_AB - Fatal error!\n";
		cerr << "  Input value of SEED = 0.\n";
		exit(1);
	}

	k = seed / 127773;
	seed = 16807 * (seed - k * 127773) - k * 2836;

	if (seed < 0)
	{
		seed = seed + i4_huge;
	}

	value = (double)(seed) * 4.656612875E-10;

	value = a + (b - a) * value;

	return value;
}

//输出每一代进化的结果
void report(int Xnration)
{
	double avg;
	double best_val;
	int i;
	double square_sum;
	double stddev;
	double sum;
	double sum_square;

	if (Xnration == 0)
	{
		cout << "\n";
		cout << "  Xnration       Best            Average       Standard \n";
		cout << "  number           value           Fitness       deviation \n";
		cout << "\n";
	}
	sum = 0.0;
	sum_square = 0.0;

	for (i = 0; i < GROUP_SCALE; i++)
	{
		sum = sum + Population[i].Fitness;
		sum_square = sum_square + Population[i].Fitness * Population[i].Fitness;
	}

	avg = sum / (double)GROUP_SCALE;
	square_sum = avg * avg * GROUP_SCALE;
	stddev = sqrt((sum_square - square_sum) / (GROUP_SCALE - 1));
	best_val = Population[GROUP_SCALE].Fitness;

	cout << "  " << setw(8) << Xnration
		<< "  " << setw(14) << best_val
		<< "  " << setw(14) << avg
		<< "  " << setw(14) << stddev << "\n";

	return;
}

//选择有交配权的父代
void selector(int& seed)
{
	struct Individual NewPopulation[GROUP_SCALE + 1];//临时存放挑选的后代个体
	const double a = 0.0;
	const double b = 1.0;
	int i;
	int j;
	int mem;
	double p;
	double sum;

	sum = 0.0;
	for (mem = 0; mem < GROUP_SCALE; mem++)
	{
		sum = sum + Population[mem].Fitness;
	}
	//计算概率密度
	for (mem = 0; mem < GROUP_SCALE; mem++)
	{
		Population[mem].ReFitness = Population[mem].Fitness / sum;
	}
	// 计算累加分布，思想是轮盘法
	Population[0].SumFitness = Population[0].ReFitness;
	for (mem = 1; mem < GROUP_SCALE; mem++)
	{
		Population[mem].SumFitness = Population[mem - 1].SumFitness +
			Population[mem].ReFitness;
	}
	// 选择个体为下一代繁殖，选择优秀的可能性大，这是轮盘法的奥秘之处
	for (i = 0; i < GROUP_SCALE; i++)
	{
		p = r8_uniform_ab(a, b, seed);
		if (p < Population[0].SumFitness)
		{
			NewPopulation[i] = Population[0];
		}
		else
		{
			for (j = 0; j < GROUP_SCALE; j++)
			{
				if (Population[j].SumFitness <= p && p < Population[j + 1].SumFitness)
				{
					NewPopulation[i] = Population[j + 1];
				}
			}
		}
	}
	//更新后代个体 
	for (i = 0; i < GROUP_SCALE; i++)
	{
		Population[i] = NewPopulation[i];
	}
	return;
}

//显示系统时间
void showTime()
{
# define TIME_SIZE 40

	static char time_buffer[TIME_SIZE];
	const struct tm* tm;
	size_t len;
	time_t now;

	now = time(NULL);
	tm = localtime(&now);

	len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

	cout << time_buffer << "\n";

	return;
# undef TIME_SIZE
}

//交叉产生子代
void Xover(int one, int two, int& seed)
{
	int i;
	int point;
	double t;
	//随机选择交叉点，这里的点是以变量的整个长度为单位
	point = randT<int>(0, N_VARS - 1);
	//point = i4_uniform_ab(0, N_VARS - 1, seed);
	//交叉
	for (i = 0; i < point; i++)
	{
		t = Population[one].Xn[i];
		Population[one].Xn[i] = Population[two].Xn[i];
		Population[two].Xn[i] = t;
	}
	return;
}
