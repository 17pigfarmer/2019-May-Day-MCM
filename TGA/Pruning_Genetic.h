//Pruning_Genetic.h
//pure genetic algorithm is from:https://blog.csdn.net/Touch_Dream/article/details/68066055


#ifndef _GENETIC_H_
#define _GENETIC_H_
#define HighThree
# include <fstream>
using namespace std;
//#define PROFIT

#ifdef HighThree
#define TIME 30
#endif // HighThree


#define  PI    3.14159265358979323846

//遗传算法参数，种群规模（0~100）、繁殖代数、函数变量个数、交叉概率、变异概率
# define GROUP_SCALE    50     
# define MAX_GENS       1000
# define N_VARS         80
# define P_MATING       0.8
# define P_MUTATION     0.20


//木板参数：木板长度，宽度，面积，四种零件的长度宽度，零件种类数,最小零件尺寸，使用的零件数
#define B_H		1500
#define B_W		3000
#define B_AREA	4500000
#define P1_L	373
#define P1_W	201
#define P2_L	406
#define P2_W	229
#define P3_L	477
#define P3_W	282
#define P4_L	311
#define P4_W	225
#define P_Num	8
#define P_MIN_WIDTH 201
#define TRY_TIME 2

//基因编码
struct GCode {
	double Area;//面积
	double Length;//长度
	double Width;//宽度
};

//每种零件有竖直，水平两种摆放
extern struct GCode Genetics[P_Num]; 


extern double Profit[P_Num];


struct Individual
{
	int Xn[N_VARS];      //存放变量值
	double Fitness;         //适应值
	double ReFitness;       //适应值概率密度
	double SumFitness;      //累加分布，为轮盘转
	
};
struct X_Range
{
	int Upper;           //变量的上界取值
	int Lower;           //变量的下界取值
};

template<typename T>
T randT(T Lower, T Upper); //产生任意类型随机数函数

void crossover(int& seed);
void elitist();        //基因保留
void evaluate();

void initGroup(int& seed);

void selectBest();
void mutate(int& seed);

double r8_uniform_ab(double a, double b, int& seed);
int i4_uniform_ab(int a, int b, int& seed);

void report(int Xnration);
void selector(int& seed);
void showTime();
void Xover(int one, int two, int& seed);
int Cutting(int* Xn);
int showCutting(int* Xn);
inline int CanBePut(int Lwidth, int Plength, int Pwidth, int Lheight);

void randperm(int N, int* a);

#endif // !_GENETIC_H_
