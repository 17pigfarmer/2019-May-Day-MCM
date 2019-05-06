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

//�Ŵ��㷨��������Ⱥ��ģ��0~100������ֳ��������������������������ʡ��������
# define GROUP_SCALE    50     
# define MAX_GENS       1000
# define N_VARS         80
# define P_MATING       0.8
# define P_MUTATION     0.20


//ľ�������ľ�峤�ȣ���ȣ��������������ĳ��ȿ�ȣ����������,��С����ߴ磬ʹ�õ������
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

//�������
struct GCode {
	double Area;//���
	double Length;//����
	double Width;//���
};

//ÿ���������ֱ��ˮƽ���ְڷ�
extern struct GCode Genetics[P_Num]; 


extern double Profit[P_Num];


struct Individual
{
	int Xn[N_VARS];      //��ű���ֵ
	double Fitness;         //��Ӧֵ
	double ReFitness;       //��Ӧֵ�����ܶ�
	double SumFitness;      //�ۼӷֲ���Ϊ����ת
	
};
struct X_Range
{
	int Upper;           //�������Ͻ�ȡֵ
	int Lower;           //�������½�ȡֵ
};

template<typename T>
T randT(T Lower, T Upper); //���������������������

void crossover(int& seed);
void elitist();        //������
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
