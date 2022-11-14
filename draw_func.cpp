#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <string>

//double sig(double factor,double a,double b,double c);
//double gauss(double factor, double a, double b, double c);
double func(double factor, double mental, int func_num);

//int func_num[FACTOR_NUM][SIGNAL_NUM] = {{0,1},{2,3},{4,5},{6,7}};

void out_data(int func_num){
	double val[100];
	for (int f = 0; f < 100; f++) {
		double factor = f / 100.0;
		for (int m = 0; m < 3; m++) {
			double mental = m * 0.45 + 0.05;
			//turn_coe(mental,coe);
			val[f] = func(factor,mental,func_num);
			std::cout<<val[f];
			if(m<2){
				std::cout<<",";
			}
		}
		std::cout<<"\n";
	}
}


int main(int argc, char* argv[]){
	int func_num = std::stoi(argv[1]);
	out_data(func_num);
	return 0;
}
