#include <math.h>
#include <iostream>
#include <fstream>

#define FACTOR_NUM 4
#define SIGNAL_NUM 2

int func_num[FACTOR_NUM][SIGNAL_NUM] = {{0,1},{2,3},{4,5},{6,7}};

double sig(double factor,double a,double b,double c){
  return 1.0/(1.0+exp(a*(b*factor-c)));
}

double gauss(double factor, double a, double b, double c){
  return exp(-1.0*pow(a*(factor-b),2.0)/c);
}

double inv_down(double factor, double a, double b){
  double A = pow(2.0,a);
  double B = A - 1.0;
  double C = b * factor + 1;
  double D = pow(C,a);
  return A/(B*D) - 1.0/B;
}

double inv_up(double factor,double a,double b){
  double A = pow(2.0,a);
  double B = A - 1.0;
  double C = b * factor + 1;
  double D = pow(C,a);
  return (D-1.0)/B;
}

//1〜10種類ある予測関数番号をもらうと計算をする関数を作ります
double func(double factor, double mental, int func_num){
  double a,b,c;
  double ret;
  switch(func_num){
    case 0:
      //timing positive, gauss
      a = 2.0;
      b = 0.5;
      c = mental;
      ret = gauss(factor,a,b,c);
      break;
    case 1:
      //-------------------TODO
      a = 2.0*(mental+5.0);
      b = 1.0/80.0;
      c = 0.5;
      ret = sig(factor,a,b,c);
      break;
    case 2:
      //oqtype positive, sig
      a = 1.5;
      b = -1.0;
      c = 3.0*mental*mental;
      ret = sig(factor,a,b,c);
      break;
    case 3:
      //-------------------TODO
      a = 1.0/80.0;
      b = 0.3;
      c = 0.05/(mental+1.0);
      ret = gauss(factor,a,b,c);
      break;
    case 4:
      //topic positive, sig
      a = -0.5*(mental-0.5);
      b = 10.0;
      c = 0.5*pow(mental,mental);
      ret = sig(factor,a,b,c);
      break;
    case 5:
      //-------------------TODO
      a = 10.0/(mental+1.0);
      b = 1.0;
      ret = inv_up(factor,a,b);
      break;
    case 6:
      //turn positive, sig
      a = -3.0;
      b = 1.0;
      c = 0.5/mental;
      ret = sig(factor,a,b,c);
      break;
    case 7:
      //-------------------TODO
      a = 2.0*(mental/2.0+4.0);
      b = 1.0;
      ret = inv_down(factor,a,b);
      break;
  return ret;
  }
}


//1〜10種類ある予測関数番号をもらうと計算をする関数を作ります
double func2(double factor,double mental, int func_num){
	switch(func_num){
		case 0:
			return 2.0/(1.0+exp((factor/80.0-0.15-mental/50.0)*30.0))-1.0;
			break;
		case 1:
			return (pow((factor+1.0),(11.0-mental))-1.0)/(pow(2.0,(11.0-mental))-1.0)*2.0-1.0;
			break;
		case 2:
			return -exp(-factor*factor/((11.0-mental)/500.0))+exp(-(factor-0.5)*(factor-0.5)/(mental/300.0));
			break;
		case 3:
			return -exp(-factor*factor/((11.0-mental)/500.0))+exp(-(factor-0.5)*(factor-0.5)/(mental/300.0));
			break;
		case 4:
			return (pow(2.0,(11.0-mental)*0.8)/(pow(2.0,(11.0-mental)*0.8)-1.0)/pow((factor/80.0+1.0),(0.8*(11.0-mental)))-1.0/(pow(2.0,(0.8*(11.0-mental)))-1))*2.0-1.0;
			break;
		case 5:
			return exp(-(factor-0.2)*(factor-0.2)/(mental/700.0))-exp(-(factor-0.7)*(factor-0.7)/((11.0-mental)/300.0));
			break;
		case 6:
			return exp(-factor*factor/(mental/600.0))-exp(-(factor-0.5)*(factor-0.5)/((11.0-mental)/200.0));
			break;
		case 7:
			return 2.0/(1.0+exp(-35.0*(factor/240.0-0.2+mental/50.0)))-1.0;
			break;
		case 8:
			return (pow((factor/80.0+1.0),(9.0*(11.0-mental)))-1.0)/(pow(2.0,(9.0*(11.0-mental)))-1.0);
			break;
		case 9:
			return 1.0/(1.0+exp(50.0*(factor/80.0-0.1-mental/50.0)))-1.0;
			break;
		case 10:
			return 0.2/(1.0+exp(10.0*mental))-0.1;//心的状態弾性
			break;
	}
}

double inv_norm(int num, double val){
  double norm_val[40] = {80,80,80,80,
            1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,
            320,320,320,320,
            80,80,80,80,
            80,80,80,80};

  if((num == 28)||(num == 29)||(num==30)||(num==31)){
    val = val * norm_val[num] - 80.0;
  }else{
    val = val*norm_val[num];
  }
  return val;
}
