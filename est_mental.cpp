//現在のMを推定する関数
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
 
#define MAX_MENTAL 10 
#define STEP -2 
  
//#pragma warning(disable:4996)

double func(double factor, double mental, int func_num);
double func2(double factor, double mental, int func_num);
double inv_norm(int num, double val);


double est_mental(double observed[4],//観測した表情認識結果の配列を引数にして下さい
  double factor[10],//状況も引数に入れます
  double para[10][4])//パラメタも必要でした
{
  double error[4] = { 0 };//4表情の誤差
  double error_all = 0;
  double tmp_error = 1e+100;
  double pro_face[4] = { 0 };//推定した表情のことでーす
  double step = pow(10.0, STEP);
  double tmp_pro_mental = 0;

  int func_num[10][4]={{0,1,2,3},{4,5,6,7},{12,13,14,15},{8,9,10,11},//もとのやつパート１
  {16,17,18,19},{20,21,22,23},//7つの状況に対する予測関数、各情動4つ分ずつある
  {24,25,26,27},{28,29,30,31},{32,33,34,35},{36,37,38,39}};//もとのやつパート２

  for (double pro_mental = 0.001; pro_mental < MAX_MENTAL; pro_mental = pro_mental + step){
    for (int facial_type = 0; facial_type < 4; facial_type++){
      for (int situation_type = 0; situation_type < 10; situation_type++){
        pro_face[facial_type]+=para[situation_type][facial_type]*func(inv_norm(func_num[situation_type][facial_type],factor[situation_type]), pro_mental, func_num[situation_type][facial_type]);
      }
      error[facial_type] = pow(observed[facial_type] - pro_face[facial_type],2);//観測値との二乗誤差を算出
      error_all+=error[facial_type];
    }//ここまでで、あるpro_mentalのときの、各表情の予測値が出ましたよと

    if (error_all < tmp_error){
      tmp_error = error_all;//tmperror更新
      tmp_pro_mental = pro_mental;//暫定メンタルも更新
    }
    error_all = 0;
    for (int i = 0; i < 4; i++){
      pro_face[i] = 0;
    }
  }
  std::cout << "estimated mental is " << tmp_pro_mental<<std::endl;
  return tmp_pro_mental;
}
