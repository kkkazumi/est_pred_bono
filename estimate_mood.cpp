#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <stdlib.h>
#include <math.h>

//#include <sys/stat.h>
//#include <sys/types.h>

#define N 30

//#define ERROR 0.01;//1e-6;

using namespace std;

double func(double factor, double mental, int func_num);
double func2(double factor, double mental, int func_num);
double inv_norm(int num, double val);

double est_mental(double observed[4], double factor[10], double para[10][4]);

int main(int argc, const char* argv[]){//(argv[0]=./estimate),argv[1]=usrname
		ifstream ifs;
		ifstream ifs2;
		string line;
		string line2;

		double factor_b[100][N];//データ数は１５個、要素数は７こ、
		double factor[100][N];//データ数は１５個、要素数は７こ、
		double signal[100][4];//Sのでーた
		string data_f[N];
		string data_s[4];//signalの読み込みを途中まで
		double hoga;
		double hage;
		int i;
		double tmp=0;
		int row=0;//全部で83ある
		stringstream ss;
		stringstream ss2;
		int j=0;//読み込んだ行数を覚えておく？じゃなくて、j=50になったらリセット

		char out_name[100],out_name2[100];
		char date_num[100];
		char factor_name[100];
		char signal_name[100];
		char face_name[100];
		char err_log[50];
	  char hap_in[100],sup_in[100],ang_in[100],sad_in[100];
	  double ment_est_para[10][4] = { 0 };//現在のM推定するときに使う(mode=1)

    int set_num;

  for(int id_num=0;id_num<9;id_num++){
    for(int num =0; num<8; num++){
      set_num = (1+num)* 5;
      //set_num = 40;
      for(int t=0;t<30;t++){
        sprintf(out_name,"./data/%d/estimated_phi.csv",id_num+1);//コマンドライン引数のmental.csvのタイトルをファイル名にする
        sprintf(factor_name,"./data/%d/factor_test.csv",id_num+1);

        ifs.open(factor_name);//状況データ読み込み
        while(getline(ifs,line)){
          stringstream ss(line);
          for(i=0;i<N;i++){
            getline(ss,data_f[i],',');
            hoga=std::stod(data_f[i]);
            hage=std::stod(data_f[i]);
            factor[j][i]=hage;
          }
          j++;
        }

        j=0;
        sprintf(signal_name,"./data/%d/signal_test.csv",id_num+1);
        ifs2.open(signal_name);//Sデータを読み込み
        while(getline(ifs2,line2)){
          stringstream ss2(line2);
          for(i=0;i<4;i++){
            getline(ss2,data_s[i],',');
            hoga=std::stod(data_s[i]);
            hage=std::stod(data_s[i]);
            signal[j][i]=hage;
          }
          j++;
        }

        //表情観測結果から現在のM推定するための表情モデルパラメタを読み込み
        double facial_para[4] = {};

        sprintf(hap_in, "./data/%d/hap_weight.csv", id_num+1);//happyのパラメタを読み込む
        sprintf(sup_in, "./data/%d/sup_weight.csv", id_num+1);//supのパラメタを読み込む
        sprintf(ang_in, "./data/%d/ang_weight.csv", id_num+1);//angのパラメタを読み込む
        sprintf(sad_in, "./data/%d/sad_weight.csv", id_num+1);//sadのパラメタを読み込む

        ifstream hap_model_in(hap_in);//happyのおもみを読み込み
        ifstream sup_model_in(sup_in);//supのおもみを読み込みf
        ifstream ang_model_in(ang_in);//angのおもみを読み込み
        ifstream sad_model_in(sad_in);//sadのおもみを読み込み

        for (int hap_para_in = 0; hap_para_in < 10; hap_para_in++){
          hap_model_in >> ment_est_para[hap_para_in][0];
          facial_para[0]+= ment_est_para[hap_para_in][0];
        }
        for (int sup_para_in = 0; sup_para_in < 10; sup_para_in++){
          sup_model_in >> ment_est_para[sup_para_in][1];
          facial_para[1]+= ment_est_para[sup_para_in][1];
        }
        for (int ang_para_in = 0; ang_para_in < 10; ang_para_in++){
          ang_model_in >> ment_est_para[ang_para_in][2];
          facial_para[2]+= ment_est_para[ang_para_in][2];
        }
        for (int sad_para_in = 0; sad_para_in < 10; sad_para_in++){
          sad_model_in >> ment_est_para[sad_para_in][3];
          facial_para[3]+= ment_est_para[sad_para_in][3];
        }

        double mental[N];
        double m_ans[N];
        ofstream mout(out_name);
        for(int i=0;i<N;i++){
          mental[i] = est_mental(signal[i], factor[i],ment_est_para);
          mout<<mental[i]<<endl;
        }

      }
    }
  }
}
