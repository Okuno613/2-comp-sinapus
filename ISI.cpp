
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
 
#define pi 3.1415926535         // 円周率

 // 離散フーリエ変換（読み込むファイル名, 書き込むファイル名）
int dft(char filename1[], char filename2[])
{  
        int k, n, N;
	int spike=0;
	int spikecnt=0;
	int count=0;
	int tmp_v=0;
        int max = 100000;  // 読み込むデータ数の上限
        double f1[max+1];
	double ans;
	
        FILE *fp1, *fp2;
    // ファイルオープン(フーリエ変換したいデータファイル, フーリエ変換後のデータ保存用ファイル)
        if((fp1=fopen(filename1,"r"))==NULL){
	  printf("FILE1 not open\n");
        return -1;
        }

        //データの読み込み
        for(N=0; N<max; N++) {
	  fscanf(fp1,"%lf", &f1[N]);
        }
	if((fp2=fopen(filename2,"w"))==NULL){
	  printf("FILE2 not open\n");
        return -1;
        }
	
        //実数部分と虚数部分に分けてフーリエ変換
        for(k=0; k<N; k++){
	  tmp_v=f1[k];

	  //	  if( -0.005<= (f1[k]-tmp_v)&& (f1[k]-tmp_v) <= 0 && count==0 ){
	  if( (f1[k-1]-f1[k-2])=>0 && (f1[k]-f1[k-1])<=0 && f1[k]>0  ){
	    spike=k-spikecnt;
	    spikecnt=k;
	    count==1;
	    if(spike != 1 && spike != 0){
	      fprintf(fp2,"%02d\n",int(spike*0.1));
	    }
	  }
	  
	  if(f1[k]<=0){
	    count=0;
	  }

	}
  
        fclose(fp1);
        fclose(fp2);
	
	
        return 0;
}

int main()
{

  dft("Vs_volt.txt","isi.txt");
  //dft("test.txt","t.txt");
}
