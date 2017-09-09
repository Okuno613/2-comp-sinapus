
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>


int main()
{  
    int k, n, N;
    double time,volt;
    int max = 100000;  // 読み込むデータ数の上限
    double f1[max+1];  // 格納配列
    double f2[max+1];
    int i=0;

    FILE *fp1;
    fp1=fopen("eodam-circle_10_10_0_0.txt","w");

    std::ifstream ifs("eodam-circle_10_10_0_0.dat");
	std::string str;

	while(getline(ifs,str)){		
    	std::cout<< str << std::endl;
    	sscanf(str.c_str(),"%lf %lf",&time,&volt);
    	f1[i]=time;
    	f2[i]=volt;
    	i++;
	}
	n=i;

	for (int i = 0; i < n; ++i){
		//fprintf(fp1,"%f \t %f \n",f1[i],f2[i]);
		fprintf(fp1,"%f\n",abs(f2[i]) );
	}

	fclose(fp1);

    return 0;
}
