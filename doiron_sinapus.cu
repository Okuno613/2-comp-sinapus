/*=======================================
  Since : May/19/2008
  Update: <2016/02/26>

float
  =======================================*/
#include "dendrite.h"
#include <omp.h>

__device__ double dv(double v ,double inp){
  return (-v+inp+V0)/TAU_recep;
}

__global__ double vrunge(double *v, double *inp, int i){
  double kv1 = DT*dv(v[i], inp[i]);
  double kv2 = DT*dv(v[i]+kv1*0.5, inp[i]);
  double kv3 = DT*dv(v[i]+kv2*0.5, inp[i]);
  double kv4 = DT*dv(v[i]+kv3,inp[i]);
  v[i] += (kv1 + 2.0*kv2 + 2.0*kv3 + kv4)/6.0;
}

__device__ double ds(double s,double o){
  return (-s+lambda*o)/TAU_sinapus;
}

__global__ double srunge(double *s, double o,int i){
  double ks1 = DT*ds(s[i], o);
  double ks2 = DT*ds(s[i]+ks1*0.5, o);
  double ks3 = DT*ds(s[i]+ks2*0.5, o);
  double ks4 = DT*ds(s[i]+ks3,o);
  s[i] += (ks1 + 2.0*ks2 + 2.0*ks3 + ks4)/6.0;
}

__device__ double dsegp(double A_egp,double Iz){ //true
  return (-A_egp + (w_egp*Iz*Iegp_sinapus) )/TAU_egp;
}

__global__ double segprunge(double *A_egp, double Iz,int tcnt){
  double ksegp1 = DT*dsegp(A_egp[tcnt], Iz);
  double ksegp2 = DT*dsegp(A_egp[tcnt]+ksegp1*0.5, Iz);
  double ksegp3 = DT*dsegp(A_egp[tcnt]+ksegp2*0.5, Iz);
  double ksegp4 = DT*dsegp(A_egp[tcnt]+ksegp3,Iz);
  A_egp[tcnt] += (ksegp1 + 2.0*ksegp2 + 2.0*ksegp3 + ksegp4)/6.0;
}

__device__ double dV_d(double V_d,double V_s,double h_d,double n_d,double p_d){
  double m_inf_d=1/(1+exp(-(V_d-V12_d)/k_m_inf_d));
  return( g_na_d* m_inf_d*m_inf_d *h_d * (V_na - V_d) + g_dr_d * n_d*n_d *p_d *(V_k -V_d)+( (g_c/(1-kappa))*(V_s-V_d)) + g_leak *(V_l - V_d));

}
      
__device__ double dh_d(double h_d, double V_d){
  double h_inf_d=1/(1+exp(-(V_d-V12_h_d)/k_h_d));
  return( (h_inf_d - h_d) /tau_h_d );
}

__device__ double dn_d(double n_d, double V_d){
  double n_inf_d=1/(1+exp(-(V_d-V12_n_d)/k_n_d));
  return( (n_inf_d - n_d) /tau_n_d );
}

__device__ double dp_d(double p_d, double V_d){
  double p_inf_d=1/(1+exp(-(V_d-V12_p_d)/k_p_d));
  return( (p_inf_d - p_d) /tau_p_d );
}
	    
__device__ double dV_s(double V_s, double inp, double V_d,double n_s){
  double m_inf_s=1/(1+exp(-(V_s-V12_s)/k_m_inf_s));
  return( inp + g_na_s * m_inf_s*m_inf_s * (1-n_s )* (V_na - V_s) +g_dr_s * n_s*n_s *(V_k -V_s)+ ((g_c/kappa)*(V_d-V_s)) +g_leak *(V_l-V_s));
}

__device__ double dn_s(double n_s, double V_s){
  double n_inf_s=1/(1+exp(-(V_s-V12_s)/k_m_inf_s));
  return( (n_inf_s - n_s) /tau_n_s );
}


__global__ double runge(double *u,double *V_s, double *n_s,double *V_d,double *h_d, double *n_d,double *p_d,int i){
  
  double kV_s1 = DT*dV_s(V_s[i], u[i], V_d[i],n_s[i]);
  double kn_s1 = DT*dn_s(n_s[i], V_s[i]);
  double kV_d1 = DT*dV_d(V_d[i], V_s[i] ,h_d[i], n_d[i], p_d[i]);
  double kh_d1 = DT*dh_d(h_d[i], V_d[i]);
  double kn_d1 = DT*dn_d(n_d[i], V_d[i]);
  double kp_d1 = DT*dp_d(p_d[i], V_d[i]);
  
  double kV_s2 = DT*dV_s(V_s[i]+kV_s1*0.5, u[i], V_d[i]+kV_d1*0.5 ,n_s[i]+kn_s1*0.5);
  double kn_s2 = DT*dn_s(n_s[i]+kn_s1*0.5, V_s[i]+kV_s1*0.5);
  double kV_d2 = DT*dV_d(V_d[i]+kV_d1*0.5, V_s[i]+kV_s1*0.5 ,h_d[i]+kh_d1*0.5, n_d[i]+kn_d1*0.5, p_d[i]+kp_d1*0.5);
  double kh_d2 = DT*dh_d(h_d[i]+kh_d1*0.5, V_d[i]+kV_d1*0.5);
  double kn_d2 = DT*dn_d(n_d[i]+kn_d1*0.5, V_d[i]+kV_d1*0.5);
  double kp_d2 = DT*dp_d(p_d[i]+kp_d1*0.5, V_d[i]+kV_d1*0.5);

  double kV_s3 = DT*dV_s(V_s[i]+kV_s2*0.5, u[i], V_d[i]+kV_d2*0.5 ,n_s[i]+kn_s2*0.5);
  double kn_s3 = DT*dn_s(n_s[i]+kn_s2*0.5, V_s[i]+kV_s2*0.5);
  double kV_d3 = DT*dV_d(V_d[i]+kV_d2*0.5, V_s[i]+kV_s2*0.5 ,h_d[i]+kh_d2*0.5, n_d[i]+kn_d2*0.5, p_d[i]+kp_d2*0.5);
  double kh_d3 = DT*dh_d(h_d[i]+kh_d2*0.5, V_d[i]+kV_d2*0.5);
  double kn_d3 = DT*dn_d(n_d[i]+kn_d2*0.5, V_d[i]+kV_d2*0.5);
  double kp_d3 = DT*dp_d(p_d[i]+kp_d2*0.5, V_d[i]+kV_d2*0.5);

  double kV_s4 = DT*dV_s(V_s[i]+kV_s3, u[i], V_d[i]+kV_d2 ,n_s[i]+kn_s2);
  double kn_s4 = DT*dn_s(n_s[i]+kn_s3, V_s[i]+kV_s3);
  double kV_d4 = DT*dV_d(V_d[i]+kV_d3, V_s[i]+kV_s3 ,h_d[i]+kh_d3, n_d[i]+kn_d3, p_d[i]+kp_d3);
  double kh_d4 = DT*dh_d(h_d[i]+kh_d3, V_d[i]+kV_d3);
  double kn_d4 = DT*dn_d(n_d[i]+kn_d3, V_d[i]+kV_d3);
  double kp_d4 = DT*dp_d(p_d[i]+kp_d3, V_d[i]+kV_d3);
  

  V_s[i] += (kV_s1 + 2.0*kV_s2 + 2.0*kV_s3 + kV_s4)/6.0;
  n_s[i] += (kn_s1 + 2.0*kn_s2 + 2.0*kn_s3 + kn_s4)/6.0;
  V_d[i] += (kV_d1 + 2.0*kV_d2 + 2.0*kV_d3 + kV_d4)/6.0;
  h_d[i] += (kh_d1 + 2.0*kh_d2 + 2.0*kh_d3 + kh_d4)/6.0;
  n_d[i] += (kn_d1 + 2.0*kn_d2 + 2.0*kn_d3 + kn_d4)/6.0;
  p_d[i] += (kp_d1 + 2.0*kp_d2 + 2.0*kp_d3 + kp_d4)/6.0;

  
}

  

__grobal__ void init(double *v,double *u,double *s,double *s_egp,double *A_egp,double *V_s, double *n_s,double *V_d,double *h_d, double *n_d,double *p_d, int *spike_s,int *spike_d, double *inp,double *THl,int *spikecnt_s,int *spikecnt_d,int *count_s,int *count_d)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  V_s[i] = -70;
  //n_s[i] = 1/(1+exp(-(-55-V12_n_d)/k_n_s));
  V_d[i] = -70.5;
  /*    h_d[i] = 1/(1+exp(-(-54.5-V12_n_d)/k_h_d));
	n_d[i] = 1/(1+exp(-(-54.5-V12_n_d)/k_n_d));
	p_d[i] = 1/(1+exp(-(-54.5-V12_n_d)/k_p_d));
  */
  inp[i] = 0;
  spike_s[i] = 0;
  spike_d[i] =0;
  spikecnt_s[i]=0;
  spikecnt_d[i]=0;
  count_s[i]=0;
  count_d[i]=0;
  THl[i]=TH;
  v[i]=0;
  u[i]=0;
  s[i]=0;
}
  
__grobal__ void init_egp(double *s_egp,double *A_egp)
{
  int tcnt = threadIdx.x + blockIdx.x * blockDim.x;
  s_egp[tcnt]=0;
  A_egp[tcnt]=0;
  }
}

__grobal__ void calv(double *v,double *u,double *s,double *s_egp,double *A_egp,double *V_s, double *n_s,double *V_d,double *h_d, double *n_d,double *p_d,  int *spike_s,int *spike_d, double *inp, double t,int *spikecnt_s,int *spikecnt_d,int *count_s,int *count_d,double *THl,FILE *fp1,FILE *fp2,FILE *fp3,FILE *fp4,FILE *fp5,FILE *fp6,FILE *fp7)
{
  int i0=(t/TEND)*NUM/2;
  int tcnt=int(t);
  //  int j0=sq/2;
  double Iz=0;
  //s_egp[tcnt]=0;
  //A_egp[tcnt]=0;

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  inp[i] = I0 * exp( -( ((i-i0)*(i-i0)) / (2*sigma*sigma)));//*sin(2.0*M_PI*t*(400/1000));
    

  vrunge(v,inp,i);
  if(v[i]>=20){
    v[i]=V0;
  }
    
  if(v[i] > TH){
    v[i] =30;
    srunge(s,10,i);
    spike_d[i] = spike_d[i]+1;
    spikecnt_d[i]=spike_d[i];
    // fprintf(fp2,"%d\t %d\t %d\n \n",i,int(t),spikecnt_d[i]);  
  }else{
    srunge(s,0,i);
  }

  if(i%2==0 or i==(NUM+1)/2){
    if (tcnt>11){
      if(i<2){
	u[i/2] = w_match*s[i] + w_plus*s[i+1] - w_minus*s[i+2]-w_egp_out*s_egp[tcnt-10];
      }else{
	u[i/2] = -w_minus*s[i-2] + w_plus*s[i-1] + w_match*s[i] + w_plus*s[i+1] - w_minus*s[i+2]-w_egp_out*s_egp[tcnt-10];
      }
      if( ((NUM/2)-2)<i){
	u[i/2] = -w_minus*s[i-2] + w_plus*s[i-1] + w_match*s[i]-w_egp_out*s_egp[tcnt-10];
      }else{
	u[i/2] = -w_minus*s[i-2] + w_plus*s[i-1] + w_match*s[i] + w_plus*s[i+1] - w_minus*s[i+2]-w_egp_out*s_egp[tcnt-10];
      }

      
    }else if(tcnt<11){
      if(i<2){
	u[i/2] = w_match*s[i] + w_plus*s[i+1] - w_minus*s[i+2];
      }else{
	u[i/2] = -w_minus*s[i-2] + w_plus*s[i-1] + w_match*s[i] + w_plus*s[i+1] - w_minus*s[i+2];
      }
	
      if( ((NUM/2)-2)<i){
	u[i/2] = -w_minus*s[i-2] + w_plus*s[i-1] + w_match*s[i];
      }else{
	u[i/2] = -w_minus*s[i-2] + w_plus*s[i-1] + w_match*s[i] + w_plus*s[i+1] - w_minus*s[i+2];
      }
    }
    	
    runge(u,V_s, n_s,V_d,h_d,n_d,p_d,i/2);
  }

    
  if(V_s[i] > 20 and count_s[i] == 0){
    spike_s[i] = spike_s[i]+1;
    //spikecnt_s[i]=spike_s[i];
    //     s[i]=s[i]-(s[i]/TAU_sinapus);
    //THl[i]=THl[i]+THup;
    count_s[i]=1;
    fprintf(fp1,"%d\t %d\t %03d\n",i,int(t),(spikecnt_d[i]-spikecnt_s[i]) );
    spikecnt_s[i]=int(t);     
  }
  spikecnt_d[i]=int(t);
      
  /*
    if(V_d[i] > 20 and count_d[i] == 0){
    spike_d[i] = spike_d[i]+1;
    spikecnt_d[i]=spike_d[i];
    count_d[i]=1;
    //THl[i]=THl[i]+THup;
    //fprintf(fp2,"%d\t %d\t %d\n \n",i,int(t),spikecnt_d[i] );
    }
  */

  if(int(t)%1==0){
    spike_s[i]=0;
    spike_d[i]=0;
  }
    
  if(count_d[i]==1 and V_d[i]<=-55){
    count_d[i]=0;
  }
  if(count_s[i]==1 and V_s[i]<=-55){
    count_s[i]=0;
  }
    

  
  fprintf(fp3,"%lf \t %lf \n",t,V_s[3]);  
  fprintf(fp4,"%lf \t %lf \n",t,V_d[3]);
  fprintf(fp6,"%lf \t %lf \n",t,v[3]);
  fprintf(fp2,"%lf \t %lf \n",t,s[3]*lambda);
}

__grobal__ void calegp(double *s_egp,double *A_egp,double t){

  int tcnt=int(t);
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  Iz = Iz+ count_s[i];
  
  
  segprunge(A_egp,Iz,int(tcnt));
  s_egp[tcnt]= 1/(1+exp((-(A_egp[tcnt])+theta_egp)/epshiron_egp));
  
  
  fprintf(fp7,"%lf \t %lf \n",t,s_egp[tcnt]);
  Iz=0;
}

void Simulation::sim()
{
    int count = 0.0;

    int size_d = sizeof(double)*NUM;
    int size_i = sizeof(int)*NUM;
    int size_t = sizeof(int)*int(TEND/DT)

    double *v =(double *)malloc(size_d);
    double *u =(double *)malloc(size_d);
    double *s =(double *)malloc(size_d);

    
    double *V_s =(double *)malloc(size_d);
    double *n_s =(double *)malloc(size_d);
    double *V_d =(double *)malloc(size_d);
    double *h_d =(double *)malloc(size_d);
    double *n_d =(double *)malloc(size_d);
    double *p_d =(double *)malloc(size_d);
    double *inp =(double *)malloc(size_d);
    cudaMalloc((void **)&d_V_s, size_d);
    cudaMalloc((void **)&d_n_s, size_d);
    cudaMalloc((void **)&d_V_d, size_d);
    cudaMalloc((void **)&d_h_d, size_d);
    cudaMalloc((void **)&d_n_d, size_d);
    cudaMalloc((void **)&d_p_d, size_d);
    cudaMalloc((void **)&d_inp, size_d);


    //    int sq = sqrt(NUM);

    int *count_s=(int *) malloc(size_i);;
    int *count_d=(int *) malloc(size_i);;
    count_s = (int *)malloc(size_i);
    count_d = (int *)malloc(size_i);
    cudaMalloc((void **)&d_count_s, size_i);
    cudaMalloc((void **)&d_count_d, size_i);

    double t = 0.0;
    double *d_t;
    cudaMalloc((void **)&d_t, size_d);
    
    int *spike_s= (int *) malloc(size_i);;
    int *spikecnt_s = (int *) malloc(size_i);;
    int *spike_d= (int *) malloc(size_i);;
    int *spikecnt_d = (int *) malloc(size_i);;
    double *THl =(double *)malloc(size_d);

    spike_s = (int *)malloc(size_i);
    spikecnt_s = (int *)malloc(size_i);
    spike_d = (int *)malloc(size_i);
    spikecnt_d = (int *)malloc(size_i);
    THl = (double *)malloc(size_d);
    cudaMalloc((void **)&d_spike_s, size_i);
    cudaMalloc((void **)&d_spikecnt_s, size_i);
    cudaMalloc((void **)&d_spike_d, size_i);
    cudaMalloc((void **)&d_spikecnt_d, size_i);
    cudaMalloc((void **)&d_THl, size_d);

    
    double *A_egp =(double *)malloc(size_t);
    double *s_egp =(double *)malloc(size_t);
    cudaMalloc((void **)&d_A_egp, size_t);
    cudaMalloc((void **)&d_s_egp, size_t);    

    FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7;
    fp1=fopen("Vs_moved.txt","w");
    fp2=fopen("s.txt","w");
    fp3=fopen("Vs_volt.txt","w");
    fp4=fopen("Vd_volt.txt","w");
    fp5=fopen("ns.txt","w");
    fp6=fopen("v_volt.txt","w");
    fp7=fopen("segp.txt","w");

    init<<<NUM/ Threads, Threads>>(v,u,s,s_egp,A_egp,V_s,n_s, V_d, h_d, n_d, p_d, spike_s,spike_d, inp, THl,spikecnt_s,spikecnt_d,count_s,count_d);
    init_egp<<<NUM/ Threads, Threads>>(s_egp,A_egp);
    for(count=0;;count++){
      
      calv<<<NUM/ Threads, Threads>>(v,u,s,s_egp,A_egp,V_s, n_s, V_d, h_d, n_d, p_d, spike_s,spike_d,inp, t, spikecnt_s,spikecnt_d,count_s,count_d,THl,fp1,fp2,fp3,fp4,fp5,fp6,fp7);

      cudaMemcpy(v,d_v, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(u,d_u, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(s,d_s, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(V_s,d_V_s, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(n_s,d_n_s, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(V_d,d_V_d, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(h_d,d_h_d, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(n_d,d_n_d, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(p_d,d_p_d, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(inp,d_inp, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(count_s,d_count_s, size_i, cudaMemcpyDeviceToHost);
      cudaMemcpy(spike_s,d_spike_s, size_i, cudaMemcpyDeviceToHost);
      cudaMemcpy(spikecnt_s,d_spikecnt_s, size_i, cudaMemcpyDeviceToHost);
      cudaMemcpy(count_d,d_count_d, size_i, cudaMemcpyDeviceToHost);
      cudaMemcpy(spike_d,d_spike_d, size_i, cudaMemcpyDeviceToHost);
      cudaMemcpy(spikecnt_d,d_spikecnt_d, size_i, cudaMemcpyDeviceToHost);
      cudaMemcpy(s_egp,d_s_egp, size_t, cudaMemcpyDeviceToHost);
      cudaMemcpy(A_egp,d_s_egp, size_t, cudaMemcpyDeviceToHost);

      
      calegp<<<NUM/ Threads, Threads>>(s_egp,A_egp,t);
      cudaMemcpy(s_egp,d_s_egp, size_t, cudaMemcpyDeviceToHost);
      cudaMemcpy(A_egp,d_s_egp, size_t, cudaMemcpyDeviceToHost);
      
                  
      t = count * DT;
      if( t > TPERIOD){
	break;
      }
    }

    cudaMemcpy(V_s,d_V_s, size_d, cudaMemcpyDeviceToHost);
    cudaMemcpy(n_s,d_n_s, size_d, cudaMemcpyDeviceToHost);
    cudaMemcpy(V_d,d_V_d, size_d, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_d,d_h_d, size_d, cudaMemcpyDeviceToHost);
    cudaMemcpy(n_d,d_n_d, size_d, cudaMemcpyDeviceToHost);
    cudaMemcpy(p_d,d_p_d, size_d, cudaMemcpyDeviceToHost);
    cudaMemcpy(inp,d_inp, size_d, cudaMemcpyDeviceToHost);
    cudaMemcpy(count_s,d_count_s, size_i, cudaMemcpyDeviceToHost);
    cudaMemcpy(spike_s,d_spike_s, size_i, cudaMemcpyDeviceToHost);
    cudaMemcpy(spikecnt_s,d_spikecnt_s, size_i, cudaMemcpyDeviceToHost);
    cudaMemcpy(count_d,d_count_d, size_i, cudaMemcpyDeviceToHost);
    cudaMemcpy(spike_d,d_spike_d, size_i, cudaMemcpyDeviceToHost);
    cudaMemcpy(spikecnt_d,d_spikecnt_d, size_i, cudaMemcpyDeviceToHost);
    
    free(V_s);
    free(n_s);
    free(V_d);
    free(h_d);
    free(n_d);
    free(p_d);
    cudaFree(d_V_s);
    cudaFree(d_n_s);
    cudaFree(d_V_d);
    cudaFree(d_h_d);
    cudaFree(d_n_d);
    cudaFree(d_p_d);
    free(count_s);
    free(count_d);
    cudaFree(d_count_s);
    cudaFree(d_count_d);
    free(spike_s);
    free(spike_d);
    free(spikecnt_s);
    free(spikecnt_d);
    cudaFree(THl);
    cudaFree(d_spike_s);
    cudaFree(d_spike_d);
    cudaFree(d_spikecnt_s);
    cudaFree(d_spikecnt_d);
    cudaFree(d_t);
    free(A_egp);
    free(s_egp);
    cudaFree(A_egp);
    cudaFree(s_egp);
    
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
    fclose(fp5);
    fclose(fp6);
    fclose(fp7);
}

 
 
int main(int argc, char* argv[]){
 Simulation sim;
 sim.sim();
 return(0);
 }
