/*=======================================
  Since : May/19/2008
  Update: <2016/02/26>

float
  =======================================*/
#include "dendrite.h"
#include <omp.h>

double dv(double v ,double inp){
  return (-v+inp+V0)/TAU_recep;
}
double vrunge(double *v, double *inp, int i){
  double kv1 = DT*dv(v[i], inp[i]);
  double kv2 = DT*dv(v[i]+kv1*0.5, inp[i]);
  double kv3 = DT*dv(v[i]+kv2*0.5, inp[i]);
  double kv4 = DT*dv(v[i]+kv3,inp[i]);
  v[i] += (kv1 + 2.0*kv2 + 2.0*kv3 + kv4)/6.0;
}

double ds(double s,double o){
  return (-s+lambda*o)/TAU_sinapus;
}

double srunge(double *s, double o,int i){
  double ks1 = DT*ds(s[i], o);
  double ks2 = DT*ds(s[i]+ks1*0.5, o);
  double ks3 = DT*ds(s[i]+ks2*0.5, o);
  double ks4 = DT*ds(s[i]+ks3,o);
  s[i] += (ks1 + 2.0*ks2 + 2.0*ks3 + ks4)/6.0;
}

double dsegp(double A_egp,double Iz){ //true
  return (-A_egp + (w_egp*Iz*Iegp_sinapus) )/TAU_egp;
}

double segprunge(double *A_egp, double Iz,int tcnt){
  double ksegp1 = DT*dsegp(A_egp[tcnt], Iz);
  double ksegp2 = DT*dsegp(A_egp[tcnt]+ksegp1*0.5, Iz);
  double ksegp3 = DT*dsegp(A_egp[tcnt]+ksegp2*0.5, Iz);
  double ksegp4 = DT*dsegp(A_egp[tcnt]+ksegp3,Iz);
  A_egp[tcnt] += (ksegp1 + 2.0*ksegp2 + 2.0*ksegp3 + ksegp4)/6.0;
}

double dV_d(double V_d,double V_s,double h_d,double n_d,double p_d){
  double m_inf_d=1/(1+exp(-(V_d-V12_d)/k_m_inf_d));
  return( g_na_d* m_inf_d*m_inf_d *h_d * (V_na - V_d) + g_dr_d * n_d*n_d *p_d *(V_k -V_d)+( (g_c/(1-kappa))*(V_s-V_d)) + g_leak *(V_l - V_d));

}

double dh_d(double h_d, double V_d){
  double h_inf_d=1/(1+exp(-(V_d-V12_h_d)/k_h_d));
  return( (h_inf_d - h_d) /tau_h_d );
}

double dn_d(double n_d, double V_d){
  double n_inf_d=1/(1+exp(-(V_d-V12_n_d)/k_n_d));
  return( (n_inf_d - n_d) /tau_n_d );
}

double dp_d(double p_d, double V_d){
  double p_inf_d=1/(1+exp(-(V_d-V12_p_d)/k_p_d));
  return( (p_inf_d - p_d) /tau_p_d );
}

double dV_s(double V_s, double inp, double V_d,double n_s){
  double m_inf_s=1/(1+exp(-(V_s-V12_s)/k_m_inf_s));
  return( inp + g_na_s * m_inf_s*m_inf_s * (1-n_s )* (V_na - V_s) +g_dr_s * n_s*n_s *(V_k -V_s)+ ((g_c/kappa)*(V_d-V_s)) +g_leak *(V_l-V_s));
}

double dn_s(double n_s, double V_s){
  double n_inf_s=1/(1+exp(-(V_s-V12_s)/k_m_inf_s));
  return( (n_inf_s - n_s) /tau_n_s );
}


double runge(double *u,double *V_s, double *n_s,double *V_d,double *h_d, double *n_d,double *p_d,int i){

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



void init(double *v,double *u,double *s,double *s_egp,double *A_egp,double *V_s, double *n_s,double *V_d,double *h_d, double *n_d,double *p_d, int *spike_s,int *spike_d, double *inp,double *THl,int *spikecnt_s,int *spikecnt_d,int *spikecnt_v,int *count_s,int *count_d)
{
  #pragma omp parallel for
  for(int i=0;i<NUM;i++){
    V_s[i] = -70;
    //n_s[i] = 1/(1+exp(-(-55-V12_n_d)/k_n_s));
    V_d[i] = -70.5;
    /*    h_d[i] = 1/(1+exp(-(-54.5-V12_n_d)/k_h_d));
    n_d[i] = 1/(1+exp(-(-54.5-V12_n_d)/k_n_d));
    p_d[i] = 1/(1+exp(-(-54.5-V12_n_d)/k_p_d));
    a
    */
    inp[i] = -10;
    spike_s[i] = 0;
    spike_d[i] =0;
    spikecnt_s[i]=0;
    spikecnt_d[i]=0;
    spikecnt_v[i]=0;
    count_s[i]=0;
    count_d[i]=0;
    THl[i]=TH;
    v[i]=0;
    u[i]=0;
    s[i]=0;
  }
#pragma omp parallel for
  for(int tcnt=TSTART;tcnt<int(TEND/DT);tcnt++){
    s_egp[tcnt]=0;
    A_egp[tcnt]=0;
  }
}

void calv(double *v,double *u,double *s,double *s_egp,double *A_egp,double *V_s, double *n_s,double *V_d,double *h_d, double *n_d,double *p_d,  int *spike_s,int *spike_d, double *inp, double t,int *spikecnt_s,int *spikecnt_d,int *spikecnt_v,int *count_s,int *count_d,double *THl,FILE *fp1,FILE *fp2,FILE *fp3,FILE *fp4,FILE *fp5,FILE *fp6,FILE *fp7)
{
  int i0=NUM/4;
  //int i0=(t/TEND)*NUM/4;
  int tcnt=int(t);
  //  int j0=sq/2;
  double Iz=0;
  //s_egp[tcnt]=0;
  //A_egp[tcnt]=0;

  #pragma omp parallel for
  for(int i=0;i<NUM;i++){

    //inp[i] = I0 * exp( -( ((i-i0)*(i-i0)) / (2*sigma*sigma)));//*sin(2.0*M_PI*t*(400/1000));
    /*
    if(t<(TEND/2)){
      inp[i]=I0* (0.6+0.4*(TEND-t)/(TEND)) *  exp( -( ((i-(i0))*(i-(i0))) / (2*sigma*sigma)));
    }
    if(t>(TEND/2)){
      inp[i]=I0* (0.6+0.4*t/(TEND)) *  exp( -( ((i-(i0))*(i-(i0))) / (2*sigma*sigma)));
    }
    if(t==500){
      fprintf(fp2,"%d \t %lf \n",i,inp[i]);
    }
    */

    if(	0	<t and t<=	100	){	inp[i]=I0* ( (	0.690165 	 +t*	0.000005 	) *  exp( -( (i-(	95.22298	 +t*	0.000107 	))*(i-(	95.22298	 +t*	0.000107 	))) / (2*(	24.75849	 +t*	-0.000333 	 )*(	24.75849	 +t*	-0.000333 	) ))); }
  	if(	100	<t and t<=	200	){	inp[i]=I0* ( (	0.737317 	 +(t-100)*	0.000003 	) *  exp( -( (i-(	96.29371	 +(t-100)*	0.000292 	))*(i-(	96.29371	 +(t-100)*	0.000292 	))) / (2*(	21.42548	 +(t-100)*	-0.000102 	 )*(	21.42548	 +(t-100)*	-0.000102 	) ))); }
  	if(	200	<t and t<=	300	){	inp[i]=I0* ( (	0.762683 	 +(t-200)*	0.000006 	) *  exp( -( (i-(	99.21003	 +(t-200)*	-0.000063 	))*(i-(	99.21003	 +(t-200)*	-0.000063 	))) / (2*(	20.40085	 +(t-200)*	0.000251 	 )*(	20.40085	 +(t-200)*	0.000251 	) ))); }
  	if(	300	<t and t<=	400	){	inp[i]=I0* ( (	0.818083 	 +(t-300)*	0.000001 	) *  exp( -( (i-(	98.57951	 +(t-300)*	0.000000 	))*(i-(	98.57951	 +(t-300)*	0.000000 	))) / (2*(	22.90907	 +(t-300)*	-0.000340 	 )*(	22.90907	 +(t-300)*	-0.000340 	) ))); }
  	if(	400	<t and t<=	500	){	inp[i]=I0* ( (	0.825708 	 +(t-400)*	0.000001 	) *  exp( -( (i-(	98.57862	 +(t-400)*	-0.000115 	))*(i-(	98.57862	 +(t-400)*	-0.000115 	))) / (2*(	19.50725	 +(t-400)*	0.000336 	 )*(	19.50725	 +(t-400)*	0.000336 	) ))); }
  	if(	500	<t and t<=	600	){	inp[i]=I0* ( (	0.835512 	 +(t-500)*	0.000002 	) *  exp( -( (i-(	97.43015	 +(t-500)*	-0.000144 	))*(i-(	97.43015	 +(t-500)*	-0.000144 	))) / (2*(	22.86506	 +(t-500)*	0.000435 	 )*(	22.86506	 +(t-500)*	0.000435 	) ))); }
  	if(	600	<t and t<=	700	){	inp[i]=I0* ( (	0.851229 	 +(t-600)*	0.000004 	) *  exp( -( (i-(	95.98898	 +(t-600)*	-0.000189 	))*(i-(	95.98898	 +(t-600)*	-0.000189 	))) / (2*(	27.21644	 +(t-600)*	-0.000448 	 )*(	27.21644	 +(t-600)*	-0.000448 	) ))); }
  	if(	700	<t and t<=	800	){	inp[i]=I0* ( (	0.887955 	 +(t-700)*	0.000005 	) *  exp( -( (i-(	91.27445	 +(t-700)*	-0.000471 	))*(i-(	91.27445	 +(t-700)*	-0.000471 	))) / (2*(	19.5221	 +(t-700)*	-0.000769 	 )*(	19.5221	 +(t-700)*	-0.000769 	) ))); }
  	if(	800	<t and t<=	900	){	inp[i]=I0* ( (	0.904451 	 +(t-800)*	0.000002 	) *  exp( -( (i-(	91.32518	 +(t-800)*	0.000005 	))*(i-(	91.32518	 +(t-800)*	0.000005 	))) / (2*(	26.1434	 +(t-800)*	0.000662 	 )*(	26.1434	 +(t-800)*	0.000662 	) ))); }
  	if(	900	<t and t<=	1000	){	inp[i]=I0* ( (	0.919546 	 +(t-900)*	-0.000001 	) *  exp( -( (i-(	89.96951	 +(t-900)*	-0.000136 	))*(i-(	89.96951	 +(t-900)*	-0.000136 	))) / (2*(	28.62637	 +(t-900)*	0.000248 	 )*(	28.62637	 +(t-900)*	0.000248 	) ))); }
  	if(	1000	<t and t<=	1100	){	inp[i]=I0* ( (	0.913321 	 +(t-1000)*	0.000004 	) *  exp( -( (i-(	86.24475	 +(t-1000)*	-0.000372 	))*(i-(	86.24475	 +(t-1000)*	-0.000372 	))) / (2*(	23.88185	 +(t-1000)*	-0.000474 	 )*(	23.88185	 +(t-1000)*	-0.000474 	) ))); }
  	if(	1100	<t and t<=	1200	){	inp[i]=I0* ( (	0.955026 	 +(t-1100)*	0.000000 	) *  exp( -( (i-(	85.67346	 +(t-1100)*	-0.000057 	))*(i-(	85.67346	 +(t-1100)*	-0.000057 	))) / (2*(	26.86954	 +(t-1100)*	0.000299 	 )*(	26.86954	 +(t-1100)*	0.000299 	) ))); }
  	if(	1200	<t and t<=	1300	){	inp[i]=I0* ( (	0.956116 	 +(t-1200)*	0.000000 	) *  exp( -( (i-(	83.19856	 +(t-1200)*	-0.000247 	))*(i-(	83.19856	 +(t-1200)*	-0.000247 	))) / (2*(	31.18158	 +(t-1200)*	0.000431 	 )*(	31.18158	 +(t-1200)*	0.000431 	) ))); }
  	if(	1300	<t and t<=	1400	){	inp[i]=I0* ( (	0.954560 	 +(t-1300)*	0.000003 	) *  exp( -( (i-(	80.06755	 +(t-1300)*	-0.000313 	))*(i-(	80.06755	 +(t-1300)*	-0.000313 	))) / (2*(	28.80237	 +(t-1300)*	-0.000238 	 )*(	28.80237	 +(t-1300)*	-0.000238 	) ))); }
  	if(	1400	<t and t<=	1500	){	inp[i]=I0* ( (	0.985528 	 +(t-1400)*	0.000000 	) *  exp( -( (i-(	75.7286	 +(t-1400)*	-0.000434 	))*(i-(	75.7286	 +(t-1400)*	-0.000434 	))) / (2*(	28.82731	 +(t-1400)*	0.000002 	 )*(	28.82731	 +(t-1400)*	0.000002 	) ))); }
  	if(	1500	<t and t<=	1600	){	inp[i]=I0* ( (	0.983193 	 +(t-1500)*	0.000000 	) *  exp( -( (i-(	72.03123	 +(t-1500)*	-0.000370 	))*(i-(	72.03123	 +(t-1500)*	-0.000370 	))) / (2*(	28.36692	 +(t-1500)*	-0.000046 	 )*(	28.36692	 +(t-1500)*	-0.000046 	) ))); }
  	if(	1600	<t and t<=	1700	){	inp[i]=I0* ( (	0.985216 	 +(t-1600)*	0.000002 	) *  exp( -( (i-(	71.78178	 +(t-1600)*	-0.000025 	))*(i-(	71.78178	 +(t-1600)*	-0.000025 	))) / (2*(	28.58217	 +(t-1600)*	0.000022 	 )*(	28.58217	 +(t-1600)*	0.000022 	) ))); }
  	if(	1700	<t and t<=	1800	){	inp[i]=I0* ( (	1.003424 	 +(t-1700)*	0.000000 	) *  exp( -( (i-(	68.32759	 +(t-1700)*	-0.000345 	))*(i-(	68.32759	 +(t-1700)*	-0.000345 	))) / (2*(	28.89617	 +(t-1700)*	0.000031 	 )*(	28.89617	 +(t-1700)*	0.000031 	) ))); }
  	if(	1800	<t and t<=	1900	){	inp[i]=I0* ( (	1.000000 	 +(t-1800)*	-0.000005 	) *  exp( -( (i-(	66.15402	 +(t-1800)*	-0.000217 	))*(i-(	66.15402	 +(t-1800)*	-0.000217 	))) / (2*(	31.90236	 +(t-1800)*	0.000301 	 )*(	31.90236	 +(t-1800)*	0.000301 	) ))); }
  	if(	1900	<t and t<=	2000	){	inp[i]=I0* ( (	0.950825 	 +(t-1900)*	-0.000004 	) *  exp( -( (i-(	62.22859	 +(t-1900)*	-0.000393 	))*(i-(	62.22859	 +(t-1900)*	-0.000393 	))) / (2*(	28.74841	 +(t-1900)*	-0.000315 	 )*(	28.74841	 +(t-1900)*	-0.000315 	) ))); }
  	if(	2000	<t and t<=	2100	){	inp[i]=I0* ( (	0.909742 	 +(t-2000)*	-0.000002 	) *  exp( -( (i-(	60.12451	 +(t-2000)*	-0.000210 	))*(i-(	60.12451	 +(t-2000)*	-0.000210 	))) / (2*(	26.91961	 +(t-2000)*	-0.000183 	 )*(	26.91961	 +(t-2000)*	-0.000183 	) ))); }
  	if(	2100	<t and t<=	2200	){	inp[i]=I0* ( (	0.886399 	 +(t-2100)*	-0.000007 	) *  exp( -( (i-(	57.40345	 +(t-2100)*	-0.000272 	))*(i-(	57.40345	 +(t-2100)*	-0.000272 	))) / (2*(	28.22486	 +(t-2100)*	0.000131 	 )*(	28.22486	 +(t-2100)*	0.000131 	) ))); }
  	if(	2200	<t and t<=	2300	){	inp[i]=I0* ( (	0.818083 	 +(t-2200)*	-0.000004 	) *  exp( -( (i-(	55.42279	 +(t-2200)*	-0.000198 	))*(i-(	55.42279	 +(t-2200)*	-0.000198 	))) / (2*(	27.4255	 +(t-2200)*	-0.000080 	 )*(	27.4255	 +(t-2200)*	-0.000080 	) ))); }
  	if(	2300	<t and t<=	2400	){	inp[i]=I0* ( (	0.779490 	 +(t-2300)*	-0.000006 	) *  exp( -( (i-(	52.05361	 +(t-2300)*	-0.000337 	))*(i-(	52.05361	 +(t-2300)*	-0.000337 	))) / (2*(	28.57791	 +(t-2300)*	0.000115 	 )*(	28.57791	 +(t-2300)*	0.000115 	) ))); }
  	if(	2400	<t and t<=	2500	){	inp[i]=I0* ( (	0.716931 	 +(t-2400)*	-0.000004 	) *  exp( -( (i-(	51.18429	 +(t-2400)*	-0.000087 	))*(i-(	51.18429	 +(t-2400)*	-0.000087 	))) / (2*(	27.32592	 +(t-2400)*	-0.000125 	 )*(	27.32592	 +(t-2400)*	-0.000125 	) ))); }
  	if(	2500	<t and t<=	2600	){	inp[i]=I0* ( (	0.679116 	 +(t-2500)*	-0.000005 	) *  exp( -( (i-(	49.85307	 +(t-2500)*	-0.000133 	))*(i-(	49.85307	 +(t-2500)*	-0.000133 	))) / (2*(	27.79217	 +(t-2500)*	0.000047 	 )*(	27.79217	 +(t-2500)*	0.000047 	) ))); }
  	if(	2600	<t and t<=	2700	){	inp[i]=I0* ( (	0.633364 	 +(t-2600)*	-0.000006 	) *  exp( -( (i-(	49.83444	 +(t-2600)*	-0.000002 	))*(i-(	49.83444	 +(t-2600)*	-0.000002 	))) / (2*(	26.70638	 +(t-2600)*	-0.000109 	 )*(	26.70638	 +(t-2600)*	-0.000109 	) ))); }
  	if(	2700	<t and t<=	2800	){	inp[i]=I0* ( (	0.568783 	 +(t-2700)*	-0.000005 	) *  exp( -( (i-(	47.88513	 +(t-2700)*	-0.000195 	))*(i-(	47.88513	 +(t-2700)*	-0.000195 	))) / (2*(	29.3359	 +(t-2700)*	0.000263 	 )*(	29.3359	 +(t-2700)*	0.000263 	) ))); }
  	if(	2800	<t and t<=	2900	){	inp[i]=I0* ( (	0.520853 	 +(t-2800)*	-0.000005 	) *  exp( -( (i-(	47.616	 +(t-2800)*	-0.000027 	))*(i-(	47.616	 +(t-2800)*	-0.000027 	))) / (2*(	29.93584	 +(t-2800)*	0.000060 	 )*(	29.93584	 +(t-2800)*	0.000060 	) ))); }
  	if(	2900	<t and t<=	3000	){	inp[i]=I0* ( (	0.466542 	 +(t-2900)*	-0.000005 	) *  exp( -( (i-(	48.02598	 +(t-2900)*	0.000041 	))*(i-(	48.02598	 +(t-2900)*	0.000041 	))) / (2*(	26.57341	 +(t-2900)*	-0.000336 	 )*(	26.57341	 +(t-2900)*	-0.000336 	) ))); }
  	if(	3000	<t and t<=	3100	){	inp[i]=I0* ( (	0.419701 	 +(t-3000)*	-0.000007 	) *  exp( -( (i-(	54.16288	 +(t-3000)*	0.000614 	))*(i-(	54.16288	 +(t-3000)*	0.000614 	))) / (2*(	32.95787	 +(t-3000)*	0.000638 	 )*(	32.95787	 +(t-3000)*	0.000638 	) ))); }
  	if(	3100	<t and t<=	3200	){	inp[i]=I0* ( (	0.353252 	 +(t-3100)*	-0.000002 	) *  exp( -( (i-(	61.20131	 +(t-3100)*	0.000704 	))*(i-(	61.20131	 +(t-3100)*	0.000704 	))) / (2*(	30.37332	 +(t-3100)*	-0.000258 	 )*(	30.37332	 +(t-3100)*	-0.000258 	) ))); }
  	if(	3200	<t and t<=	3300	){	inp[i]=I0* ( (	0.330843 	 +(t-3200)*	0.000000 	) *  exp( -( (i-(	75.28752	 +(t-3200)*	0.001409 	))*(i-(	75.28752	 +(t-3200)*	0.001409 	))) / (2*(	44.24161	 +(t-3200)*	0.001387 	 )*(	44.24161	 +(t-3200)*	0.001387 	) ))); }
  	if(	3300	<t and t<=	3400	){	inp[i]=I0* ( (	0.329132 	 +(t-3300)*	0.000007 	) *  exp( -( (i-(	85.78446	 +(t-3300)*	0.001050 	))*(i-(	85.78446	 +(t-3300)*	0.001050 	))) / (2*(	24.30534	 +(t-3300)*	-0.001994 	 )*(	24.30534	 +(t-3300)*	-0.001994 	) ))); }
  	if(	3400	<t and t<=	3500	){	inp[i]=I0* ( (	0.397915 	 +(t-3400)*	0.000011 	) *  exp( -( (i-(	91.28252	 +(t-3400)*	0.000550 	))*(i-(	91.28252	 +(t-3400)*	0.000550 	))) / (2*(	26.89889	 +(t-3400)*	0.000259 	 )*(	26.89889	 +(t-3400)*	0.000259 	) ))); }
  	if(	3500	<t and t<=	3600	){	inp[i]=I0* ( (	0.508092 	 +(t-3500)*	0.000018 	) *  exp( -( (i-(	95.22298	 +(t-3500)*	0.000394 	))*(i-(	95.22298	 +(t-3500)*	0.000394 	))) / (2*(	24.75849	 +(t-3500)*	-0.000214 	 )*(	24.75849	 +(t-3500)*	-0.000214 	) ))); }


    vrunge(v,inp,i);
    if(v[i]>=20){
	     v[i]=V0;
    }
    if(v[i] > TH){
      v[i] =30;
      srunge(s,10,i);
      //spike_d[i] = spike_d[i]+1;
      //spikecnt_v[i]=spike_d[i];
      fprintf(fp4,"%d\t %d\t %03d\n",i,int(t),(spikecnt_d[i]-spikecnt_v[i]) );
      spikecnt_v[i]=int(t);
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

  }

  fprintf(fp3,"%lf \t %lf \n",t,V_s[3]);
  //fprintf(fp4,"%lf \t %lf \n",t,u[20]);
  fprintf(fp6,"%lf \t %lf \n",t,v[3]);
  //fprintf(fp2,"%lf \t %lf \n",t,inp[3]);
  //  fprintf(fp2,"%lf \t %lf \n",t, (g_c/kappa)*(V_d[0]-V_s[0]) );


#pragma omp parallel for
  for(int i=0;i<NUM;i++){
    Iz = Iz+ count_s[i];
    //spike_s[i]=0;
  }


  segprunge(A_egp,Iz,int(tcnt));


  s_egp[tcnt]= 1/(1+exp((-(A_egp[tcnt])+theta_egp)/epshiron_egp));


  fprintf(fp7,"%lf \t %lf \n",t,s_egp[tcnt]);
  Iz=0;
}

void Simulation::sim()
{
    int count = 0.0;


    double *v =new double[NUM];
    double *u =new double[NUM];
    double *s =new double[NUM];

    double *V_s =new double[NUM];
    double *n_s =new double[NUM];
    double *V_d =new double[NUM];
    double *h_d =new double[NUM];
    double *n_d =new double[NUM];
    double *p_d =new double[NUM];

    double *inp =new double[NUM];
    //    int sq = sqrt(NUM);

    int *count_s=new int[NUM];
    int *count_d=new int[NUM];

    double t = 0.0;
    int *spike_s= new int[NUM];
    int *spikecnt_s = new int[NUM];
    int *spike_d= new int[NUM];
    int *spikecnt_d = new int[NUM];
    int *spikecnt_v = new int[NUM];
    double *THl =new double[NUM];

    double *A_egp =new double[int(TEND/DT)];
    double *s_egp =new double[int(TEND/DT)];

    FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7;
    fp1=fopen("Vs_moved.txt","w");
    fp2=fopen("inp.txt","w");
    fp3=fopen("Vs_volt.txt","w");
    fp4=fopen("V_moved.txt","w");
    fp5=fopen("ns.txt","w");
    fp6=fopen("v_volt.txt","w");
    fp7=fopen("segp.txt","w");

    init(v,u,s,s_egp,A_egp,V_s,n_s, V_d, h_d, n_d, p_d, spike_s,spike_d, inp, THl,spikecnt_s,spikecnt_d,spikecnt_v,count_s,count_d);
    for(count=0;;count++){


      calv(v,u,s,s_egp,A_egp,V_s, n_s, V_d, h_d, n_d, p_d, spike_s,spike_d,inp, t, spikecnt_s,spikecnt_d,spikecnt_v,count_s,count_d,THl,fp1,fp2,fp3,fp4,fp5,fp6,fp7);

      t = count * DT;
      if( t > TPERIOD){
	break;
      }
    }
    free(V_s);
    free(n_s);
    free(V_d);
    free(h_d);
    free(n_d);
    free(p_d);
    free(count_s);
    free(count_d);
    free(spike_s);
    free(spike_d);
    free(spikecnt_s);
    free(spikecnt_d);

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
    fclose(fp5);
}



int main(int argc, char* argv[]){
 Simulation sim;
 sim.sim();
 return(0);
 }
