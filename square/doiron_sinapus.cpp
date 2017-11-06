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
    if(	0	<t and t<=	100	){	inp[i]=I0* ( (	0.891284 	 +t*	0.000000 	) *  exp( -( (i-(	59.27471	 +t*	0.000189 	))*(i-(	59.27471	 +t*	0.000189 	))) / (2*(	24.38624	 +t*	0.000424 	 )*(	24.38624	 +t*	0.000424 	) ))); }
  	if(	100	<t and t<=	200	){	inp[i]=I0* ( (	0.891100 	 +(t-100)*	0.000001 	) *  exp( -( (i-(	61.16131	 +(t-100)*	-0.000069 	))*(i-(	61.16131	 +(t-100)*	-0.000069 	))) / (2*(	28.62596	 +(t-100)*	-0.000200 	 )*(	28.62596	 +(t-100)*	-0.000200 	) ))); }
  	if(	200	<t and t<=	300	){	inp[i]=I0* ( (	0.901419 	 +(t-200)*	0.000003 	) *  exp( -( (i-(	60.47362	 +(t-200)*	0.000142 	))*(i-(	60.47362	 +(t-200)*	0.000142 	))) / (2*(	26.62462	 +(t-200)*	-0.000071 	 )*(	26.62462	 +(t-200)*	-0.000071 	) ))); }
  	if(	300	<t and t<=	400	){	inp[i]=I0* ( (	0.927584 	 +(t-300)*	0.000004 	) *  exp( -( (i-(	61.89209	 +(t-300)*	-0.000345 	))*(i-(	61.89209	 +(t-300)*	-0.000345 	))) / (2*(	25.91791	 +(t-300)*	-0.000452 	 )*(	25.91791	 +(t-300)*	-0.000452 	) ))); }
  	if(	400	<t and t<=	500	){	inp[i]=I0* ( (	0.969965 	 +(t-400)*	0.000001 	) *  exp( -( (i-(	58.43922	 +(t-400)*	0.000134 	))*(i-(	58.43922	 +(t-400)*	0.000134 	))) / (2*(	21.40212	 +(t-400)*	0.000418 	 )*(	21.40212	 +(t-400)*	0.000418 	) ))); }
  	if(	500	<t and t<=	600	){	inp[i]=I0* ( (	0.977520 	 +(t-500)*	-0.000001 	) *  exp( -( (i-(	59.77896	 +(t-500)*	-0.000121 	))*(i-(	59.77896	 +(t-500)*	-0.000121 	))) / (2*(	25.5784	 +(t-500)*	0.000007 	 )*(	25.5784	 +(t-500)*	0.000007 	) ))); }
  	if(	600	<t and t<=	700	){	inp[i]=I0* ( (	0.963884 	 +(t-600)*	0.000000 	) *  exp( -( (i-(	58.57392	 +(t-600)*	-0.000223 	))*(i-(	58.57392	 +(t-600)*	-0.000223 	))) / (2*(	25.65056	 +(t-600)*	-0.000053 	 )*(	25.65056	 +(t-600)*	-0.000053 	) ))); }
  	if(	700	<t and t<=	800	){	inp[i]=I0* ( (	0.967201 	 +(t-700)*	0.000002 	) *  exp( -( (i-(	54.94731	 +(t-700)*	-0.000363 	))*(i-(	54.94731	 +(t-700)*	-0.000363 	))) / (2*(	25.5331	 +(t-700)*	-0.000012 	 )*(	25.5331	 +(t-700)*	-0.000012 	) ))); }
  	if(	800	<t and t<=	900	){	inp[i]=I0* ( (	0.982863 	 +(t-800)*	-0.000001 	) *  exp( -( (i-(	53.73758	 +(t-800)*	-0.000121 	))*(i-(	53.73758	 +(t-800)*	-0.000121 	))) / (2*(	26.70461	 +(t-800)*	0.000117 	 )*(	26.70461	 +(t-800)*	0.000117 	) ))); }
  	if(	900	<t and t<=	1000	){	inp[i]=I0* ( (	0.972545 	 +(t-900)*	0.000002 	) *  exp( -( (i-(	52.61918	 +(t-900)*	-0.000112 	))*(i-(	52.61918	 +(t-900)*	-0.000112 	))) / (2*(	27.3978	 +(t-900)*	0.000069 	 )*(	27.3978	 +(t-900)*	0.000069 	) ))); }
  	if(	1000	<t and t<=	1100	){	inp[i]=I0* ( (	0.992261 	 +(t-1000)*	-0.000001 	) *  exp( -( (i-(	48.90762	 +(t-1000)*	-0.000371 	))*(i-(	48.90762	 +(t-1000)*	-0.000371 	))) / (2*(	27.32849	 +(t-1000)*	-0.000007 	 )*(	27.32849	 +(t-1000)*	-0.000007 	) ))); }
  	if(	1100	<t and t<=	1200	){	inp[i]=I0* ( (	0.986364 	 +(t-1100)*	0.000001 	) *  exp( -( (i-(	46.60458	 +(t-1100)*	-0.000230 	))*(i-(	46.60458	 +(t-1100)*	-0.000230 	))) / (2*(	30.67393	 +(t-1100)*	0.000335 	 )*(	30.67393	 +(t-1100)*	0.000335 	) ))); }
  	if(	1200	<t and t<=	1300	){	inp[i]=I0* ( (	0.992629 	 +(t-1200)*	-0.000001 	) *  exp( -( (i-(	43.61977	 +(t-1200)*	-0.000298 	))*(i-(	43.61977	 +(t-1200)*	-0.000298 	))) / (2*(	26.67837	 +(t-1200)*	-0.000400 	 )*(	26.67837	 +(t-1200)*	-0.000400 	) ))); }
  	if(	1300	<t and t<=	1400	){	inp[i]=I0* ( (	0.978257 	 +(t-1300)*	0.000000 	) *  exp( -( (i-(	42.37567	 +(t-1300)*	-0.000124 	))*(i-(	42.37567	 +(t-1300)*	-0.000124 	))) / (2*(	29.13253	 +(t-1300)*	0.000245 	 )*(	29.13253	 +(t-1300)*	0.000245 	) ))); }
  	if(	1400	<t and t<=	1500	){	inp[i]=I0* ( (	0.974940 	 +(t-1400)*	0.000002 	) *  exp( -( (i-(	38.69343	 +(t-1400)*	-0.000368 	))*(i-(	38.69343	 +(t-1400)*	-0.000368 	))) / (2*(	27.28079	 +(t-1400)*	-0.000185 	 )*(	27.28079	 +(t-1400)*	-0.000185 	) ))); }
  	if(	1500	<t and t<=	1600	){	inp[i]=I0* ( (	0.999631 	 +(t-1500)*	-0.000001 	) *  exp( -( (i-(	34.63195	 +(t-1500)*	-0.000406 	))*(i-(	34.63195	 +(t-1500)*	-0.000406 	))) / (2*(	30.39242	 +(t-1500)*	0.000311 	 )*(	30.39242	 +(t-1500)*	0.000311 	) ))); }
  	if(	1600	<t and t<=	1700	){	inp[i]=I0* ( (	0.990234 	 +(t-1600)*	-0.000001 	) *  exp( -( (i-(	31.37259	 +(t-1600)*	-0.000326 	))*(i-(	31.37259	 +(t-1600)*	-0.000326 	))) / (2*(	29.42133	 +(t-1600)*	-0.000097 	 )*(	29.42133	 +(t-1600)*	-0.000097 	) ))); }
  	if(	1700	<t and t<=	1800	){	inp[i]=I0* ( (	0.977520 	 +(t-1700)*	0.000002 	) *  exp( -( (i-(	29.99838	 +(t-1700)*	-0.000137 	))*(i-(	29.99838	 +(t-1700)*	-0.000137 	))) / (2*(	26.01006	 +(t-1700)*	-0.000341 	 )*(	26.01006	 +(t-1700)*	-0.000341 	) ))); }
  	if(	1800	<t and t<=	1900	){	inp[i]=I0* ( (	1.000000 	 +(t-1800)*	-0.000006 	) *  exp( -( (i-(	26.74168	 +(t-1800)*	-0.000326 	))*(i-(	26.74168	 +(t-1800)*	-0.000326 	))) / (2*(	27.39393	 +(t-1800)*	0.000138 	 )*(	27.39393	 +(t-1800)*	0.000138 	) ))); }
  	if(	1900	<t and t<=	2000	){	inp[i]=I0* ( (	0.942141 	 +(t-1900)*	-0.000006 	) *  exp( -( (i-(	25.98319	 +(t-1900)*	-0.000076 	))*(i-(	25.98319	 +(t-1900)*	-0.000076 	))) / (2*(	26.29298	 +(t-1900)*	-0.000110 	 )*(	26.29298	 +(t-1900)*	-0.000110 	) ))); }
  	if(	2000	<t and t<=	2100	){	inp[i]=I0* ( (	0.878939 	 +(t-2000)*	-0.000004 	) *  exp( -( (i-(	23.16812	 +(t-2000)*	-0.000282 	))*(i-(	23.16812	 +(t-2000)*	-0.000282 	))) / (2*(	26.00406	 +(t-2000)*	-0.000029 	 )*(	26.00406	 +(t-2000)*	-0.000029 	) ))); }
  	if(	2100	<t and t<=	2200	){	inp[i]=I0* ( (	0.836926 	 +(t-2100)*	-0.000007 	) *  exp( -( (i-(	20.79276	 +(t-2100)*	-0.000238 	))*(i-(	20.79276	 +(t-2100)*	-0.000238 	))) / (2*(	27.94151	 +(t-2100)*	0.000194 	 )*(	27.94151	 +(t-2100)*	0.000194 	) ))); }
  	if(	2200	<t and t<=	2300	){	inp[i]=I0* ( (	0.768196 	 +(t-2200)*	-0.000007 	) *  exp( -( (i-(	18.15388	 +(t-2200)*	-0.000264 	))*(i-(	18.15388	 +(t-2200)*	-0.000264 	))) / (2*(	27.66358	 +(t-2200)*	-0.000028 	 )*(	27.66358	 +(t-2200)*	-0.000028 	) ))); }
  	if(	2300	<t and t<=	2400	){	inp[i]=I0* ( (	0.701493 	 +(t-2300)*	-0.000007 	) *  exp( -( (i-(	15.97034	 +(t-2300)*	-0.000218 	))*(i-(	15.97034	 +(t-2300)*	-0.000218 	))) / (2*(	26.98634	 +(t-2300)*	-0.000068 	 )*(	26.98634	 +(t-2300)*	-0.000068 	) ))); }
  	if(	2400	<t and t<=	2500	){	inp[i]=I0* ( (	0.636079 	 +(t-2400)*	-0.000004 	) *  exp( -( (i-(	16.83395	 +(t-2400)*	0.000086 	))*(i-(	16.83395	 +(t-2400)*	0.000086 	))) / (2*(	24.8848	 +(t-2400)*	-0.000210 	 )*(	24.8848	 +(t-2400)*	-0.000210 	) ))); }
  	if(	2500	<t and t<=	2600	){	inp[i]=I0* ( (	0.592593 	 +(t-2500)*	-0.000007 	) *  exp( -( (i-(	12.54824	 +(t-2500)*	-0.000429 	))*(i-(	12.54824	 +(t-2500)*	-0.000429 	))) / (2*(	28.21471	 +(t-2500)*	0.000333 	 )*(	28.21471	 +(t-2500)*	0.000333 	) ))); }
  	if(	2600	<t and t<=	2700	){	inp[i]=I0* ( (	0.522020 	 +(t-2600)*	-0.000006 	) *  exp( -( (i-(	14.58176	 +(t-2600)*	0.000203 	))*(i-(	14.58176	 +(t-2600)*	0.000203 	))) / (2*(	28.58086	 +(t-2600)*	0.000037 	 )*(	28.58086	 +(t-2600)*	0.000037 	) ))); }
  	if(	2700	<t and t<=	2800	){	inp[i]=I0* ( (	0.461397 	 +(t-2700)*	-0.000005 	) *  exp( -( (i-(	13.63782	 +(t-2700)*	-0.000094 	))*(i-(	13.63782	 +(t-2700)*	-0.000094 	))) / (2*(	31.80619	 +(t-2700)*	0.000323 	 )*(	31.80619	 +(t-2700)*	0.000323 	) ))); }
  	if(	2800	<t and t<=	2900	){	inp[i]=I0* ( (	0.409619 	 +(t-2800)*	-0.000007 	) *  exp( -( (i-(	17.25287	 +(t-2800)*	0.000362 	))*(i-(	17.25287	 +(t-2800)*	0.000362 	))) / (2*(	36.48495	 +(t-2800)*	0.000468 	 )*(	36.48495	 +(t-2800)*	0.000468 	) ))); }
  	if(	2900	<t and t<=	3000	){	inp[i]=I0* ( (	0.339783 	 +(t-2900)*	-0.000004 	) *  exp( -( (i-(	26.6096	 +(t-2900)*	0.000936 	))*(i-(	26.6096	 +(t-2900)*	0.000936 	))) / (2*(	34.73353	 +(t-2900)*	-0.000175 	 )*(	34.73353	 +(t-2900)*	-0.000175 	) ))); }
  	if(	3000	<t and t<=	3100	){	inp[i]=I0* ( (	0.299982 	 +(t-3000)*	-0.000001 	) *  exp( -( (i-(	34.2642	 +(t-3000)*	0.000765 	))*(i-(	34.2642	 +(t-3000)*	0.000765 	))) / (2*(	36.98869	 +(t-3000)*	0.000226 	 )*(	36.98869	 +(t-3000)*	0.000226 	) ))); }
  	if(	3100	<t and t<=	3200	){	inp[i]=I0* ( (	0.291505 	 +(t-3100)*	0.000003 	) *  exp( -( (i-(	46.0694	 +(t-3100)*	0.001181 	))*(i-(	46.0694	 +(t-3100)*	0.001181 	))) / (2*(	31.39383	 +(t-3100)*	-0.000559 	 )*(	31.39383	 +(t-3100)*	-0.000559 	) ))); }
  	if(	3200	<t and t<=	3300	){	inp[i]=I0* ( (	0.324120 	 +(t-3200)*	0.000006 	) *  exp( -( (i-(	53.19747	 +(t-3200)*	0.000713 	))*(i-(	53.19747	 +(t-3200)*	0.000713 	))) / (2*(	33.40203	 +(t-3200)*	0.000201 	 )*(	33.40203	 +(t-3200)*	0.000201 	) ))); }
  	if(	3300	<t and t<=	3400	){	inp[i]=I0* ( (	0.383453 	 +(t-3300)*	0.000013 	) *  exp( -( (i-(	55.73193	 +(t-3300)*	0.000253 	))*(i-(	55.73193	 +(t-3300)*	0.000253 	))) / (2*(	27.04269	 +(t-3300)*	-0.000636 	 )*(	27.04269	 +(t-3300)*	-0.000636 	) ))); }
  	if(	3400	<t and t<=	3500	){	inp[i]=I0* ( (	0.514096 	 +(t-3400)*	0.000016 	) *  exp( -( (i-(	59.2389	 +(t-3400)*	0.000351 	))*(i-(	59.2389	 +(t-3400)*	0.000351 	))) / (2*(	27.77806	 +(t-3400)*	0.000074 	 )*(	27.77806	 +(t-3400)*	0.000074 	) ))); }
  	if(	3500	<t and t<=	3600	){	inp[i]=I0* ( (	0.676064 	 +(t-3500)*	0.000022 	) *  exp( -( (i-(	59.27471	 +(t-3500)*	0.000004 	))*(i-(	59.27471	 +(t-3500)*	0.000004 	))) / (2*(	24.38624	 +(t-3500)*	-0.000339 	 )*(	24.38624	 +(t-3500)*	-0.000339 	) ))); }


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
