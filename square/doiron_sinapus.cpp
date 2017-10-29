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



void init(double *v,double *u,double *s,double *s_egp,double *A_egp,double *V_s, double *n_s,double *V_d,double *h_d, double *n_d,double *p_d, int *spike_s,int *spike_d, double *inp,double *THl,int *spikecnt_s,int *spikecnt_d,int *count_s,int *count_d)
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
#pragma omp parallel for
  for(int tcnt=TSTART;tcnt<int(TEND/DT);tcnt++){
    s_egp[tcnt]=0;
    A_egp[tcnt]=0;
  }
}

void calv(double *v,double *u,double *s,double *s_egp,double *A_egp,double *V_s, double *n_s,double *V_d,double *h_d, double *n_d,double *p_d,  int *spike_s,int *spike_d, double *inp, double t,int *spikecnt_s,int *spikecnt_d,int *count_s,int *count_d,double *THl,FILE *fp1,FILE *fp2,FILE *fp3,FILE *fp4,FILE *fp5,FILE *fp6,FILE *fp7)
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
 if(	0	<t<	100	){	inp[i]=I0* ( (	1	 +t*	-0.000001615686 	) *  exp( -( (i-(	45.71217	 +t*	-0.000023926000 	))*(i-(	45.71217	 +t*	-0.000023926000 	))) / (2*(	30.07002	 +t*	-0.000061631000 	 )*(	30.07002	 +t*	-0.000061631000 	) ))); }
 if(	100	<t<	200	){	inp[i]=I0* ( (	0.983843137	 +t*	-0.000004956863 	) *  exp( -( (i-(	45.47291	 +t*	-0.000055396000 	))*(i-(	45.47291	 +t*	-0.000055396000 	))) / (2*(	29.45371	 +t*	0.000029343000 	 )*(	29.45371	 +t*	0.000029343000 	) ))); }
 if(	200	<t<	300	){	inp[i]=I0* ( (	0.93427451	 +t*	-0.000002196078 	) *  exp( -( (i-(	44.91895	 +t*	0.000068512000 	))*(i-(	44.91895	 +t*	0.000068512000 	))) / (2*(	29.74714	 +t*	0.000008062000 	 )*(	29.74714	 +t*	0.000008062000 	) ))); }
 if(	300	<t<	400	){	inp[i]=I0* ( (	0.912313725	 +t*	-0.000001317647 	) *  exp( -( (i-(	45.60407	 +t*	-0.000157442000 	))*(i-(	45.60407	 +t*	-0.000157442000 	))) / (2*(	29.82776	 +t*	-0.000084933000 	 )*(	29.82776	 +t*	-0.000084933000 	) ))); }
 if(	400	<t<	500	){	inp[i]=I0* ( (	0.899137255	 +t*	-0.000001364706 	) *  exp( -( (i-(	44.02965	 +t*	-0.000195686000 	))*(i-(	44.02965	 +t*	-0.000195686000 	))) / (2*(	28.97843	 +t*	-0.000165377000 	 )*(	28.97843	 +t*	-0.000165377000 	) ))); }
 if(	500	<t<	600	){	inp[i]=I0* ( (	0.885490196	 +t*	-0.000003498039 	) *  exp( -( (i-(	42.07279	 +t*	-0.000078013000 	))*(i-(	42.07279	 +t*	-0.000078013000 	))) / (2*(	27.32466	 +t*	0.000019226000 	 )*(	27.32466	 +t*	0.000019226000 	) ))); }
 if(	600	<t<	700	){	inp[i]=I0* ( (	0.850509804	 +t*	-0.000002400000 	) *  exp( -( (i-(	41.29266	 +t*	-0.000136690000 	))*(i-(	41.29266	 +t*	-0.000136690000 	))) / (2*(	27.51692	 +t*	0.000190447000 	 )*(	27.51692	 +t*	0.000190447000 	) ))); }
 if(	700	<t<	800	){	inp[i]=I0* ( (	0.826509804	 +t*	-0.000003717647 	) *  exp( -( (i-(	39.92576	 +t*	-0.000330413000 	))*(i-(	39.92576	 +t*	-0.000330413000 	))) / (2*(	29.42139	 +t*	0.000114294000 	 )*(	29.42139	 +t*	0.000114294000 	) ))); }
 if(	800	<t<	900	){	inp[i]=I0* ( (	0.789333333	 +t*	-0.000001239216 	) *  exp( -( (i-(	36.62163	 +t*	-0.000123526000 	))*(i-(	36.62163	 +t*	-0.000123526000 	))) / (2*(	30.56433	 +t*	0.000234202000 	 )*(	30.56433	 +t*	0.000234202000 	) ))); }
 if(	900	<t<	1000	){	inp[i]=I0* ( (	0.776941176	 +t*	-0.000002274510 	) *  exp( -( (i-(	35.38637	 +t*	-0.000014315000 	))*(i-(	35.38637	 +t*	-0.000014315000 	))) / (2*(	32.90635	 +t*	-0.000334828000 	 )*(	32.90635	 +t*	-0.000334828000 	) ))); }
 if(	1000	<t<	1100	){	inp[i]=I0* ( (	0.754196078	 +t*	-0.000003560784 	) *  exp( -( (i-(	35.24322	 +t*	-0.000140959000 	))*(i-(	35.24322	 +t*	-0.000140959000 	))) / (2*(	29.55807	 +t*	0.000501315000 	 )*(	29.55807	 +t*	0.000501315000 	) ))); }
 if(	1100	<t<	1200	){	inp[i]=I0* ( (	0.718588235	 +t*	-0.000001443137 	) *  exp( -( (i-(	33.83363	 +t*	0.000006808000 	))*(i-(	33.83363	 +t*	0.000006808000 	))) / (2*(	34.57122	 +t*	0.000216583000 	 )*(	34.57122	 +t*	0.000216583000 	) ))); }
 if(	1200	<t<	1300	){	inp[i]=I0* ( (	0.704156863	 +t*	0.000000156863 	) *  exp( -( (i-(	33.90171	 +t*	0.000316350000 	))*(i-(	33.90171	 +t*	0.000316350000 	))) / (2*(	36.73705	 +t*	-0.000296699000 	 )*(	36.73705	 +t*	-0.000296699000 	) ))); }
 if(	1300	<t<	1400	){	inp[i]=I0* ( (	0.70572549	 +t*	0.000000329412 	) *  exp( -( (i-(	37.06521	 +t*	0.000130278000 	))*(i-(	37.06521	 +t*	0.000130278000 	))) / (2*(	33.77006	 +t*	0.000327440000 	 )*(	33.77006	 +t*	0.000327440000 	) ))); }
 if(	1400	<t<	1500	){	inp[i]=I0* ( (	0.709019608	 +t*	0.000001788235 	) *  exp( -( (i-(	38.36799	 +t*	0.000237319000 	))*(i-(	38.36799	 +t*	0.000237319000 	))) / (2*(	37.04446	 +t*	0.000112600000 	 )*(	37.04446	 +t*	0.000112600000 	) ))); }
 if(	1500	<t<	1600	){	inp[i]=I0* ( (	0.726901961	 +t*	0.000006541176 	) *  exp( -( (i-(	40.74118	 +t*	0.000369612000 	))*(i-(	40.74118	 +t*	0.000369612000 	))) / (2*(	38.17046	 +t*	-0.000427560000 	 )*(	38.17046	 +t*	-0.000427560000 	) ))); }
 if(	1600	<t<	1700	){	inp[i]=I0* ( (	0.792313725	 +t*	0.000007764706 	) *  exp( -( (i-(	44.4373	 +t*	0.000024745000 	))*(i-(	44.4373	 +t*	0.000024745000 	))) / (2*(	33.89486	 +t*	-0.000101505000 	 )*(	33.89486	 +t*	-0.000101505000 	) ))); }
 if(	1700	<t<	1800	){	inp[i]=I0* ( (	0.869960784	 +t*	0.000013003922 	) *  exp( -( (i-(	44.68475	 +t*	0.000102742000 	))*(i-(	44.68475	 +t*	0.000102742000 	))) / (2*(	32.87981	 +t*	-0.000280979000 	 )*(	32.87981	 +t*	-0.000280979000 	) ))); }
 if(	1800	<t<	1900	){	inp[i]=I0* ( (	1	 +t*	0.000000000000 	) *  exp( -( (i-(	45.71217	 +t*	0.000000000000 	))*(i-(	45.71217	 +t*	0.000000000000 	))) / (2*(	30.07002	 +t*	0.000000000000 	 )*(	30.07002	 +t*	0.000000000000 	) ))); }



    else{inp[i]=0;}





    vrunge(v,inp,i);
    if(v[i]>=20){
	v[i]=V0;
    }

    if(v[i] > TH){
      v[i] =30;
      srunge(s,10,i);
      spike_d[i] = spike_d[i]+1;
      spikecnt_d[i]=spike_d[i];
      fprintf(fp4,"%d\t %d\t %d\n \n",i,int(t),spikecnt_d[i]);
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

    init(v,u,s,s_egp,A_egp,V_s,n_s, V_d, h_d, n_d, p_d, spike_s,spike_d, inp, THl,spikecnt_s,spikecnt_d,count_s,count_d);
    for(count=0;;count++){


      calv(v,u,s,s_egp,A_egp,V_s, n_s, V_d, h_d, n_d, p_d, spike_s,spike_d,inp, t, spikecnt_s,spikecnt_d,count_s,count_d,THl,fp1,fp2,fp3,fp4,fp5,fp6,fp7);

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
