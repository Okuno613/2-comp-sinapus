/*=======================================
  Since : May/19/2008
  Update: <2016/02/26>
  =======================================*/
#include "receptor.h"


__device__ double du(double u, double inp, double vrest){
    return (inp + vrest);
}

__device__ double runge(double u, double inp, double vrest){
    double k1 = DT*du(u, inp, vrest);
    double k2 = DT*du(u+k1*0.5, inp, vrest);
    double k3 = DT*du(u+k2*0.5, inp, vrest);
    double k4 = DT*du(u+k3, inp, vrest);
    u += (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
    return u;
}

__global__ void init(double *th, int *spike, double *inp)
{
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    th[index] = TH0;
    inp[index] = 0.0;
    spike[index] = 0;
}

__global__ void calth(double *th,  int *spike, double *inp, double *t)
{
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    inp[index] = sin(2.0*M_PI*t[0]/100.0);
    th[index] = runge(th[index], 0, - b/a);
    spike[index] = 0;

    if(inp[index] > th[index]){
        th[index] += b;
        spike[index] = 1;
    }
}

void Simulation::sim()
{
    int count = 0;

    double t = 0.0;
    double *d_t;
    cudaMalloc((void **)&d_t, sizeof(double));

    double *th, * inp;
    double *d_th, *d_inp;

    int size_d = sizeof(double)*NUM;
    th = (double *)malloc(size_d);
    inp = (double *)malloc(size_d);
    cudaMalloc((void **)&d_th, size_d);
    cudaMalloc((void **)&d_inp, size_d);

    int *spike;
    int *d_spike;
    int size_i = sizeof(int)*NUM;
    spike  = (int *)malloc(size_i);
    cudaMalloc((void **)&d_spike, size_i);

    init<<<NUM/Threads, Threads>>>(d_th, d_spike, d_inp);
    cudaMemcpy(d_t, &t, sizeof(int), cudaMemcpyHostToDevice);
    for(;;){
        calth<<<NUM/Threads, Threads>>>(d_th, d_spike, d_inp, d_t);

//        std::cout<<inp[0]<<" "<<th[0]<<"\n";

        count++;
        t = count * DT;
        if( t > TEND){
            break;
        }
//        cudaMemcpy(th, d_th, sizeof(double)*NUM, cudaMemcpyDeviceToHost);
//        cudaMemcpy(inp, d_inp, sizeof(double)*NUM, cudaMemcpyDeviceToHost);
        cudaMemcpy(d_t, &t, sizeof(double), cudaMemcpyHostToDevice);
    }
    cudaMemcpy(th, d_th, sizeof(double)*NUM, cudaMemcpyDeviceToHost);
    free(th);
    free(inp);
    free(spike);
    cudaFree(d_th);
    cudaFree(d_inp);
    cudaFree(d_spike);
}

int main(int argc, char* argv[]){
    Simulation sim;

    sim.sim();
    return(0);
}
