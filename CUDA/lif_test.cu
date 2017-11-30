/*=======================================
  Since : May/19/2008
  Update: <2015/11/24>
  =======================================*/
#include "lif_test.h"

__device__ double du(double u, double inp, double vrest){
    return (- u + inp + vrest)/TAU;
}

__device__ double runge(double u, double inp, double vrest){
    double k1 = DT*du(u, inp, vrest);
    double k2 = DT*du(u+k1*0.5, inp, vrest);
    double k3 = DT*du(u+k2*0.5, inp, vrest);
    double k4 = DT*du(u+k3, inp, vrest);
    u += (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
    return u;
}

__global__ void init(double *v, int *spike, double *inp)
{
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    v[index] = V0;
    inp[index] = 20.0;
    spike[index] = 0;
}

__global__ void calV(double *v,  int *spike, double *inp)
{
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    v[index] = runge(v[index], inp[index], V0);
    spike[index] = 0;

    if(v[index] > TH){
        v[index] = V0;
        spike[index] = 1;
    }
}

void Simulation::sim()
{
    int count = 0;
    double t = 0.0;

    double *v, * inp;
    double *d_v, *d_inp;

    int size_d = sizeof(double)*NUM;
    v = (double *)malloc(size_d);
    inp = (double *)malloc(size_d);
    cudaMalloc((void **)&d_v, size_d);
    cudaMalloc((void **)&d_inp, size_d);

    int *spike;
    int *d_spike;
    int size_i = sizeof(int)*NUM;
    spike  = (int *)malloc(size_i);
    cudaMalloc((void **)&d_spike, size_i);

    init<<<NUM/ 64, 64>>>(d_v, d_spike, d_inp);

    for(;;){
        calV<<<NUM/64, 64>>>(d_v, d_spike, d_inp);

//        std::cout<<v[0]<<"\n";

        count++;
        t = count * DT;
        if( t > TEND){
            break;
        }
    }
    cudaMemcpy(v, d_v, sizeof(double)*NUM, cudaMemcpyDeviceToHost);

    free(v);
    free(inp);
    free(spike);
    cudaFree(d_v);
    cudaFree(d_inp);
    cudaFree(d_spike);
}

int main(int argc, char* argv[]){
    Simulation sim;

    sim.sim();
    return(0);
}
