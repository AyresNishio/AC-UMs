#include <stdio.h>
#include <stdlib.h>

__device__ int d_c[3] = {3,2,1};
__global__ void hello_const(){
    printf("%i %i %i, ",d_c[0],d_c[1],d_c[2]);
}

int main(){
    int c[3] = {5,4,2};
    cudaMemcpyToSymbol(d_c,&c,3*sizeof(int),0,cudaMemcpyHostToDevice);
    hello_const<<<2,4>>>();
}