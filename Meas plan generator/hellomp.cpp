#include <stdio.h>
#include <iostream>

#include<omp.h>


int main(){

    int thread_id;

    #pragma omp parallel
    {
        printf("Hello from process: %d\n", omp_get_thread_num());
    }
    return 0;
}