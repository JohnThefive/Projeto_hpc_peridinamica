// wrapper.cu
#include <thrust/scan.h>
#include <thrust/device_ptr.h>

extern "C" {
    // Função que o Fortran vai chamar
    void thrust_scan_wrapper(int *d_input, int *d_output, int n) {
        // Transforma ponteiros brutos em ponteiros do Thrust
        thrust::device_ptr<int> dev_in(d_input);
        thrust::device_ptr<int> dev_out(d_output);

        // Executa a soma acumulada (Exclusive Scan)
        // O '1' no final é o valor inicial (init) que seu código usa (pointfam(1)=1)
        thrust::exclusive_scan(dev_in, dev_in + n, dev_out, 1);
    }
}