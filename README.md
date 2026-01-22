Este repositório contém implementações de um solver de **Peridinâmica baseada em ligas (Bond-based Peridynamics)** para simulação de mecânica da fratura em uma placa 2D. 

O objetivo principal deste projeto é realizar um estudo comparativo de desempenho (**benchmarking**) e validação numérica entre diferentes paradigmas de programação paralela utilizando Fortran moderno.

## Objetivos do Projeto

O código simula uma placa discretizada com aproximadamente 250.000 nós materiais (`ndivx=500`, `ndivy=500`), onde o desafio computacional principal reside na **Busca de Vizinhança (Neighbor Search)** e no cálculo de forças interpartículas dentro de um horizonte $\delta$.

O projeto consiste em 4 versões distintas do mesmo algoritmo:

1.  **Serial:** Versão de referência (Single-core) para validação da física e `baseline` de tempo.
2.  **OpenMP (CPU Multi-threading):** Paralelização de loops para processadores multi-core utilizando `!$omp parallel do`.
3.  **OpenMP (GPU Offloading):** Utilização de diretivas `!$omp target` para descarregar o processamento para a GPU (agnóstico de vendor).
4.  **CUDA Fortran:** Implementação nativa para GPUs NVIDIA utilizando kernels explícitos e gerenciamento de memória na GPU.