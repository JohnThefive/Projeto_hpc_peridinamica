# Estratégias de Paralelização OpenMP Offloading (GPU) no Código Peridinâmico (`openmp_offloading.f90`)

Este documento apresenta uma análise técnica da versão com aceleração heterogênea em **GPU via OpenMP Offloading**.
Esta análise complementa a documentação da versão CPU ([`openmp_onnode.f90`])
---

## Sumário Comparativo: Versão CPU (`onnode`) vs Versão GPU (`offloading`)

| Aspecto | Versão CPU (`openmp_onnode.f90`) | Versão GPU Offloading (`openmp_offloading.f90`) |
| :--- | :--- | :--- |
| **Ambiente de Execução** | Memória compartilhada unificada em núcleos de CPU (Host). | Arquitetura heterogênea Host (CPU) + Device (GPU dedicada). |
| **Diretivas Principais** | `!$OMP PARALLEL DO`, `!$OMP BARRIER`, `!$OMP SINGLE`. | `!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO`, `TARGET DATA`, `TARGET ENTER/EXIT DATA`. |
| **Residência de Memória** | Arrays sempre residem na RAM da CPU. | Mapeamentos explícitos na VRAM da GPU (`MAP(TO)`, `MAP(ALLOC)`, `UPDATE FROM`). |
| **Topologia - Fase 1 (Contagem)**| Paralelo na CPU (`PARALLEL DO`). | Paralelo na GPU via kernel offload (`TARGET TEAMS DISTRIBUTE PARALLEL DO`). |
| **Topologia - Fase 2 (Prefix Sum)**| **Paralelo na CPU** (Prefix sum com 2 fases + barreiras). | **Serial na CPU** (com transferência PCIe intermediária `UPDATE FROM`). |
| **Topologia - Fase 3 (Montagem)**| Paralelo na CPU (`PARALLEL DO`). | Paralelo na GPU (`TARGET TEAMS DISTRIBUTE PARALLEL DO`). |
| **Condições de Contorno (BC)** | Sequencial na CPU. | **Paralelizado na GPU** (`TARGET TEAMS DISTRIBUTE PARALLEL DO`). |
| **Cálculo de Forças e Dano** | Paralelo na CPU (`SCHEDULE(GUIDED)`). | **Paralelizado na GPU** (distribuído massivamente em *Teams* de threads). |
| **Integração Cinemática** | Sequencial na CPU (3 laços seriais). | **Paralelizado na GPU** (3 kernels distintos para evitar tráfego de `pforce` pela PCIe). |
| **Escrita em Disco (Snapshots)**| Acesso direto à memória RAM da CPU. | Sincronização seletiva sob demanda (`TARGET UPDATE FROM`) nos passos 750, 1000 e 1250. |

---

## 1. O que é Similar ou Idêntico à Versão CPU

Para evitar redundância com a documentação da versão CPU, destacam-se os pontos cuja formulação matemática e lógica algorítmica permanecem equivalentes:

1. **Formulação Centrada no Nó (*Gather Pattern*)**:
   - Assim como na CPU, cada iteração do laço de forças calcula estritamente as interações exercidas sobre o ponto material $i$.
   - **Motivação na GPU:** Se fosse utilizada a 3ª Lei de Newton com dispersão (*scatter*), seria obrigatório o uso de operações atômicas na memória global da GPU (`!$OMP ATOMIC`), o que provocaria serialização nos controladores de memória e colapso de desempenho (*memory pipeline stalls*). A abordagem *Gather* mantém a execução 100% livre de bloqueios.
2. **Representação CSR da Vizinhança**:
   - Utilização da trinca `numfam(i)`, `pointfam(i)` e do vetor unidimensional `nodefam(:)`.
3. **Escopo de Variáveis Locais (`PRIVATE`)**:
   - As variáveis auxiliares (`j, cnode, idist, nlength, fac, theta, scx, scy, scr, dforce1, dforce2, dmgpar1, dmgpar2`) continuam marcadas como `PRIVATE`, sendo alocadas nos registradores rápidos (*registers*) de cada thread da GPU.

---

## 2. Inovações e Estratégias Exclusivas da Versão GPU

### 2.1 Modelo de Execução Heterogênea (`TEAMS DISTRIBUTE PARALLEL DO`)

Na CPU, uma região paralela cria um conjunto plano de threads que compartilham a memória cache. Na GPU, a arquitetura de hardware é hierárquica (composta por blocos de multiprocessadores, ex.: *Streaming Multiprocessors* / *Compute Units*).

O OpenMP traduz isso através da diretiva combinada:

```fortran
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
!$OMP& PRIVATE(j, cnode, idist, nlength, fac, theta, scx, scy, scr, &
!$OMP&         dforce1, dforce2, dmgpar1, dmgpar2)
do i = 1, totnode
    ...
enddo
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
```

#### Decomposição dos Modificadores:
- **`TARGET`**: Transfere o fluxo de controle e o ambiente de execução da CPU (*Host*) para a GPU (*Device*).
- **`TEAMS`**: Cria uma liga de equipes de threads independentes (análogo aos *Thread Blocks* em CUDA ou *Work-Groups* em OpenCL/SYCL).
- **`DISTRIBUTE`**: Fatia o laço externo (`i = 1 .. totnode`) em grandes blocos e distribui esses blocos entre as equipes de threads (*teams*).
- **`PARALLEL DO`**: As threads dentro de cada equipe executam em paralelo as iterações atribuídas àquela equipe (execução em *SIMD/Warps/Sub-groups*).

---

### 2.2 Gerenciamento do Ciclo de Vida dos Dados (Data Lifetime & PCIe)

O custo de tráfego de dados pelo barramento PCIe é a principal fonte de gargalo em aceleração por GPU. O código emprega duas estratégias distintas para contornar esse problema:

```
                            FLUXO DE DADOS HOST <--> DEVICE
  
  [ CPU (Host) ]                                               [ GPU (Device VRAM) ]
        │                                                              │
        ├── TARGET ENTER DATA MAP(TO: coord, delta) MAP(ALLOC: numfam) ──► coord, delta, numfam
        │                                                              │
        │                        (Kernel Fase 1)                       │
        │   ◄── TARGET UPDATE FROM(numfam) ────────────────────────────┤ (Calcula numfam)
        │                                                              │
  (Aloca nodefam)                                                      │
        │                                                              │
        ├── TARGET ENTER DATA MAP(TO: pointfam) MAP(ALLOC: nodefam) ───► pointfam, nodefam
        │                                                              │
        │                        (Kernel Fase 3)                       │
        │   ◄── TARGET UPDATE FROM(nodefam) ───────────────────────────┤ (Preenche nodefam)
        │                                                              │
  (Setup Inicial)                                                      │
        │                                                              │
        ├── TARGET DATA MAP(TO: fail, disp...) MAP(ALLOC: pforce...) ──► Mantém tudo residente!
        │   │                                                          │
        │   ├── [ Loop Temporal 1..1250 ] (SEM CÓPIAS DE DADOS) ───────┤ Executa BC, Forças e
        │   │                                                          │ Integração 100% na GPU
        │   └── Snapshots (750/1000/1250):                             │
        │         ◄── TARGET UPDATE FROM(disp, dmg) ───────────────────┤ Envia só para salvar txt
        │                                                              │
        └── TARGET EXIT DATA MAP(RELEASE: ...) ────────────────────────► Desaloca VRAM
```

---

### 2.3 Pré-Processamento Híbrido com Diretivas Não-Estruturadas

Por que o código não usou um único bloco estruturado `!$OMP TARGET DATA` no pré-processamento?
- Porque o vetor `nodefam` **ainda não teve seu tamanho definitivo determinado** na inicialização.
- Se fosse alocado estaticamente com superestimativa (ex: 10.000.000 de inteiros), haveria desperdício massivo de VRAM e tempo de cópia inútil.

#### Solução: Diretivas Não-Estruturadas (`TARGET ENTER DATA`)
1. **Passo 1 (Carga Inicial)**:
   ```fortran
   !$OMP TARGET ENTER DATA MAP(TO: coord(1:totnode,:), delta) &
   !$OMP& MAP(ALLOC: numfam(1:totnode))
   ```
   - `MAP(TO: ...)`: Transfere `coord` e `delta` para a GPU apenas uma vez.
   - `MAP(ALLOC: ...)`: Aloca espaço para `numfam` na VRAM sem gastar largura de banda transferindo lixo da memória RAM.
2. **Passo 2 (Execução e Retorno da Contagem)**:
   - A GPU executa a Fase 1 da topologia.
   - O comando `!$OMP TARGET UPDATE FROM(numfam(1:totnode))` transfere **apenas** o vetor de contagem de volta para a CPU.
3. **Passo 3 (Intervenção Serial na CPU)**:
   - A CPU executa o prefix sum simples de `numfam` para gerar `pointfam`.
   - A CPU calcula `total_family_size = pointfam(totnode) + numfam(totnode) - 1` e aloca exatamente `nodefam(:)`.
4. **Passo 4 (Anexação Tardia na GPU)**:
   ```fortran
   !$OMP TARGET ENTER DATA MAP(TO: pointfam(1:totnode)) &
   !$OMP& MAP(ALLOC: nodefam(1:total_family_size))
   ```
   - Os novos arrays são mapeados para a GPU sem reinicializar ou perder `coord` e `delta`, que permaneceram residentes na memória da placa gráfica.
5. **Passo 5 (Preenchimento e Retorno Parcial)**:
   - A GPU executa a Fase 3 e preenche `nodefam`.
   - `TARGET UPDATE FROM(nodefam)` traz os dados para a CPU configurar a trinca e calibrações de superfície que rodam no Host.

---

### 2.4 Residência Permanente no Laço Temporal (`TARGET DATA` Global)

A maior otimização de desempenho da versão de offloading reside nas linhas 436–610:

```fortran
!$OMP TARGET DATA &
!$OMP& MAP(TO: fail(1:totnode,1:maxfam), fncst(1:totnode,:), bforce(1:totnode,:), &
!$OMP&         disp(1:totnode,:), vel(1:totnode,:)) &
!$OMP& MAP(ALLOC: pforce(1:totnode,:), acc(1:totnode,:), dmg(1:totnode))

do tt = 1, nt
    ...
enddo

!$OMP END TARGET DATA
```

#### Por que isso é determinante?
- Se a região de dados fosse aberta e fechada a cada iteração de `tt`, seriam realizadas **$1.250 \times 2 = 2.500$ transferências completas de malha** através do barramento PCIe.
- Mantendo uma **região global aberta sobre todo o laço temporal**, as variáveis persistem nos chips de memória GDDR/HBM da GPU entre um passo e outro.
- **Diferenciação Semântica dos Mapas:**
  - `MAP(TO: ...)`: Arrays cujos valores iniciais (definidos no pré-processamento na CPU) são necessários para o primeiro passo (`fail`, `fncst`, `bforce`, `disp`, `vel`).
  - `MAP(ALLOC: ...)`: Arrays que são integralmente sobrescritos a cada iteração (`pforce`, `acc`, `dmg`). O modificador `ALLOC` evita transferir gigabytes de dados desnecessários na inicialização.
  - **Escalares são `FIRSTPRIVATE` por Padrão**: Constantes escalares como `bc`, `vol`, `dt`, `radij`, `dens` são passadas por valor no cabeçalho do kernel, sem custo de alocação de mapa.

---

### 2.5 Paralelização das Etapas de Contorno e Cinemática na GPU

Ao contrário da versão CPU (onde as condições de contorno e a integração cinemática rodavam sequencialmente na thread mestre), a versão GPU **moveu todos esses laços para a GPU**:

#### 1. Condições de Contorno na GPU (`tempo_total_bc_gpu_s`):
```fortran
!$omp target teams distribute parallel do
do i = (totint+1), totbottom
    vel(i,2)  = -20.0d0
    disp(i,2) = -20.0d0 * tt * dt
enddo
!$omp end target teams distribute parallel do
```
- Embora o número de nós de borda seja pequeno ($\approx 3.000$ nós), rodar na GPU é vantajoso para manter os dados de `vel` e `disp` sincronizados no mesmo dispositivo de memória.

#### 2. Integração Cinemática Explícita na GPU:
```fortran
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
do i = 1, totint
    acc(i,1)  = (pforce(i,1) + bforce(i,1)) / dens
    ...
    disp(i,2) = disp(i,2) + vel(i,2) * dt
enddo
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
```

> **Justificativa Arquitetural Fundamental:**  
> Se a integração cinemática fosse mantida na CPU, a cada passo de tempo $tt$ seria obrigatório trazer `pforce` da GPU para a CPU e, em seguida, reenviar `disp` e `vel` da CPU para a GPU. Esse tráfego contínuo anularia completamente qualquer ganho de processamento obtido pela placa de vídeo. Ao transferir a cinemática para a GPU, o laço temporal executa **completamente fechado no dispositivo**.

---

### 2.6 Extração Sob Demanda de Resultados (`TARGET UPDATE FROM`)

A devolução de dados para a CPU ocorre apenas nos passos em que um arquivo de saída precisa ser gravado em disco:

```fortran
if (tt.eq.750 .or. tt.eq.1000 .or. tt.eq.1250) then
    ! Copia da GPU para a CPU apenas as variáveis e nós necessários para o snapshot
    !$OMP TARGET UPDATE FROM(disp(1:totint,:), dmg(1:totint))
    open(26, file = ...)
    do i = 1, totint
        write(26, ...) coord(i,1), coord(i,2), disp(i,1), disp(i,2), dmg(i)
    enddo
    close(26)
endif
```

---

## 3. Análise Crítica de Oportunidades de Otimização Específicas da GPU

Durante o desenvolvimento da versão para GPU, foram identificados alguns pontos que poderiam ser melhorados para reduzir o tempo de execução. As principais oportunidades encontradas foram:

### 3.1 Sobrecarga de lançamento de múltiplos kernels

A cada iteração do laço temporal `tt`, são executados seis kernels diferentes, relacionados às condições de contorno, ao cálculo das forças e às etapas de cinemática. Como são realizadas 1.250 iterações, isso resulta em 7.500 lançamentos de kernels.

Esse número elevado de lançamentos pode gerar uma sobrecarga considerável, principalmente nos kernels responsáveis pelas condições de contorno, que possuem poucas iterações. Uma possibilidade identificada foi juntar algumas dessas etapas em menos kernels, reduzindo a quantidade de lançamentos realizados durante a simulação.

### 3.2 Prefix Sum executado na CPU

Na implementação desenvolvida, a segunda etapa do pré-processamento, responsável pelo cálculo do `pointfam` a partir do `numfam`, ainda é realizada na CPU. Isso exige a transferência dos dados entre CPU e GPU antes que o processamento possa continuar.

Uma melhoria possível seria realizar também essa etapa na GPU. Dessa forma, o pré-processamento poderia ser executado inteiramente no dispositivo, evitando as transferências de dados utilizadas na implementação atual.

### 3.3 Acesso à memória

Outro ponto observado foi o padrão de acesso à memória durante o cálculo das interações entre os pontos. Alguns acessos são sequenciais e aproveitam melhor a organização dos dados, enquanto outros dependem da lista `nodefam` e acabam acessando posições diferentes da memória.

Uma possível melhoria seria reorganizar a forma como os pontos são armazenados e acessados, buscando manter pontos próximos fisicamente também próximos na memória. Essa abordagem poderia reduzir o custo dos acessos à memória durante o cálculo das interações.


## 4. Guia de Compilação para OpenMP Offloading

### 4.1 Compilação com Intel Fortran Compiler (`ifx`) para GPUs Intel
```bash
ifx -O3 -fiopenmp -fopenmp-targets=spir64 openmp_offloading.f90 -o perid_offload.exe
```

### 4.2 Compilação com NVIDIA HPC SDK (`nvfortran`) para GPUs NVIDIA
```bash
nvfortran -O3 -mp=gpu -gpu=cc80,cuda12.0 -Minfo=mp openmp_offloading.f90 -o perid_offload.exe
```

### 4.3 Variáveis de Diagnóstico em Tempo de Execução
Para inspecionar se os kernels e transferências de memória estão realmente sendo despachados para o dispositivo gráfico correto:

```bash
# Para compiladores baseados em LLVM / Intel
export LIBOMPTARGET_DEBUG=1
export LIBOMPTARGET_INFO=4

# Para verificar a GPU padrão ativa
export OMP_DEFAULT_DEVICE=0
```
