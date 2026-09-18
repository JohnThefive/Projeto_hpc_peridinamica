# Estratégias de Paralelização em CUDA Fortran (`kernels.cuf` e `main.cuf`)

## Sumário Comparativo entre as Três Versões do Projeto

| Dimensão | OpenMP CPU (`onnode`) | OpenMP Offloading (`offloading`) | CUDA Fortran (`cuda`) |
| :--- | :--- | :--- | :--- |
| **Organização do Código** | Monolítico (1 arquivo `.f90`). | Monolítico (1 arquivo `.f90`). | **Modular**: separação estrita entre módulo de kernels (`kernels.cuf`) e orquestrador host (`main.cuf`). |
| **Paralelismo** | Diretivas de alto nível (`!$OMP PARALLEL DO`). | Diretivas de alto nível (`!$OMP TARGET TEAMS DISTRIBUTE`). | **Paralelismo Nativo GPU**: subrotinas `attributes(global)` com sintaxe de chevron `<<<grid, block>>>`. |
| **Gerenciamento de Memória** | Memória RAM única (Host). | Mapeamento declarativo abstrato (`MAP(TO/ALLOC)`). | **Ponteiros Explícitos de VRAM**: qualificadores `attributes(device)` e alocações diretas na GPU (`coord_d`, `d_disp`). |
| **Transferências de Dados** | Inexistente (memória compartilhada). | Diretivas abstratas de runtime (`TARGET UPDATE FROM/TO`). | **Atribuições Diretas Fortran**: `coord_d = coord` e `disp = d_disp` mapeadas para `cudaMemcpy` síncronas. |
| **Controle de Threads e Blocos** | Automático pelo runtime OpenMP. | Estimado heuristicamente pelo compilador (`teams`/`threads`). | **Configuração Determinística**: `dim3(256, 1, 1)` ajustado para múltiplos de *warps* de 32 threads. |
| **Kernels por Passo de Tempo** | 1 laço paralelo (Forças). | Originalmente 6 laços paralelizados na GPU. | **3 Kernels Unificados**: BC (1 kernel sob medida), Forças (1 kernel) e Cinemática (1 kernel fundido). |

---

## 1. O que Permanece Equivalente (Invariantes Físicos)

Para evitar repetições desnecessárias:
1. **Modelo Físico e Discretização**: Permanece a malha bidimensional de $N = 253.000$ pontos materiais, horizonte $\delta = 3{,}015\,\Delta x$, $N_t = 1.250$ passos de tempo e critério de falha de vínculo por estiramento crítico.
2. **Formulação Centrada no Nó (*Gather Pattern*)**: Cada thread global calcula exclusivamente as forças exercidas sobre o seu nó $i$ (`p_force(i,:)`), lendo os vizinhos `cnode` a partir da lista CSR. Não há necessidade de primitivas atômicas (`atomicAdd`) para a 3ª Lei de Newton.
3. **Representação CSR da Vizinhança**: Mantém-se o uso de `numfam` (contagem de vizinhos), `pointfam` (deslocamento inicial) e `nodefam` (lista plana de vizinhos).

---

## 2. Inovações e Diferenças Arquiteturais Exclusivas do CUDA Fortran

### 2.1 Arquitetura Modular e Subrotinas `attributes(global)`

Ao contrário das diretivas inseridas diretamente em laços das versões OpenMP, o CUDA Fortran isola a lógica executada na placa de vídeo dentro do módulo `kernels_peridinamica`:

```
                  ARQUITETURA DO SISTEMA CUDA FORTRAN
  
  ┌─────────────────────────────────────────────────────────────┐
  │                        main.cuf                             │
  │  • Alocação Host (RAM) e Device (VRAM)                      │
  │  • Dimensionamento dos Grids: dim3(256, 1, 1)               │
  │  • Laço Temporal Explícito (nt = 1250)                      │
  └───────────────┬─────────────────────────────┬───────────────┘
                  │                             │
    Cópia H2D     │ Chamadas <<<grid, block>>>  │ Cópia D2H (Snapshots)
    (Atribuição)  ▼                             ▼ (Atribuição)
  ┌─────────────────────────────────────────────────────────────┐
  │                   kernels.cuf (Device)                      │
  │  1. count_neighbors_gpu      (Fase 1 Topologia)             │
  │  2. fill_nodefam_gpu         (Fase 3 Topologia)             │
  │  3. apply_bc_gpu             (Kernel 3 - BCs Borda)         │
  │  4. compute_force_damage_gpu (Kernel 4 - Forças e Dano)     │
  │  5. integrate_kinematics_gpu (Kernel 5 - Cinemática Geral)  │
  └─────────────────────────────────────────────────────────────┘
```

Todas as subrotinas executadas na GPU são qualificadas com `attributes(global)` e recebem ponteiros com `attributes(device)`.

---

### 2.2 Gerenciamento Explícito de VRAM e Sintaxe de Cópia Fortran

Em vez de mapeamentos semânticos como `MAP(TO: ...)` do OpenMP, o CUDA Fortran declara variáveis fisicamente separadas para a memória da GPU:

```fortran
! Arrays no Host (CPU)
real(8) :: coord(totnode,2), disp(totnode,2), vel(totnode,2)

! Arrays no Device (GPU) - declarados com atributo 'device'
real(8), device, allocatable :: coord_d(:,:), d_disp(:,:), d_vel(:,:)
allocate(coord_d(totnode, 2))
allocate(d_disp(totnode, 2))
allocate(d_vel(totnode, 2))
```

#### Transferências Host $\longleftrightarrow$ Device Simplificadas:
No CUDA Fortran, o compilador sobrecarrega o operador de atribuição (`=`):
- **Host para Device:** `coord_d(1:totnode, 1:2) = coord(1:totnode, 1:2)`  
  *(Compilado internamente como uma chamada assíncrona/síncrona de `cudaMemcpyHostToDevice`)*.
- **Device para Host:** `disp(1:totnode, 1:2) = d_disp(1:totnode, 1:2)`  
  *(Compilado como `cudaMemcpyDeviceToHost`)*.

Isso elimina ambiguidades de ambiente de dados presentes no OpenMP Offloading, permitindo saber com exatidão onde cada variável está alocada.

---

### 2.3 Mapeamento Explícito de Threads e Proteção de Borda

Em OpenMP, a divisão de iterações entre threads é tratada pelo compilador. No CUDA Fortran, o programador controla o cálculo manual do índice global em 1-based index (padrão Fortran):

```fortran
! Identificação da thread global
i = (blockIdx%x - 1) * blockDim%x + threadIdx%x

! Cláusula de Guarda (Boundary Check)
if (i <= totnode) then
    ! Trabalho do ponto material i
end if
```

#### Dimensionamento Geométrico do Grid:
```fortran
block        = dim3(256, 1, 1)
grid_familia = dim3((totnode + 255) / 256, 1, 1)
```
- O bloco de **256 threads** é múltiplo exato do tamanho do *warp* da NVIDIA (32 threads), garantindo alta taxa de ocupação dos multiprocessadores (*Streaming Multiprocessors - SMs*).
- A divisão arredondada para cima `(totnode + 255) / 256` resulta em $\lceil 253.000 / 256 \rceil = 989$ blocos ($253.184$ threads totais). A condição de guarda `if (i <= totnode)` desativa as 184 threads excedentes no último bloco, evitando acessos ilegais à memória.

---

## 3. Análise Detalhada dos 5 Kernels CUDA

### 3.1 Kernel 1 e Kernel 2: Topologia Híbrida (Busca de Família)
* **Kernel 1 (`count_neighbors_gpu`)**: Cada thread $i$ varre os $N$ nós do domínio calculando $\sqrt{\Delta x^2 + \Delta y^2} \le \delta$. O uso da instrução `cycle` quando `i == j` substitui o aninhamento condicional, otimizando o fluxo de instruções no *warp*.
* **Intervenção Serial na CPU (Scan)**:
  ```fortran
  numfam(1:totnode, 1) = numfam_d(1:totnode) ! D2H
  pointfam(1,1) = 1
  do i = 2, totnode
      pointfam(i,1) = pointfam(i-1,1) + numfam(i-1,1)
  enddo
  pointfam_d(1:totnode) = pointfam(1:totnode, 1) ! H2D
  ```
  Assim como no OpenMP Offloading, o cálculo de soma prefixada é realizado no Host, exigindo transferências intermediárias.
* **Kernel 2 (`fill_nodefam_gpu`)**: Preenche `family_nodes` utilizando ponteiro de array dimensional indefinido `family_nodes(*)`, onde cada thread escreve a partir de `start_idx = family_pointer(i)`.

---

### 3.2 Kernel 3: Condições de Contorno com Grid Customizado (`apply_bc_gpu`)

Diferente do OpenMP Offloading (que executava dois laços independentes para borda inferior e superior), o CUDA Fortran implementa uma solução engenhosa:

```fortran
! Configuração no main.cuf
n_bordas = tottop - totint  ! Apenas 3.000 nós de borda!
grid     = dim3((n_bordas + block%x - 1) / block%x, 1, 1)

call aplica_bc_gpu<<<grid, block>>>(totint, totbottom, tottop, tt, dt, d_vel, d_disp)
```

No corpo do kernel (`kernels.cuf`):
```fortran
i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
idx = total_internal + i

if (idx <= total_top) then
    if (idx <= total_bottom) then
        velocity(idx, 2)     = -20.0_dp
        displacement(idx, 2) = -20.0_dp * tt * dt
    else
        velocity(idx, 2)     = 20.0_dp
        displacement(idx, 2) = 20.0_dp * tt * dt
    end if
end if
```

#### Vantagens desta Abordagem:
1. **Grid Sob Medida**: Em vez de lançar um grid sobre todos os $253.000$ nós para processar apenas $3.000$, o grid lança estritamente $\lceil 3.000 / 256 \rceil = 12$ blocos.
2. **Kernel Único**: Unifica as bordas inferior e superior em uma única chamada, reduzindo o overhead de despacho do driver.

---

### 3.3 Kernel 4: Forças Peridinâmicas e Dano (`compute_force_damage_gpu`)

Representa o núcleo computacional da peridinâmica.
* **Passagem de Parâmetros**: Escalares constantes são passados com o modificador `value` (ex: `totnode`, `delta`, `pi`, `bond_constant`, `vol`, `scr0`), sendo transferidos diretamente pela memória constante / registradores de parâmetros do kernel.
* **Alocação de Registradores**: Variáveis intermediárias (`dmg_par1`, `idist`, `nlength`, `fac`, `theta`, `scx`, `scy`, `scr`, `dforce1`, `dforce2`) são declaradas localmente dentro da subroutine, sendo alocadas pelo compilador nos registradores de alta velocidade do núcleo de hardware (*hardware registers*).

---

### 3.4 Kernel 5: Integrador Cinemático Unificado (`integrate_kinematics_gpu`)

Esta é uma das principais otimizações do código CUDA em relação ao OpenMP Offloading:

* **No OpenMP Offloading**: Eram necessários **3 laços separados** (nós internos, borda inferior e borda superior), resultando em 3 lançamentos de kernel.
* **No CUDA Fortran**: Foi unificado em **um único kernel geral**:

```fortran
i = (blockIdx%x - 1) * blockDim%x + threadIdx%x

if (i <= totnode) then
    acc1 = (p_force(i,1) + b_force(i,1)) / dens
    acc2 = (p_force(i,2) + b_force(i,2)) / dens

    if (i <= total_internal) then
        ! Nós Internos: Integração livre em X e Y
        velocity(i,1)  = velocity(i,1) + acc1 * dt
        velocity(i,2)  = velocity(i,2) + acc2 * dt
        displacement(i,1) = displacement(i,1) + velocity(i,1) * dt
        displacement(i,2) = displacement(i,2) + velocity(i,2) * dt
    else
        ! Nós de Contorno: Integra apenas X (Y é prescrito pelas BCs)
        velocity(i,1)  = velocity(i,1) + acc1 * dt
        displacement(i,1) = displacement(i,1) + velocity(i,1) * dt
    end if
end if
```

> **Redução Massiva de Lançamentos de Kernel:**  
> A cada passo temporal, o CUDA executa apenas **3 kernels** (`apply_bc_gpu`, `compute_force_damage_gpu` e `integrate_kinematics_gpu`), em comparação com os 6 lançamentos da versão OpenMP Offloading inicial, minimizando o impacto de latência de despacho da API.

---

## 4. Análise Crítica e Oportunidades de Melhoria

Para desenvolvedores que continuarem o projeto:

### 4.1 Alocação Estática de `nodefam` vs Alocação Dinâmica
* No `main.cuf.txt`, o vetor na GPU foi alocado estaticamente com tamanho arbitrário de 10 milhões de posições:
  ```fortran
  allocate(nodefam_d(10000000))
  ```
  Isso aloca $\approx 40$ MB fixos na GPU e transfere 10 milhões de inteiros na linha 298.
* **Oportunidade**: Calcular `total_family_size` no Host logo após o Scan da Parte B e alocar `nodefam_d(total_family_size)` com o tamanho exato, alinhando-se com a estratégia de economia de VRAM desenvolvida na versão de offloading.

### 4.2 Divergência de Branches no Kernel 5
* O teste `if (i <= total_internal)` no Kernel 5 faz com que o último *warp* da região interna execute os dois ramos do condicional (*warp divergence*).
* Como apenas 1 *warp* em 989 sofre com isso na fronteira entre interno e borda, o impacto prático é desprezível ($< 0{,}1\%$), confirmando que a fusão dos kernels de cinemática foi vantajosa.


## 5. Guia de Compilação com NVIDIA HPC SDK

O código CUDA Fortran requer o compilador **`nvfortran`** com suporte nativo a CUDA:

```bash
# Compilação otimizada com arquitetura de GPU específica (ex: Ampere CC 8.0)
nvfortran -O3 -cuda -gpu=cc80,cuda12.0 -mp -Minfo=accel kernels.cuf main.cuf -o perid_cuda.exe
```

### Flags Explicadas:
- `-cuda`: Habilita a extensão CUDA Fortran e reconhecimento de sufixos `.cuf`.
- `-gpu=cc80`: Gera código SASS otimizado para a microarquitetura do chip (ex: A100/RTX 30xx).
- `-mp`: Habilita suporte ao OpenMP utilizado para a telemetria (`omp_get_wtime`).
- `-Minfo=accel`: Emite relatório detalhado da alocação de registradores e geração de código para os kernels.
