module kernels_peridinamica
    use cudafor
    implicit none

contains

    ! --- KERNEL 1: Conta quantos vizinhos cada nó tem ---
    attributes(global) subroutine conta_vizinhos_gpu(n, coord, delta, numfam)
        integer, value :: n
        real(8), value :: delta
        real(8), device :: coord(n, 2)
        integer, device :: numfam(n)
        
        integer :: i, j
        real(8) :: dist, dx, dy

        ! Identifica qual thread é responsável por qual nó 'i'
        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x

        if (i <= n) then
            numfam(i) = 0
            do j = 1, n
                if (i == j) cycle
                
                dx = coord(j,1) - coord(i,1)
                dy = coord(j,2) - coord(i,2)
                dist = sqrt(dx*dx + dy*dy)

                if (dist <= delta) then
                    numfam(i) = numfam(i) + 1
                end if
            end do
        end if
    end subroutine conta_vizinhos_gpu

    ! --- KERNEL 2: Preenche a lista nodefam ---
    attributes(global) subroutine preenche_nodefam_gpu(n, coord, delta, pointfam, nodefam)
        integer, value :: n
        real(8), value :: delta
        real(8), device :: coord(n, 2)
        integer, device :: pointfam(n)
        integer, device :: nodefam(*) ! Tamanho desconhecido aqui
        
        integer :: i, j, k, inicio
        real(8) :: dist, dx, dy

        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x

        if (i <= n) then
            k = 0 ! Contador local
            inicio = pointfam(i) ! Onde começa a lista deste nó
            
            do j = 1, n
                if (i == j) cycle

                dx = coord(j,1) - coord(i,1)
                dy = coord(j,2) - coord(i,2)
                dist = sqrt(dx*dx + dy*dy)

                if (dist <= delta) then
                    ! Escreve na memória global da GPU
                    nodefam(inicio + k) = j
                    k = k + 1
                end if
            end do
        end if
    end subroutine preenche_nodefam_gpu

end module kernels_peridinamica
    
    
    
    
    
program node_fortran_cuda

use cudafor 
use kerlnels_peridinamica 
implicit none

integer ndivx, ndivy, totnode, nt, maxfam, nnum, cnode, i, j, tt, nbnd, totint, totbottom, tottop
!ndivx: Number of divisions in x direction - except boundary region
parameter(ndivx = 500)
!ndivy: Number of divisions in y direction - except boundary region
parameter(ndivy = 500)
!nbnd: Number of divisions in the boundary region
parameter(nbnd = 3)
!totnode: Total number of material points
parameter (totnode = ndivx*(ndivy + 2 * nbnd)) 
!nt: Total number of time step
parameter(nt = 1250)
!maxfam: Maximum number of material points inside a horizon of a material point
parameter(maxfam = 100)

real *8 length,width, dx, delta, thick, dens, emod, pratio, area, vol, bc  
real *8 sedload1, sedload2, dt, totime, ctime, idist, fac, radij, nlength, dforce1, dforce2 
real *8 crlength, scr0, pi, tmpdx, tmpvol, tmpcx, tmpcy, tmpux, tmpuy, dmgpar1, dmgpar2, theta 
real *8 scx, scy, scr

real *8 coord(totnode,2), pforce(totnode,2), pforceold(totnode,2), bforce(totnode,2), stendens(totnode,2)
real *8 fncst(totnode,2), disp(totnode,2), vel(totnode,2), velhalfold(totnode,2), velhalf(totnode,2)
real *8 acc(totnode,2), massvec(totnode,2), enddisp(nt,1), endtime(nt,1), dmg(totnode,1)
integer numfam(totnode,1), pointfam(totnode,1), nodefam(10000000,1), fail(totnode,maxfam)

! -- variaveis da GPU 
real(8), device, allocatable :: coord_d(:,:)
integer, device, allocatable :: numfam_d(:)
integer, device, allocatable :: pointfam_d(:)
integer, device, allocatable :: nodefam_d(:

type(dim3) :: grid, block
integer :: threadsPorBloco, numBlocos
integer :: err_cuda

! --- VARIÁVEIS PARA TIMER ---
integer :: t_inicio_sim, t_final_sim, taxa_clock, kount
real *8 :: tempo_total_sim_s
integer :: total_vizinhos_usados

pi = dacos(-1.0d0)

do i = 1, totnode 
    !coord: Material point locations, 1:x-coord, 2:y-coord
	coord(i,1) = 0.0d0
	coord(i,2) = 0.0d0
    !numfam: Number of family members of each material point
	numfam(i,1) = 0
    !pointfam: index array to find the family members in nodefam array
	pointfam(i,1) = 0
    !pforce: total peridynamic force acting on a material point, 1:x-coord, 2:y-coord
	pforce(i,1) = 0.0d0
	pforce(i,2) = 0.0d0
    !pforceold: total peridynamic force acting on a material point in the previous time step
    !1:x-coord, 2:y-coord
	pforceold(i,1) = 0.0d0
	pforceold(i,2) = 0.0d0
    !bforce: body load acting on a material point, 1:x-coord, 2:y-coord
	bforce(i,1) = 0.0d0
	bforce(i,2) = 0.0d0
    !stendens: strain energy of a material point, 1:loading 1, 2:loading 2
	stendens(i,1) = 0.0d0
	stendens(i,2) = 0.0d0
    !fncst: surface correction factors of a material point, 1:loading 1, 2:loading 2
	fncst(i,1) = 1.0d0 
	fncst(i,2) = 1.0d0  
    !disp: displacement of a material point, 1:x-coord, 2:y-coord
	disp(i,1) = 0.0d0
	disp(i,2) = 0.0d0
    !vel: velocity of a material point, 1:x-coord, 2:y-coord
	vel(i,1) = 0.0d0
	vel(i,2) = 0.0d0
	velhalfold(i,1) = 0.0d0
	velhalfold(i,2) = 0.0d0
	velhalf(i,1) = 0.0d0
	velhalf(i,2) = 0.0d0
    !acc: acceleration of a material point, 1:x-coord, 2:y-coord 
	acc(i,1) = 0.0d0
	acc(i,2) = 0.0d0
    !massvec: massvector for adaptive dynamic relaxation, 1:x-coord, 2:y-coord
	massvec(i,1) = 0.0d0
	massvec(i,2) = 0.0d0
    !fail: Failure array
	do j = 1, maxfam
		fail(i,j) = 0
    enddo
    !dmg: Damage of a material point
	dmg(i,1) = 0.0d0
enddo

do i = 1, 1000000
    !nodefam: array containing family members of all material points
	nodefam(i,1) = 0
enddo

!length: Total length of the plate
length = 0.05d0
!width: Total width of the plate
width = 0.05d0
!dx: Spacing between material points
dx = length / ndivx
!delta: Horizon
delta = 3.015 * dx
!thick: Thickness of the plate
thick = dx
!dens: Density
dens = 8000.0d0
!emod: Elastic modulus
emod = 192.0d9
!pratio12 = Poisson's ratio
pratio = 1.0d0 / 3.0d0
!area: Cross-sectional area
area = dx * dx
!vol: Volume of a material point
vol = area * dx
!bc: Bond constant 
bc = 9.0d0 * emod / (pi * thick * (delta**3))
!sedload1: Strain energy density for the first loading
sedload1 = 9.0d0 / 16.0d0 * emod * 1.0d-6   
!sedload2: Strain energy density for the second loading
sedload2 = 9.0d0 / 16.0d0 * emod * 1.0d-6
!dt: Time interval
dt = 0.8d0 * dsqrt(2.0d0*dens*dx/(pi*delta**2*dx*bc))
!totime: Total time
totime = nt * dt
!ctime: Current time
ctime = 0.0d0
!idist: Initial distance
idist = 0.0d0
do i = 1, nt
	enddisp(i,1) = 0.0d0
	endtime(i,1) = 0.0d0
enddo
!fac: Volume correction factor
fac = 0.0d0
!radij: Material point radius
radij = dx / 2.0d0
!nnum: Material point number
nnum = 0
!cnode: Current material point
cnode = 0
!Length of deformed bond
nlength  = 0.0d0
!dforce1: x component of the PD force between two material points
dforce1 = 0.0d0
!dforce1: y component of the PD force between two material points
dforce2 = 0.0d0
!crlength: Crack length
crlength = 0.01d0
!scr0: Critical stretch
scr0 = 0.04472d0

!Initialization of fail flag array
!1 means no failure, 0 means failure of the PD bond
do i = 1,totnode
	do j = 1,maxfam
		fail(i,j) = 1
    enddo
enddo

!Specification of the locations of material points
!Material points of the internal region
do i = 1,ndivy
    do j = 1,ndivx
        nnum = nnum + 1
        coord(nnum,1) = (-1.0d0 * length / 2.0d0) + (dx / 2.0d0) + (j-1) * dx
        coord(nnum,2) = (-1.0d0 * width / 2.0d0) + (dx / 2.0d0) + (i-1) * dx
    enddo
enddo

totint = nnum

!Material points of the boundary region - bottom
do i = 1,nbnd
    do j = 1,ndivx
        nnum = nnum + 1
        coord(nnum,1) = -1.0d0 /2.0d0 * length + (dx / 2.0d0) + (j - 1) * dx
        coord(nnum,2) = -1.0d0 /2.0d0 * width - (dx / 2.0d0) - (i - 1) * dx
    enddo
enddo

totbottom = nnum

!Material points of the boundary region - top
do i = 1,nbnd
    do j = 1,ndivx
        nnum = nnum + 1
        coord(nnum,1) = -1.0d0 /2.0d0 * length + (dx / 2.0d0) + (j - 1) * dx
        coord(nnum,2) = 1.0d0 /2.0d0 * width + (dx / 2.0d0) + (i - 1) * dx
    enddo
enddo

tottop = nnum

! preparação da GPU 

allocate(coord_d(totnode, 2))
allocate(numfam_d(totnode))
allocate(pointfam_d(totnode))
allocate(nodefam_d(10000000)) 

! Copia as coordenadas calculadas na CPU para a GPU
coord_d = coord(1:totnode, :)
    
! Configura a execução: blocos de 256 threads
threadsPorBloco = 256
numBlocos = (totnode + threadsPorBloco - 1) / threadsPorBloco
block = dim3(threadsPorBloco, 1, 1)
grid = dim3(numBlocos, 1, 1)

!T_Inicio = Ler o tempo aqui
call SYSTEM_CLOCK(count=t_inicio_sim, count_rate=taxa_clock) 

!Determination of material points inside the horizon of each material poinT

! gpu CONTA quantps vizinhos cada um tem 
call conta_vizinhos_gpu<<<grid, block>>>(totnode, coord_d, delta, numfam_d)

! Espera a GPU terminar e traz o resultado 'numfam' de volta pra CPU
err_cuda = cudaDeviceSynchronize()
numfam(1:totnode, 1) = numfam_d(1:totnode)

! --- PARTE 2: CPU faz a soma acumulada (Pointfam) ---
    ! Isso precisa ser serial, mas é muito rápido
if (totnode > 0) pointfam(1,1) = 1
do i = 2, totnode
    pointfam(i,1) = pointfam(i-1,1) + numfam(i-1,1)
 enddo  

! Manda o 'pointfam' calculado de volta pra GPU para o próximo passo
pointfam_d(1:totnode) = pointfam(1:totnode, 1)

! --- PARTE 3: GPU preenche a lista de vizinhos (Nodefam) ---
call preenche_nodefam_gpu<<<grid, block>>>(totnode, coord_d, delta, pointfam_d, nodefam_d)
    
! Espera terminar e traz a lista gigante de volta pra CPU
err_cuda = cudaDeviceSynchronize()
nodefam(1:10000000, 1) = nodefam_d(1:10000000)

!T_Final = Ler o tempo aqui
call system_clock(count=t_final_sim)

        
!T_Total = T_Final - T_Inicio
tempo_total_sim_s = real(t_final_sim - t_inicio_sim) / real(taxa_clock)


! Limpeza da memória da GPU 
deallocate(coord_d, numfam_d, pointfam_d, nodefam_d)


! 1) Escrever um arquivo de saida com os seguintes dados:
!    Indice dos pontos de cada familia (pointfam)
!    Pontos que compoem a familia de cada ponto (nodefam)
!    Tempo que a simulação levou 

! 2) Estudar o OpenMP para CPU e para GPU-offloading
!    Gerar uma versão do codigo para CPU e outra para GPU-offloading

! 3) Estudar o CUDAFortran 
!    Gerar uma versão do código para comparação. 

!Debug no teu PC, resultados no LabSin.

!Fazer 4 comparações: 1) Seriado, 2) OpenMP(CPU), 3) OpenMP(GPU-offloading), 4)FortranCUDA  
!Comparações 1, 2, 3 -> Compilador IntelFortran(ifx) 4 -> Compilador nvfortran
!IDE Visual Studio Community 
!Intel Fortran: https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html?packages=hpc-toolkit&hpc-toolkit-os=windows&hpc-toolkit-win=offline#collapseCollapsible1761228516178
! RESULTADOS: 1) Comparar se as familias nos 4 testes são iguais e qual o tempo que cada teste levou. 

! - ESCREVER ARQUIVOS NO Fortran

!open(26,file = 'coord_disp_pd_750_pwc_v20.txt')
!do i = 1, totint
!	write(26,111) coord(i,1), coord(i,2), disp(i,1), disp(i,2), dmg(i,1)
!enddo
!111 format(e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5)
!222 format(i3,3x,i3,3x,....)
! 1234567.89012   1234567.89012   1234567.89012   1234567.89012 ...
!close(26)

!------------------------------- RODAR ATÉ AQUI! --------------------------------------------!

! Calcula o limite real de dados no nodefam para não imprimir lixo/zeros
! O último índice usado é onde começa a família do último nó + quantos vizinhos ele tem - 1
print *, "Iniciando a escrita do arquivo 'familia_resultados.txt'..."

!bre o arquivo (Unit 26)
open(unit=26, file='familia_resultados_cuda.txt', status='replace', action='write')
! --- 1. CABEÇALHO E TEMPO ---
write(26, *) "===================================================="
write(26, *) "RELATORIO DE EXECUCAO - PERIDINAMICA"
write(26, *) "===================================================="
write(26, '(A, F15.6, A)') "Tempo Total de Simulacao: ", tempo_total_sim_s, " segundos"
write(26, *) ""
! --- 2. DADOS DE POINTFAM (Índices) ---
write(26, *) "===================================================="
write(26, *) "TABELA DE INDICES (POINTFAM)"
write(26, *) "Significado: (No Global) -> (Onde comeca a lista em Nodefam)"
write(26, *) "===================================================="

do i = 1, totnode
    ! Formato: Inteiro(10 digitos) | Texto | Inteiro(10 digitos)
    write(26, '(I10, A, I10)') i, " -> ", pointfam(i,1)
enddo
write(26, *) ""
! --- 3. DADOS DE NODEFAM (Listagem detalhada das Famílias) ---
write(26, *) "===================================================="
write(26, *) "CONTEUDO DAS FAMILIAS (NODEFAM)"
write(26, *) "Lista de vizinhos para cada ponto material"
write(26, *) "===================================================="
do i = 1, totnode
    ! Cabeçalho de cada família: "Familia do Ponto 1 (Total: 50)"
    write(26, '(A, I8, A, I6, A)') "Familia do Ponto ", i, " (Total: ", numfam(i,1), ")"
    
    if (numfam(i,1) > 0) then
        ! Loop interno para escrever os vizinhos na mesma linha
        do j = 1, numfam(i,1)
            ! O calculo do índice é: Onde começa a familia + o contador atual - 1
            ! advance='no' impede que pule de linha a cada número
            write(26, '(I8, 1x)', advance='no') nodefam(pointfam(i,1) + j - 1, 1)
        enddo
        ! Pula uma linha depois de escrever todos os vizinhos deste ponto
        write(26, *) 
    else
        write(26, '(A)') "      (Sem vizinhos / Familia vazia)"
    endif
    
    ! Linha em branco extra para separar visualmente os blocos (opcional)
    write(26, *) 
enddo
close(26)
print *, "Escrita do arquivo concluida com sucesso."


end program node_fortran_cuda
