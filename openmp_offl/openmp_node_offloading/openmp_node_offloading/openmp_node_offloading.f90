program openmp_node_offloading
use omp_lib
implicit none

! !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO ... para fazer offloading de um loop

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
real *8 acc(totnode,2), massvec(totnode,2)

! - viram vetores 1D ---
real *8 enddisp(nt), endtime(nt), dmg(totnode)

!     10.000.000 posicoes e passa a ser "allocatable", com o tamanho exato    ---
!     calculado na CPU (total_family_size) antes de ir para a GPU             ---
integer numfam(totnode), pointfam(totnode), fail(totnode,maxfam)
integer, allocatable :: nodefam(:)
integer :: total_family_size

! declarações para registro de tempo 
integer :: kount
real *8 :: start_time, end_time,t_start, t_end
real *8 :: t_parte1, t_parte2, t_parte3, tempo_total_sim_s1

!declarações para registro de tempo (boundary conditions)
real*8  :: t_inicio_bc, t_final_bc, tempo_total_bc_gpu_s

!declarações para registro de tempo (calculo de forcas/dano)
real*8  :: t_inicio_forca, t_final_forca, tempo_total_forca_dano_gpu_s

!declarações para registro de tempo (escrita dos arquivos de saida - I/O)
real*8  :: t_inicio_escrita, t_final_escrita, tempo_total_escrita_gpu_s

!declaração para registro do tempo total do programa (wall-clock)
real*8  :: t_prog_inicio, t_prog_fim, tempo_total_programa_gpu_s

pi = dacos(-1.0d0)

t_prog_inicio = omp_get_wtime() ! marca o inicio do wall-clock total do programa

do i = 1, totnode 
    !coord: Material point locations, 1:x-coord, 2:y-coord
	coord(i,1) = 0.0d0
	coord(i,2) = 0.0d0
    !numfam: Number of family members of each material point
	numfam(i) = 0
    !pointfam: index array to find the family members in nodefam array
	pointfam(i) = 0
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
	dmg(i) = 0.0d0
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
	enddisp(i) = 0.0d0
	endtime(i) = 0.0d0
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

! /// trabalho começa aqui ///

t_start = omp_get_wtime()

!Determination of material points inside the horizon of each material point

!    - MAP(TO:...) copia dados da CPU para a GPU (só ida)
!    - MAP(ALLOC:...) aloca espaço na GPU, mas não copia nada
!    - MAP(FROM:...) copia dados da GPU para a CPU (só volta)

! --- usamos TARGET ENTER DATA (nao-estruturado)
!     em vez de TARGET DATA aqui, porque o tamanho de nodefam só é conhecido
!     DEPOIS de contar numfam na GPU. Assim coord/delta ficam residentes na
!     GPU e, mais adiante, "anexamos" pointfam/nodefam sem ter que fechar
!     esta região e reenviar coord de novo.
!$OMP TARGET ENTER DATA MAP(TO: coord(1:totnode,:), delta) &
!$OMP& MAP(ALLOC: numfam(1:totnode))

! PROCESSAMENTO NA GPU
t_start = omp_get_wtime() !maracador de tempo no incio da parte paralela.

   !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(j, idist)
    do i = 1, totnode
        numfam(i) = 0
        do j = 1, totnode
            if (i .ne. j) then
                idist = dsqrt((coord(j,1) - coord(i,1))**2 + (coord(j,2) - coord(i,2))**2)
                if (idist <= delta) then
                    numfam(i) = numfam(i) + 1
                endif
            endif
        enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
    
    ! sincronizando GPU -> CPU para parte serial 
    !'UPDATE FROM' copia o array da gpu para a cpu .
    
    ! 1.2 Trazer resultado para CPU (Gargalo de transferência)
    
    !$OMP TARGET UPDATE FROM(numfam(1:totnode))
    
    t_end = omp_get_wtime() !marcador de tempo do fim da parte paralela
    t_parte1 = t_end - t_start 
    
    ! PARTE 2 SERIAL 
    
    t_start = omp_get_wtime() ! marcador de tempo no inicio da parte serial.
    
    pointfam(1) = 1
    do i = 2, totnode
        pointfam(i) = pointfam(i-1) + numfam(i-1)
    enddo

    !     tamanho exato calculado na CPU e
    !     alocação de nodefam SÓ AGORA, com o tamanho definitivo
    total_family_size = pointfam(totnode) + numfam(totnode) - 1
    allocate(nodefam(total_family_size))
    
    t_end = omp_get_wtime()
    t_parte2 = t_end - t_start
    
    ! agora vamos fazer o caminho CPU -> GPU 
    !'UPDATE TO' copia o array da CPU para a GPU.
    
    t_start = omp_get_wtime() !marcador de tempo no inicio da segunda parte paralela 
    
    ! pointfam e nodefam só entram no ambiente de dados da GPU agora, já com
    ! o tamanho exato de nodefam (evitando superalocar 10.000.000 posições)
    
    !$OMP TARGET ENTER DATA MAP(TO: pointfam(1:totnode)) &
    !$OMP& MAP(ALLOC: nodefam(1:total_family_size))
    
    ! Preencher Nodefam
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(j, idist, kount)
    do i = 1, totnode
        kount = 0
        do j = 1, totnode
            if (i .ne. j) then
                idist = dsqrt((coord(j,1) - coord(i,1))**2 + (coord(j,2) - coord(i,2))**2)
                if (idist <= delta) then
                    kount = kount + 1
                    nodefam(pointfam(i) + kount - 1) = j
                endif
            endif
        enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
    
    ! nodefam precisa voltar para a CPU porque os laços seriais de
    ! pré-processamento (fenda e Loading 1 e 2), logo abaixo, ainda rodam na CPU
    !$OMP TARGET UPDATE FROM(nodefam(1:total_family_size))
    
t_end = omp_get_wtime() !marcador de tempo do fim da segunda parte paralela 
t_parte3 = t_end - t_start


! Somando tudo
tempo_total_sim_s1 = t_parte1 + t_parte2 + t_parte3

! coord, delta, numfam, pointfam e nodefam permanecem residentes na GPU
! (entraram via TARGET ENTER DATA e não foram liberados) para serem
! reaproveitados sem nova cópia dentro do loop temporal, mais abaixo.

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

!Definition of the crack surface
!PD bonds penetrating through the crack surface are broken
do i = 1,totnode
    do j = 1,numfam(i)
        cnode = nodefam(pointfam(i)+j-1)
        if ((coord(cnode,2) > 0.0d0).and.(coord(i,2) < 0.0d0)) then
            if ((dabs(coord(i,1)) - (crlength / 2.0d0)).le.1.0d-10) then
				fail(i,j) = 0
            elseif ((dabs(coord(cnode,1)) - (crlength / 2.0d0)).le.1.0d-10) then
				fail(i,j) = 0
            endif
        elseif ((coord(i,2) > 0.0d0).and.(coord(cnode,2) < 0.0d0)) then
            if((dabs(coord(i,1)) - (crlength / 2.0d0)).le.1.0d-10) then 
				fail(i,j) = 0
			elseif((dabs(coord(cnode,1)) - (crlength / 2.0d0)).le.1.0e-10) then
				fail(i,j) = 0
            endif
        endif        
    enddo
enddo

!Loading 1
do i = 1,totnode
    disp(i,1) = 0.001d0 * coord(i,1)
    disp(i,2) = 0.0d0
enddo

do i = 1,totnode
    stendens(i,1) = 0.0d0
    do j = 1,numfam(i)
        cnode = nodefam(pointfam(i)+j-1)
        idist = dsqrt((coord(cnode,1) - coord(i,1))**2 + (coord(cnode,2) - coord(i,2))**2)
        nlength = dsqrt((coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1))**2 + (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2))**2)
        if (idist.le.delta-radij) then
            fac = 1.0d0
        elseif (idist.le.delta+radij) then
            fac = (delta+radij-idist)/(2.0d0*radij)
        else
            fac = 0.0d0
        endif
                       
        stendens(i,1) = stendens(i,1) + 0.5d0 * 0.5d0 * bc * ((nlength - idist) / idist)**2 * idist * vol * fac  
    enddo
    !Calculation of surface correction factor in x direction 
    !by finding the ratio of the analytical strain energy density value
    !to the strain energy density value obtained from PD Theory
    fncst(i,1) = sedload1 / stendens(i,1)
enddo
    
!Loading 2
do i = 1,totnode
    disp(i,1) = 0.0d0
    disp(i,2) = 0.001d0 * coord(i,2)
enddo

do i = 1,totnode
    stendens(i,2) = 0.0d0
    do j = 1,numfam(i)
        cnode = nodefam(pointfam(i)+j-1)
        idist = dsqrt((coord(cnode,1) - coord(i,1))**2 + (coord(cnode,2) - coord(i,2))**2)
        nlength = dsqrt((coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1))**2 + (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2))**2)
        if (idist.le.delta-radij) then
            fac = 1.0d0
        elseif (idist.le.delta+radij) then
            fac = (delta+radij-idist)/(2.0d0*radij)
        else
            fac = 0.0d0
        endif   
                      
        stendens(i,2) = stendens(i,2) + 0.5d0 * 0.5d0 * bc * ((nlength - idist) / idist)**2 * idist * vol * fac 
    enddo
    !Calculation of surface correction factor in y direction 
    !by finding the ratio of the analytical strain energy density value
    !to the strain energy density value obtained from PD Theory
    fncst(i,2) = sedload2 / stendens(i,2)
enddo

!Initialization of displacements and velocities
do i = 1,totnode
    vel(i,1) = 0.0d0
    vel(i,2) = 0.0d0
    disp(i,1) = 0.0d0
    disp(i,2) = 0.0d0         
enddo

! Inicializa o acumulador de tempo das BCs
tempo_total_bc_gpu_s = 0.0d0

! Inicializa o acumulador de tempo do calculo de forcas/dano
tempo_total_forca_dano_gpu_s = 0.0d0

! Inicializa o acumulador de tempo de escrita dos arquivos de saida (I/O)
tempo_total_escrita_gpu_s = 0.0d0

! --------------------------------------------
! coord, delta, numfam, pointfam e nodefam já estão residentes na GPU desde a
! fase acima (TARGET ENTER DATA) e não precisam ser remapeados aqui. Esta
! regiao cobre tudo que o loop temporal le/escreve a cada passo de tempo, em
! UMA ÚNICA TARGET DATA que envolve o "do tt = 1, nt" inteiro - sem abrir e
! fechar regiões de dados a cada iteração.
! disp e vel entram só com MAP(TO) (valor inicial = 0) porque a devolução
! para a CPU é feita sob demanda, via TARGET UPDATE FROM, só nos passos de
! tempo em que o resultado precisa ser escrito em arquivo (750/1000/1250).
! pforce, acc e dmg entram com MAP(ALLOC) porque são reescritos por completo
! a cada iteração, então copiar o valor inicial da CPU seria desperdício.
! Escalares como bc, vol, radij, scr0, dens, dt e length são tratados
! automaticamente como firstprivate por cada TARGET, sem necessidade de MAP.


!$OMP TARGET DATA &
!$OMP& MAP(TO: fail(1:totnode,1:maxfam), fncst(1:totnode,:), bforce(1:totnode,:), &
!$OMP&         disp(1:totnode,:), vel(1:totnode,:)) &
!$OMP& MAP(ALLOC: pforce(1:totnode,:), acc(1:totnode,:), dmg(1:totnode))

!Time integration
do tt = 1,nt
    write(*,*) 'tt = ', tt
	ctime = tt * dt

!______________________________DAQUI____________________________________________

    !Application of boundary conditions at the top and bottom edges

    ! captura de tempo inicial 
    t_inicio_bc = omp_get_wtime()
    
    ! --- item 3: nada de "target data" por iteração aqui - vel/disp já estão
    !     residentes na GPU desde a região ampla aberta acima do loop tt
      !$omp target teams distribute parallel do
      do i = (totint+1), totbottom
          vel(i,2) = -20.0d0
          disp(i,2) = -20.0d0 * tt * dt
      enddo
      !$omp end target teams distribute parallel do
    
      !$omp target teams distribute parallel do
      do i = (totbottom+1), tottop
          vel(i,2) = 20.0d0
          disp(i,2) = 20.0d0 * tt * dt
      enddo   
      !$omp end target teams distribute parallel do
    
    ! Agora capturamos o tempo final
    t_final_bc = omp_get_wtime()
    tempo_total_bc_gpu_s = tempo_total_bc_gpu_s + (t_final_bc - t_inicio_bc)
      
    
!______________________________ATÉ AQUI ____________________________________________ 

!______________________________________AQUI__________________________________________
    ! ---  -------------------------------
    ! Diretiva imediatamente acima do cálculo de pforce. Todas as variáveis
    ! de uso local a cada iteração de i (thread) vão em PRIVATE, para que
    ! cada ponto material tenha sua própria cópia e não haja corrida de dados
    ! entre pontos processados em paralelo.
    ! captura de tempo inicial do calculo de forcas/dano
    t_inicio_forca = omp_get_wtime()

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
    !$OMP& PRIVATE(j, cnode, idist, nlength, fac, theta, scx, scy, scr, &
    !$OMP&         dforce1, dforce2, dmgpar1, dmgpar2)
    do i = 1,totnode
        dmgpar1 = 0.0d0
        dmgpar2 = 0.0d0
        pforce(i,1) = 0.0d0
        pforce(i,2) = 0.0d0
        
        do j = 1,numfam(i)
                cnode = nodefam(pointfam(i)+j-1)
                idist = dsqrt((coord(cnode,1) - coord(i,1))**2 + (coord(cnode,2) - coord(i,2))**2)
                nlength = dsqrt((coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1))**2 + (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2))**2)
                
                !Volume correction
                if (idist.le.delta-radij) then
                    fac = 1.0d0
                elseif (idist.le.delta+radij) then
                    fac = (delta+radij-idist)/(2.0d0*radij)
                else
                    fac = 0.0d0
                endif
                
                if (dabs(coord(cnode,2) - coord(i,2)).le.1.0d-10) then 
                    theta = 0.0d0
                elseif (dabs(coord(cnode,1) - coord(i,1)).le.1.0d-10) then
                    theta = 90.0d0 * pi / 180.0d0
                else
                    theta = datan(dabs(coord(cnode,2) - coord(i,2)) / dabs(coord(cnode,1) - coord(i,1)))
                endif
                
                !Determination of the surface correction between two material points
                scx = (fncst(i,1) + fncst(cnode,1)) / 2.0d0
                scy = (fncst(i,2) + fncst(cnode,2)) / 2.0d0
                scr = 1.0d0 / (((dcos(theta))**2 / (scx)**2) + ((dsin(theta))**2 / (scy)**2))
                scr = dsqrt(scr)
                
                if (fail(i,j).eq.1) then
                    !Calculation of the peridynamic force in x and y directions 
                    !acting on a material point i due to a material point j
                    dforce1 = bc * (nlength - idist) / idist * vol * scr * fac * (coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1)) / nlength             
                    dforce2 = bc * (nlength - idist) / idist * vol * scr * fac * (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2)) / nlength             
                else
                    dforce1 = 0.0d0
                    dforce2 = 0.0d0
                endif 
                pforce(i,1) = pforce(i,1) + dforce1             
                pforce(i,2) = pforce(i,2) + dforce2             
                
                !Definition of a no-fail zone             
                if (dabs((nlength - idist) / idist) > scr0) then
                    if (dabs(coord(i,2)).le.(length/4.0d0)) then
                        fail(i,j) = 0 
                    endif
                endif                      
                            
                dmgpar1 = dmgpar1 + fail(i,j) * vol * fac
                dmgpar2 = dmgpar2 + vol * fac             
        enddo
        !Calculation of the damage parameter
        dmg(i) = 1.0d0 - dmgpar1 / dmgpar2
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    ! Agora capturamos o tempo final do calculo de forcas/dano
    t_final_forca = omp_get_wtime()
    tempo_total_forca_dano_gpu_s = tempo_total_forca_dano_gpu_s + (t_final_forca - t_inicio_forca)
    ! _________________________________________________ATÉ AQUI_____________________________________________________________ 
    
    ! A integração (acc/vel/disp) continua na GPU, dentro da mesma região de
    ! dados ampla: se ficasse na CPU, seria preciso trazer pforce de volta e
    ! reenviar vel/disp a cada passo de tempo, anulando a otimização das familias 
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
    do i = 1,totint
        !Calculation of acceleration of material point i
        acc(i,1) = (pforce(i,1) + bforce(i,1)) / dens
        acc(i,2) = (pforce(i,2) + bforce(i,2)) / dens
        !Calculation of velocity of material point i
        !by integrating the acceleration of material point i
        vel(i,1) = vel(i,1) + acc(i,1) * dt
        vel(i,2) = vel(i,2) + acc(i,2) * dt
        !Calculation of displacement of material point i
        !by integrating the velocity of material point i
        disp(i,1) = disp(i,1) + vel(i,1) * dt
        disp(i,2) = disp(i,2) + vel(i,2) * dt
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
    
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
    do i = (totint+1), totbottom
        acc(i,1) = (pforce(i,1) + bforce(i,1)) / dens
        vel(i,1) = vel(i,1) + acc(i,1) * dt
        disp(i,1) = disp(i,1) + vel(i,1) * dt
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
    do i = (totbottom+1), tottop
        acc(i,1) = (pforce(i,1) + bforce(i,1)) / dens
        vel(i,1) = vel(i,1) + acc(i,1) * dt
        disp(i,1) = disp(i,1) + vel(i,1) * dt
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
               
    endtime(tt) = ctime

    ! captura de tempo inicial da escrita (I/O) - so ha custo real quando
    ! tt bate com um dos passos abaixo; nos demais passos o bloco e vazio
    t_inicio_escrita = omp_get_wtime()

    if (tt.eq.750) then
        ! Só busca da GPU exatamente o que vai ser escrito em arquivo
        !$OMP TARGET UPDATE FROM(disp(1:totint,:), dmg(1:totint))
        open(26,file = 'coord_disp_pd_750_pwc_v20.txt')
        do i = 1, totint
            write(26,111) coord(i,1), coord(i,2), disp(i,1), disp(i,2), dmg(i)
        enddo
        close(26)
    elseif (tt.eq.1000) then
        !$OMP TARGET UPDATE FROM(disp(1:totint,:), dmg(1:totint))
        open(26,file = 'coord_disp_pd_1000_pwc_v20.txt')
        do i = 1, totint
            write(26,111) coord(i,1), coord(i,2), disp(i,1), disp(i,2), dmg(i)
        enddo
        close(26)
    elseif (tt.eq.1250) then
        !$OMP TARGET UPDATE FROM(disp(1:totint,:), dmg(1:totint))
        open(26,file = 'coord_disp_pd_1250_pwc_v20.txt')
        do i = 1, totint
            write(26,111) coord(i,1), coord(i,2), disp(i,1), disp(i,2), dmg(i)
        enddo
        close(26)
    endif

    ! captura de tempo final da escrita (I/O) e acumula
    t_final_escrita = omp_get_wtime()
    tempo_total_escrita_gpu_s = tempo_total_escrita_gpu_s + (t_final_escrita - t_inicio_escrita)

enddo

!$OMP END TARGET DATA

! marca o fim do wall-clock total do programa (equivalente ao "Tempo Total
! do Programa" da versao CPU) - captura logo apos o fim da simulacao,
! antes da escrita do arquivo de depuracao familia_resultados_offloading.txt
t_prog_fim = omp_get_wtime()
tempo_total_programa_gpu_s = t_prog_fim - t_prog_inicio

111 format(e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5)

! Libera o que ficou residente na GPU desde a fase de construção das famílias
!$OMP TARGET EXIT DATA MAP(RELEASE: coord(1:totnode,:), delta, numfam(1:totnode), &
!$OMP& pointfam(1:totnode), nodefam(1:total_family_size))

    ! ====================================================================================================
    ! ESCRITA DO ARQUIVO DE RESULTADOS GERAIS
    ! ====================================================================================================
    print *, "Iniciando a escrita do arquivo 'familia_resultados_offloading.txt'..."
    open(unit=26, file='familia_resultados_offloading.txt', status='replace')

    ! Escreve o tempo
    write(26, '(A, F12.6, A)') "Tempo de execucao (OMP_GET_WTIME): ", tempo_total_sim_s1, " segundos"
    write(26, *)
    write(26, '(A, F12.6, A)') "tempo de execucao parte 1 (contagem numfam): ", t_parte1, " segundos"
    write(26, '(A, F12.6, A)') "tempo de execucao parte 2 (calculo pointfam - scan paralelo ): ", t_parte2, " segundos"
    write(26, '(A, F12.6, A)') "tempo de execucao parte 3 (preenchimento nodefam): ", t_parte3, " segundos"
    write(26, *)
    write(26, *) "===================================================="
    write(26, *) "TEMPO DAS CONDICOES DE CONTORNO (BC)"
    write(26, *) "===================================================="
    write(26, '(A, F12.6, A)') "Tempo acumulado das condicoes de contorno: ", tempo_total_bc_gpu_s, " segundos"
    write(26, *)
    write(26, *) "===================================================="
    write(26, *) "TEMPO DE CALCULO DE FORCAS/DANO"
    write(26, *) "===================================================="
    write(26, '(A, F12.6, A)') "Tempo acumulado do calculo de forcas/dano: ", tempo_total_forca_dano_gpu_s, " segundos"
    write(26, *)
    write(26, *) "===================================================="
    write(26, *) "TEMPO DE ESCRITA (I/O)"
    write(26, *) "===================================================="
    write(26, '(A, F12.6, A)') "Tempo acumulado de escrita dos arquivos de saida: ", tempo_total_escrita_gpu_s, " segundos"
    write(26, *)
    write(26, *) "===================================================="
    write(26, *) "TEMPO TOTAL DO PROGRAMA (WALL-CLOCK)"
    write(26, *) "===================================================="
    write(26, '(A, F12.6, A)') "Tempo total do programa: ", tempo_total_programa_gpu_s, " segundos"
    write(26, *)
    write(26, *) "===================================================="
    write(26, *) "INDICE DOS PONTOS DE CADA FAMILIA (POINTFAM)"
    write(26, *) "Formato: (Indice do Ponto, Indice de Inicio em NODEFAM)"
    write(26, *) "===================================================="
    do i = 1, totnode
        write(26, '(I10, A, I10)') i, " , ", pointfam(i)
    enddo
    write(26, *) ""
    write(26, *) "===================================================="
    write(26, *) "PONTOS QUE COMPOEM A FAMILIA (NODEFAM)"
    write(26, *) "Formato: Familia do Ponto X (Total: Y)"
    write(26, *) "         [lista de pontos j]"
    write(26, *) "===================================================="
    do i = 1, totnode
        write(26, '(A, I10, A, I6, A)') "Familia do Ponto ", i, " (Total: ", numfam(i), ")"

        if (numfam(i) > 0) then
            ! Escreve todos os 'j' membros da família 'i'
            ! Usamos advance='no' para tentar colocar vários na mesma linha
            do j = 1, numfam(i)
                write(26, '(I10, 1x)', advance='no') nodefam(pointfam(i) + j - 1)
            enddo
            write(26, *) ! Quebra de linha para a proxima familia
        else
            write(26, '(A)') "  (Familia vazia)"
        endif
        write(26, *) ! Linha em branco extra para separar
    enddo

    close(26)
    print *, "Escrita do arquivo concluida."

    deallocate(nodefam)

end program openmp_node_offloading