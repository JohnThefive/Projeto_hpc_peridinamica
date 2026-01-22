program openmp_node

  use OMP_LIB       
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

! mudando as variaveis para 1d 
integer numfam(totnode), pointfam(totnode)
integer nodefam(10000000,1), fail(totnode,maxfam)


! TESTANDO A FUNÇÃO DO OPENMP
integer :: kount, scan_sum  ! 'scan_sum' é a variável de acúmulo
real *8 :: start_time, end_time
real *8 :: tempo_total_sim_s1

! variaveis para ajudar com o scan paralelo 
integer :: prefix_offsets(256) ! Vetor auxiliar para soma das threads
integer :: tid, nthreads, i_start, i_end, my_sum, my_offset, global_accum

pi = dacos(-1.0d0)

do i = 1, totnode
    ! numfam, pointfam agora são acessados como 1D
    numfam(i) = 0
    pointfam(i) = 0
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
    dmg(i,1) = 0.0d0
enddo

do i = 1, 10000000 ! <- Corrigido de 1.000.000 para 10.000.000
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



! inicio do codigo paralelizado (parte 1)
start_time = omp_get_wtime()
!$OMP PARALLEL DO PRIVATE(j, idist) SHARED(numfam, coord, delta)
do i = 1, totnode
    numfam(i) = 0  ! Acesso 1D
    do j = 1, totnode
        if (i .ne. j) then
            idist = dsqrt((coord(j,1) - coord(i,1))**2 + (coord(j,2) - coord(i,2))**2)
            if (idist <= delta) then
                numfam(i) = numfam(i) + 1 ! Acesso 1D
            endif
        endif
    enddo
enddo
!$OMP END PARALLEL DO

! segunda parte - calcula pointfam usando scan  (prefix sum) paralelo
! Existe uma função pra fazer essa operção menos verbosa, mas o compilador ainda não suporta

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, tid, nthreads, i_start, i_end, my_sum, my_offset)
   tid = omp_get_thread_num() 
   nthreads = omp_get_num_threads()
   
   ! divisão para distribuir o trabalho 
   i_start = (totnode * tid) / nthreads + 1
   i_end   = (totnode * (tid + 1)) / nthreads
   
   !scan local manual 
   my_sum = 0
   do i = i_start, i_end
       pointfam(i) = my_sum
       my_sum = my_sum + numfam(i)
   end do 
   
   !salva o total parcial de cada therad 
   prefix_offsets(tid + 1) = my_sum
   
   ! Espera todas as threads terminarem a Fase 1
   !$OMP BARRIER
   
   !  Cálculo dos Offsets Globais (Feito por uma só thread)
   !$OMP SINGLE
     global_accum = 1
     do i = 1, nthreads
         my_offset = prefix_offsets(i)   ! Pega o total da thread 'i'
         prefix_offsets(i) = global_accum ! Define onde a thread 'i' deve começar
         global_accum = global_accum + my_offset ! Incrementa o acumulador global
     end do 
   !$OMP END SINGLE
   !$OMP BARRIER  
     
   my_offset = prefix_offsets(tid + 1)
   do i = i_start, i_end
       pointfam(i) = pointfam(i) + my_offset
   end do
!$OMP END PARALLEL
       

! terceira parte: Preencher a lista de famílias (nodefam)
!$OMP PARALLEL DO PRIVATE(j, idist, kount) SHARED(pointfam, nodefam, coord, delta)
do i = 1, totnode
    kount = 0
    do j = 1, totnode
        if (i .ne. j) then
            idist = dsqrt((coord(j,1) - coord(i,1))**2 + (coord(j,2) - coord(i,2))**2)
            if (idist <= delta) then
                kount = kount + 1
                ! pointfam(i) agora é 1D
                nodefam(pointfam(i) + kount - 1, 1) = j
            endif
        endif
    enddo
enddo
!$OMP END PARALLEL DO


end_time = omp_get_wtime()

tempo_total_sim_s1 = end_time - start_time

print *, "Tempo total de execucao (OMP_GET_WTIME): ", tempo_total_sim_s1, " segundos"


! 1) Escrever um arquivo de saida com os seguintes dados:
!    Indice dos pontos de cada familia (pointfam)
!    Pontos que compoem a familia de cada ponto (nodefam)
!    Tempo que a simulação levou

print *, "Iniciando a escrita do arquivo 'familia_resultados.txt'..."

open(unit=26, file='familia_resultados.txt', status='replace')

! Escreve o tempo
write(26, '(A, F12.6, A)') "Tempo de execucao (OMP_GET_WTIME): ", tempo_total_sim_s1, " segundos"
write(26, *) ""
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
            write(26, '(I10, 1x)', advance='no') nodefam(pointfam(i) + j - 1, 1)
        enddo
        write(26, *) ! Quebra de linha para a proxima familia
    else
        write(26, '(A)') "  (Familia vazia)"
    endif
    write(26, *) ! Linha em branco extra para separar
enddo

close(26)

print *, "Escrita do arquivo concluida."

end program openmp_node