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

integer numfam(totnode), pointfam(totnode)
integer nodefam(10000000,1), fail(totnode,maxfam)

! TESTANDO A FUNÇÃO DO OPENMP
integer :: kount, scan_sum  ! 'scan_sum' é a variável de acúmulo
real *8 :: start_time, end_time
real *8 :: t_parte1, t_parte2, t_parte3, tempo_total_sim_s1

! para marcar o tempo do loop de boundary conditions 
real*8 :: tempo_total_bc_s
real*8 :: t_inicio_bc, t_final_bc

! Variáveis para marcar o tempo do loop de forças e dano
real*8 :: tempo_total_forca_s
real*8 :: t_inicio_forca, t_final_forca

! Variáveis para o tempo de I/O e Wall-clock total
real*8 :: t_inicio_io, t_final_io, tempo_total_io_s
real*8 :: t_inicio_total, t_final_total, tempo_wall_clock

! variaveis para ajudar com o scan paralelo 
integer :: prefix_offsets(256) ! Vetor auxiliar para soma das threads
integer :: tid, nthreads, i_start, i_end, my_sum, my_offset, global_accum

! Inicializa o cronômetro total do programa
t_inicio_total = omp_get_wtime()

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

do i = 1, 10000000 !
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


! ------------------ Inicio do codigo paralelizado (parte 1) ------------------
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
end_time = omp_get_wtime()   
t_parte1 = end_time - start_time

! ------------------ segunda parte - calcula pointfam usando scan ------------------
start_time = omp_get_wtime() 

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, tid, nthreads, i_start, i_end, my_sum, my_offset)
   tid = omp_get_thread_num() 
   nthreads = omp_get_num_threads()
   
   i_start = (totnode * tid) / nthreads + 1
   i_end   = (totnode * (tid + 1)) / nthreads
   
   my_sum = 0
   do i = i_start, i_end
       pointfam(i) = my_sum
       my_sum = my_sum + numfam(i)
   end do 
   
   prefix_offsets(tid + 1) = my_sum
   
   !$OMP BARRIER
   
   !$OMP SINGLE
     global_accum = 1
     do i = 1, nthreads
         my_offset = prefix_offsets(i)   
         prefix_offsets(i) = global_accum 
         global_accum = global_accum + my_offset 
     end do 
   !$OMP END SINGLE
   !$OMP BARRIER  
     
   my_offset = prefix_offsets(tid + 1)
   do i = i_start, i_end
       pointfam(i) = pointfam(i) + my_offset
   end do
!$OMP END PARALLEL

end_time = omp_get_wtime()   
t_parte2 = end_time - start_time 
        

! ------------------ terceira parte: Preencher a lista de famílias ------------------
start_time = omp_get_wtime()

!$OMP PARALLEL DO PRIVATE(j, idist, kount) SHARED(pointfam, nodefam, coord, delta)
do i = 1, totnode
    kount = 0
    do j = 1, totnode
        if (i .ne. j) then
            idist = dsqrt((coord(j,1) - coord(i,1))**2 + (coord(j,2) - coord(i,2))**2)
            if (idist <= delta) then
                kount = kount + 1
                nodefam(pointfam(i) + kount - 1, 1) = j
            endif
        endif
    enddo
enddo
!$OMP END PARALLEL DO

end_time = omp_get_wtime()   
t_parte3 = end_time - start_time

tempo_total_sim_s1 = t_parte1 + t_parte2 + t_parte3
! ------------------ Fim Do trecho paralelizado da topologia ------------------


    
!Definition of the crack surface
do i = 1,totnode
    do j = 1,numfam(i) 
        cnode = nodefam(pointfam(i)+j-1,1)
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
        cnode = nodefam(pointfam(i)+j-1,1) 
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
        cnode = nodefam(pointfam(i)+j-1,1) 
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
    fncst(i,2) = sedload2 / stendens(i,2)
enddo
    
!Initialization of displacements and velocities
do i = 1,totnode
    vel(i,1) = 0.0d0
    vel(i,2) = 0.0d0
    disp(i,1) = 0.0d0
    disp(i,2) = 0.0d0         
enddo


! Zera os acumuladores de tempo 
tempo_total_bc_s = 0.0d0
tempo_total_forca_s = 0.0d0
    
!Time integration
do tt = 1,nt
    write(*,*) 'tt = ', tt
    ctime = tt * dt
    
    !------------------------------------------------------------------- Daqui ---------------------------------------------------------------------

    t_inicio_bc = omp_get_wtime() 
    
    !$omp parallel default(shared) private(i)
    !$omp do
    do i = (totint+1), totbottom
        vel(i,2) = -20.0d0
        disp(i,2) = -20.0d0 * tt * dt
    enddo
    !$omp end do

    !$omp do
    do i = (totbottom+1), tottop
        vel(i,2) = 20.0d0
        disp(i,2) = 20.0d0 * tt * dt
    enddo   
    !$omp end do
    !$omp end parallel
    
    t_final_bc = omp_get_wtime()   
    tempo_total_bc_s = tempo_total_bc_s + (t_final_bc - t_inicio_bc) 
    
    !------------------------------------------------------------------- Até aqui ---------------------------------------------------------------------
    
    
    !------------------------------------------------------------------- Daqui ---------------------------------------------------------------------

    t_inicio_forca = omp_get_wtime() 
    
    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP PRIVATE(i, j, dmgpar1, dmgpar2, cnode, idist, nlength, fac, &
    !$OMP         theta, scx, scy, scr, dforce1, dforce2) &
    !$OMP SCHEDULE(GUIDED)
    do i = 1,totnode
        dmgpar1 = 0.0d0
        dmgpar2 = 0.0d0
        pforce(i,1) = 0.0d0
        pforce(i,2) = 0.0d0
        do j = 1,numfam(i)          
                cnode = nodefam(pointfam(i)+j-1,1) 
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
        dmg(i,1) = 1.0d0 - dmgpar1 / dmgpar2
    enddo
    !$OMP END PARALLEL DO
    
    t_final_forca = omp_get_wtime() 
    tempo_total_forca_s = tempo_total_forca_s + (t_final_forca - t_inicio_forca) 
    
    !------------------------------------------------------------------- Até aqui ---------------------------------------------------------------------
   
    do i = 1,totint
        acc(i,1) = (pforce(i,1) + bforce(i,1)) / dens
        acc(i,2) = (pforce(i,2) + bforce(i,2)) / dens
        vel(i,1) = vel(i,1) + acc(i,1) * dt
        vel(i,2) = vel(i,2) + acc(i,2) * dt
        disp(i,1) = disp(i,1) + vel(i,1) * dt
        disp(i,2) = disp(i,2) + vel(i,2) * dt
    enddo
    
    do i = (totint+1), totbottom
        acc(i,1) = (pforce(i,1) + bforce(i,1)) / dens
        vel(i,1) = vel(i,1) + acc(i,1) * dt
        disp(i,1) = disp(i,1) + vel(i,1) * dt
    enddo

    do i = (totbottom+1), tottop
            acc(i,1) = (pforce(i,1) + bforce(i,1)) / dens
            vel(i,1) = vel(i,1) + acc(i,1) * dt
            disp(i,1) = disp(i,1) + vel(i,1) * dt
    enddo
               
    endtime(tt,1) = ctime

    if (tt.eq.750) then
        open(26,file = 'coord_disp_pd_750_pwc_v20.txt')
        do i = 1, totint
            write(26,111) coord(i,1), coord(i,2), disp(i,1), disp(i,2), dmg(i,1)
        enddo
        close(26)
    elseif (tt.eq.1000) then
        open(26,file = 'coord_disp_pd_1000_pwc_v20.txt')
        do i = 1, totint
            write(26,111) coord(i,1), coord(i,2), disp(i,1), disp(i,2), dmg(i,1)
        enddo
        close(26)
    elseif (tt.eq.1250) then
        open(26,file = 'coord_disp_pd_1250_pwc_v20.txt')
        do i = 1, totint
            write(26,111) coord(i,1), coord(i,2), disp(i,1), disp(i,2), dmg(i,1)
        enddo
        close(26)
    endif

enddo

! --------------------------- Inicio da escrita dos resultados -------------------------------   
t_inicio_io = omp_get_wtime()

open(unit=26, file='familia_resultados_onloading.txt', status='replace')

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
        do j = 1, numfam(i)
            write(26, '(I10, 1x)', advance='no') nodefam(pointfam(i) + j - 1, 1)
        enddo
        write(26, *) 
    else
        write(26, '(A)') "  (Familia vazia)"
    endif
    write(26, *) 
enddo

! Fechamos os relógios de I/O e Total antes do encerramento do arquivo
t_final_io = omp_get_wtime()
tempo_total_io_s = t_final_io - t_inicio_io

t_final_total = omp_get_wtime()
tempo_wall_clock = t_final_total - t_inicio_total

! Escrevemos o bloco consolidado de performance no fim do log
write(26, *) ""
write(26, *) "===== RESUMO DE PERFORMANCE (PARALELO) ====="
write(26, '(A, F12.6, A)') "Tempo Total da Topologia (Fase 1+2+3): ", tempo_total_sim_s1, " s"
write(26, '(A, F12.6, A)') "  -> Parte 1 (Contagem):                 ", t_parte1, " s"
write(26, '(A, F12.6, A)') "  -> Parte 2 (Scan/Prefix Sum):          ", t_parte2, " s"
write(26, '(A, F12.6, A)') "  -> Parte 3 (Preenchimento):            ", t_parte3, " s"
write(26, '(A, F12.6, A)') "Tempo das Condicoes (BC):              ", tempo_total_bc_s, " s"
write(26, '(A, F12.6, A)') "Tempo de Forcas/Dano:                  ", tempo_total_forca_s, " s"
write(26, '(A, F12.6, A)') "Tempo de Escrita (I/O):                ", tempo_total_io_s, " s"
write(26, '(A, F12.6, A)') "Tempo Total do Programa (Wall-clock):  ", tempo_wall_clock, " s"
write(26, *) "============================================"

close(26)

print *, "Simulacao finalizada com sucesso!"
print *, "O relatorio de performance foi salvo no final de 'familia_resultados_onloading.txt'."

111 format(e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5)

end program openmp_node