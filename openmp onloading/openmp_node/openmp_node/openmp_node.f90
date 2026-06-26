program openmp_node

    use omp_lib       
    implicit none

    ! =========================================================================
    ! PARÂMETROS E CONSTANTES
    ! =========================================================================
    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: ndivx = 500
    integer, parameter :: ndivy = 500
    integer, parameter :: nbnd = 3
    integer, parameter :: totnode = ndivx * (ndivy + 2 * nbnd)
    integer, parameter :: nt = 1250
    integer, parameter :: maxfam = 100

    ! =========================================================================
    ! DECLARAÇÕES DE VARIÁVEIS
    ! =========================================================================
    real(dp) :: length, width, dx, delta, thick, dens, emod
    real(dp) :: area, vol, bc
    real(dp) :: sedload1, sedload2, dt, totime, ctime, idist, fac
    real(dp) :: radij, nlength, dforce1, dforce2
    real(dp) :: crlength, pi, dmgpar1, dmgpar2, theta, scx, scy, scr, scr0
    

    integer  :: nnum, cnode, i, j, tt, totint, totbottom, tottop
    integer  :: kount, total_family_size

    ! Arrays 
    integer  :: numfam(totnode)
    integer  :: pointfam(totnode)
    real(dp) :: dmg(totnode)

    ! Arrays dinâmicos alocados em tempo de execução para evitar estouro de memória
    integer, allocatable :: nodefam(:)
    integer, allocatable :: prefix_offsets(:)

    ! Matrizes padrão
    real(dp) :: coord(totnode,2), pforce(totnode,2), bforce(totnode,2), stendens(totnode,2)
    real(dp) :: fncst(totnode,2), disp(totnode,2), vel(totnode,2), acc(totnode,2)
    integer  :: fail(totnode, maxfam)

    ! Variáveis de tempo
    real(dp) :: start_time, end_time
    real(dp) :: t_parte1, t_parte2, t_parte3, tempo_total_sim_s1
    real(dp) :: tempo_total_bc_s = 0.0_dp
    real(dp) :: tempo_total_forca_s = 0.0_dp
    real(dp) :: t_inicio_io, t_final_io, tempo_total_io_s = 0.0_dp
    real(dp) :: t_inicio_total, t_final_total, tempo_wall_clock = 0.0_dp

    ! Variáveis para o Scan Paralelo
    integer :: tid, nthreads, i_start, i_end, my_sum, my_offset, global_accum

    ! ======================== INÍCIO DO CÓDIGO EXECUTÁVEL ========================

    t_inicio_total = omp_get_wtime()
    
    
    pi = acos(-1.0_dp)

    ! Inicialização das variaveis
    numfam   = 0
    pointfam = 0
    coord    = 0.0_dp
    pforce   = 0.0_dp
    bforce   = 0.0_dp
    stendens = 0.0_dp
    fncst    = 1.0_dp
    disp     = 0.0_dp
    vel      = 0.0_dp
    acc      = 0.0_dp
    dmg      = 0.0_dp
    fail     = 1 

    ! Propriedades Físicas e da Malha
    length   = 0.05_dp
    width    = 0.05_dp
    dx       = length / ndivx
    delta    = 3.015_dp * dx
    thick    = dx
    dens     = 8000.0_dp
    emod     = 192.0d9
    area     = dx * dx
    vol      = area * dx
    bc       = 9.0_dp * emod / (pi * thick * (delta**3))
    sedload1 = 9.0_dp / 16.0_dp * emod * 1.0d-6
    sedload2 = 9.0_dp / 16.0_dp * emod * 1.0d-6
    dt       = 0.8_dp * sqrt(2.0_dp * dens * dx / (pi * delta**2 * dx * bc))
    totime   = nt * dt
    ctime    = 0.0_dp
    radij    = dx / 2.0_dp
    crlength = 0.01_dp
    scr0     = 0.04472_dp
    nnum     = 0

    ! Geração dos nós - Região Interna
    do i = 1, ndivy
        do j = 1, ndivx
            nnum = nnum + 1
            coord(nnum,1) = (-1.0_dp * length / 2.0_dp) + (dx / 2.0_dp) + (j-1) * dx
            coord(nnum,2) = (-1.0_dp * width / 2.0_dp) + (dx / 2.0_dp) + (i-1) * dx
        enddo
    enddo
    totint = nnum

    ! Geração dos nós - Borda Inferior
    do i = 1, nbnd
        do j = 1, ndivx
            nnum = nnum + 1
            coord(nnum,1) = -1.0_dp / 2.0_dp * length + (dx / 2.0_dp) + (j - 1) * dx
            coord(nnum,2) = -1.0_dp / 2.0_dp * width - (dx / 2.0_dp) - (i - 1) * dx
        enddo
    enddo
    totbottom = nnum

    ! Geração dos nós - Borda Superior
    do i = 1, nbnd
        do j = 1, ndivx
            nnum = nnum + 1
            coord(nnum,1) = -1.0_dp / 2.0_dp * length + (dx / 2.0_dp) + (j - 1) * dx
            coord(nnum,2) =  1.0_dp / 2.0_dp * width + (dx / 2.0_dp) + (i - 1) * dx
        enddo
    enddo
    tottop = nnum

    ! ------------------ Início do código paralelizado (Parte 1: Contagem) ------------------
    ! Usa OpenMP para descobrir quantos vizinhos existem dentro do horizonte (delta)
    start_time = omp_get_wtime()
    
    !$OMP PARALLEL DO PRIVATE(j, idist) SHARED(numfam, coord, delta)
    do i = 1, totnode
        do j = 1, totnode
            if (i /= j) then
                idist = sqrt((coord(j,1) - coord(i,1))**2 + (coord(j,2) - coord(i,2))**2)
                if (idist <= delta) then
                    numfam(i) = numfam(i) + 1
                endif
            endif
        enddo
    enddo
    !$OMP END PARALLEL DO
    
    end_time = omp_get_wtime()   
    t_parte1 = end_time - start_time

    ! ------------------ Segunda parte: Calcula pointfam usando Scan ------------------
    ! Utiliza o Prefix Sum paralelo para determinar onde a vizinhança de cada nó deve começar.
    start_time = omp_get_wtime() 

    ! Alocação do offset dinâmico,
    allocate(prefix_offsets(omp_get_max_threads() + 1))
    prefix_offsets = 0

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
            
    ! Calcula o tamanho total que o array de vizinhos precisará ter e aloca
    total_family_size = pointfam(totnode) + numfam(totnode) - 1
    allocate(nodefam(total_family_size))
    nodefam = 0

    ! ------------------ Terceira parte: Preencher a lista de famílias ------------------
    start_time = omp_get_wtime()

    !$OMP PARALLEL DO PRIVATE(j, idist, kount) SHARED(pointfam, nodefam, coord, delta)
    do i = 1, totnode
        kount = 0
        do j = 1, totnode
            if (i /= j) then
                idist = sqrt((coord(j,1) - coord(i,1))**2 + (coord(j,2) - coord(i,2))**2)
                if (idist <= delta) then
                    nodefam(pointfam(i) + kount) = j
                    kount = kount + 1
                endif
            endif
        enddo
    enddo
    !$OMP END PARALLEL DO

    end_time = omp_get_wtime()   
    t_parte3 = end_time - start_time
    tempo_total_sim_s1 = t_parte1 + t_parte2 + t_parte3
    ! ------------------ Fim Do trecho paralelizado da topologia ------------------
        
    ! Trinca (Definition of the crack surface)
    do i = 1, totnode
        do j = 1, numfam(i) 
            cnode = nodefam(pointfam(i) + j - 1)
            if ((coord(cnode,2) > 0.0_dp) .and. (coord(i,2) < 0.0_dp)) then
                if ((abs(coord(i,1)) - (crlength / 2.0_dp)) <= 1.0d-10) then
                    fail(i,j) = 0
                elseif ((abs(coord(cnode,1)) - (crlength / 2.0_dp)) <= 1.0d-10) then
                    fail(i,j) = 0
                endif
            elseif ((coord(i,2) > 0.0_dp) .and. (coord(cnode,2) < 0.0_dp)) then
                if((abs(coord(i,1)) - (crlength / 2.0_dp)) <= 1.0d-10) then 
                    fail(i,j) = 0
                elseif((abs(coord(cnode,1)) - (crlength / 2.0_dp)) <= 1.0d-10) then
                    fail(i,j) = 0
                endif
            endif        
        enddo
    enddo

    ! Loading 1
    disp(:,1) = 0.001_dp * coord(:,1)
    disp(:,2) = 0.0_dp

    do i = 1, totnode
        stendens(i,1) = 0.0_dp
        do j = 1, numfam(i) 
            cnode = nodefam(pointfam(i) + j - 1) 
            idist = sqrt((coord(cnode,1) - coord(i,1))**2 + (coord(cnode,2) - coord(i,2))**2)
            nlength = sqrt((coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1))**2 + &
                           (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2))**2)
            if (idist <= delta - radij) then
                fac = 1.0_dp
            elseif (idist <= delta + radij) then
                fac = (delta + radij - idist) / (2.0_dp * radij)
            else
                fac = 0.0_dp
            endif
            stendens(i,1) = stendens(i,1) + 0.5_dp * 0.5_dp * bc * ((nlength - idist) / idist)**2 * idist * vol * fac  
        enddo
        fncst(i,1) = sedload1 / stendens(i,1)
    enddo
        
    ! Loading 2
    disp(:,1) = 0.0_dp
    disp(:,2) = 0.001_dp * coord(:,2)

    do i = 1, totnode
        stendens(i,2) = 0.0_dp
        do j = 1, numfam(i)
            cnode = nodefam(pointfam(i) + j - 1) 
            idist = sqrt((coord(cnode,1) - coord(i,1))**2 + (coord(cnode,2) - coord(i,2))**2)
            nlength = sqrt((coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1))**2 + &
                           (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2))**2)
            if (idist <= delta - radij) then
                fac = 1.0_dp
            elseif (idist <= delta + radij) then
                fac = (delta + radij - idist) / (2.0_dp * radij)
            else
                fac = 0.0_dp
            endif   
            stendens(i,2) = stendens(i,2) + 0.5_dp * 0.5_dp * bc * ((nlength - idist) / idist)**2 * idist * vol * fac 
        enddo
        fncst(i,2) = sedload2 / stendens(i,2)
    enddo
        
    ! Reset inicial antes do laço temporal
    vel = 0.0_dp
    disp = 0.0_dp         

    ! ======================== LOOP DE INTEGRAÇÃO TEMPORAL ========================
    do tt = 1, nt
        ctime = tt * dt
        
        ! ------------------ Boundary Conditions (BC) Paralelo ------------------
        start_time = omp_get_wtime() 
        
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
        !$OMP DO
        do i = (totint + 1), totbottom
          vel(i,2) = -20.0_dp
          disp(i,2) = -20.0_dp * tt * dt
        enddo
        !$OMP END DO

        !$OMP DO
        do i = (totbottom + 1), tottop
          vel(i,2) = 20.0_dp
          disp(i,2) = 20.0_dp * tt * dt
        enddo   
        !$OMP END DO
        !$OMP END PARALLEL
        
        end_time = omp_get_wtime()   
        tempo_total_bc_s = tempo_total_bc_s + (end_time - start_time) 
        
        ! ------------------ Cálculo de Forças e Dano Paralelo ------------------
        start_time = omp_get_wtime() 
        
        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP PRIVATE(i, j, dmgpar1, dmgpar2, cnode, idist, nlength, fac, &
        !$OMP         theta, scx, scy, scr, dforce1, dforce2) &
        !$OMP SCHEDULE(GUIDED)
        do i = 1, totnode
            dmgpar1 = 0.0_dp
            dmgpar2 = 0.0_dp
            pforce(i,1) = 0.0_dp
            pforce(i,2) = 0.0_dp

            do j = 1, numfam(i)          
                cnode = nodefam(pointfam(i) + j - 1) 
                idist = sqrt((coord(cnode,1) - coord(i,1))**2 + (coord(cnode,2) - coord(i,2))**2)
                nlength = sqrt((coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1))**2 + &
                               (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2))**2)
                
                ! Volume correction
                if (idist <= delta - radij) then
                    fac = 1.0_dp
                elseif (idist <= delta + radij) then
                    fac = (delta + radij - idist) / (2.0_dp * radij)
                else
                    fac = 0.0_dp
                endif

                if (abs(coord(cnode,2) - coord(i,2)) <= 1.0d-10) then 
                    theta = 0.0_dp
                elseif (abs(coord(cnode,1) - coord(i,1)) <= 1.0d-10) then
                    theta = 90.0_dp * pi / 180.0_dp
                else
                    theta = atan(abs(coord(cnode,2) - coord(i,2)) / abs(coord(cnode,1) - coord(i,1)))
                endif
                
                ! Surface correction
                scx = (fncst(i,1) + fncst(cnode,1)) / 2.0_dp
                scy = (fncst(i,2) + fncst(cnode,2)) / 2.0_dp
                scr = 1.0_dp / (((cos(theta))**2 / (scx)**2) + ((sin(theta))**2 / (scy)**2))
                scr = sqrt(scr)
                
                if (fail(i,j) == 1) then
                    dforce1 = bc * (nlength - idist) / idist * vol * scr * fac * &
                              (coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1)) / nlength              
                    dforce2 = bc * (nlength - idist) / idist * vol * scr * fac * &
                              (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2)) / nlength              
                else
                    dforce1 = 0.0_dp
                    dforce2 = 0.0_dp
                endif 
                pforce(i,1) = pforce(i,1) + dforce1              
                pforce(i,2) = pforce(i,2) + dforce2              
                
                ! No-fail zone              
                if (abs((nlength - idist) / idist) > scr0) then
                    if (abs(coord(i,2)) <= (length / 4.0_dp)) then
                        fail(i,j) = 0 
                    endif
                endif                      
                            
                dmgpar1 = dmgpar1 + fail(i,j) * vol * fac
                dmgpar2 = dmgpar2 + vol * fac              
            enddo

            ! Calcula dano e previne divisão por zero em pontos sem família válida
            if (dmgpar2 > 0.0_dp) then
                dmg(i) = 1.0_dp - dmgpar1 / dmgpar2
            else
                dmg(i) = 0.0_dp
            end if
        enddo
        !$OMP END PARALLEL DO
        
        end_time = omp_get_wtime() 
        tempo_total_forca_s = tempo_total_forca_s + (end_time - start_time) 
        
        ! --- Integração Cinemática Explícita ---
        do i = 1, totint
            acc(i,1) = (pforce(i,1) + bforce(i,1)) / dens
            acc(i,2) = (pforce(i,2) + bforce(i,2)) / dens
            vel(i,1) = vel(i,1) + acc(i,1) * dt
            vel(i,2) = vel(i,2) + acc(i,2) * dt
            disp(i,1) = disp(i,1) + vel(i,1) * dt
            disp(i,2) = disp(i,2) + vel(i,2) * dt
        enddo
        
        do i = (totint + 1), totbottom
            acc(i,1) = (pforce(i,1) + bforce(i,1)) / dens
            vel(i,1) = vel(i,1) + acc(i,1) * dt
            disp(i,1) = disp(i,1) + vel(i,1) * dt
        enddo

        do i = (totbottom + 1), tottop
            acc(i,1) = (pforce(i,1) + bforce(i,1)) / dens
            vel(i,1) = vel(i,1) + acc(i,1) * dt
            disp(i,1) = disp(i,1) + vel(i,1) * dt
        enddo
                
        ! Snapshots
        if (tt == 750) then
            open(26, file = 'coord_disp_pd_750_pwc_v20.txt')
            do i = 1, totint
                write(26, '(5(E12.5, 3X))') coord(i,1), coord(i,2), disp(i,1), disp(i,2), dmg(i)
            enddo
            close(26)
        elseif (tt == 1000) then
            open(26, file = 'coord_disp_pd_1000_pwc_v20.txt')
            do i = 1, totint
                write(26, '(5(E12.5, 3X))') coord(i,1), coord(i,2), disp(i,1), disp(i,2), dmg(i)
            enddo
            close(26)
        elseif (tt == 1250) then
            open(26, file = 'coord_disp_pd_1250_pwc_v20.txt')
            do i = 1, totint
                write(26, '(5(E12.5, 3X))') coord(i,1), coord(i,2), disp(i,1), disp(i,2), dmg(i)
            enddo
            close(26)
        endif
    enddo

    ! --------------------------- Escrita de Resultados e Performance -------------------------------   
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
                write(26, '(I10, 1X)', advance='no') nodefam(pointfam(i) + j - 1)
            enddo
            write(26, *) 
        else
            write(26, '(A)') "  (Familia vazia)"
        endif
        write(26, *) 
    enddo

    t_final_io = omp_get_wtime()
    tempo_total_io_s = t_final_io - t_inicio_io
    t_final_total = omp_get_wtime()
    tempo_wall_clock = t_final_total - t_inicio_total

    ! Relatório impresso direto no arquivo txt (Sem prints no terminal)
    write(26, *) ""
    write(26, *) "===== RESUMO DE PERFORMANCE (OPENMP CPU) ====="
    write(26, '(A, F12.6, A)') "Tempo Total da Topologia (Fase 1+2+3): ", tempo_total_sim_s1, " s"
    write(26, '(A, F12.6, A)') "  -> Parte 1 (Contagem):                 ", t_parte1, " s"
    write(26, '(A, F12.6, A)') "  -> Parte 2 (Scan/Prefix Sum):          ", t_parte2, " s"
    write(26, '(A, F12.6, A)') "  -> Parte 3 (Preenchimento):            ", t_parte3, " s"
    write(26, '(A, F12.6, A)') "Tempo das Condicoes (BC):              ", tempo_total_bc_s, " s"
    write(26, '(A, F12.6, A)') "Tempo de Forcas/Dano:                  ", tempo_total_forca_s, " s"
    write(26, '(A, F12.6, A)') "Tempo de Escrita (I/O):                ", tempo_total_io_s, " s"
    write(26, '(A, F12.6, A)') "Tempo Total do Programa (Wall-clock):  ", tempo_wall_clock, " s"
    write(26, *) "=============================================="

    close(26)

end program openmp_node