program openmp_node_offloading 
    
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
    real(dp) :: crlength, pi, dmgpar1, dmgpar2, theta, scx, scy, scr
    real(dp) :: scr0

    integer  :: nnum, cnode, i, j, tt, totint, totbottom, tottop
    integer  :: kount, total_family_size

    ! Arrays 1D limpos (Removido o antipadrão Nx1)
    integer  :: numfam(totnode)
    integer  :: pointfam(totnode)
    real(dp) :: dmg(totnode)

    ! Arrays dinâmicos alocados em tempo de execução
    integer, allocatable :: nodefam(:)

    ! Matrizes padrão
    real(dp) :: coord(totnode,2), pforce(totnode,2), bforce(totnode,2), stendens(totnode,2)
    real(dp) :: fncst(totnode,2), disp(totnode,2), vel(totnode,2), acc(totnode,2)
    integer  :: fail(totnode,maxfam)

    ! Declarações para registro de tempo 
    real(dp) :: start_time, end_time, t_start, t_end
    real(dp) :: t_parte1, t_parte2, t_parte3, tempo_total_sim_s1
    real(dp) :: t_inicio_bc, t_final_bc, tempo_total_bc_gpu_s = 0.0_dp
    real(dp) :: t_inicio_forca, t_final_forca, tempo_total_forca_s = 0.0_dp
    real(dp) :: t_inicio_io, t_final_io, tempo_total_io_s = 0.0_dp
    real(dp) :: t_inicio_total, t_final_total, tempo_wall_clock = 0.0_dp

    ! ======================== INÍCIO DO CÓDIGO EXECUTÁVEL ========================

    t_inicio_total = omp_get_wtime()
    pi = acos(-1.0_dp)

    ! Inicialização Limpa
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
    fail     = 1 ! Inicializa a falha em 1 diretamente

    ! Propriedades Físicas e Geométricas
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


    ! ======================== TOPOLOGIA (GPU OFFLOADING) ========================
    ! Mapeia apenas coord e delta para a GPU nesta fase
    !$OMP TARGET DATA MAP(TO: coord, delta)
    
    t_start = omp_get_wtime() 

    ! PARTE 1 (GPU): Conta os vizinhos
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(j, idist) MAP(FROM: numfam)
    do i = 1, totnode
        numfam(i) = 0
        do j = 1, totnode
            if (i /= j) then
                idist = sqrt((coord(j,1) - coord(i,1))**2 + (coord(j,2) - coord(i,2))**2)
                if (idist <= delta) then
                    numfam(i) = numfam(i) + 1
                endif
            endif
        enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
    
    t_end = omp_get_wtime() 
    t_parte1 = t_end - t_start 
    
    ! PARTE 2 (CPU - Serial): Scan / Prefix Sum
    t_start = omp_get_wtime() 
    
    pointfam(1) = 1
    do i = 2, totnode
        pointfam(i) = pointfam(i-1) + numfam(i-1)
    enddo
    
    t_end = omp_get_wtime()
    t_parte2 = t_end - t_start
    
    ! Alocação dinâmica de nodefam
    total_family_size = pointfam(totnode) + numfam(totnode) - 1
    allocate(nodefam(total_family_size))
    nodefam = 0
    
    ! PARTE 3 (GPU): Preenchimento
    t_start = omp_get_wtime() 
    
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(j, idist, kount) MAP(TO: pointfam) MAP(FROM: nodefam)
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
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
    
    !$OMP END TARGET DATA

    t_end = omp_get_wtime() 
    t_parte3 = t_end - t_start

    tempo_total_sim_s1 = t_parte1 + t_parte2 + t_parte3


    ! ======================== CÁLCULOS INICIAIS DA CPU ========================
    
    !Definition of the crack surface
    !PD bonds penetrating through the crack surface are broken
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
        
    !  Loading 2
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
        
    ! Reseta velocidade e deslocamento antes do loop
    vel = 0.0_dp
    disp = 0.0_dp         

    ! ======================== LOOP TEMPORAL ========================
    ! Mapeamos todos os dados estáticos necessários para a GPU
    !$OMP TARGET DATA MAP(TO: coord, nodefam, pointfam, numfam, fncst, fail) &
    !$OMP MAP(TOFROM: vel, disp, dmg)

    do tt = 1, nt
        write(*,*) 'tt = ', tt
        ctime = tt * dt

        ! ------------------ Boundary Conditions na GPU ------------------
        t_inicio_bc = omp_get_wtime()
        
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
        do i = (totint+1), totbottom
            vel(i,2) = -20.0_dp
            disp(i,2) = -20.0_dp * tt * dt
        enddo
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
        
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
        do i = (totbottom+1), tottop
            vel(i,2) = 20.0_dp
            disp(i,2) = 20.0_dp * tt * dt
        enddo   
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
        
        t_final_bc = omp_get_wtime()
        tempo_total_bc_gpu_s = tempo_total_bc_gpu_s + (t_final_bc - t_inicio_bc)
        
        ! ------------------ Forças Peridinâmicas na GPU ------------------
        t_inicio_forca = omp_get_wtime()
        
        ! NOTA: Mapeamos pforce como FROM, para que a GPU possa enviar 
        ! as forças de volta para a CPU calcular a cinemática depois.
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(SHARED) &
        !$OMP PRIVATE(j, dmgpar1, dmgpar2, cnode, idist, nlength, fac, &
        !$OMP         theta, scx, scy, scr, dforce1, dforce2) MAP(FROM: pforce)
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
                    dforce1 = bc * (nlength - idist) / idist * vol * scr * fac * (coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1)) / nlength              
                    dforce2 = bc * (nlength - idist) / idist * vol * scr * fac * (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2)) / nlength              
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
            
            ! Evita divisão por zero silenciosa
            if (dmgpar2 > 0.0_dp) then
                dmg(i) = 1.0_dp - dmgpar1 / dmgpar2
            else
                dmg(i) = 0.0_dp
            end if
        enddo
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
        
        t_final_forca = omp_get_wtime() 
        tempo_total_forca_s = tempo_total_forca_s + (t_final_forca - t_inicio_forca)
        
        ! ------------------ Integração Cinemática (CPU) ------------------
        ! A CPU precisa saber as novas velocidades e deslocamentos que a GPU calculou no BC
        !$OMP TARGET UPDATE FROM(vel, disp)

        do i = 1, totint
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

        ! A CPU atualizou vel e disp. Temos que enviar de volta para a GPU para a proxima iteração.
        !$OMP TARGET UPDATE TO(vel, disp)
                   
        if (tt == 750 .or. tt == 1000 .or. tt == 1250) then
            ! A GPU calcula o dano, atualizamos na CPU para salvar no TXT
            !$OMP TARGET UPDATE FROM(dmg)
            
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
        endif
    enddo
    
    !$OMP END TARGET DATA


    ! ======================== I/O E RELATÓRIO ========================
    t_inicio_io = omp_get_wtime()

    open(unit=26, file='familia_resultados_offloading.txt', status='replace')

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
                write(26, '(I10, 1x)', advance='no') nodefam(pointfam(i) + j - 1)
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

    ! Relatório impresso direto no TXT
    write(26, *) ""
    write(26, *) "===== RESUMO DE PERFORMANCE (OFFLOADING GPU) ====="
    write(26, '(A, F12.6, A)') "Tempo Total da Topologia (Fase 1+2+3): ", tempo_total_sim_s1, " s"
    write(26, '(A, F12.6, A)') "  -> Parte 1 (GPU numfam + transfer):    ", t_parte1, " s"
    write(26, '(A, F12.6, A)') "  -> Parte 2 (Serial CPU pointfam):      ", t_parte2, " s"
    write(26, '(A, F12.6, A)') "  -> Parte 3 (GPU nodefam + transfer):   ", t_parte3, " s"
    write(26, '(A, F12.6, A)') "Tempo das Condicoes (BC na GPU):       ", tempo_total_bc_gpu_s, " s"
    write(26, '(A, F12.6, A)') "Tempo de Forcas/Dano (GPU):            ", tempo_total_forca_s, " s"
    write(26, '(A, F12.6, A)') "Tempo de Escrita (I/O):                ", tempo_total_io_s, " s"
    write(26, '(A, F12.6, A)') "Tempo Total do Programa (Wall-clock):  ", tempo_wall_clock, " s"
    write(26, *) "--------------------------------------------------"
    if (tempo_wall_clock > 0.0_dp) then
        write(26, '(A, F12.2, A)') "Impacto das BCs no Tempo Total:        ", (tempo_total_bc_gpu_s / tempo_wall_clock) * 100.0_dp, " %"
        write(26, '(A, F12.2, A)') "Impacto das Forcas no Tempo Total:     ", (tempo_total_forca_s / tempo_wall_clock) * 100.0_dp, " %"
    endif
    write(26, *) "=================================================="

    close(26)

end program openmp_node_offloading