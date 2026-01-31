module kernels_peridinamica
    use cudafor
    implicit none

contains

    ! KERNEL 1: Conta quantos vizinhos cada nó tem
    attributes(global) subroutine conta_vizinhos_gpu(n, coord, delta, numfam)
        integer, value :: n
        real(8), value :: delta
        real(8), device :: coord(n, 2)
        integer, device :: numfam(n)
        
        integer :: i, j
        real(8) :: dist, dx, dy

        ! Identifica qual thread é responsável por qual nó i
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

    ! KERNEL 2: Preenche a lista nodefam 
    attributes(global) subroutine preenche_nodefam_gpu(n, coord, delta, pointfam, nodefam)
        integer, value :: n
        real(8), value :: delta
        real(8), device :: coord(n, 2)
        integer, device :: pointfam(n)
        integer, device :: nodefam(*) ! Ponteiro para array gigante
        
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
                    nodefam(inicio + k) = j
                    k = k + 1
                end if
            end do
        end if
    end subroutine preenche_nodefam_gpu

  attributes(global) subroutine calcula_forca_gpu(n, coord, disp, pointfam, nodefam, &
                                                    numfam, bc, scr0, pforce, dmg)
        integer, value :: n
        real(8), value :: bc, scr0
        real(8), device :: coord(n, 2)    ! Coordenada Inicial (Referência)
        real(8), device :: disp(n, 2)     ! Deslocamento atual
        integer, device :: pointfam(n)    ! Onde começa a lista de vizinhos
        integer, device :: nodefam(*)     ! Lista de vizinhos
        integer, device :: numfam(n)      ! Quantos vizinhos
        real(8), device :: pforce(n, 2)   ! SAÍDA: Força resultante (x, y)
        real(8), device :: dmg(n)         ! SAÍDA: Índice de dano
        
        integer :: i, j, k, start_idx, num_neigh
        real(8) :: idist, ndist, stretch, force_mod, fx, fy
        real(8) :: xi, yi, xj, yj, uxi, uyi, uxj, uyj
        real(8) :: dx, dy, cdx, cdy
        integer :: intact_bonds

        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x

        if (i <= n) then
            ! Inicializa acumuladores para este nó
            fx = 0.0d0
            fy = 0.0d0
            intact_bonds = 0
            
            ! Carrega dados do nó i para registradores 
            xi = coord(i, 1)
            yi = coord(i, 2)
            uxi = disp(i, 1)
            uyi = disp(i, 2)
            
            start_idx = pointfam(i)
            num_neigh = numfam(i)
            
            ! Se não tiver vizinhos, sai
            if (num_neigh == 0) return

            ! Loop pelos vizinhos do nó i
            do k = 0, num_neigh - 1
                j = nodefam(start_idx + k)
                
                ! Coordenadas do vizinho j
                xj = coord(j, 1)
                yj = coord(j, 2)
                uxj = disp(j, 1)
                uyj = disp(j, 2)

                ! 1. Distância Inicial 
                dx = xj - xi
                dy = yj - yi
                idist = sqrt(dx*dx + dy*dy)
                
                ! Evita divisão por zero (não deve ocorrer se i!=j)
                if (idist < 1.0d-12) cycle

                ! 2. Distância Deformada 
                cdx = (xj + uxj) - (xi + uxi)
                cdy = (yj + uyj) - (yi + uyi)
                ndist = sqrt(cdx*cdx + cdy*cdy)

                ! 3. Cálculo do Stretch 
                stretch = (ndist - idist) / idist

                ! 4. Critério de Falha 
                ! Se stretch > scr0, a ligação quebra (não soma força)
                if (stretch < scr0) then
                    ! Ligação intacta
                    intact_bonds = intact_bonds + 1
                    
                    ! Força Escalar (Lei constitutiva linear elástica)
                    ! F = c * s * (direção) * vol_correction (simplificado aqui)
                    force_mod = bc * stretch 
                    
                    ! Projeta componentes (cos theta, sin theta)
                    fx = fx + force_mod * (cdx / ndist)
                    fy = fy + force_mod * (cdy / ndist)
                end if
            end do
            
            ! Grava resultado na memória global
            pforce(i, 1) = fx
            pforce(i, 2) = fy
            
            ! Atualiza Dano: 1 - (ligações atuais / ligações iniciais.
            dmg(i) = 1.0d0 - (real(intact_bonds, 8) / real(num_neigh, 8))
        end if
    end subroutine calcula_forca_gpu



    end module 
    
program node_fortran_cuda
    use cudafor 
    use kernels_peridinamica 
    implicit none

    !  PARAMETROS DE MALHA 
    integer, parameter :: ndivx = 500
    integer, parameter :: ndivy = 500
    integer, parameter :: nbnd = 3
    integer, parameter :: totnode = ndivx*(ndivy + 2 * nbnd) 
    integer, parameter :: nt = 2000 
    integer, parameter :: size_nodefam = 40000000 

    ! VARIAVEIS DA FÍSICA 
    real(8) :: length, width, dx, delta, thick, dens, emod, vol, bc, scr0
    real(8) :: pi, dt, v_pull, limit_y_bot, limit_y_top
    
    ! ARRAYS
    real(8), allocatable :: coord(:,:)
    integer, allocatable :: numfam(:), pointfam(:)
    real(8), allocatable :: disp_final(:,:), dmg_final(:)

    real(8), device, allocatable :: coord_d(:,:)
    integer, device, allocatable :: numfam_d(:), pointfam_d(:), nodefam_d(:)
    real(8), device, allocatable :: pforce_d(:,:), vel_d(:,:), disp_d(:,:), dmg_d(:)

    ! Controle
    type(dim3) :: grid, block
    integer :: threadsPorBloco, numBlocos, err_cuda, i, j, nnum, totint
    integer :: t_inicio, t_final, taxa_clock
    real(8) :: tempo_total

    pi = 3.141592653589793d0

    ! CONFIGURAÇÃO
    length = 0.05d0
    width = 0.05d0
    dx = length / ndivx
    delta = 3.015d0 * dx
    thick = dx
    dens = 8000.0d0
    emod = 192.0d9
    bc = 9.0d0 * emod / (pi * thick * (delta**3)) 
    scr0 = 0.0005d0 
    vol = dx * dx * thick
    dt = 0.8d0 * sqrt(2.0d0 * dens * dx / (pi * delta**2 * dx * bc))
    v_pull = 1.0d0 

    limit_y_bot = (-width/2.0d0) + (nbnd * dx)
    limit_y_top = (width/2.0d0) - (nbnd * dx)

    print *, "== SIMULACAO PERIDINAMICA CUDA =="
    print *, "Nos:", totnode

    !  ALOCAÇÃO CPU
    allocate(coord(totnode, 2))
    allocate(numfam(totnode))
    allocate(pointfam(totnode))
    
    allocate(disp_final(totnode, 2))
    allocate(dmg_final(totnode))
 
    
    ! Gera Malha
    nnum = 0
    do i = 1,ndivy
        do j = 1,ndivx
            nnum = nnum + 1
            coord(nnum,1) = (-1.0d0 * length / 2.0d0) + (dx / 2.0d0) + (j-1) * dx
            coord(nnum,2) = (-1.0d0 * width / 2.0d0) + (dx / 2.0d0) + (i-1) * dx
        enddo
    enddo
    do i = 1,nbnd
        do j = 1,ndivx
            nnum = nnum + 1
            coord(nnum,1) = -1.0d0 /2.0d0 * length + (dx / 2.0d0) + (j - 1) * dx
            coord(nnum,2) = -1.0d0 /2.0d0 * width - (dx / 2.0d0) - (i - 1) * dx
        enddo
    enddo
    do i = 1,nbnd
        do j = 1,ndivx
            nnum = nnum + 1
            coord(nnum,1) = -1.0d0 /2.0d0 * length + (dx / 2.0d0) + (j - 1) * dx
            coord(nnum,2) = 1.0d0 /2.0d0 * width + (dx / 2.0d0) + (i - 1) * dx
        enddo
    enddo

    !  PREPARAÇÃO GPU 
    print *, "Alocando memoria na GPU..."
    
    ! Separe estas alocações também para evitar novos erros F-0000
    allocate(coord_d(totnode, 2))
    
    allocate(numfam_d(totnode))
    allocate(pointfam_d(totnode))
    allocate(nodefam_d(size_nodefam))
    
    allocate(pforce_d(totnode, 2))
    allocate(vel_d(totnode, 2))
    allocate(disp_d(totnode, 2))
    allocate(dmg_d(totnode))

    ! índices explícitos (1:totnode, :) para ajudar o compilador
    coord_d(1:totnode, 1:2) = coord(1:totnode, 1:2)
    
    ! Inicialização explícita
    vel_d(1:totnode, 1:2) = 0.0d0
    disp_d(1:totnode, 1:2) = 0.0d0
    pforce_d(1:totnode, 1:2) = 0.0d0
    dmg_d(1:totnode) = 0.0d0

    threadsPorBloco = 256
    numBlocos = (totnode + threadsPorBloco - 1) / threadsPorBloco
    block = dim3(threadsPorBloco, 1, 1)
    grid = dim3(numBlocos, 1, 1)

    !  BUSCA DE VIZINHOS 
    print *, "Gerando lista de vizinhos (GPU)..."
    call SYSTEM_CLOCK(count=t_inicio, count_rate=taxa_clock) 

    call conta_vizinhos_gpu<<<grid, block>>>(totnode, coord_d, delta, numfam_d)
    
    ! Trazendo dados de volta 
    numfam(1:totnode) = numfam_d(1:totnode) 
    
    pointfam(1) = 1
    do i = 2, totnode
        pointfam(i) = pointfam(i-1) + numfam(i-1)
    enddo
    
    if ((pointfam(totnode)+numfam(totnode)) > size_nodefam) then
        print *, "ERRO: Aumente size_nodefam"
        stop
    endif
    
    ! Enviando dados 
    pointfam_d(1:totnode) = pointfam(1:totnode)

    call preenche_nodefam_gpu<<<grid, block>>>(totnode, coord_d, delta, pointfam_d, nodefam_d)
    err_cuda = cudaDeviceSynchronize()
    
    print *, "Vizinhanca concluida. Total interacoes:", pointfam(totnode)

    !  LOOP DE TEMPO
    print *, "Iniciando Loop de Tempo..."
    
    do i = 1, nt
        call calcula_forca_gpu<<<grid, block>>>(totnode, coord_d, disp_d, pointfam_d, &
                                                nodefam_d, numfam_d, bc, scr0, pforce_d, dmg_d)
        
        call aplica_bc_gpu<<<grid, block>>>(totnode, coord_d, vel_d, disp_d, &
                                            limit_y_bot, limit_y_top, v_pull)

        call integracao_tempo_gpu<<<grid, block>>>(totnode, pforce_d, vel_d, disp_d, &
                                                   dt, dens, vol)
        
        if (mod(i, 500) == 0) print *, "Passo:", i, "/", nt
    end do
    
    err_cuda = cudaDeviceSynchronize()

    call SYSTEM_CLOCK(count=t_final)
    tempo_total = real(t_final - t_inicio) / real(taxa_clock)
    print *, "Simulacao concluida em:", tempo_total, "s"

    ! RESULTADOS 
    print *, "Baixando resultados da GPU..."
    
    ! Copia de volta usando índices explícitos
    disp_final(1:totnode, 1:2) = disp_d(1:totnode, 1:2)
    dmg_final(1:totnode) = dmg_d(1:totnode)
    
    open(unit=26, file='resultado_final_cuda.txt', status='replace')
    write(26, *) "X Y DispX DispY Dano"
    do i = 1, totnode
        write(26, '(5E15.6)') coord(i,1), coord(i,2), disp_final(i,1), disp_final(i,2), dmg_final(i)
    enddo
    close(26)
    print *, "Arquivo 'resultado_final_cuda.txt' gerado."

    deallocate(coord_d, numfam_d, pointfam_d, nodefam_d, pforce_d, vel_d, disp_d, dmg_d)

end program node_fortran_cuda
    

    
