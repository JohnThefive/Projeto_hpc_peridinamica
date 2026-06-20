program main
  use omp_lib
  implicit none

    integer ndivx, ndivy, totnode, nt, maxfam, nnum, cnode, i, j, tt, nbnd, totint, totbottom, tottop
    parameter(ndivx = 500)
    parameter(ndivy = 500)
    parameter(nbnd = 3)
    parameter (totnode = ndivx*(ndivy + 2 * nbnd)) 
    parameter(nt = 1250)
    parameter(maxfam = 100)

    real*8 length, width, dx, delta, thick, dens, emod, pratio, area, vol, bc  
    real*8 sedload1, sedload2, dt, totime, ctime, idist, fac, radij, nlength, dforce1, dforce2 
    real*8 crlength, scr0, pi, tmpdx, tmpvol, tmpcx, tmpcy, tmpux, tmpuy, dmgpar1, dmgpar2, theta 
    real*8 scx, scy, scr

    real*8 coord(totnode,2), pforce(totnode,2), pforceold(totnode,2), bforce(totnode,2), stendens(totnode,2)
    real*8 fncst(totnode,2), disp(totnode,2), vel(totnode,2), velhalfold(totnode,2), velhalf(totnode,2)
    real*8 acc(totnode,2), massvec(totnode,2), enddisp(nt,1), endtime(nt,1), dmg(totnode,1)
    integer numfam(totnode,1), pointfam(totnode,1), nodefam(10000000,1), fail(totnode,maxfam)
    
    ! --- VARIÁVEIS PARA TIMER (OpenMP wtime) ---
    real*8 :: t_inicio_sim, t_final_sim, tempo_total_sim_s ! tempo das distancias das familias de contagem (numfam, nodefam e pointfam) 
    real*8 :: t_inicio_bc, t_final_bc, tempo_total_bc_s
    real*8 :: t_inicio_forca, t_final_forca, tempo_total_forca_s
    real*8 :: t_inicio_io, t_final_io, tempo_total_io_s
    real*8 :: t_inicio_total, t_final_total, tempo_wall_clock  ! Tempo total do programa 
    
    ! ----------------------------
    
    ! Inicia o cronômetro absoluto do programa
    t_inicio_total = omp_get_wtime()

    pi = dacos(-1.0d0)

    do i = 1, totnode 
        coord(i,1) = 0.0d0
        coord(i,2) = 0.0d0
        numfam(i,1) = 0
        pointfam(i,1) = 0
        pforce(i,1) = 0.0d0
        pforce(i,2) = 0.0d0
        pforceold(i,1) = 0.0d0
        pforceold(i,2) = 0.0d0
        bforce(i,1) = 0.0d0
        bforce(i,2) = 0.0d0
        stendens(i,1) = 0.0d0
        stendens(i,2) = 0.0d0
        fncst(i,1) = 1.0d0 
        fncst(i,2) = 1.0d0  
        disp(i,1) = 0.0d0
        disp(i,2) = 0.0d0
        vel(i,1) = 0.0d0
        vel(i,2) = 0.0d0
        velhalfold(i,1) = 0.0d0
        velhalfold(i,2) = 0.0d0
        velhalf(i,1) = 0.0d0
        velhalf(i,2) = 0.0d0
        acc(i,1) = 0.0d0
        acc(i,2) = 0.0d0
        massvec(i,1) = 0.0d0
        massvec(i,2) = 0.0d0
        do j = 1, maxfam
            fail(i,j) = 0
        enddo
        dmg(i,1) = 0.0d0
    enddo

    do i = 1, 10000000
        nodefam(i,1) = 0
    enddo

    length = 0.05d0
    width = 0.05d0
    dx = length / ndivx
    delta = 3.015 * dx
    thick = dx
    dens = 8000.0d0
    emod = 192.0d9
    pratio = 1.0d0 / 3.0d0
    area = dx * dx
    vol = area * dx
    bc = 9.0d0 * emod / (pi * thick * (delta**3))
    sedload1 = 9.0d0 / 16.0d0 * emod * 1.0d-6   
    sedload2 = 9.0d0 / 16.0d0 * emod * 1.0d-6
    dt = 0.8d0 * dsqrt(2.0d0*dens*dx/(pi*delta**2*dx*bc))
    totime = nt * dt
    ctime = 0.0d0
    idist = 0.0d0
    
    do i = 1, nt
        enddisp(i,1) = 0.0d0
        endtime(i,1) = 0.0d0
    enddo
    
    fac = 0.0d0
    radij = dx / 2.0d0
    nnum = 0
    cnode = 0
    nlength  = 0.0d0
    dforce1 = 0.0d0
    dforce2 = 0.0d0
    crlength = 0.01d0
    scr0 = 0.04472d0

    do i = 1,totnode
        do j = 1,maxfam
            fail(i,j) = 1
        enddo
    enddo

    ! Material points of the internal region
    do i = 1,ndivy
        do j = 1,ndivx
            nnum = nnum + 1
            coord(nnum,1) = (-1.0d0 * length / 2.0d0) + (dx / 2.0d0) + (j-1) * dx
            coord(nnum,2) = (-1.0d0 * width / 2.0d0) + (dx / 2.0d0) + (i-1) * dx
        enddo
    enddo
    totint = nnum

    ! Bottom boundary
    do i = 1,nbnd
        do j = 1,ndivx
            nnum = nnum + 1
            coord(nnum,1) = -1.0d0 /2.0d0 * length + (dx / 2.0d0) + (j - 1) * dx
            coord(nnum,2) = -1.0d0 /2.0d0 * width - (dx / 2.0d0) - (i - 1) * dx
        enddo
    enddo
    totbottom = nnum

    ! Top boundary
    do i = 1,nbnd
        do j = 1,ndivx
            nnum = nnum + 1
            coord(nnum,1) = -1.0d0 /2.0d0 * length + (dx / 2.0d0) + (j - 1) * dx
            coord(nnum,2) = 1.0d0 /2.0d0 * width + (dx / 2.0d0) + (i - 1) * dx
        enddo
    enddo
    tottop = nnum
    
# ___________________________________________________ Inicio do calculo DA familia de cada ponto _________________________________________________________________________    

    t_inicio_sim = omp_get_wtime()
    
    do i = 1,totnode
        if (i.eq.1) then 
            pointfam(i,1) = 1
        else
            pointfam(i,1) = pointfam(i-1,1) + numfam(i-1,1)
        endif
        do j = 1,totnode
            idist = dsqrt((coord(j,1) - coord(i,1))**2 + (coord(j,2) - coord(i,2))**2)
            if (i.ne.j) then
                if(idist.le.delta) then
                    numfam(i,1) = numfam(i,1) + 1
                    nodefam(pointfam(i,1)+numfam(i,1)-1,1) = j
                endif
            endif
        enddo
    enddo
    
    t_final_sim = omp_get_wtime()
    tempo_total_sim_s = t_final_sim - t_inicio_sim
# ___________________________________________________ final do calculo DA familia de cada ponto _________________________________________________________________________  
    

    
    !Definition of the crack surface
!PD bonds penetrating through the crack surface are broken
do i = 1,totnode
    do j = 1,numfam(i,1)
        cnode = nodefam(pointfam(i,1)+j-1,1)
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
    do j = 1,numfam(i,1)
        cnode = nodefam(pointfam(i,1)+j-1,1)
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
    do j = 1,numfam(i,1)
        cnode = nodefam(pointfam(i,1)+j-1,1)
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

! contagem de tempo para aplicãção das condições de contorno
tempo_total_bc_s = 0.0d0
!Time integration
do tt = 1,nt
    write(*,*) 'tt = ', tt
	ctime = tt * dt

! ______________________________DAQUI____________________________________________
    t_inicio_bc = omp_get_wtime()
    
    !Application of boundary conditions at the top and bottom edges
    do i = (totint+1), totbottom
        vel(i,2) = -20.0d0
        disp(i,2) = -20.0d0 * tt * dt
    enddo

    do i = (totbottom+1), tottop
        vel(i,2) = 20.0d0
        disp(i,2) = 20.0d0 * tt * dt
    enddo   
    
    t_final_bc = omp_get_wtime()
    tempo_total_bc_s = tempo_total_bc_s + (t_final_bc - t_inicio_bc)
    
!______________________________ATÉ AQUI ____________________________________________ 

!______________________________________AQUI__________________________________________
    t_inicio_forca = omp_get_wtime()
    
    do i = 1,totnode
        dmgpar1 = 0.0d0
        dmgpar2 = 0.0d0
        pforce(i,1) = 0.0d0
        pforce(i,2) = 0.0d0
        do j = 1,numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1)
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
        dmg(i,1) = 1.0d0 - dmgpar1 / dmgpar2
    enddo
    
    !contagem de tempo para o caulculo das forças internas e do dano
    t_final_forca = omp_get_wtime()
    tempo_total_forca_s = tempo_total_forca_s + (t_final_forca - t_inicio_forca)
   ! _________________________________________________ATÉ AQUI_____________________________________________________________ 
    
    
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
        !printing results to an output file
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

! --------------------------- Inicio DA escrita dos resultados -------------------------------

    t_inicio_io = omp_get_wtime()
    
    open(unit=26, file='familia_resultados.txt', status='replace')

    write(26, *) "===================================================="
    write(26, *) "INDICE DOS PONTOS DE CADA FAMILIA (POINTFAM)"
    write(26, *) "===================================================="
    do i = 1, totnode
        write(26, '(I10, A, I10)') i, " , ", pointfam(i, 1)
    enddo

    write(26, *) ""
    write(26, *) "===================================================="
    write(26, *) "PONTOS QUE COMPOEM A FAMILIA (NODEFAM)"
    write(26, *) "===================================================="
    do i = 1, totnode
        write(26, '(A, I10, A, I6, A)') "Familia do Ponto ", i, " (Total: ", numfam(i,1), ")"
        if (numfam(i,1) > 0) then
            do j = 1, numfam(i,1)
                write(26, '(I10, 1x)', advance='no') nodefam(pointfam(i,1) + j - 1, 1)
            enddo
            write(26, *)
        else
            write(26, '(A)') "  (Familia vazia)"
        endif
        write(26, *)
    enddo

    ! Fechamos os cronometros de I/O e Total ANTES de fechar o arquivo
    t_final_io = omp_get_wtime()
    tempo_total_io_s = t_final_io - t_inicio_io
    
    t_final_total = omp_get_wtime()
    tempo_wall_clock = t_final_total - t_inicio_total

    !  Escrevemos o bloco de performance direto no final do arquivo
    write(26, *) ""
    write(26, *) "===== RESUMO DE PERFORMANCE (SEQUENCIAL) ====="
    write(26, '(A, F12.6, A)') "Tempo de Pre-processamento: ", tempo_total_sim_s, " s"
    write(26, '(A, F12.6, A)') "Tempo das Condicoes (BC)  : ", tempo_total_bc_s, " s"
    write(26, '(A, F12.6, A)') "Tempo de Forcas/Dano      : ", tempo_total_forca_s, " s"
    write(26, '(A, F12.6, A)') "Tempo de Escrita (I/O)    : ", tempo_total_io_s, " s"
    write(26, '(A, F12.6, A)') "Tempo Total do Programa   : ", tempo_wall_clock, " s"
    write(26, *) "=============================================="

    !Fechamos o arquivo 
    close(26)

    print *, "O relatorio de performance foi salvo no final de 'familia_resultados.txt'."

111 format(e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5)

end program main