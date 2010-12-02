! program tmp
! 
! 	implicit none
! 	integer 			:: it, ip0, ip1, i, info
! 	double precision	:: dR0(292, 4), dR1(292,5), dR2(292,5), dvectors(5,4), dvalues(4), &
! 							tau(5, 3), rho(3, 1), psi(5, 1), kappa(3, 3), beta(5, 1), dHmat(5,4), &
! 							dX(292, 4)
! 	
! ! 	open(unit=11, file='/Users/andreasnoackjensen/R0.txt', recl=292*4)
! ! 	do i = 1, 292
! ! 		read(unit=11, fmt=*) dR0(i,:)		
! ! 	end do
! ! 	close(11)
! ! 	open(unit=11, file='/Users/andreasnoackjensen/R1.txt')
! ! 	do i = 1, 292
! ! 		read(unit=11, fmt=*) dR1(i,:)
! ! 	end do
! ! 	close(11)
! ! 	open(unit=11, file='/Users/andreasnoackjensen/R2.txt')
! ! 	do i = 1, 292
! ! 		read(unit=11, fmt=*) dR2(i,:)
! ! 	end do
! ! 	close(11)
! 	open(unit=11, file='/Users/andreasnoackjensen/X.txt')
! 	do i = 1, 292
! 		read(unit=11, fmt=*) dX(i,:)
! 	end do
! 	close(11)
! 	it = 292
! 	ip0 = 4
! 	ip1 = 5
! 
! !	call rrr(dR1(:,1:4), dR2, 292, 4, 5, dvalues, dvectors)
! !	tau	= reshape((/1,0,0,0,0,0,1,0,0,0,0,0,1,0,0/), shape(tau))
! !	dHmat	= reshape((/1,0,0,-1,0,0,1,0,-1,0,0,0,1,-1,0,0,0,0,0,1/), shape(dHmat))
! !	write(*,*) tau
! !	call tauswitch(dR0, dR1, dR2, dHmat, it, ip0, ip1, 1, 1, 4, tau, rho, psi, kappa, 10000, 1.0d-15, info)
! !	call mlm(dR2, dR0, 292, 5, 4, beta)
! !	call dgemm('N', 'N', 5, 1, 3, 1.0d+0, tau, 5, rho, 3, 0.0d+0, beta, 5)
! 	dR0 = dX
! 	call fraclag(dR0, (/292, 4/), 0.4d+0)
! ! 	dX = dX - dR0
! ! 	do i = 1, 100
! ! 		write(*,*) dX(i,:)		
! ! 	end do
! !	write(*,*) info
! 
! end program tmp

subroutine rrr(dR0, dR1, it, ip0, ip1, dvalues, dvectors)

	implicit none

	integer, intent(in)				:: it, ip0, ip1
	double precision, intent(in) 	:: dR0(it, ip0), dR1(it, ip1)
	double precision, intent(out)	:: dvectors(ip1, ip0), dvalues(ip0)
	double precision				:: tau0(ip0), tau1(ip1), work(10*(2*ip1+ (ip1 + 1)*it)), qrR1(ip1, ip1), &
										dZ(ip1, ip0), vt(ip0, ip1), dR0cp(it, ip0), dR1cp(it, ip1)
	integer							:: jpvt0(ip0), jpvt1(ip1), lwork, info, iwork(8*ip0), i, j
	
	lwork	= 10*(2*ip1+ (ip1 + 1)*it)
	dR0cp = dR0
	dR1cp = dR1

	call dgeqp3(it, ip0, dR0cp, it, jpvt0, tau0, work, lwork, info)
	call dormqr('L', 'T', it, ip1, ip0, dR0cp, it, tau0, dR1cp, it, work, lwork, info)
	call dgeqp3(it, ip1, dR1cp, it, jpvt1, tau1, work, lwork, info)
	qrR1 = dR1cp(1:ip1,:)
	call dorgqr(it, ip1, ip1, dR1cp, it, tau1, work, lwork, info)
	dZ	= transpose(dR1cp(1:ip0,:))
	call dgesdd('S', ip1, ip0, dZ, ip1, dvalues, dvectors, ip1, vt, ip1, work, lwork, iwork, info)
	
	dvalues		= dvalues**2
	
	call dtrsm('L', 'U', 'N', 'N', ip1, ip0, sqrt(dble(it)), qrR1, ip1, dvectors, ip1)
	do i = 1, ip1
		dZ(jpvt1(i),:) = dvectors(i,:)
	end do
	dvectors = dZ
	
	return
		
end subroutine rrr

subroutine mlm(dX, dY, it, ipx, ipy, dbeta)
	
	implicit none

	integer, intent(in) 			:: it, ipx, ipy
	double precision, intent(inout)	:: dY(it,ipy), dbeta(ipx, ipy)
	double precision, intent(in)	:: dX(ipx,ipx)
	double precision				:: dXX(ipx, ipx)
	integer							:: info
	
	dbeta = 0.0d+0
	call dgemm('T', 'N', ipx, ipy, it, 1.0d+0, dX, it, dY, it, 0.0d+0, dbeta, ipx)
	dXX = 0.0d+0
	call dsyrk('U', 'T', ipx, it, 1.0d+0, dX, it, 0.0d+0, dXX, ipx)
	call dposv('U', ipx, ipy, dXX, ipx, dbeta, ipx, info)
	dy = 0.0d+0
	call dgemm('N', 'N', it, ipy, ipx, 1.0d+0, dX, it, dbeta, ipx, 0.0d+0, dY, it)

	return
	
end subroutine mlm

subroutine tauswitch(dR0, dR1, dR2, dH, it, ip0, ip1, ir, is2, im, dtau, drho, dpsi, dkappa, &
	imaxiter, dtol, info)

	implicit none

	integer, intent(in)				:: ip0, ip1, ir, is2, im, imaxiter
	integer, intent(inout)			:: it
	integer, intent(out)			:: info
	double precision, intent(in)	:: dR0(it, ip0), dR1(it, ip1), dR2(it,ip1), dH(ip1, im), dtol
	double precision, intent(inout)	:: dtau(ip1, ip0 - is2)
	double precision, intent(out)	:: drho(ip0 - is2, ir), dkappa(ip0 - is2, ip0 - ir), dpsi(ip1,ir)
	double precision				:: dtmp11(it**2), dtmp12(it**2), dtmp21(it, it), dtmp22(it, it), &
		dtmp23(it, it), dA(ip0, ip0), dB(ip0, ip0), dR1tau(it, ip0 - is2), dR2tauR1(it, ip0 + ip1 - is2), &
		dalph(ip0, ir), dalphort(ip0, ip0 - ir), dalphbar(ip0, ir), &
		dOmega1(ir, ir), dOmega2(ip0-ir,ip0-ir), dmA(ip0 - is2, ip0 - is2), dmB(im, im), &
		dmC(ip0 - is2, ip0 - is2), dmD(im, im), dmE(ip0 - is2, ip1), &
		reldif, dtinv, dresR2R0(it, ip1), dresR0R0(it, ip0), dll
	integer							:: is1, iwork(it**2), i, j, k
	
	reldif	= 1.0d+0
	dll		= -1.0d+100
	
	is1		= ip0 - ir - is2
	dtinv	= 1 / dble(it)

	do i = 1, imaxiter
		dR1tau		= 0.0d+0
		dR2tauR1	= 0.0d+0
		call dgemm('N', 'N', it, ir + is1, ip1, 1.0d+0, dR1, it, dtau, ip1, 0.0d+0, dR1tau, it)
		call dgemm('N', 'N', it, ir + is1, ip1, 1.0d+0, dR2, it, dtau, ip1, 0.0d+0, dR2tauR1, it)
		dR2tauR1(:,(ir + is1 + 1):(ir + is1 + ip1)) = dR1

		dtmp21(:,1:ip0) = dR0
		call mlm(dR1tau, dtmp21, it, ir + is1, ip0, dtmp22(1:(ir + is1),1:ip0))
		dB = 0.0d+0
		call dsyrk('U', 'T', ip0, it, dtinv, dR0 - dtmp21(1:it,1:ip0), it, 0.0d+0, dB, ip0)

		dtmp21(:,1:ip0) = dR0		
		call mlm(dR2tauR1, dtmp21(:, 1:ip0), it, ir + is1 + ip1, ip0, dtmp22(1:(ir + is1 + ip1),1:ip0))
		dA = 0.0d+0
		call dsyrk('U', 'T', ip0, it, dtinv, dR0 - dtmp21(1:it,1:ip0), it, 0.0d+0, dA, ip0)
		
		dtmp23(1:ip0,1:ip0) = dB
		call dsygvd(1, 'V', 'U', ip0, dA, ip0, dB, ip0, dtmp11, dtmp12, it**2, iwork, it**2, info)
		
		do j = 1, ip0
			if ( j <= ip0 - ir ) then
				dalphort(:,j) 			= dA(:,ip0 - j + 1)
			else
				dtmp21(1:ip0,j-ip0+ir)	= dA(:,ip0 - j + 1)
			endif
		end do
		
		call dsymm('L', 'U', ip0, ir, 1.0d+0, dtmp23(1:ip0,1:ip0), ip0, dtmp21(1:ip0, 1:ir), ip0, 0.0d+0, &
			dalph, ip0)
		
		dalphbar = dalph
		call mbar(dalphbar, ip0, ir)

! Step Two		
		call dgemm('N', 'N', it, ip0 - ir, ip0, 1.0d+0, dR0, it, dalphort, ip0, 0.0d+0, dtmp21(:, 1:(ip0 - ir)), it)
		dtmp21(1:it,(ip0 - ir + 1):(ip0 + ip1 + is1)) = dR2tauR1
		call dgemm('N', 'N', it, ir, ip0, 1.0d+0, dR0, it, dalphbar, ip0, 0.0d+0, dtmp22(:, 1:ir), it)
		
! Caldulation of A
		dtmp22(1:it,(ir + 1):(2*ir)) = dtmp22(1:it,1:ir)
		call mlm(dtmp21(1:it,1:(ip0 + ip1 + is1)), dtmp22(1:it,1:ir), it, ip0 + ip1 + is1, ir, dtmp23(1:(ip0 + ip1 + is1),1:ir))
		drho	= dtmp23((ip0 - ir + 1):(ip0 + is1), 1:ir)
		dpsi	= dtmp23((ip0 + is1 + 1):(ip1 + ip1 + is1), 1:ir)
		dtmp21(:,1:ir)	= dtmp22(:,(ir+1):(2*ir)) - dtmp22(:,1:ir)
	! The next is for use in matrix E
		call dsyrk('U', 'T', ir, it, 1.0d+0, dtmp21(:, 1:ir), it, 0.0d+0, dOmega1, ir)
	! Back to A
		call dgeqrf(it, ir, dtmp21(:,1:ir), it, dtmp11(1:ir), dtmp12, it**2, info)
		call dtrsm('R', 'U', 'N', 'N', ir + is1, ir, 1.0d+0, dtmp21(1:ir, 1:ir), ir, &
				dtmp23((ip0 - ir + 1):(ip0 + is1),1:ir), ir + is1)
		call dsyrk('U', 'N', ir + is1, ir, 1.0d+0, dtmp23((ip0 - ir + 1):(ip0 + is1),1:ir), ir + is1, &
				0.0d+0, dmA, ir + is1)

! Calculation of B
		call dgemm('N', 'N', it, ip0 - ir, ip0, 1.0d+0, dR0, it, dalphort, ip0, 0.0d+0, dtmp21(:, 1:(ip0 - ir)), it)
		dtmp21(:,(ip0 - ir + 1):(ip0 - ir + ip1)) 			= dR1
		dtmp21(:,(ip0 - ir + ip1 + 1):(ip0 - ir + 2*ip1))	= dR2
		call mlm(dtmp21(:,1:(ip0 - ir + ip1)), dtmp21(:, (ip0 - ir + ip1 + 1):(ip0 - ir + 2*ip1)), it, &
				ip0 - ir + ip1, ip1, dtmp22(1:(ip0 - ir + ip1),1:ip1))
		dresR2R0	= dR2 - dtmp21(:, (ip0 - ir + ip1 + 1):(ip0 - ir + 2*ip1))
	! The next if for matrix E
		dtmp21(:,(ip0 - ir + ip1 + 1):(2*ip0 - ir + ip1))	= dR0
		call mlm(dtmp21(:,1:(ip0 - ir + ip1)), dtmp21(:, (ip0 - ir + ip1 + 1):(2*ip0 - ir + ip1)), it, &
				ip0 - ir + ip1, ip0, dtmp22(1:(ip0 - ir + ip1),1:ip0))
		dresR0R0	= dR0 - dtmp21(:, (ip0 - ir + ip1 + 1):(2*ip0 - ir + ip1))
	! Back to B
		call dgemm('N', 'N', it, im, ip1, 1.0d+0, dresR2R0, it, dH, ip1, 0.0d+0, dtmp21(:, 1:im), it)
		call dsyrk('U', 'T', im, it, 1.0d+0, dtmp21(:, 1:im), it, 0.0d+0, dmB, im)
		
! Calculation of C
		call dgemm('N', 'N', it, ip0 - ir, ip0, 1.0d+0, dR0, it, dalphort, ip0, 0.0d+0, dtmp21(:, 1:(ip0 - ir)), it)
		dtmp21(1:it, (ip0 - ir + 1):(2*(ip0 - ir))) = dtmp21(1:it, 1:(ip0 - ir))
		call mlm(dR1tau, dtmp21(:,1:(ip0 - ir)), it, ir + is1, ip0 - ir, dtmp22(1:(ir + is1),1:(ip0 - ir)))
		dkappa	= dtmp22(1:(ir + is1),1:(ip0 - ir))
		dtmp21(:, 1:(ip0 - ir)) = dtmp21(1:it, (ip0 - ir + 1):(2*(ip0 - ir))) - dtmp21(:, 1:(ip0 - ir))
	! The next is for use in matrix E
		call dsyrk('U', 'T', ip0 - ir, it, 1.0d+0, dtmp21(:, 1:(ip0 - ir)), it, &
				0.0d+0, dOmega2, ip0 - ir)
	! Back to C
		call dgeqrf(it, ip0 - ir, dtmp21(:,1:(ip0 - ir)), it, dtmp11(1:(ip0 - ir)), dtmp12, it**2, info)
		call dtrsm('R', 'U', 'N', 'N', ir + is1, ip0 - ir, 1.0d+0, dtmp21(1:(ip0 - ir), 1:(ip0 - ir)), ip0 - ir, &
				dtmp22(1:(ir + is1),1:(ip0 - ir)), ir + is1)
		call dsyrk('U', 'N', ir + is1, ip0 - ir, 1.0d+0, dtmp22(1:(ir + is1),1:(ip0 - ir)), ir + is1, &
				0.0d+0, dmC, ir + is1)

! Calculation of D
		call dgemm('N', 'N', it, im, ip1, 1.0d+0, dR1, it, dH, ip1, 0.0d+0, dtmp21(:, 1:im), it)		
		call dsyrk('U', 'T', im, it, 1.0d+0, dtmp21(:,1:im), it, 0.0d+0, dmD, im)
! Calculation of E
	! First part
		dtmp21(1:ir, 1:ip0) = transpose(dalphbar)
		call dposv('U', ir, ip0, dOmega1, ir, dtmp21(1:ir,1:ip0), ir, info)
		call dgemm('T', 'N', ip0, ip1, it, 1.0d+0, dresR0R0, it, dresR2R0, it, 0.0d+0, dtmp22(1:ip0, 1:ip1), ip0)
		call dgemm('N', 'N', ir, ip1, ip0, 1.0d+0, dtmp21(1:ir,1:ip0), ir, dtmp22(1:ip0,1:ip1), ip0, &
				0.0d+0, dtmp23(1:ir, 1:ip1), ir)
		call dgemm('N', 'N', ir + is1, ip1, ir, 1.0d+0, drho, ir + is1, dtmp23(1:ir,1:ip1), ir, &
				0.0d+0, dtmp21(1:(ir + is1), 1:ip1), ir+is1)
		dmE = dtmp21(1:(ir+is1),1:ip1)

	! Second part
		dtmp21(1:(ip0 - ir), 1:ip0) = transpose(dalphort)
		call dposv('U', ip0 - ir, ip0, dOmega2, ip0 - ir, dtmp21(1:(ip0 - ir),1:ip0), ip0 - ir, info)
		call dgemm('T', 'N', ip0, ip1, it, 1.0d+0, dR0, it, dR1, it, 0.0d+0, dtmp22(1:ip0, 1:ip1), ip0)
		call dgemm('N', 'N', ip0 - ir, ip1, ip0, 1.0d+0, dtmp21(1:(ip0 - ir),1:ip0), ip0 - ir, &
				dtmp22(1:ip0,1:ip1), ip0, 0.0d+0, dtmp23(1:(ip0 - ir), 1:ip1), ip0 - ir)
		call dgemm('N', 'N', ir + is1, ip1, ip0 - ir, 1.0d+0, dkappa, ir + is1, &
				dtmp23(1:(ip0 - ir),1:ip1), ip0 - ir, 0.0d+0, dtmp21(1:(ir + is1), 1:ip1), ir+is1)
		dtmp22(1:(ir + is1), 1:ip1) = dmE + dtmp21(1:(ir+is1),1:ip1)
		call dgemm('N', 'N', ir + is1, im, ip1, 1.0d+0, dtmp22(1:(ir + is1), 1:ip1), ir + is1, &
				dH, ip1, 0.0d+0, dmE(:,1:im), ir + is1)

		info = 1

 		dtmp21(1:(ir + is1), 1:(ir + is1)) = dmA + dmC
		call dsygvd(1, 'V', 'U', ir + is1, dmC, ir + is1, dtmp21(1:(ir + is1), 1:(ir + is1)), ir + is1, &
 				dtmp11(1:(ir + is1)), dtmp12, it**2, iwork, it**2, info)
		dtmp21(1:im, 1:im) = dmB + dmD
		call dsygvd(1, 'V', 'U', im, dmD, im, dtmp21(1:im, 1:im), im, &
 				dtmp11((ir + is1 + 1):(ir + is1 + im)), dtmp12, it**2, iwork, it**2, info)
		call dgemm('T', 'N', ir + is1, im, ir + is1, 1.0d+0, dmC, ir + is1, dmE(:,1:im), ir + is1, 0.0d+0, &
 				dtmp21(1:(ir + is1), 1:im), ir + is1)
 		call dgemm('N', 'N', ir + is1, im, im, 1.0d+0, dtmp21(1:(ir + is1), 1:im), ir + is1, &
 				dmD, im, 0.0d+0, dtmp22(1:(ir + is1), 1:im), ir + is1)
 		do j = 1, ir + is1
 			do k = 1, im
 				dtmp21(j, k) = dtmp22(j, k) / ((1.0d+0 - dtmp11(j)) * (1.0d+0 - dtmp11(ir + is1 + k)) + &
 					dtmp11(j) * dtmp11(ir + is1 + k))
 			end do
 		end do
		
 		call dgemm('N', 'N', ir + is1, im, ir + is1, 1.0d+0, dmC, ir + is1, &
 			dtmp21(1:(ir + is1), 1:im), ir + is1, 0.0d+0, dtmp22(1:(ir + is1), 1:im), ir + is1)
 		call dgemm('N', 'T', ir + is1, im, im, 1.0d+0, dtmp22(1:(ir + is1), 1:im), ir + is1, &
 			dmD, im, 0.0d+0, dtmp21(1:(ir + is1), 1:im), ir + is1)
! 		write(*,*) dtmp21(1:(ir + is1), 1:ip1)
		
! 		dtmp22(1:im, 1:(ir + is1)) = transpose(dtmp21(1:(ir + is1), 1:im))
		
		call dgemm('N', 'T', ip1, ir + is1, im, 1.0d+0, dH, ip1, dtmp21(1:(ir + is1), 1:im), ir + is1, &
				0.0d+0, dtau, ip1)

	! Calculation fo the likelihood
		call dgemm('N', 'N', it, ir + is1, ip1, 1.0d+0, dR2, it, dtau, ip1, 0.0d+0, &
			dtmp21(:,1:(ir + is1)), it)
		call dgemm('N', 'N', it, ir, ir + is1, 1.0d+0, dtmp21(:,1:(ir + is1)), it, drho, ir + is1, 0.0d+0, &
			dtmp22(:, 1:ir), it)
		call dgemm('N', 'N', it, ir, ip1, 1.0d+0, dR1, it, dpsi, ip1, 0.0d+0, dtmp22(:,(ir + 1):(2*ir)), it)
		dtmp21(:,1:ir) = dtmp22(:,1:ir) + dtmp22(:,(ir + 1):(2*ir))
		call dgemm('N', 'N', it, ir + is1, ip1, 1.0d+0, dR1, it, dtau, ip1, 0.0d+0, dtmp21(:,(ir + 1):(2*ir + is1)), it)
		dtmp22(:,1:ip0) = dR0
		call mlm(dtmp21(:, 1:(2*ir + is1)), dtmp22(:,1:ip0), it, 2*ir + is1, ip0, dtmp23(1:(2*ir + is1),1:ip0))
		dtmp21(:,1:ip0) = (dR0 - dtmp22(:,1:ip0)) * sqrt(dtinv)
		call dgeqp3(it, ip0, dtmp21(:,1:ip0), it, iwork(1:ip0), dtmp11(1:ip0), dtmp12, it**2, info)
		dtmp11(1) = 1
		do j = 1, ip0
			dtmp11(1) = dtmp11(1) * dtmp21(j,j)
		end do
		dtmp11(1) = -0.5d+0 * it * log(dtmp11(1)**2)

 		reldif = abs((dll - dtmp11(1)) / dll)
		dll = dtmp11(1)

		if (reldif < dtol) then
			it = i
			exit
		endif
	end do
	if ( i == imaxiter ) it = 0
	return
	
end subroutine tauswitch

subroutine ortcomp(dA, dAort, m, n)

	implicit none

	integer, intent(in)				:: m, n
	double precision, intent(in) 	:: dA(m, n)
	double precision, intent(out)	:: dAort(m, m-n)
	integer							:: i, jpvt(n), lwork, info
	double precision				:: tau(n), dQ(m, m), work(4*(2*n+(n+1)*10))
	
	lwork = 4*(2*n+(n+1)*10)
	
	do i = 1, n
		jpvt(i) = i
	end do

	call dgeqp3(m,n,dA,m,jpvt,tau,work,lwork,info)
	dQ = 0.0d+0
	dQ(:,1:n) = dA
!	tau((n+1):m) = 0.0d+0
	call dorgqr(m,m,n,dQ,m,tau,work,lwork,info)
!	write(*,*) info
	
	dAort = dQ(:,(n+1):m)
	return
	
end subroutine ortcomp

subroutine mbar(dA, m, n)
	integer, intent(in)				:: m, n
	double precision, intent(inout) :: dA(m, n)
	double precision				:: dAtmp(m, n), tau(n), dtmp(m**2*n**2)
	integer							:: itmp(m**2*n**2), info

	call dgeqp3(m, n, dA, m, itmp(1:n), tau, dtmp, m**2*n**2, info)
	dAtmp	= dA
	call dorgqr(m, n, n, dA, m, tau, dtmp, m**2*n**2, info)
	call dtrsm('R', 'U', 'C', 'N', m, n, 1.0d+0, dAtmp(1:n,:), n, dA, m)
	
	return
	
end subroutine mbar