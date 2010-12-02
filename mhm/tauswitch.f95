!program tmp

!	implicit none
!	integer 			:: it, ip0, ip1
!	double precision	:: dR0(5,2), dR1(5,3), dvectors(3,3), dvalues(3)
	
!	dR0 = reshape((/1,3,5,7,9,2,4,6,8,10/), shape(dR0))
!	dR1 = reshape((/1,2,3,4,5,6,7,8,9,11,1,1,1,1,1/), shape(dR1))
!	it = 5
!	ip0 = 2
!	ip1 = 3
	
!	write(*,*) dR0, shape(dR0)
	
!	call rrr(dR0, dR1, it, ip0, ip1, dvalues, dvectors)

!	write(*,*) dvectors
!	write(*,*) dvalues
!end program tmp

subroutine rrr(dR0, dR1, it, ip0, ip1, dvalues, dvectors)

	implicit none

	integer, intent(in)				:: it, ip0, ip1
	double precision, intent(in) 	:: dR0(it, ip0), dR1(it, ip1)
	double precision, intent(out)	:: dvectors(ip1, ip1), dvalues(ip1)
	double precision				:: dS00(ip0, ip0), dS01(ip0, ip1), dS01c(ip0, ip1), dA(ip1, ip1), work(10*(1 + 6*ip1 + 2*ip1**2))
	integer							:: lwork, liwork, iwork(5*(3 + 5*ip1)), info, i
	
	lwork	= 10*(1 + 6*ip1 + 2*ip1**2)
	liwork	= 5*(3 + 5*ip1)

	dS00		= 0.0d+0
	dS01		= 0.0d+0
	dS01c		= 0.0d+0
	dvectors	= 0.0d+0
	dA			= 0.0d+0

	call dsyrk('U', 'T', ip0, it, 1.0d+0, dR0, it, 1.0d+0, dS00, ip0)
	call dsyrk('U', 'T', ip1, it, 1.0d+0, dR1, it, 1.0d+0, dvectors, ip1)
	call dgemm('T', 'N', ip0, ip1, it, 1.0d+0, dR0, it, dR1, it, 0.0d+0, dS01, ip0)
	dS01c = dS01
	call dposv('U', ip0, ip1, dS00, ip0, dS01c, ip0, info)
	call dgemm('T', 'N', ip1, ip1, ip0, 1.0d+0, dS01, ip0, dS01c, ip0, 0.0d+0, dA, ip1)
	call dsygvd(1, 'V', 'U', ip1, dA, ip1, dvectors, ip1, dvalues, work, lwork, iwork, liwork, info)

	do i = 1, ip1
		dvectors(:,ip1 - i + 1) = dA(:,i)
		work(i)	= dvalues(ip1 - i + 1)
	end do
	dvalues		= work(1:ip1)
	
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
	call dsyrk('U', 'T', ipx, it, 1.0d+0, dX, ipx, 0.0d+0, dXX, ipx)
	call dposv('U', ipx, ipy, dXX, ipx, dbeta, ipx, info)
	dy = 0.0d+0
	call dgemm('N', 'N', it, ipy, ipx, 1.0d+0, dX, it, dbeta, ipx, 0.0d+0, dY, it)

	return
	
end subroutine mlm

subroutine tauswitch(dR0, dR1, dR2, it, ip0, ip1, ir, is2, dbeta, dtau, ddelta)

	implicit none

	integer, intent(in)				:: it, ip0, ip1, ir, is2
	double precision, intent(in)	:: dR0(it, ip0), dR1(it, ip1), dR2(it,ip1)
	double precision, intent(out)	:: dbeta(ip1, ir), dtau(ip1, ip0-is2), ddelta(ip1,ir)
	double precision				:: dtmp11(it**2), dtmp12(it**2), dtmp21(it, it), dtmp22(it, it), dtmp23(it, it), &
		dA(ip0, ip0), dB(ip0, ip0), dSSR0R1tau(ip0, ip0), dR1tau(it, ip0 - is2), dR2tauR1(it, ip0 + ip1 - is2), &
		dalph(ip0, ir), dalphort(ip0, ip0 - ir), dalphbar(ip0, ir), drhot(ir,ip0 - is2), dOmega1(ir, ir), &
		dOmega2(ip0 - ir, ip0 - ir), dkappa(ip0 - is2, ip0 - ir), dmA(ip0 - is2, ip0 - is2), dmB(ip1, ip1), &
		dmC(ip0 - is2, ip0 - is2), dmD(ip1, ip1), dmE(ip0 - is2, ip1)
	integer							:: is1, iwork(it**2), i, info
	logical 						:: repeat
	
	is1	= ip0 - ir - is2
	
	call rrr(dR1(2:it,1:ip0), dR2(1:(it-1),:), it-1, ip0, ip1, dtmp11, dtmp21)
	dtau = dtmp21(1:ip1,1:ir)
	if ( is1 > 0 ) then
		call ortcomp(dtau(:,1:ir), dtmp21, ip1, ir)
		dtau(:, (ir + 1):(ir + is1)) = dtmp21(1:ip1,1:is1)
	end if
	
	repeat = .true.
	do while(repeat)
		dR1tau		= 0.0d+0
		dR2tauR1	= 0.0d+0
		call dgemm('N', 'N', it, ir + is1, ip1, 1.0d+0, dR1, it, dtau, ip1, 0.0d+0, dR1tau, it)
		call dgemm('N', 'N', it, ir + is1, ip1, 1.0d+0, dR2, it, dtau, ip1, 0.0d+0, dR2tauR1, it)
		dR2tauR1(:,(ir + is1 + 1):(ir + is1 + ip1)) = dR1

		dtmp21(:,1:ip0) = dR0
		call mlm(dR1tau, dtmp21, it, ir + is1, ip0, dtmp22(1:(ir + is1),1:ip0))
		dSSR0R1tau = 0.0d+0
		call dsyrk('U', 'T', ip0, it, 1.0d+0, dR0 - dtmp21(1:it,1:ip0), it, 0.0d+0, dSSR0R1tau, ip0)
		dB = dSSR0R1tau

		dtmp21(:,1:ip0) = dR0		
		call mlm(dR2tauR1, dtmp21, it, ir + is1 + ip1, ip0, dtmp22(1:(ir + is1),1:ip0))
		dA = 0.0d+0
		call dsyrk('U', 'T', ip0, it, 1.0d+0, dR0 - dtmp21(1:it,1:ip0), it, 0.0d+0, dA, ip0)
		
		call dsygvd(1, 'V', 'U', ip0, dA, ip0, dB, ip0, dtmp11, dtmp12, it**2, iwork, it**2, info)
		
		do i = 1, ip0 - ir
			dalphort(:,i) 	= dA(:,ip0 - i + 1)
			dtmp21(1:ip0,i)	= dA(:,ip0 - i + 1)
		end do
		
		dalph = 0.0d+0
		call dgemm('N', 'N', ip0, ir, ip0, 1.0d+0, dSSR0R1tau, ip0, dtmp21(1:ip0,(ip0 - ir + 1):ip0), ip0, &
			0.0d+0, dalph, ip0)
		
		dtmp23(1:ir,1:ip0) = transpose(dalph)
		call dsyrk('U', 'T', ir, ip0, 1.0d+0, dalph, ip0, 0.0d+0, dtmp21(1:ir,1:ir), ir)
		call dposv('U', ir, ip0, dtmp21(1:ir,1:ir), ir, dtmp23(1:ir,1:ip0), ir, info)
		dalphbar = transpose(dtmp23(1:ir,1:ip0))
		
		call dgemm('N', 'N', it, ip0 - ir, ip0, 1.0d+0, dR0, it, dalphort, ip0, 0.0d+0, dtmp21(1:it, 1:(ip0 - ir)), it)
		dtmp21(1:it,(ip0 - ir + 1):(ip0 + ip1 + is1)) = dR2tauR1
		call dgemm('N', 'N', it, ir, ip0, 1.0d+0, dR0, it, dalphbar, ip0, 0.0d+0, dtmp22(1:it, 1:ir), it)
		
		dtmp22(1:it,(ir + 1):(2*ir)) = dtmp22(1:it,1:ir)
		call mlm(dtmp21(1:it,1:(ip0 + ip1 + is1)), dtmp22(1:it,1:ir), it, ip0 + ip1 + is1, ir, dtmp23(1:(ip0 + ip1 + is1),1:ir))
		drhot	= transpose(dtmp23((ip0 - ir + 1):(ip0 + is1),1:ir))
		call dsyrk('U', 'T', ir, it, 1.0d+0 / it, dtmp22(1:it,(ir + 1):(2*ir)) - dtmp22(1:it,1:ir), it, 0.0d+0, dOmega1, ir)
		
		dtmp21(1:it,(ip0 - ir + 1):(2*(ip0 - ir))) = dtmp21(1:it,1:(ip0 - ir))
		call mlm(dR1tau, dtmp21(1:it,1:(ip0 - ir)), it, ir + is1, ip0 - ir, dtmp23(1:(ir + is1),1:(ip0 - ir)))
		dkappa	= dtmp23(1:(ir + is1),1:(ip0 - ir))
		call dsyrk('U', 'T', ip0 - ir, it, 1.0d+0 / it, dtmp21(1:it,(ip0 - ir + 1):(2*(ip0 - ir))) - dtmp21(1:it,1:(ip0 - ir)), &
			it, 0.0d+0, dOmega2, ip0 - ir)
		
!		Calculation of mA
		
		
	end do 
	
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
	write(*,*) info
	
	dAort = dQ(:,(n+1):m)
	return
	
end subroutine ortcomp