! program tmp
! 
! 	implicit none
! 	double precision	:: tmp1(5, 2), x, gamma
! 	
! 	tmp1 = reshape((/1,3,5,7,9,2,4,6,8,10/), shape(tmp1))
! 
! !	write(*,*) tmp1
! !	call fracdiff(tmp1, (/5,2/), 2.0d+0)
! !	write(*,*) tmp1
! 
! 	write(*,*) gamma(5.0d+0)
! 
! end program tmp

subroutine rmfilter(x, nl, coef, ndim, init)
!
!	Author: Andreas Noack Jensen
!	Date: April 2010
!
!	A subroutine for calculating the multivariate recursive filter of the series x. 
!	
	implicit none
	integer, intent(in) :: nl(2), ndim(3)
	double precision, intent(inout) :: x(nl(1),nl(2))
	double precision, intent(in) :: coef(ndim(1),ndim(2),ndim(3)), init(ndim(3),nl(2))
	double precision :: y(nl(1),nl(2))
	integer :: i, j
	
	y = 0.0d+0      		
	
	do i = 1, nl(1), 1
		do j = 1, ndim(3), 1
			if ( i <= j ) then
				call dgemm('N', 'T', 1, ndim(1), nl(2), 1.0d+0, init(ndim(3)-j+i,:), 1, &
				coef(:,:,j), ndim(1), 1.0d+0, y(i,:), 1)
	  		else
				call dgemm('N', 'T', 1, ndim(1), nl(2), 1.0d+0, y(i-j,:), 1,coef(:,:,j), &
					ndim(1), 1.0d+0, y(i,:), 1)
			end if
		end do
		y(i,:) = y(i,:) + x(i,:)
	end do
	x = y
	return
end subroutine rmfilter

subroutine fracdiff(dx, idim, dd)
!
!
!
	implicit none
	integer, intent(in) 			:: idim(2)
	double precision, intent(inout)	:: dx(idim(1), idim(2))
	double precision, intent(in)	:: dd
	double precision				:: dxtmp(idim(1), idim(2)), dtmp
	integer							:: i, j, k
	
	dxtmp = 0.0d+0
	do i = 2, idim(1)
		do j = 1, i - 1
			if ( abs(mod(dd, 1.0d+0)) < 1e-15 .AND. j > dd) then
				exit				
			else
				dtmp = 1.0d+0
				do k = 1, j
					dtmp		= dtmp * (dd - j + k) / k
				end do
				dxtmp(i,:)	= dxtmp(i,:) + (-1)**(j+1) * dtmp * dx(i - j,:)					
			end if
		end do
	end do
	
	dx(2:idim(1),:) = dx(2:idim(1),:) - dxtmp(2:idim(1),:)
	return
	
end subroutine fracdiff
