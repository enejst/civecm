subroutine rmfilter(x, nl, coef, ndim, init)

	implicit none
	integer, intent(in)				:: nl(2), ndim(3)
	double precision, intent(inout) :: x(nl(1),nl(2))
	double precision, intent(in)	:: coef(ndim(1),ndim(2),ndim(3)),init(ndim(3),ndim(nl(2)))
	double precision				:: y(nl(1),nl(2))
	integer							:: i, j

	y = 0.0d+0
	do i = 1, nl(1), 1
		do j = 1, ndim(3), 1
			if ( i <= j ) then
				call dgemm('N', 'T', 1, ndim(1), nl(2), 1.0d+0, init(ndim(3)-j+1,:), 1, coef(:,:,j), ndim(1), 1.0d+0, y(i,:), 1)
			else
				call dgemm('N', 'T', 1, ndim(1), nl(2), 1.0d+0, y(i-j,:), 1, coef(:,:,j), ndim(1), 1.0d+0, y(i,:), 1)
			end if
		end do
		y(i,:) = y(i,:) + x(i,:)
	end do
	x = y
	return
end subroutine rmfilter