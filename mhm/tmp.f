        program tmp

        	implicit none
        	double precision :: x, dgamma

			write(*,*) x
			x = dgamma(5.0d+0)
			write(*,*) x

        end program tmp


C       subroutine mygamma(x)
C 
C       implicit none
C       double precision, intent(inout) :: x
C       double precision :: dgamma
C 
C       x = dgamma(x)
C 
C       end subroutine mygamma
