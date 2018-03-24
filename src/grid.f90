subroutine grid
use lib_array
use para
implicit none
integer*8 i,j,k
real*8:: x(ndim),y(ndim),z(ndim)
call linspace(-dim/2.d0,dim/2.d0,x)
call linspace(-dim/2.d0,dim/2.d0,y)
call linspace(-dim/2.d0,dim/2.d0,z)
!x=x+x0
!y=y+y0
!z=z+z0
!!$OMP PARALLEL DO default(shared) PRIVATE(I,j,k),SCHEDULE(dynamic,8)
!$OMP PARALLEL DO default(shared) PRIVATE(I,j,k),SCHEDULE(auto)
do i=1,ndim
	do j=1,ndim
		do k=1,ndim
			grids(1,i,j,k)=x(i)
			grids(2,i,j,k)=y(j)
			grids(3,i,j,k)=z(k)
			cyl(1,i,j,k)=dsqrt(y(j)**2d0+x(i)**2d0)
			cyl(2,i,j,k)=datan2(y(j),x(i))
			cyl(3,i,j,k)=z(k)
			call integu(grids(4,i,j,k),i,j,k)
!			write(*,'(4e)') grids(:,i,j,k)
		enddo
	enddo
enddo
!$OMP END PARALLEL DO
do i=1,ndim
    do j=1,ndim
        do k=1,ndim
            write(*,'(4e)') grids(:,i,j,k)
        enddo
    enddo
enddo

endsubroutine grid
