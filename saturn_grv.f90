program saturn_grv
use param
implicit none
CALL GETARG(1,arg1)

select case (arg1)
case('1')
	call test_sph
case('2')
	call test_orbit
case default
	stop
end select

end program saturn_grv

!真实轨道仿真
subroutine test_orbit

end subroutine test_orbit

!土星环的引力场与土星带谐项相关分析
subroutine test_sph
use param
use grav
use lib_array
implicit none

integer :: i,j,nr,np
!real(dp),parameter:: rs=60356000.d0 !环半径(m)
real(dp)::x(3),xbf(3),fs(3),ufs(3),fr(3),ufr(3),r1,r2,p1,p2,dtemp,dtemp1,vdot
real(dp),allocatable:: rho(:),phi(:)
!rho 	:径向
!phi	:纬向
CALL GETARG(2,arg2)
CALL GETARG(3,arg3)
CALL GETARG(4,arg4)
CALL GETARG(5,arg5)
CALL GETARG(6,arg6)
CALL GETARG(7,arg7)
CALL GETARG(8,arg8)

read(arg3,*)nr
read(arg4,*)np
read(arg5,*)r1
read(arg6,*)r2
read(arg7,*)p1
read(arg8,*)p2

r1=r1*rsat
r2=r2*rsat
p1=p1*pi/180.d0
p2=p2*pi/180.d0

allocate(rho(nr),phi(np))
!xbf(1)=-47174.707491d3
!xbf(2)=-145456.356920d3
!xbf(3)=152.790709d3
grv_file=trim(arg2)
grv_file_id=20
grv_body_id=699
grv_order=8
call getgrvp

!call linspace(66900000.d0,136775000.d0,rho)
call linspace(r1,r2,rho)
!call linspace(-pi/2,pi/2,phi)
call linspace(p1,p2,phi)
!rho=rs+rho

do i=1,nr
	do j=1,np
		xbf(1)=0.d0
		xbf(2)=rho(i)*dcos(phi(j))
		xbf(3)=rho(i)*dsin(phi(j))
!		write(*,*)rho(i)
		xbf(1)=0.d0
		xbf(2)=116775000.d0
		xbf(3)=1000.d0
		call gforcebf(xbf,fs)
		call UNORM ( fs, ufs, dtemp )
		call ringbf(xbf,fr)
		write(*,*)fr
		stop
		call UNORM ( fr, ufr, dtemp1 )
		write(*,'(5e)')rho(i)/rsat,phi(j)*180.d0/pi,dtemp,dtemp1,dabs(vdot(ufs,ufr))
	enddo
enddo

end subroutine test_sph


