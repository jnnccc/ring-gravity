subroutine ringbf(xbf,f)
use param
implicit none
real(dp):: xbf(3),f(3),r(2),ftemp(3)
xbf=xbf/rsat
cyl(1)=dsqrt(xbf(1)*xbf(1)+xbf(2)*xbf(2))
cyl(2)=datan2(xbf(2),xbf(1))
cyl(3)=xbf(3)
cyl(4)=xbf(1)
cyl(5)=xbf(2)
cyl(6)=xbf(3)
!write(*,*)cyl

f=0.d0
!D ring
call ring_force(ringd/rsat,ftemp)
!write(*,*)ftemp
f=f+ftemp

!C ring
call ring_force(ringc/rsat,ftemp)
!write(*,*)ftemp
f=f+ftemp

!B ring
call ring_force(ringb/rsat,ftemp)
!write(*,*)ftemp
f=f+ftemp

!A ring
call ring_force(ringa/rsat,ftemp)
!write(*,*)ftemp
f=f+ftemp

end subroutine ringbf

subroutine ring_force(rthrea,f)
!完成单个环的引力积分
!xbf    土星固定系坐标
!rthrea 环的内径和外径(归一化)
!f		固定系坐标下的受力

use param
implicit none
!real*8,parameter:: pi=dacos(-1.d0)
integer,parameter::MAXTRI=500,MEVALS=60000
real*8:: f(3),rthrea(2)
REAL*8:: X(3,2),Y(3,2),DATA(10*MAXTRI),RES,ERR,theta1,theta2
INTEGER IWORK(3*MAXTRI),NU,ND,NEVALS,IFLAG,i,j,k
EXTERNAL intx,inty,intz
!id(1)=i
!id(2)=j
!id(3)=k
!print*,id
theta1=-pi
theta2=pi
X(1,1)=rthrea(1)
Y(1,1)=theta1
X(2,1)=rthrea(2)
Y(2,1)=theta1
X(3,1)=rthrea(2)
Y(3,1)=theta2
X(1,2)=rthrea(1)
Y(1,2)=theta1
X(2,2)=rthrea(2)
Y(2,2)=theta2
X(3,2)=rthrea(1)
Y(3,2)=theta2

!NU=0
!ND=0
!IFLAG=1
!CALL TWODQ(intx,2,X,Y,1.E-08,1,MAXTRI,MEVALS,f(1),ERR,NU,ND, NEVALS,IFLAG,DATA,IWORK)
!PRINT*,ERR,NEVALS,IFLAG
f(1)=0.d0
NU=0
ND=0
IFLAG=1
CALL TWODQ(inty,2,X,Y,1.E-08,1,MAXTRI,MEVALS,f(2),ERR,NU,ND, NEVALS,IFLAG,DATA,IWORK)
!PRINT*,ERR,NEVALS,IFLAG
NU=0
ND=0
IFLAG=1
CALL TWODQ(intz,2,X,Y,1.E-08,1,MAXTRI,MEVALS,f(3),ERR,NU,ND, NEVALS,IFLAG,DATA,IWORK)
!PRINT*,ERR,NEVALS,IFLAG
end subroutine ring_force

real*8 FUNCTION intx(t,theta)
use param
implicit none
real*8 t,theta,r(3),nr(3),dtemp
r(1)=t*dcos(theta)-cyl(4)
r(2)=t*dsin(theta)-cyl(5)
r(3)=-cyl(6)
call UNORM ( r, nr, dtemp )
!write(*,*)nr
intx=t*nr(1)/(cyl(3)**2+t**2+cyl(1)**2-2.d0*cyl(1)*t*dcos(cyl(2)-theta))
RETURN
END FUNCTION intx

real*8 FUNCTION inty(t,theta)
use param
implicit none
real*8 t,theta,r(3),nr(3),dtemp
r(1)=t*dcos(theta)-cyl(4)
r(2)=t*dsin(theta)-cyl(5)
r(3)=-cyl(6)
call UNORM ( r, nr, dtemp )
inty=t*nr(2)/(cyl(3)**2+t**2+cyl(1)**2-2.d0*cyl(1)*t*dcos(cyl(2)-theta))
RETURN
END FUNCTION inty

real*8 FUNCTION intz(t,theta)
use param
implicit none
real*8 t,theta,r(3),nr(3),dtemp
r(1)=t*dcos(theta)-cyl(4)
r(2)=t*dsin(theta)-cyl(5)
r(3)=-cyl(6)
call UNORM ( r, nr, dtemp )
intz=t*nr(3)/(cyl(3)**2+t**2+cyl(1)**2-2.d0*cyl(1)*t*dcos(cyl(2)-theta))
RETURN
END FUNCTION intz


