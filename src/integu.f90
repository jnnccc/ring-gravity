subroutine integf(xbf,rthrea,f)
!完成单个环的引力积分
!xbf    土星固定系坐标
!rthrea 环的内径和外径(m)
!f		固定系坐标下的受力

!use para
implicit none
real*8,parameter:: pi=dacos(-1.d0)
real*8:: xbf(3),f(3),rthrea(2)
REAL*8:: X(3,2),Y(3,2),DATA(4500),RES,ERR,theta1,theta2
INTEGER IWORK(1000),NU,ND,NEVALS,IFLAG,ID(3),i,j,k
EXTERNAL Fx,fy,fz
!id(1)=i
!id(2)=j
!id(3)=k
!print*,id
theta1=0.d0
theta2=2.d0*pi
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
NU=0
ND=0
IFLAG=1

CALL TWODQ(intx,2,X,Y,1.E-08,1,200,32000,fx,ERR,NU,ND, NEVALS,IFLAG,DATA,IWORK,id)
CALL TWODQ(inty,2,X,Y,1.E-08,1,200,32000,fy,ERR,NU,ND, NEVALS,IFLAG,DATA,IWORK,id)
CALL TWODQ(intz,2,X,Y,1.E-08,1,200,32000,fz,ERR,NU,ND, NEVALS,IFLAG,DATA,IWORK,id)
!PRINT*,u,ERR,NEVALS,IFLAG
end subroutine integf

real*8 FUNCTION intx(t,theta,id)
use para
implicit none
real*8 t,theta
intx=t/dsqrt(cyl(3,id(1),id(2),id(3))**2+t**2+cyl(1,id(1),id(2),id(3))**2-2.d0*cyl(1,id(1),id(2),id(3))*t*dcos(cyl(2,id(1),id(2),id(3))-theta))
RETURN
END FUNCTION intx

real*8 FUNCTION inty(t,theta,id)
use para
implicit none
real*8 t,theta
inty=t/dsqrt(cyl(3,id(1),id(2),id(3))**2+t**2+cyl(1,id(1),id(2),id(3))**2-2.d0*cyl(1,id(1),id(2),id(3))*t*dcos(cyl(2,id(1),id(2),id(3))-theta))
RETURN
END FUNCTION inty

real*8 FUNCTION intz(t,theta,id)
use para
implicit none
real*8 t,theta
intz=t/dsqrt(cyl(3,id(1),id(2),id(3))**2+t**2+cyl(1,id(1),id(2),id(3))**2-2.d0*cyl(1,id(1),id(2),id(3))*t*dcos(cyl(2,id(1),id(2),id(3))-theta))
RETURN
END FUNCTION intz


