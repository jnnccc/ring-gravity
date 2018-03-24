module par
IMPLICIT NONE
real*8 r,z,phi
real*8,parameter::pi=dacos(-1.d0)
endmodule par

program twodintx
use  par
implicit none

REAL*8 X(3,2),Y(3,2),DATA(4500),RES,ERR,r1,r2,theta1,theta2
INTEGER IWORK(1000),NU,ND,NEVALS,IFLAG
EXTERNAL F
r=2.82842712474619d0
z=-2.d0
phi=-2.35619449019234d0

r1=1.52651490011283d0
r2=1.95095241255724d0
theta1=0.d0
theta2=2.d0*pi
X(1,1)=r1
Y(1,1)=theta1
X(2,1)=r2
Y(2,1)=theta1
X(3,1)=r2
Y(3,1)=theta2
X(1,2)=r1
Y(1,2)=theta1
X(2,2)=r2
Y(2,2)=theta2
X(3,2)=r1
Y(3,2)=theta2
NU=0
ND=0
IFLAG=1
CALL TWODQ(F,2,X,Y,1.E-08,1,200,32000,RES,ERR,NU,ND, NEVALS,IFLAG,DATA,IWORK)
PRINT*,RES,ERR,NEVALS,IFLAG

end program twodintx

real*8 FUNCTION F(t,theta)
use par
implicit none
real*8 t,theta
F=t/dsqrt(z**2+t**2+r**2-2.d0*r*t*dcos(phi-theta))
RETURN
END FUNCTION F

