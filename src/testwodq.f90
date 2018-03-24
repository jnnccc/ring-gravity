program twodintx
use twod_quad
REAL*8 X(3,2),Y(3,2),DATA(450),RES,ERR,s(3)
INTEGER IWORK(100),NU,ND,NEVALS,IFLAG
EXTERNAL F
X(1,1)=0.
Y(1,1)=0.
X(2,1)=1.
Y(2,1)=0.
X(3,1)=1.
Y(3,1)=1.
X(1,2)=0.
Y(1,2)=0.
X(2,2)=1.
Y(2,2)=1.
X(3,2)=0.
Y(3,2)=1.
NU=0
ND=0
IFLAG=1
CALL TWODQ(F,2,X,Y,1.E-08,1,50,4000,RES,ERR,NU,ND, NEVALS,IFLAG,DATA,IWORK)
PRINT*,RES,ERR,NEVALS,IFLAG
!CALL twodqd(f, 2, x, y, 1.e-04, 1, 50, 4000, result,  &
!                  error, nu, nd, nevals, iflag)
!PRINT*,RES,ERR,NEVALS,IFLAG
end program twodintx

real*8 FUNCTION F(X,Y)
real*8 x,y
F=dCOS(X+Y)
RETURN
END

