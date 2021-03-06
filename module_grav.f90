module grav
use param
!use utilmod
!use coordsys
!integer,parameter::grv_order=100
real(dp),allocatable::c(:,:),s(:,:)
real(dp)::gm ,a
character*150::grv_file
integer:: grv_order,grv_file_id,grv_body_id
private 
public getgrvp,getgrvp1,gravity,gforcebf,grv_order,grv_file,grv_file_id,grv_body_id	   	
contains
!public
!getgrvp		read grv file
!getgrvp1
!gravity		compute gravity force in J2000 frame
!gforcebf		compute gravity force in body-fixed frame
!grv_order		order of gravity
!grv_file		name of grv file
!grv_file_id	file number of 
!grv_body_id	central body id

subroutine gravity(et,r,f)                                
! et:time wrt. spacecraft                                
!  r:inertial coord m                                   
!  f:force wrt. inertial      m/s^2                     
implicit none                                            
real(dp):: et,f(3),accbf(3),r(3),xbf(3),mf2i(3,3),mi2f(3,3)           
select case (grv_body_id)
	case(499)                                            
    	call f2i_mars(et,mf2i)    !body fixed to J2000
	case(699)
  	 	call pXFORM ( 'IAU_SATURN', 'J2000', ET, Mf2i ) 
case default
		write(*,*)'wrong body id'
		stop
endselect                     
!    call pXFORM ( 'IAU_MARS', 'J2000', ET, Mf2i )       
    call INVert ( mf2i, mi2f )                           
    call MXV ( Mi2f, r, xbf )                            
!    xbf=1000.d0*xbf                                     
   call gforcebf(xbf,accbf)                             
!   write(*,*)accbf                                      
    call MXV ( Mf2i, accbf, f )                          
endsubroutine gravity                                    


subroutine gforcebf(xbf,accbf)
implicit none
real(dp)::xbf(3),k,x,y,z,r,rx,ry,rz,px,py,pz,lx,ly,lz,x2y,&
	&ur,up,ul,fnm,cosp,sinp,normp(0:grv_order,0:grv_order),cltemp,sltemp,dpnm(0:grv_order,0:grv_order),&
	&cl(0:grv_order),sl(0:grv_order),q,h,kk,g,accbf(3)
integer:: n,m,i

x=xbf(1)
y=xbf(2)
z=xbf(3)
r=dsqrt(x*x+y*y+z*z)
x2y=x*x+y*y
rx=x/r
ry=y/r
rz=z/r
px=-x*z/(dsqrt(x2y)*r**2)
py=-y*z/(dsqrt(x2y)*r**2)
pz=dsqrt(x2y)/r**2
lx=-y/x2y
ly=x/x2y
lz=0.d0
sinp=z/r
cosp=dsqrt(x2y)/r
cltemp=x/dsqrt(x2y)
sltemp=y/dsqrt(x2y)
!initial
cl=0.d0
sl=0.d0
normp=0.d0
dpnm=0.d0
ur=0.d0
up=0.d0
ul=0.d0
accbf=0.d0

cl(0)=1.d0
sl(0)=0.d0
do i=1,grv_order
!	call legendre_associated ( grv_order, i, sinp,normp(:,i) )
	if (i.ne.0) then
		cl(i)=cl(i-1)*cltemp-sl(i-1)*sltemp
		sl(i)=sl(i-1)*cltemp+cl(i-1)*sltemp
	endif
enddo
!write(*,'(7f10.3)')(normp(m,:),m=0,grv_order)
!write(*,'(2f)')(cl(i),sl(i),i=1,grv_order)
!stop
dpnm(0,0)=0.d0
dpnm(1,0)= dsqrt(3.d0)*cosp
dpnm(1,1)=-dsqrt(3.d0)*sinp
normp(0,0)=1.d0
normp(1,0)= dsqrt(3.d0)*sinp
normp(1,1)= dsqrt(3.d0)*cosp
!write(*,*)cosp,sinp

do n=2,grv_order
	q=dsqrt((2.d0*n+1.d0)/(2.d0*n))
	g=dsqrt(2.d0*n+1.d0)
	do m=0,n
		h=dsqrt((2.d0*n-1.d0)*(2.d0*n+1.d0)/((n+m)*(n-m)))
		kk=dsqrt((n+m-1.0)*(n-m-1.0)*(2.d0*n+1.d0)/((n+m)*(n-m)*(2.d0*n-3.d0)))
		if(m.eq.n) then
			dpnm(n,m)=-q*sinp*normp(n-1,n-1)+q*cosp*dpnm(n-1,n-1)
			normp(n,m)=normp(n-1,n-1)*q*cosp
		elseif(m.eq.n-1) then
			 dpnm(n,m)=g*cosp*normp(n-1,n-1)+g*sinp*dpnm(n-1,n-1)
			normp(n,m)=normp(n-1,n-1)*g*sinp
		else
			dpnm(n,m)=h*cosp*normp(n-1,m)+h*sinp*dpnm(n-1,m)-kk*dpnm(n-2,m)
			normp(n,m)=normp(n-1,m)*h*sinp-normp(n-2,m)*kk
		endif

		ur=-(n+1.d0)*(gm/r**2)*(a/r)**n*normp(n,m)*(c(n,m)*cl(m)+s(n,m)*sl(m))
		up= (gm/r)*(a/r)**n*dpnm(n,m)*(c(n,m)*cl(m)+s(n,m)*sl(m))
		ul=(gm/r)*m*(a/r)**n*normp(n,m)*(-c(n,m)*sl(m)+s(n,m)*cl(m))
		accbf(1)=accbf(1)+ur*rx+up*px+ul*lx
		accbf(2)=accbf(2)+ur*ry+up*py+ul*ly
		accbf(3)=accbf(3)+ur*rz+up*pz+ul*lz
	enddo
enddo
endsubroutine gforcebf

!gravity file of jpl
subroutine getgrvp
implicit none
integer::n,m,rn,rm
allocate(c(2:grv_order,0:grv_order),s(2:grv_order,0:grv_order))
OPEN(grv_file_id,FILE=grv_file,STATUS='old')
read(grv_file_id,*)a,gm
a=a*1.d3
gm=gm*1.d9
do n=2,grv_order
	do m=0,n
	read(grv_file_id,'(I5,1x,I5,1x,d23.16,1x,d23.16)',end=112)rn,rm,c(n,m),s(n,m)
	if (rn.ne.n.or.rm.ne.m) then
		write(*,*)'read grv error!'
		stop
	endif	
	enddo
enddo
112 continue
close(grv_file_id)
endsubroutine getgrvp

!gravity file of gsfc
subroutine getgrvp1
implicit none
integer::n,m,rn,rm
allocate(c(2:grv_order,0:grv_order),s(2:grv_order,0:grv_order))
OPEN(22,FILE='/home/jnc/work/spice/src/trsim/gravtest/geok',STATUS='old')
read(22,*)ctemp
read(22,'(A20,2E20.10)')ctemp,gm,a
read(22,*)ctemp
write(*,*)gm,a
!a=a*1.d3
!gm=gm*1.d9
do n=2,grv_order
    do m=0,n
    read(22,'(A6,2I3,2D21.14)',end=112)ctemp,rn,rm,c(n,m),s(n,m)
	if (rn.ne.n.or.rm.ne.m) then
        write(*,*)'read grv error!'
        stop
    endif
!   write(*,*)c(n,m),s(n,m)
    enddo
enddo
112 continue
close(22)
endsubroutine getgrvp1



subroutine f2i_mars(et,m)       
!fixed system to J2000 inertial system matrix                     
real(dp)::et,m(3,3),spd,rpd,day,N,J,psi,I,PHI
M=0.d0
M(1,1)=1.d0       
M(2,2)=1.d0                                           
M(3,3)=1.d0                                           
day=et/spd()                                          
N = 3.37919183d0 !deg                                 
J = 24.67682669d0 !deg                                
psi=81.9684072598d0+day*(-0.0000057901d0)             
I=25.1893863984d0+day*(-0.000000003d0)                
phi=133.38462d0+day*350.891985305d0                   
!call ROTMAT ( M, -N*rpd(), 3, Mtemp1 ) !R_z(-N)      
!call ROTMAT ( Mtemp1, -J*rpd(), 1, Mtemp2 )          
!call ROTMAT ( Mtemp2, -psi*rpd(), 3, Mtemp1 )        
!call ROTMAT ( Mtemp1, -I*rpd(), 1, Mtemp2 )          
!call ROTMAT ( Mtemp2, -phi*rpd(), 3, M )             
call ROTMAT ( M, -phi*rpd(), 3, Mtemp1 ) !R_z(-N)     
call ROTMAT ( Mtemp1, -I*rpd(), 1, Mtemp2 )           
call ROTMAT ( Mtemp2, -psi*rpd(), 3, Mtemp1 )         
call ROTMAT ( Mtemp1, -J*rpd(), 1, Mtemp2 )           
call ROTMAT ( Mtemp2, -N*rpd(), 3, M )                
end subroutine f2i_mars   

endmodule grav
