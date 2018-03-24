subroutine readriv
use para
implicit none
OPEN(driv_id,FILE='input/driv',STATUS='old')
READ(driv_id,*)ndim,dim
allocate(grids(4,ndim,ndim,ndim),cyl(3,ndim,ndim,ndim))
READ(driv_id,*)r1,r2
r1=r1/r_saturn
r2=r2/r_saturn
close(driv_id)
endsubroutine readriv
