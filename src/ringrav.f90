program ringrav
use para
implicit none
!r1=92000.d0
!r2=117580.d0
!r1=r1/r_saturn
!r2=r2/r_saturn
!dim=4.d0
!r1=0.5d0
!r2=1.5d0

call readriv
print*, "内径=",r1,"外径=",r2
print*,"1单位长度(土星半径)=60268Km"
print*,"格网边长=",dim

call grid 

end program ringrav

