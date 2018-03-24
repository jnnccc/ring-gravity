module para
implicit none
integer*8 ::ndim
integer,parameter::driv_id=20
real*8,parameter:: r_saturn=60268.d0!,Grho=
real*8,allocatable :: grids(:,:,:,:)!存储网格坐标,dim为网格边长
real*8,allocatable :: cyl(:,:,:,:) !网格柱坐标，积分时用到
real*8 :: dim,r1,r2 !环尺寸
!real*8 :: cyl_r,cyl_phi,cyl_z !被积函数中的柱坐标
real*8,parameter::pi=dacos(-1.d0)
end module para
