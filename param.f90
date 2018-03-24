module param
IMPLICIT NONE

!public read_driv,eph_gen,occultation_gen,GLOBAL_PAR!,tracking_gen,doppler_gen,occul_gen
!private EPH_PAR, TRACKING_PAR,OCCULTATION_PAR,DOPPLER_PAR,driv_id,driv_length,itemp,dtemp,dtemp1,cte
!private FRAME, ABCORR,KER,out_file, ntarget, obsrvr, utc0,utc1,spk_id,et0,&
!&step
integer,parameter :: sp = selected_real_kind(p=6,r=37)
integer,parameter :: dp = selected_real_kind(p=15,r=307)
integer,parameter :: qp = selected_real_kind(p=33,r=4931)
real(dp),parameter :: pi=dacos(-1.d0),rsat=60356000d0,ringd(2)=(/66900000d0,74510000d0/),ringc(2)=(/74658000d0,92000000d0/),ringb(2)=(/92800000.d0,117580000.d0/),ringa(2)=(/122170000d0,136775000.d0/)
real(dp)::cyl(6)
CHARACTER*80::ctemp
CHARACTER *1 :: arg1 !运行模式
CHARACTER *20 :: arg2 !土星引力场文件
CHARACTER *3 :: arg3 !径向格点数
CHARACTER *3 :: arg4 !纬向格点数	
CHARACTER *4 :: arg5 !归一化内径
CHARACTER *4 :: arg6 !归一化外径
CHARACTER *4 :: arg7 !纬度下界(度)
CHARACTER *4 :: arg8 !纬度上界(度)
endmodule param
