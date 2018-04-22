;********************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;********************************************************
; open file
;********************************************************
f1=addfile("./wrfout_d01_2005-08-21_00:00:00.nc","r")
a=asciiread("list", (/4/), "string")
;printVarSummary(a)
ntimes=dimsizes(a)
do nt=0,3

print("******************************Processing:"+(nt+1)+"/"+ntimes+"********************************")  
f=addfile(a(nt)+".nc","r")

UU = wrf_user_getvar(f, "A_U", 0)      ; u-wind compt.
;printVarSummary(UU)
VV = wrf_user_getvar(f, "A_V", 0)      ; v-wind compt.
;U=tofloat(wrf_user_unstagger(UU,UU@stagger))
;V=tofloat(wrf_user_unstagger(VV,VV@stagger))
vr = wrf_user_getvar(f1, "avo", 0)
;printVarSummary(U)
lat = f1->XLAT(0,:,0)
lon = f1->XLONG(0,0,:)

MSFU = f1->MAPFAC_U(0,:,:)
;printVarSummary(MSFU)
MSFV = f1->MAPFAC_V(0,:,:)
MSFM = f1->MAPFAC_M(0,:,:)
COR  = f1->F(0,:,:)
;print(COR)
COR  = 0
DX = f@DX
DY = f@DY
;print(DX)
DD = wrf_avo( UU, VV, MSFU, MSFV, MSFM, COR, DX, DY, 0)
;printVarSummary(DD)
ff = new(dimsizes(DD),typeof(DD))
do n=0,33
  do j = 0,149
	do i= 0,164
ff(n,j,i) = DD(n,j,i)*MSFM(j,i)*MSFM(j,i) 
    end do
  end do
end do
ff = ff(:,:,:)
;printVarSummary(DD)

copyatt(ff,vr)
ff@description = ""

ff = ff*10^-5
;system("/bin/rm -f fff1.nc")
ncdf     = addfile(a(nt)+"_fff.nc" ,"c")  ; open output netCDF file
ncdf->ff=ff
;printVarSummary(ff)

end do
end
