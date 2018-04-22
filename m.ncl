load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"

begin

f=addfile("./wrfout_d01_2005-08-20_12:00:00.nc","r")
;u00 = wrf_user_getvar(f, "ua", 0)
;printVarSummary(u00)
;v00 = wrf_user_getvar(f, "va", 0)
;slp = wrf_user_getvar(f, "slp", 0)
lat2d = f->XLAT(0,:,:)
lon2d = f->XLONG(0,:,:)
gh = wrf_user_getvar(f, "geopt", 0)
p =  wrf_user_getvar(f, "pressure", 0)
g1=gh/9.8/10
plev = 700.
gh1 = wrf_user_intrp3d(g1,p,"h",plev,0.,False)
wrf_smooth_2d( gh1(:,:), 2 )

mu = addfile("mawa_adssv1", 0)
u = mu->adssv1 

mv = addfile("mawa_adssv2", 0)
v = mv->adssv2 


u = u*10^-5
v = v*10^-5







type="x11"
;tc location
a=59
b=73

;adj domain
w=8
lon11=f->XLONG(0,a-w,b-w)
lat11=f->XLAT(0,a-w,b-w)
lon12=f->XLONG(0,a+w,b-w)
lat12=f->XLAT(0,a+w,b-w)
lon13=f->XLONG(0,a+w,b+w)
lat13=f->XLAT(0,a+w,b+w)
lon14=f->XLONG(0,a-w,b+w)
lat14=f->XLAT(0,a-w,b+w)

lon_d1=(/lon11,lon12,lon13,lon14,lon11/)
lat_d1=(/lat11,lat12,lat13,lat14,lat11/)




lat = lat2d(:,0)
lon = lon2d(0,:)
lat@long_name = "latitude"
lat@units  = "degrees_north"
u!1      = "lat"
u&lat    =  lat
v!1      = "lat"
v&lat    =  lat

lon@long_name = "longitude"
lon@units  = "degrees_east"
u!2      = "lon"
u&lon    =  lon
v!2      = "lon"
v&lon    =  lon



gh1!0 = "lat"
gh1&lat    =  lat
gh1!1      = "lon"
gh1&lon    =  lon



V=sqrt(u*u+v*v)
maxv = max(V)
print(maxv)
do k=0,33
 do i =0,149
     do j=0,164   ;74:167,120:234
         if (V(k,i,j).lt.1) then                      ;(maxv*0.2)) then
             u(k,i,j) = 0
             v(k,i,j) = 0; write branch
         end if
    end do
 end do
end do


wks = gsn_open_wks(type,"0") ; send graphics to PNG file
;       gsn_define_colormap(wks,"BkBlAqGrYeOrReViWh200")
        res                = True
        res@gsnDraw        = False
        res@gsnFrame       = False
        res@gsnMaximize    = True
        res@tmXTOn         = False
        res@tmYROn         = False
        res@gsnLeftString  = "ADSSV(10^-5) at 700hPa"
        res@gsnRightString = ""
        res@gsnAddCyclic       = False
  ;;set map;;       
        mpres                             = res
        mpres@mpDataSetName               = "Earth..1"
        mpres@mpDataBaseVersion           = "MediumRes"
        mpres@mpOutlineOn                 = True
        mpres@mpGeophysicalLineThicknessF = 2
        mpres@mpNationalLineThicknessF    = 2
        mpres@mpFillDrawOrder             = "PostDraw"
;        mpres@mpFillOn                    = True
      ;  mpres@mpFillAreaSpecifiers        = (/"water",       "land" /)
      ;  mpres@mpSpecifiedFillColors       = (/"deepskyblue2","white"/)
;       mpres@mpSpecifiedFillColors      = (/100,0/)
;        mpres@mpMaskAreaSpecifiers        = (/"China:states","Taiwan"/)

;;set area;;
        mpres@mpMinLatF                   = 10
        mpres@mpMaxLatF                   = 35
        mpres@mpMinLonF                   = 130
        mpres@mpMaxLonF                   = 160
;;set vector;;
        res_vc                            = res
        res_vc@vcGlyphStyle               = "LineArrow"
        res_vc@vcLineArrowThicknessF      = 2
        res_vc@vcMinDistanceF             = 0.015
        res_vc@vcRefLengthF               = 0.04
        res_vc@gsnLeftString              =""
;        res_vc@gsnLeftStringOrthogonalPosF = -1.010 
;        res_vc@gsnRightStringOrthogonalPosF = -1.010 

;;wind barb resources don't apply
;;      res_vc@vcGlyphStyle               = "WindBarb"
;;      res_vc@vcWindBarbLineThicknessF   = 5
;;      res_vc@vcWindBarbColor            = "Gray40"

  ;      res_vc@vcRefAnnoOn               = True
        res_vc@vcRefMagnitudeF           = 2
        res_vc@vcRefAnnoString1          = "2"
        res_vc@vcRefLengthF              = 0.035
        res_vc@vcRefAnnoSide             = "Top"
 ;       res_vc@vcRefAnnoString2On        = False
 ;       res_vc@vcRefAnnoPerimOn          = False
        res_vc@vcRefAnnoOrthogonalPosF   = -0.143
;        res_vc@vcRefAnnoParallelPosF     = 0.999
  ;      res_vc@vcRefAnnoBackgroundColor  = "Purple"
        res_vc@vcVectorDrawOrder         = "PostDraw"
  ;      res_vc@gsnRightString            = "Wind"
        res_vc@vcGlyphStyle     = "CurlyVector"



        cnres=mpres
        cnres@gsnLeftStringFontHeightF = 0.012
        cnres@gsnLeftStringOrthogonalPosF = .010
        cnres@gsnRightStringOrthogonalPosF = .010
        cnres@cnInfoLabelOn  = False

        cnres@tmXBLabelFontHeightF = 0.01                ; resize tick labels
        cnres@tmYLLabelFontHeightF = 0.01
        cnres@cnLineThicknessF  = 2
;;plot;;
 ;       map     = gsn_csm_map_ce(wks,mpres)
        
 

do k=0,33
             cnres@gsnLeftString  = "ADSSV(10^-5) at "+k
     contour = gsn_csm_contour_map_ce(wks,gh1,cnres)
        vector  = gsn_csm_vector(wks,u(k,:,:),v(k,:,:),res_vc)

;;overlay filled contours and vectors on the map;;
;        overlay(map,contour)
;        overlay(map,vector)
         overlay(contour,vector)
    mkres = True
    mkres@gsMarkerColor = "red"
    mkres@gsMarkerIndex = 16

    xbox_out1 = new(5,float)
    ybox_out1 = new(5,float)
;    xbox_out2 = new(5,float)
;    ybox_out2 = new(5,float)


;;drawing "map" will draw everything: map, contours, vectors, and text;;
        draw(contour)
        datatondc(contour,lon_d1,lat_d1,xbox_out1,ybox_out1)
        lnres = True
   lnres@gsLineThicknessF  = 6
;   lnres@gsLineDashPattern = 1
   lnres@gsLineColor     = "gray"
   gsn_polyline_ndc(wks, xbox_out1, ybox_out1, lnres)


    frame(wks)
delete(vector)
delete(contour)
end do
end
             
