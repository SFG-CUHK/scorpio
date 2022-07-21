module limiterModule
implicit none

abstract interface
   function limiter(aa,bb,cc)
      double precision::limiter,aa,bb,cc
   end function limiter

   function limiter3(aa,bb,cc,ri,rm,ro,dr)
      double precision::limiter3,aa,bb,cc,ri,rm,ro,dr
   end function limiter3

end interface

contains
function zslop(aa,bb,cc)
double precision::zslop,aa,bb,cc
zslop = 0.d0
end function

function fslop(aa,bb,cc)
double precision::fslop,aa,bb,cc,minmod,ftemp
 fslop=dmax1(dsign(1.d0,(bb-aa)*(cc-bb)),0.d0)*  &
       dsign(1.d0,(cc-aa))*dmin1(dabs(0.5d0*(cc-aa)),  &
       dmin1(2.d0*dabs((bb-aa)),2.d0*dabs((cc-bb))))
end function

function vslop(aa,bb,cc)
double precision::vslop,aa,bb,cc
   vslop=(dsign(1.d0,(bb-aa))+dsign(1.d0,(cc-bb))) &
      *(dabs((bb-aa))*dabs((cc-bb)))/         &
       (dabs((bb-aa))+dabs((cc-bb))+1.0d-7)*1.d0
end function

function minmod(aa,bb,cc)
double precision::minmod,aa,bb,cc
double precision::a,b
    a=cc-bb
    b=bb-aa
    if(a*b .lt. 0.d0) then
       minmod = 0.d0
    elseif(abs(a) .lt. abs(b)) then
       minmod = a
    else
       minmod = b
    endif
end function

function zslop3(aa,bb,cc,ri,rm,ro,dr)
double precision zslop3,aa,bb,cc,ri,rm,ro,dr
zslop3=0.d0
end function

function minmod3(aa,bb,cc,ri,rm,ro,dr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Skinner & Ostriker, 2010, ApJS, 188, 290
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision::minmod3,aa,bb,cc,ri,rm,ro,dr
double precision::a,b,gam,expdel,mexpdel

    a=cc-bb
    b=bb-aa
    gam=dr/(6.d0*rm)
    expdel=2.d0*(ro-rm)/dr-1.d0
    mexpdel=1.d0/expdel

    if(a*b .lt. 0.d0) then
      minmod3 = 0.d0
    elseif(abs(a) .lt. abs(b)) then
      !minmod3 = a/(1.d0-dr**2/(12.d0*ro*rm))*(1.d0+gam)
      minmod3=a/((expdel+1.d0)/2.d0+dr*(rm*expdel**2-ro)/(12.d0*ro*rm))*(1.d0+gam)
    else
      !minmod3 = b/(1.d0-dr**2/(12.d0*rm*ri))*(1.d0+gam)
      minmod3=b/((1.d0+mexpdel)/2.d0+dr*(ri-rm*mexpdel**2)/(12.d0*rm*ri))*(1.d0+gam)
    endif
end function
end module limiterModule
