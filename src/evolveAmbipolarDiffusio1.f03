!!!! D. A. Tilly, D. S. Balsara, C. Meyer, 2012, New Astronomy, 17, 368 !!!!


subroutine evolveAD1D(gn,gi,qn,qi)
use gridModule
implicit none
type(grid)::gn,gi
double precision,dimension(1-gn%nbuf:gn%nMesh(1)+gn%nbuf,gn%nvar)::qn
double precision,dimension(1-gi%nbuf:gi%nMesh(1)+gi%nbuf,gi%nvar)::qi
double precision,dimension(1-gn%nbuf:gn%nMesh(1)+gn%nbuf)::nden,nvx,nene_internal,nvx1,nvx2,nene_internal1,nene_internal2
double precision,dimension(1-gi%nbuf:gi%nMesh(1)+gi%nbuf)::iden,ivx,ivy,ivz,iene_internal,ivx1,ivx2
double precision,dimension(1-gi%nbuf:gi%nMesh(1)+gi%nbuf)::iene_internal1,iene_internal2,bxl,byl
double precision,dimension(1-gi%nbuf:gi%nMesh(1)+gi%nbuf)::bzl,bxr,byr,bzr,bxc,byc,bzc
double precision::dt,hdt,qdt
integer::i
double precision::alpha,mu_n,mu_i
double precision::delta,m11,m12,m21,m22,vec1,vec2
double precision::eta_n,eta_i,w,w1,w2
double precision::gam_n,gam_i


alpha=gn%alpha_ad
mu_n=gn%mu_ad
mu_i=gi%mu_ad
gam_n=gn%adiGamma
gam_i=gi%adiGamma

dt=gn%dt
hdt=dt/3.d0
qdt=0.25d0*dt

nden=qn(:,1)
 nvx=qn(:,2)/nden
 nene_internal=qn(:,3)-0.5d0*nden*nvx**2

iden=qi(:,1)
 ivx=qi(:,2)/iden
 ivy=qi(:,3)/iden
 ivz=qi(:,4)/iden
 bxl=qi(:,5)
 byl=qi(:,6)
 bzl=qi(:,7)
 bxr=qi(:,9)
 byr=qi(:,10)
 bzr=qi(:,11)
 bxc=0.5d0*(bxl+bxr)
 byc=0.5d0*(byl+byr)
 bzc=0.5d0*(bzl+bzr)
 iene_internal=qi(:,8)-0.5d0*iden*(ivx**2+ivy**2+ivz**2)-0.5d0*(bxc**2+byc**2+bzc**2) 

do i=1,gn%nMesh(1)
  delta=1.d0+alpha*(iden(i)+nden(i))*qdt
    m11=1.d0+alpha*nden(i)*qdt
    m12=alpha*iden(i)*qdt
    m21=alpha*nden(i)*qdt
    m22=1.d0+alpha*iden(i)*qdt
   vec1=(1.d0-alpha*iden(i)*qdt)*nvx(i)+alpha*iden(i)*qdt*ivx(i)
   vec2=(1.d0-alpha*nden(i)*qdt)*ivx(i)+alpha*nden(i)*qdt*nvx(i)
  nvx1(i)=1.d0/delta*(m11*vec1+m12*vec2)
  ivx1(i)=1.d0/delta*(m21*vec1+m22*vec2)

  eta_n=3.d0*alpha*mu_n*iden(i)*(gam_n-1.d0)/(mu_i+mu_n)
  eta_i=3.d0*alpha*mu_i*nden(i)*(gam_i-1.d0)/(mu_i+mu_n)
       w=alpha*nden(i)*iden(i)*(nvx(i)-ivx(i))**2/(mu_i+mu_n)
      w1=alpha*nden(i)*iden(i)*(nvx1(i)-ivx1(i))**2/(mu_i+mu_n)
  delta=1.d0+qdt*(eta_n+eta_i)
    m11=1.d0+eta_i*qdt
    m12=eta_i*qdt
    m21=eta_n*qdt
    m22=1.d0+eta_n*qdt
   vec1=(1.d0-eta_n*qdt)*nene_internal(i)+(eta_i*qdt)*iene_internal(i)+mu_i*w*qdt+mu_i*w1*qdt
   vec2=(1.d0-eta_i*qdt)*iene_internal(i)+(eta_n*qdt)*nene_internal(i)+mu_n*w*qdt+mu_n*w1*qdt
  nene_internal1(i)=1.d0/delta*(m11*vec1+m12*vec2)
  iene_internal1(i)=1.d0/delta*(m21*vec1+m22*vec2)

  delta=1.d0+alpha*(nden(i)+iden(i))*hdt
    m11=1.d0+alpha*nden(i)*hdt
    m12=alpha*iden(i)*hdt
    m21=alpha*nden(i)*hdt
    m22=1.d0+alpha*iden(i)*hdt
   vec1=(4.d0*nvx1(i)-nvx(i))/3.d0
   vec2=(4.d0*ivx1(i)-ivx(i))/3.d0
  
 nvx2(i)=1.d0/delta*(m11*vec1+m12*vec2)
 ivx2(i)=1.d0/delta*(m21*vec1+m22*vec2)

  delta=1.d0+(eta_i+eta_n)*hdt
     w2=alpha*nden(i)*iden(i)*(nvx2(i)-ivx2(i))**2/(mu_i+mu_n)
    m11=1.d0+eta_i*hdt
    m12=eta_i*hdt
    m21=eta_n*hdt
    m22=1.d0+eta_n*hdt
   vec1=(4.d0*nene_internal1(i)-nene_internal(i)+mu_i*w2*dt)/3.d0
   vec2=(4.d0*iene_internal1(i)-iene_internal(i)+mu_n*w2*dt)/3.d0
  nene_internal2(i)=1.d0/delta*(m11*vec1+m12*vec2)
  iene_internal2(i)=1.d0/delta*(m21*vec1+m22*vec2)

enddo

  qn(:,2)=nden*nvx2
  qi(:,2)=iden*ivx2
  qn(:,3)=nene_internal2+0.5d0*nden*nvx2**2
  qi(:,8)=iene_internal2+0.5d0*iden*(ivx2**2+ivy**2+ivz**2)+0.5d0*(bxc**2+byc**2+bzc**2)
 
end subroutine evolveAD1D


subroutine evolveAD2D(gn,gi,qn,qi)
use gridModule
implicit none
type(grid)::gn,gi
double precision,dimension(1-gn%nbuf:gn%nMesh(1)+gn%nbuf,1-gn%nbuf:gn%nMesh(2)+gn%nbuf,gn%nvar)::qn
double precision,dimension(1-gi%nbuf:gi%nMesh(1)+gi%nbuf,1-gi%nbuf:gi%nMesh(2)+gi%nbuf,gi%nvar)::qi
double precision,dimension(1-gn%nbuf:gn%nMesh(1)+gn%nbuf,1-gn%nbuf:gn%nMesh(2)+gn%nbuf)::nden,nvx,nvy,nene_internal, &
nvx1,nvx2,nvy1,nvy2,nene_internal1,nene_internal2
double precision,dimension(1-gi%nbuf:gi%nMesh(1)+gi%nbuf,1-gi%nbuf:gi%nMesh(2)+gi%nbuf)::iden,ivx,ivy,ivz,&
iene_internal,ivx1,ivx2,ivy1,ivy2
double precision,dimension(1-gi%nbuf:gi%nMesh(1)+gi%nbuf,1-gi%nbuf:gi%nMesh(2)+gi%nbuf)::iene_internal1,iene_internal2,bxl,byl
double precision,dimension(1-gi%nbuf:gi%nMesh(1)+gi%nbuf,1-gi%nbuf:gi%nMesh(2)+gi%nbuf)::bzl,bxr,byr,bzr,bxc,byc,bzc
double precision::dt,hdt,qdt
integer::i,j,nx,ny
double precision::alpha,mu_n,mu_i
double precision::delta,m11,m12,m21,m22,vec1,vec2
double precision::eta_n,eta_i,w,w1,w2
double precision::gam_n,gam_i


nx=gn%nMesh(1)
ny=gn%nMesh(2)
alpha=gn%alpha_ad
mu_n=gn%mu_ad
mu_i=gi%mu_ad
gam_n=gn%adiGamma
gam_i=gi%adiGamma

dt=gn%dt
hdt=dt/3.d0
qdt=0.25d0*dt

nden=qn(:,:,1)
 nvx=qn(:,:,2)/nden
 nvy=qn(:,:,3)/nden
 nene_internal=qn(:,:,4)-0.5d0*nden*(nvx**2+nvy**2)

iden=qi(:,:,1)
 ivx=qi(:,:,2)/iden
 ivy=qi(:,:,3)/iden
 ivz=qi(:,:,4)/iden
 bxl=qi(:,:,5)
 byl=qi(:,:,6)
 bzl=qi(:,:,7)
 bxr=qi(:,:,9)
 byr=qi(:,:,10)
 bzr=qi(:,:,11)
 bxc=0.5d0*(bxl+bxr)
 byc=0.5d0*(byl+byr)
 bzc=0.5d0*(bzl+bzr)
 iene_internal=qi(:,:,8)-0.5d0*iden*(ivx**2+ivy**2+ivz**2)-0.5d0*(bxc**2+byc**2+bzc**2) 

do j=1,gn%nMesh(2)
  do i=1,gn%nMesh(1)
    delta=1.d0+alpha*(iden(i,j)+nden(i,j))*qdt
      m11=1.d0+alpha*nden(i,j)*qdt
      m12=alpha*iden(i,j)*qdt
      m21=alpha*nden(i,j)*qdt
      m22=1.d0+alpha*iden(i,j)*qdt
      vec1=(1.d0-alpha*iden(i,j)*qdt)*nvx(i,j)+alpha*iden(i,j)*qdt*ivx(i,j)
      vec2=(1.d0-alpha*nden(i,j)*qdt)*ivx(i,j)+alpha*nden(i,j)*qdt*nvx(i,j)
      nvx1(i,j)=1.d0/delta*(m11*vec1+m12*vec2)
      ivx1(i,j)=1.d0/delta*(m21*vec1+m22*vec2)
      
      vec1=(1.d0-alpha*iden(i,j)*qdt)*nvy(i,j)+alpha*iden(i,j)*qdt*ivy(i,j)
      vec2=(1.d0-alpha*nden(i,j)*qdt)*ivy(i,j)+alpha*nden(i,j)*qdt*nvy(i,j)
      nvy1(i,j)=1.d0/delta*(m11*vec1+m12*vec2)
      ivy1(i,j)=1.d0/delta*(m21*vec1+m22*vec2)

    eta_n=3.d0*alpha*mu_n*iden(i,j)*(gam_n-1.d0)/(mu_i+mu_n)
    eta_i=3.d0*alpha*mu_i*nden(i,j)*(gam_i-1.d0)/(mu_i+mu_n)
         w=alpha*nden(i,j)*iden(i,j)*(( nvx(i,j)- ivx(i,j))**2+( nvy(i,j)- ivy(i,j))**2)/(mu_i+mu_n)
        w1=alpha*nden(i,j)*iden(i,j)*((nvx1(i,j)-ivx1(i,j))**2+(nvy1(i,j)-ivy1(i,j))**2)/(mu_i+mu_n)
    delta=1.d0+qdt*(eta_n+eta_i)
      m11=1.d0+eta_i*qdt
      m12=eta_i*qdt
      m21=eta_n*qdt
      m22=1.d0+eta_n*qdt
     vec1=(1.d0-eta_n*qdt)*nene_internal(i,j)+(eta_i*qdt)*iene_internal(i,j)+mu_i*w*qdt+mu_i*w1*qdt
     vec2=(1.d0-eta_i*qdt)*iene_internal(i,j)+(eta_n*qdt)*nene_internal(i,j)+mu_n*w*qdt+mu_n*w1*qdt
    nene_internal1(i,j)=1.d0/delta*(m11*vec1+m12*vec2)
    iene_internal1(i,j)=1.d0/delta*(m21*vec1+m22*vec2)

    delta=1.d0+alpha*(nden(i,j)+iden(i,j))*hdt
      m11=1.d0+alpha*nden(i,j)*hdt
      m12=alpha*iden(i,j)*hdt
      m21=alpha*nden(i,j)*hdt
      m22=1.d0+alpha*iden(i,j)*hdt
     vec1=(4.d0*nvx1(i,j)-nvx(i,j))/3.d0
     vec2=(4.d0*ivx1(i,j)-ivx(i,j))/3.d0  
   nvx2(i,j)=1.d0/delta*(m11*vec1+m12*vec2)
   ivx2(i,j)=1.d0/delta*(m21*vec1+m22*vec2)

     vec1=(4.d0*nvy1(i,j)-nvy(i,j))/3.d0
     vec2=(4.d0*ivy1(i,j)-ivy(i,j))/3.d0
   nvy2(i,j)=1.d0/delta*(m11*vec1+m12*vec2)
   ivy2(i,j)=1.d0/delta*(m21*vec1+m22*vec2)


    delta=1.d0+(eta_i+eta_n)*hdt
       w2=alpha*nden(i,j)*iden(i,j)*((nvx2(i,j)-ivx2(i,j))**2+(nvy2(i,j)-ivy2(i,j))**2)/(mu_i+mu_n)
      m11=1.d0+eta_i*hdt
      m12=eta_i*hdt
      m21=eta_n*hdt
      m22=1.d0+eta_n*hdt
     vec1=(4.d0*nene_internal1(i,j)-nene_internal(i,j)+mu_i*w2*dt)/3.d0
     vec2=(4.d0*iene_internal1(i,j)-iene_internal(i,j)+mu_n*w2*dt)/3.d0
    nene_internal2(i,j)=1.d0/delta*(m11*vec1+m12*vec2)
    iene_internal2(i,j)=1.d0/delta*(m21*vec1+m22*vec2)

  enddo
enddo

  qn(:,:,2)=nden*nvx2
  qi(:,:,2)=iden*ivx2
  qn(:,:,3)=nden*nvy2
  qi(:,:,3)=iden*ivy2
  qn(:,:,4)=nene_internal2+0.5d0*nden*(nvx2**2+nvy2**2)
  qi(:,:,8)=iene_internal2+0.5d0*iden*(ivx2**2+ivy2**2+ivz**2)+0.5d0*(bxc**2+byc**2+bzc**2)
 
end subroutine evolveAD2D
