module m_datuak
use mcf_tipos

real(kind=dp),public,parameter	::Lx=10.0_dp,Ly=10.0_dp,Lz=10.0_dp,nu=180.0_dp,P0=7.0_dp,c=331.3_dp,kappa=1.418E5_dp
real(kind=dp),public,parameter	::iturrix=5.0_dp,iturriy=9.9_dp,iturriz=5.0_dp,detektore1=4.9_dp
real(kind=dp),public,parameter	::tmax=10.0_dp/nu
character(len=1),public,parameter	::konst1="z",dim="3"

!oztopoak
integer,public,parameter	::zenbat_oztopo=3
real(kind=dp),public,parameter,dimension(3,3)::obloke=reshape((/0.0_dp,5.0_dp,0.0_dp,3.5_dp,5.0_dp,0.0_dp,7.0_dp,5.0_dp,0.0_dp/),shape=(/3,3/),order=(/2,1/))
real(kind=dp),public,parameter,dimension(3)::xbloke=(/3.0_dp,6.5_dp,Lx/),ybloke=5.5_dp,zbloke=Lz					
end module m_datuak
