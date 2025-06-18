module m_datuak
use mcf_tipos

real(kind=dp),public,parameter	::Lx=10.0_dp,Ly=10.0_dp,Lz=10.0_dp,nu=180.0_dp,P0=7.0_dp,c=333.1_dp,kappa=1.418E5_dp
real(kind=dp),public,parameter	::iturrix=5.0_dp,iturriy=9.9_dp,iturriz=5.0_dp,detektore1=0.5_dp
real(kind=dp),public,parameter	::tmax=15.0_dp/nu
character(len=1),public,parameter	::konst1="z",dim="3"

real(kind=dp),public,parameter,dimension(3)::obloke=(/4.0_dp,5.0_dp,0.0_dp/)
real(kind=dp),public,parameter			::xbloke=6.0_dp,ybloke=6.0_dp,zbloke=Lz

end module m_datuak
