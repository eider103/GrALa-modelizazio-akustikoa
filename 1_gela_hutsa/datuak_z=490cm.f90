module m_datuak
use mcf_tipos

real(kind=dp),public,parameter	::Lx=10.0_dp,Ly=10.0_dp,Lz=10.0_dp,nu=180.0_dp,P0=7.0_dp,c=333.1_dp,kappa=1.418E5_dp
real(kind=dp),public,parameter	::iturrix=10.0_dp,iturriy=5.0_dp,iturriz=5.0_dp,detektore1=4.9_dp
real(kind=dp),public,parameter	::tmax=10.0_dp/nu
character(len=1),public,parameter	::konst1="z",dim="3"

end module m_datuak
