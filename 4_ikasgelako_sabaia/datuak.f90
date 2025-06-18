module m_datuak
use mcf_tipos

real(kind=dp),public,parameter	::Lx=14.22_dp,Ly=7.1_dp,Lz=3.4_dp,nu=180.0_dp,P0=7.0_dp,c=331.3_dp,kappa=1.418E5_dp
real(kind=dp),public,parameter	::iturrix=1.0_dp,iturriy=3.5_dp,iturriz=1.60_dp,detektore1=3.2_dp
real(kind=dp),public,parameter	::tmax=20.0_dp/nu
character(len=1),public,parameter	::konst1="z",dim="3"
real(kind=dp),public,parameter	::alpha_zura=4.677964325925926E-002,&
				&alpha_ikasle=0.585280316752767,&				
				&alpha_horma=4.473591321033211E-002,&
				&alpha_lurra=1.003234000000000E-002,&
				&alpha_arbela=0.175655359409594,&
				&alpha_sabaia=1.003273210332103E-002,&
				&alpha_leiho=0.175655359409594

end module m_datuak
