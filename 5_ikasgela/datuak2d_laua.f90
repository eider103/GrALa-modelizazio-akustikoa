module m_datuak
use mcf_tipos

real(kind=dp),public,parameter	::Lx=14.22_dp,Ly=7.1_dp,Lz=3.0_dp,nu=180.0_dp,P0=7.0_dp,c=331.3_dp,kappa=1.418E5_dp
real(kind=dp),public,parameter	::iturrix=1.0_dp,iturriy=3.5_dp,iturriz=1.60_dp,detektore1=2.0_dp,detektore2=1.15_dp
real(kind=dp),public,parameter	::tmax=15.0_dp/nu
character(len=1),public,parameter	::konst1="y",dim="2",konst2="z"
real(kind=dp),public,parameter	::alpha_zura=4.677964325925926E-002,&
					&alpha_ikasle=0.585280316752767,&								
					&alpha_horma=4.473591321033211E-002,&
					&alpha_lurra=1.003234000000000E-002,&
					&alpha_arbela=0.175655359409594,&
					&alpha_sabaia=0.179020492693727,&
					&alpha_leiho=0.175655359409594
!oztopoak
integer,public,parameter::zenbat_oztopo=6 !1-tarima, 2-it kutxa, 3-irakasle mahaia, 4-5-6-iakslemahaiak
real(kind=dp),public,parameter,dimension(zenbat_oztopo,3)::obloke=reshape((/0.0_dp,0.0_dp,0.0_dp,&
								& 1.75_dp,0.9_dp,0.3_dp,&
								& 1.65_dp,2.35_dp,0.3_dp,&
								& 3.75_dp,0.0_dp,0.0_dp,&
								& 3.75_dp,2.6_dp,0.0_dp, &
								& 3.75_dp,5.4_dp,0.0_dp/),shape=(/6,3/),order=(/2,1/))
real(kind=dp),public,parameter,dimension(zenbat_oztopo)::xbloke=(/2.44_dp,2.44_dp,2.44_dp,12.7_dp,12.7_dp,12.7_dp/)
real(kind=dp),public,parameter,dimension(zenbat_oztopo)::ybloke=(/6.0_dp,1.5_dp,3.8_dp,1.7_dp,4.4_dp,7.1_dp/)
real(kind=dp),public,parameter,dimension(zenbat_oztopo)::zbloke=(/0.3_dp,1.14_dp,1.35_dp,0.77_dp,0.77_dp,0.77_dp/)
end module m_datuak
