module m_sabai_sarea
use mcf_tipos
use m_oinarria
use m_datuak

integer,public,parameter		::habe_px=13,habe_py=7
real(kind=dp),public,parameter	::loderdia=0.1_dp,sak=0.4_dp

public::sabaia_eraiki

contains

subroutine sabaia_eraiki(oix,oiy,xix,yix,xiy,yiy,iLz)
	integer,intent(inout)	::iLz
	integer,dimension(habe_px,3),intent(out)	::oix
	integer,dimension(habe_py,3),intent(out)	::oiy
	integer,dimension(habe_px),intent(out)	::xix,yix
	integer,dimension(habe_py),intent(out)	::xiy,yiy

	real(kind=dp),dimension(habe_px,3)::ozx
	real(kind=dp),dimension(habe_py,3)::ozy
	real(kind=dp),dimension(habe_px)::xzx,yzx	!erdiko karakterea z: koordenatu errealak
	real(kind=dp),dimension(habe_py)::xzy,yzy	!erdiko karakterea i: koordenatu osoak
	real(kind=dp)				::xspace,yspace
	integer	::i,j,k
	
	!hasieratu
	ozx(:,:)=0.0_dp
	ozy(:,:)=0.0_dp
	xzx(:)=0.0_dp
	yzx(:)=0.0_dp
	xzy(:)=0.0_dp
	yzy(:)=0.0_dp
	oix(:,:)=0
	oiy(:,:)=0
	xix(:)=0
	yix(:)=0
	xiy(:)=0
	yiy(:)=0
	
	xspace=(Lx-2.0_dp*loderdia*(real(habe_px-1,kind=dp)))/real(habe_px-1,kind=dp)	!habeen arteko epsazio hutsa
	yspace=(Ly-2.0_dp*loderdia*(real(habe_py-1,kind=dp)))/real(habe_py-1,kind=dp)
	do i=2,size(xzx)-1 !x ardatza ebakitzen duten habeak
		ozx(i,:)=(/loderdia+real(i-1,kind=dp)*xspace+2.0_dp*loderdia*real(i-2,kind=dp),0.0_dp,Lz-sak/)	!habeen lodiera kontuan hartuz 
		xzx(i)=ozx(i,1)+loderdia*2.0_dp
	end do
	yzx(:)=Ly	!habeek y luzera osoa hartzen dute
	do i=2,size(xzy)-1	!y ardatza ebakitzen duten
		ozy(i,:)=(/0.0_dp,loderdia+real(i-1,kind=dp)*yspace+2.0_dp*loderdia*real(i-2,kind=dp),Lz-sak/)
		yzy(i)=ozy(i,2)+loderdia*2.0_dp
	end do
	xzy(:)=Lx	!habeek x luzera osoa hartzen dute
	!ertzeetakoak eskuz
	ozx(1,:)=(/0.0_dp,0.0_dp,Lz-sak/)
	ozx(habe_px,:)=(/Lx-loderdia,0.0_dp,Lz-sak/)
	xzx(1)=loderdia
	xzx(habe_px)=Lx
	ozy(1,:)=(/0.0_dp,0.0_dp,Lz-sak/)	!ozx-ren berdina
	ozy(habe_py,:)=(/0.0_dp,Ly-loderdia,Lz-sak/)
	yzy(1)=loderdia
	yzy(habe_py)=Ly
	iLz=ind(Lz)

	!koordenatu osoak dituzten matrizeak eraiki
	do i=1,habe_px
		do j=1,3
			oix(i,j)=ind(ozx(i,j))
		end do
	end do
	do i=1,habe_py
		do j=1,3
			oiy(i,j)=ind(ozy(i,j))
		end do
	end do
	do i=1,habe_px
		xix(i)=ind(xzx(i))
		yix(i)=ind(yzx(i))
	end do
	do j=1,habe_py
		xiy(j)=ind(xzy(j))
		yiy(j)=ind(yzy(j))
	end do
end subroutine sabaia_eraiki
end module m_sabai_sarea
