program fdtd
use mcf_tipos
use m_sabai_sarea
use m_mb
use m_oinarria
use m_datuak

integer,parameter				::n=400,unit=1
integer						::i,j,k,iz0,kontagailua,itx,ity,itz
real(kind=dp),dimension(:,:,:),allocatable	::p_orain,p_berria,u_orain,u_berria,v_orain,v_berria,w_orain,w_berria
real(kind=dp)					::rho,rv,rp,x,y,tau,t,pt,L,kont_max
real(kind=dp)					::zn_leiho,zn_arbela,zn_zura,zn_horma,zn_lurra,zn_sabaia,zn_ikasle
integer,dimension(habe_px,3)	::oix
integer,dimension(habe_py,3)	::oiy
integer,dimension(habe_px)	::xix,yix
integer,dimension(habe_py)	::xiy,yiy
integer				::iLz

call L_topatu(L,Lx,Ly,Lz) !luzera handiena topatu

!konstanteak definitu
rho=kappa/c/c
h=L/real(n,kind=dp)
tau=h/sqrt(3.0_dp)/c
rv=tau/rho/h
rp=tau*kappa/h

print*,"kont_max=",floor(tmax/tau)
print*,"tau=",tau
print*,"h=",h

call neurketak(c,nu,h,tau)	!neurketa nahikoa egiten ari dela egiaztatu

zn_leiho=rho*c*(1.0_dp+sqrt(1.0_dp-alpha_leiho))/(1.0_dp-sqrt(1.0_dp-alpha_leiho))
zn_arbela=rho*c*(1.0_dp+sqrt(1.0_dp-alpha_arbela))/(1.0_dp-sqrt(1.0_dp-alpha_arbela))
zn_zura=rho*c*(1.0_dp+sqrt(1.0_dp-alpha_zura))/(1.0_dp-sqrt(1.0_dp-alpha_zura))
zn_horma=rho*c*(1.0_dp+sqrt(1.0_dp-alpha_horma))/(1.0_dp-sqrt(1.0_dp-alpha_horma))
zn_lurra=rho*c*(1.0_dp+sqrt(1.0_dp-alpha_lurra))/(1.0_dp-sqrt(1.0_dp-alpha_lurra))
zn_sabaia=rho*c*(1.0_dp+sqrt(1.0_dp-alpha_sabaia))/(1.0_dp-sqrt(1.0_dp-alpha_sabaia))
zn_ikasle=rho*c*(1.0_dp+sqrt(1.0_dp-alpha_ikasle))/(1.0_dp-sqrt(1.0_dp-alpha_ikasle))

!gela laukizuzena da
nx=nint(Lx/h)
ny=nint(Ly/h)
nz=nint(Lz/h)

!matrizeak eraiki:
allocate(p_orain(0:nx,0:ny,0:nz))
allocate(p_berria(0:nx,0:ny,0:nz))

allocate(u_orain(0:nx+1,0:ny,0:nz))
allocate(u_berria(0:nx+1,0:ny,0:nz))

allocate(v_orain(0:nx,0:ny+1,0:nz))
allocate(v_berria(0:nx,0:ny+1,0:nz))

allocate(w_orain(0:nx,0:ny,0:nz+1))
allocate(w_berria(0:nx,0:ny,0:nz+1))

allocate(oztoporik(-1:nx+1,-1:ny+1,-1:nz+1))

print*,"matrizeen tamainak"
print*,"p:",size(p_orain,1),size(p_orain,2),size(p_orain,3)

call sabaia_eraiki(oix,oiy,xix,yix,xiy,yiy,iLz)

!OHARRA: abidaura (3D) eta presio matrizeak despazatuta daude elkarrekiko
!abiadura neurketa puntua bi presio neurketa punturen arteko puntuan egiten da: Yeeren algoritmoa

!hasieratu
p_orain=0.0_dp
p_berria=0.0_dp

u_orain=0.0_dp
u_berria=0.0_dp

v_orain=0.0_dp
v_berria=0.0_dp

w_orain=0.0_dp
w_berria=0.0_dp

oztoporik(:,:,:)=1 !airea, oraingoz
!hormak
oztoporik(-1,:,:)=0
oztoporik(nx+1,:,:)=0
oztoporik(:,-1,:)=0
oztoporik(:,ny+1,:)=0
oztoporik(:,:,-1)=0
oztoporik(:,:,nz+1)=0

kontagailua=0

!maskara binarioa eraiki
do i=1,habe_px
	call oztopo(oix(i,1),oix(i,2),oix(i,3),xix(i),yix(i),iLz,oztoporik)
end do
do i=1,habe_py
	call oztopo(oiy(i,1),oiy(i,2),oiy(i,3),xiy(i),yiy(i),iLz,oztoporik)
end do
do i=1,zenbat_oztopo
	call oztopo(ind(obloke(i,1)),ind(obloke(i,2)),ind(obloke(i,3)),ind(xbloke(i)),ind(ybloke(i)),ind(zbloke(i)),oztoporik)
end do

open(unit=2,file="habeak.dat",status="replace",action="write")
print*,"maskara binarioa"
do i=-1,nx+1
	x=h*i
	do j=-1,ny+1
		y=h*j
		write(unit=2,fmt=*)x,y,oztoporik(i,j,ind(detektore1))
	end do
end do
close(unit=2)

!hasierako baldintzak
t=0.0_dp
!iturriaren posizioa:
itx=ind(iturrix)
ity=ind(iturriy)
itz=ind(iturriz)
p_orain(itx,ity,itz)=p(t)	!eta hasierako abiadura=0

!ekuazioak ebatzi
do
	if (t>tmax) then
		print*,"denbora maximora heldu da"
		stop
	end if
	
do i=0,nx
	do j=0,ny
		do k=0,nz
			if (oztoporik(i,j,k)==0) then
				p_orain(i,j,k)=0.0_dp
				u_orain(i,j,k)=0.0_dp
				v_orain(i,j,k)=0.0_dp
				w_orain(i,j,k)=0.0_dp
			end if
		end do
	end do
end do

	!p(t=0.5*tau) lortu. p(t=-0.5*tau) moduan (p=0*tau) hurbilketa hartu da.
	do i=1,nx
		do j=1,ny	!ertzetako puntuak MBtan
			do k=1,nz
				if ((i==itx).and.(j==ity).and.(k==itz)) then
					p_berria(itx,ity,itz)=pt
				else if (oztoporik(i,j,k)==1) then
					p_berria(i,j,k)=p_orain(i,j,k)-rp*(u_orain(i+1,j,k)-u_orain(i,j,k)+&
									v_orain(i,j+1,k)-v_orain(i,j,k)+w_orain(i,j,k+1)-w_orain(i,j,k))
				end if
			end do
		end do	
	end do
	
	!MB presioarentzat 
	p_berria(0,:,:)=u_orain(0,:,:)*zn_arbela
	p_berria(nx,:,:)=u_orain(nx+1,:,:)*zn_horma
	p_berria(:,0,:)=v_orain(:,0,:)*zn_leiho
	p_berria(:,ny,:)=v_orain(:,ny+1,:)*zn_horma
	p_berria(:,:,0)=w_orain(:,:,0)*zn_lurra
	p_berria(:,:,nz)=w_orain(:,:,nz+1)*zn_sabaia
	
	!oztopoak
	do i=1,habe_px !sabaia
		call mb_p_oztopo(p_berria,u_orain,v_orain,w_orain,zn_sabaia,oix(i,1),oix(i,2),oix(i,3),xix(i),yix(i),iLz)
	end do
	do i=1,habe_py
		call mb_p_oztopo(p_berria,u_orain,v_orain,w_orain,zn_sabaia,oiy(i,1),oiy(i,2),oiy(i,3),xiy(i),yiy(i),iLz)
	end do
	do i=1,3	!oholtza, it kutxa, irakaslearen mahaia
			call mb_p_oztopo(p_berria,u_orain,v_orain,w_orain,zn_zura,ind(obloke(i,1)),ind(obloke(i,2)),ind(obloke(i,3)),ind(xbloke(i)),ind(ybloke(i)),ind(zbloke(i)))
	end do
	do i=4,6	!ikasleen mahaiak
		call mb_p_oztopo(p_berria,u_orain,v_orain,w_orain,zn_ikasle,ind(obloke(i,1)),ind(obloke(i,2)),ind(obloke(i,3)),ind(xbloke(i)),ind(ybloke(i)),ind(zbloke(i)))
	end do
	
	do i=0,nx
		do j=0,ny
			do k=0,nz
				if (oztoporik(i,j,k)==0) then
					p_berria(i,j,k)=0.0_dp
				end if
			end do
		end do
	end do
	
	!datuak gorde
	call idatzi(dim,p_berria,kontagailua,unit,ind(detektore1),konst1)
	
	kontagailua=kontagailua+1
	!lortu u(t=1*tau)
	do i=0,nx
		do j=0,ny
			do k=0,nz
				if (oztoporik(i,j,k)==0) then
					u_orain(i,j,k)=0.0_dp
					v_orain(i,j,k)=0.0_dp
					w_orain(i,j,k)=0.0_dp
				end if
			end do
		end do
	end do

	do i=1,nx
		do j=1,ny
			do k=1,nz	
				if (oztoporik(i,j,k)==1) then
					u_berria(i,j,k)=u_orain(i,j,k)-rv*(p_berria(i,j,k)-p_berria(i-1,j,k))		 !benetako i koordenatua i+0.5
					v_berria(i,j,k)=v_orain(i,j,k)-rv*(p_berria(i,j,k)-p_berria(i,j-1,k))
					w_berria(i,j,k)=w_orain(i,j,k)-rv*(p_berria(i,j,k)-p_berria(i,j,k-1))
				end if
			end do
		end do	
	end do
	
	!p-n falta diren puntuak ezartzeko: onartu p(n+2)~p(n+1)
	!X norabideko MB
	u_berria(0,:,:)=p_berria(0,:,:)/zn_arbela
	u_berria(nx+1,:,:)=p_berria(nx,:,:)/zn_horma

	!Y norabideko MB
	v_berria(:,0,:)=p_berria(:,0,:)/zn_leiho
	v_berria(:,ny+1,:)=p_berria(:,ny,:)/zn_horma	

	!Z norabideko MB
	w_berria(:,:,0)=p_berria(:,:,0)/zn_lurra
	w_berria(:,:,nz+1)=p_berria(:,:,nz)/zn_sabaia

	!oztopoak
	do i=1,habe_px !sabaia
		call mb_u_oztopo(p_berria,u_berria,v_berria,w_berria,zn_sabaia,oix(i,1),oix(i,2),oix(i,3),xix(i),yix(i),iLz)
	end do
	do i=1,habe_py
		call mb_u_oztopo(p_berria,u_berria,v_berria,w_berria,zn_sabaia,oiy(i,1),oiy(i,2),oiy(i,3),xiy(i),yiy(i),iLz)
	end do
	do i=1,3	!oholtza, it kutxa, irakaslearen mahaia
		call mb_u_oztopo(p_berria,u_berria,v_berria,w_berria,zn_zura,ind(obloke(i,1)),ind(obloke(i,2)),ind(obloke(i,3)),ind(xbloke(i)),ind(ybloke(i)),ind(zbloke(i)))
	end do
	do i=4,6	!ikasleen mahaiak
		call mb_u_oztopo(p_berria,u_berria,v_berria,w_berria,zn_ikasle,ind(obloke(i,1)),ind(obloke(i,2)),ind(obloke(i,3)),ind(xbloke(i)),ind(ybloke(i)),ind(zbloke(i)))
	end do
	
do i=0,nx
	do j=0,ny
		do k=0,nz
			if (oztoporik(i,j,k)==0) then
				p_berria(i,j,k)=0.0_dp
				u_berria(i,j,k)=0.0_dp
				v_berria(i,j,k)=0.0_dp
				w_berria(i,j,k)=0.0_dp
			end if
		end do
	end do
end do

	!hurrengo iterazioa prestatu
	u_orain=u_berria
	v_orain=v_berria
	w_orain=w_berria
	p_orain=p_berria
	t=t+tau
	pt=p(t)
	p_orain(itx,ity,itz)=pt
	close(unit=unit)
end do

contains

function p(t) result(pres)
	real(kind=dp),intent(in)::t
	real(kind=dp)		::pres
	real(kind=dp),parameter	::pi=acos(-1.0_dp)
	pres=P0*cos(2.0_dp*pi*nu*t)
end function

subroutine L_topatu(L,Lx,Ly,Lz)
	real(kind=dp),intent(in)	::Lx,Ly,Lz
	real(kind=dp),intent(out)	::L
	L=Lx
	if (L<Ly) then
		L=Ly
	end if
	if (L<Lz) then
		L=Lz
	end if
end subroutine L_topatu

subroutine neurketak(c,nu,h,tau)
	real(kind=dp),intent(in)	::c,nu,h,tau
	if (c/nu/h<10.0_dp) then
		print*,"espazioan uhin luzerako 10 neurketa baino gutxiago daude"
	end if
	if (1.0_dp/nu/tau<10.0_dp) then
		print*,"denboran periodoko 10 neurketa baino gutxiago daude"
	end if	
end subroutine

end program fdtd
