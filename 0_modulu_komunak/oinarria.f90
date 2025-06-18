module m_oinarria
use mcf_tipos
use m_datuak

real(kind=dp),public :: h
integer,public	::nx,ny,nz
integer,public,dimension(:,:,:),allocatable::oztoporik

public::idatzi,ind,oztopo !oztoporik=1 airea, oztoporik=0 oztopo

contains

subroutine oztopo(x1,y1,z1,x2,y2,z2,oztoporik)
	integer,intent(in)		::x2,y2,z2	!sabaiaren kasuanezin dira gelaren luzerak aldatu.
	integer,intent(in)	::x1,y1,z1	!loditu behar badira, beherantz.
	integer, dimension(-1:nx+1,-1:ny+1,-1:nz+1),intent(inout)::oztoporik
!~ 	if ((x2-x1<=1).or.(y2-y1<=1).or.(z2-z1<=1)) then
!~              !hau beharrezkoa balitz koordenatuak intent(inout) izan behar dira, ezin dira parametroak izan
!~ 		print*,"oztopoak artifizialki loditu dira"
!~ 		x1=min(x1,x2-2)
!~ 		y1=min(y1,y2-2)
!~ 		z1=min(z1,z2-2)
!~ 	end if
	oztoporik(x1+1:x2-1,y1+1:y2-1,z1+1:z2-1)=0	!oztopoa, ez dira ekuazioak eguneratuko
	if (x1==0) then
		oztoporik(x1,y1+1:y2-1,z1+1:z2-1)=0
	end if
	if (x2==ind(Lx)) then
		oztoporik(x2,y1+1:y2-1,z1+1:z2-1)=0
	end if
	if (y1==0) then
		oztoporik(x1+1:x2-1,y1,z1+1:z2-1)=0
	end if
	if (y2==ind(Ly)) then
		oztoporik(x1+1:x2-1,y2,z1+1:z2-1)=0
	end if
	if (z1==0) then
		oztoporik(x1+1:x2-1,y1+1:y2-1,z1)=0
	end if
	if (z2==ind(Lz)) then
		oztoporik(x1+1:x2-1,y1+1:y2-1,z2)=0
	end if
end subroutine oztopo

function ind(x) result (indize)
	real(kind=dp),intent(in)::x
	integer			::indize
	indize=nint(x/h)
	if ((indize>nx).and.(indize>ny).and.(indize>nz)) then
		print*,"indizeak saretik irten dira"
		print*,x
	end if
end function ind

subroutine idatzi(dim,p,kontagailua,unit,idetektore1,konst1,idetektore2,konst2)
	real(kind=dp),dimension(0:,0:,0:),intent(in)	::p
	integer,intent(in)				::kontagailua,unit
	character(len=10)				::izena
	character(len=*),intent(in)			::dim
	integer,optional,intent(in)			::idetektore1	
	character(len=*),optional,intent(in)		::konst1
	integer,optional,intent(in)			::idetektore2
	character(len=*),optional,intent(in)		::konst2
	
	write(izena,"(i4.4,a)")kontagailua,".dat"
	open(unit=unit,file=izena,action="write",status="replace")
	if (dim=="3") then
		if (konst1=="z") then
			do i=0,nx
				x=h*i
				do j=0,ny
					y=j*h
					write(unit=unit,fmt=*)x,y,p(i,j,idetektore1)
				end do
			end do
		else if (konst1=="y") then	!grafikatzeko ordenatuta
			do k=0,nz
				z=h*k
				do i=0,nx
					x=i*h
					write(unit=unit,fmt=*)z,x,p(i,idetektore1,z)
				end do
			end do	
		else if (konst1=="x") then
			do j=0,ny
				y=h*j
				do k=0,nz
					z=k*h
					write(unit=1,fmt=*)y,z,p(idetektore1,j,z)
				end do
			end do
		else
			print*,"(1)Grafikatzeko planoa zehaztu gabe"
			stop
		end if
	else if (dim=="2") then
		if ((konst1=="z").or.(konst2=="z")) then
			if ((konst1=="y").or.(konst2=="y")) then	!aldagaia x
				do i=0,nx
					x=i*h
					write(unit=unit,fmt=*)x,p(i,idetektore1,idetektore2)	!1,2, triangelu zuzenaren ordenean
				end do
			else if ((konst1=="x").or.(konst2=="x")) then	!aldagaia y
				do j=0,ny
					y=j*h
					write(unit=unit,fmt=*)y,p(idetektore1,j,idetektore2)
				end do
			else
				print*,"(2)Grafikatzeko ardatza zehaztu gabe"
				stop
			end if
		else if ((konst1=="y").or.(konst2=="y")) then
			if	((konst1=="x").or.(konst2=="x")) then !aldagaia z
				do k=0,nz
					z=k*h
					write(unit=unit,fmt=*)z,p(idetektore1,idetektore2,k)
				end do
			else
				print*,"(3)Grafikatzeko ardatza zehaztu gabe"
				stop
			end if
		else
			print*,"(4)Grafikatzeko dimentsioak zehaztu gabe"
			stop
		end if
	else
		print*,"(5)Grafikatzeko dimentsioak zehaztu gabe"
		stop
	end if

end subroutine idatzi
end module m_oinarria
