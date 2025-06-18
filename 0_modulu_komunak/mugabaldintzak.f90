module m_mb
use mcf_tipos
use m_oinarria

public::mb_p_oztopo,mb_u_oztopo

contains

subroutine mb_p_oztopo(p,u,v,w,zn,x1,y1,z1,x2,y2,z2)
	real(kind=dp),dimension(0:,0:,0:),intent(inout)::p
	real(kind=dp),dimension(0:,0:,0:),intent(in)::u
	real(kind=dp),dimension(0:,0:,0:),intent(in)::v
	real(kind=dp),dimension(0:,0:,0:),intent(in)::w
	integer,intent(in)		::x1,y1,z1,x2,y2,z2
	real(kind=dp),intent(in)	::zn
	integer			::i,j,k
	
if (x1 < 0 .or. x2 > nx .or. y1 < 0 .or. y2 > ny .or. z1 < 0 .or. z2 > nz) then
	print*, "Oztopoa eremutik kanpo dago:"
	print*, "x1=",x1," x2=",x2," nx=",nx
	print*, "y1=",y1," y2=",y2," ny=",ny
	print*, "z1=",z1," z2=",z2," nz=",nz
	stop
end if

!islapena
if (x1/=0) then
	do j=y1,y2
		do k=z1,z2
			if ((oztoporik(x1,j,k)==1).and.(oztoporik(x1+1,j,k)==0)) then
				p(x1,j,k)=u(x1,j,k)*zn
			end if
		end do
	end do
end if
if (x2/=nx) then
	do j=y1,y2
		do k=z1,z2
			if ((oztoporik(x2,j,k)==1).and.(oztoporik(x2-1,j,k)==0)) then
				p(x2,j,k)=u(x2,j,k)*zn
			end if
		end do
	end do
end if
if (y1/=0) then
	do i=x1,x2
		do k=z1,z2
			if ((oztoporik(i,y1,k)==1).and.(oztoporik(i,y1+1,k)==0)) then
				p(i,y1,k)=v(i,y1,k)*zn
			end if
		end do
	end do
end if
if (y2/=ny) then				
	do i=x1,x2
		do k=z1,z2
			if ((oztoporik(i,y2,k)==1).and.(oztoporik(i,y2-1,k)==0)) then
				p(i,y2,k)=v(i,y2,k)*zn
			end if
		end do
	end do
end if
if (z1/=0) then
	do i=x1,x2
		do j=y1,y2
			if ((oztoporik(i,j,z1)==1).and.(oztoporik(i,j,z1+1)==0)) then
				p(i,j,z1)=w(i,j,z1)*zn
			end if
		end do
	end do
end if
if (z2/=nz) then
	do i=x1,x2
		do j=y1,y2
			if ((oztoporik(i,j,z2)==1).and.(oztoporik(i,j,z2-1)==0)) then
				p(i,j,z2)=w(i,j,z2)*zn
			end if
		end do
	end do
end if
end subroutine mb_p_oztopo

subroutine mb_u_oztopo(p,u,v,w,zn,x1,y1,z1,x2,y2,z2)
	real(kind=dp),dimension(0:,0:,0:),intent(inout)::u
	real(kind=dp),dimension(0:,0:,0:),intent(inout)::v
	real(kind=dp),dimension(0:,0:,0:),intent(inout)::w
	real(kind=dp),dimension(0:,0:,0:),intent(inout)::p
	real(kind=dp),intent(in)	::zn
	integer,intent(in)	::x1,x2,y1,y2,z1,z2
	integer ::i,j,k
	
if (x1 < 0 .or. x2 > nx .or. y1 < 0 .or. y2 > ny .or. z1 < 0 .or. z2 > nz) then
	print*, "Oztopoa eremutik kanpo dago:"
	print*, "x1=",x1," x2=",x2," nx=",nx
	print*, "y1=",y1," y2=",y2," ny=",ny
	print*, "z1=",z1," z2=",z2," nz=",nz
	stop
end if

!islapena
!X noranzkoa
if (x1/=0) then
	do j=y1,y2
		do k=z1,z2
			if ((oztoporik(x1,j,k)==1).and.(oztoporik(x1+1,j,k)==0)) then
				u(x1,j,k)=p(x1,j,k)/zn
			end if
		end do
	end do
end if
if (x2/=nx) then
	do j=y1,y2
		do k=z1,z2
			if ((oztoporik(x2,j,k)==1).and.(oztoporik(x2-1,j,k)==0)) then
				u(x2,j,k)=p(x2,j,k)/zn
			end if
		end do
	end do
end if
!Y noranzkoa
if (y1/=0) then
	do i=x1,x2
		do k=z1,z2
			if ((oztoporik(i,y1,k)==1).and.(oztoporik(i,y1+1,k)==0)) then
				v(i,y1,k)=p(i,y1,k)/zn
			end if
		end do
	end do
end if
if (y2/=ny) then
	do i=x1,x2
		do k=z1,z2
			if ((oztoporik(i,y2,k)==1).and.(oztoporik(i,y2-1,k)==0)) then
				v(i,y2,k)=p(i,y2,k)/zn
			end if
		end do
	end do
end if
!Z noranzkoa
if (z1/=0) then
	do i=x1,x2
		do j=y1,y2
			if ((oztoporik(i,j,z1)==1).and.(oztoporik(i,j,z1+1)==0)) then
				w(i,j,z1)=p(i,j,z1)/zn
			end if
		end do
	end do
end if
if (z2/=nz) then
	do i=x1,x2
		do j=y1,y2
			if ((oztoporik(i,j,z2)==1).and.(oztoporik(i,j,z2-1)==0)) then
				w(i,j,z2)=p(i,j,z2)/zn
			end if
		end do
	end do
end if
	
end subroutine mb_u_oztopo
end module m_mb
