module subprogs
	implicit none
contains

	subroutine tridag(a,b,c,r,u,n)
		integer n, NMAX
  	real(8) a(n), b(n), c(n), r(n), u(n)
		parameter (NMAX=50000)
		integer j
		real(8) bet, gam(NMAX)
		if (b(1).eq.0.) stop 'tridag: rewrite equations'
		bet=b(1)
		u(1)=r(1)/bet
		do j = 2, n
			gam(j)=c(j-1)/bet
			bet=b(j)-a(j)*gam(j)
			if (bet.eq.0.) stop 'tridag failed'
			u(j)=(r(j)-a(j)*u(j-1))/bet
		enddo
		do j = n-1, 1, -1
			u(j)=u(j)-gam(j+1)*u(j+1)
		enddo
	end subroutine tridag
end module subprogs

program main
	use subprogs
	implicit none
	integer i, j, n
	real(8) dz, dt, Ri, ub, zb, ks, ds, vs, cb, maxerr
	real(8), allocatable :: z(:), zm(:), nu_t(:), nu_tm(:)
	real(8), allocatable :: nu_tcm(:), cf(:)
	real(8), allocatable :: u(:), uo(:), um(:), D_u(:)
	real(8), allocatable :: dudzm(:), dudz(:)
	real(8), allocatable :: a_u(:), b_u(:), c_u(:), r_u(:) 
	real(8), allocatable :: c(:), co(:), cm(:), zeta(:), Xm(:), fum(:)
	real(8), allocatable :: dcdzm(:), dcdz(:)
	real(8), allocatable :: err_u(:), err_c(:)
	real(8), parameter :: kappa=0.4d0, sgm_c=1.3d0

	open(10,file='zu.d')
	open(20,file='zc.d')
	open(30,file='r_vs_Ri0p01.d')

	n=5000

	allocate (z(1:n), zm(0:n), nu_t(1:n), nu_tm(0:n))
	allocate (nu_tcm(0:n), cf(0:n))
	allocate (u(1:n), uo(1:n), um(0:n), D_u(1:n))
	allocate (dudzm(0:n), dudz(1:n))
	allocate (a_u(1:n), b_u(1:n), c_u(1:n), r_u(1:n))
	allocate (c(1:n), co(1:n), cm(0:n), zeta(0:n), Xm(0:n), fum(0:n))
	allocate (dcdzm(0:n), dcdz(1:n))
	allocate (err_u(1:n), err_c(1:n))
!
	Ri=0.01d0
	ds=1.0d0/1.0d4
	ks=2.5d0*ds
	zb=1.0d0*ds
	ub=1.0d0/kappa*dlog(30.0d0*zb/ks)
!
	dt=1.0d-7
	dz=(1.0d0-zb)/dble(n)
	do i=0,n
		zm(i)=zb+dz*dble(i)
	enddo
	z(1:n)=zm(1:n)-0.5d0*dz

	!-------- Initial velocity	-----------------------

	do j=1,n
		read(10,*) z(j), u(j)
		read(20,*) z(j), c(j)
		write(*,*) z(j), u(j), c(j)
	enddo

!	u(1:n)=1.0d0/kappa*dlog(30.0d0*z(1:n)/ks)
!	c(1:n)=1.0d0-z(1:n)

	do j=0,50
		vs=0.001d0*500d0**(0.02d0*dble(50-j))
		maxerr=1.0d0
		do while(maxerr>1d-6)

		uo(1:n)=u(1:n)
		co(1:n)=c(1:n)

		!------------  gradient	-------------

		dudzm(0)=2.0d0*(uo(1)-ub)/dz
		dudzm(1:n-1)=(uo(2:n)-uo(1:n-1))/dz
		dudzm(n)=0.0d0
		dudz(1:n)=(dudzm(0:n-1)+dudzm(1:n))/2.0d0

		dcdzm(0)=2.0d0*(co(1)-cb)/dz
		dcdzm(1:n-1)=(co(2:n)-co(1:n-1))/dz
		dcdzm(n)=2.0d0*(0d0-co(n))/dz
		dcdz(1:n)=(dcdzm(0:n-1)+dcdzm(1:n))/2.0d0

		!-------	eddy viscosity	----------------

		cf(0:n)=(dudzm(0:n)**2d0-1.350d0*Ri*dcdzm(0:n))&
					 /(dudzm(0:n)**2d0-14.85d0*Ri*dcdzm(0:n))
		nu_tm(0:n)=kappa**2d0*zm(0:n)**2d0*dudzm(0:n)*cf(0:n)
		nu_t(1:n)=(nu_tm(0:n-1)+nu_tm(1:n))/2.0d0
		nu_tcm(0:n)=nu_tm(0:n)/sgm_c

		!------- 	velocity	----------------------

		D_u(0:n)=nu_tm(0:n)/dz

		a_u(1)=0.0d0
		a_u(2:n)=-D_u(1:n-1)
		b_u(1)=dz/dt+2.0d0*D_u(0)+D_u(1)
		b_u(2:n-1)=dz/dt+D_u(1:n-2)+D_u(2:n-1)
		b_u(n)=dz/dt+D_u(n-1)
		c_u(1:n-1)=-D_u(1:n-1)
		c_u(n)=0.0d0

		r_u(1)=(uo(1)/dt+c(1))*dz+2.0d0*D_u(0)*ub
		r_u(2:n)=(uo(2:n)/dt+c(2:n))*dz

		call tridag(a_u,b_u,c_u,r_u,u,n)

		!-------------- concentration --------------------------

		zeta(0)=0d0
		do i=1,n
			zeta(i)=zeta(i-1)+vs/nu_t(i)*dz
!			write(*,*) i, zeta(i), nu_t(i)
		enddo

		fum(0:n)=dexp(-zeta(0:n))

		Xm(0)=0.0d0
		do i=1,n
			Xm(i)=Xm(i-1)+(fum(i-1)+fum(i))*dz/2.0d0
!			write(*,*) i, fum(i)
		enddo

		cb=1d0/Xm(n)
		cm(0:n)=cb*fum(0:n)
		c(1:n)=(cm(0:n-1)+cm(1:n))/2.0d0

		!---------  error	-----------------

		err_u(1:n)=dabs(uo(1:n)-u(1:n))
		err_c(1:n)=dabs(co(1:n)-c(1:n))
		maxerr=max(maxval(err_u(1:n)),maxval(err_c(1:n)))

		write(*,*) maxerr, cb

		enddo

		write(30,'(2f15.5)') vs, cm(0)

	enddo

	close(10)
	close(20)
	close(30)

end program main
