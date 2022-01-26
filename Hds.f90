program main
	implicit none
	integer i, j, n, m, k
	real(8) dz, dt, ub, zb, ks, ds, vs, cb, maxerr, S, Rit, dy
	real(8) a1,a2,a3,a4,a5,Rf,Rp,p,q,ua,ca,za,ips,h,sum_all,ru,rc
	real(8) A,E,zetu,HDs,alpha,Fr
	real(8), allocatable :: z(:), zm(:), nu_t(:), nu_tm(:)
	real(8), allocatable :: nu_tcm(:), cof(:)
	real(8), allocatable :: u(:), uo(:), um(:), D_u(:)
	real(8), allocatable :: dudzm(:), dudz(:)
	real(8), allocatable :: a_u(:), b_u(:), c_u(:), r_u(:), ak(:), bk(:), ck(:) 
	real(8), allocatable :: c(:), co(:), cm(:), zeta(:), Xm(:), X(:), fum(:)
	real(8), allocatable :: dcdzm(:), dcdz(:)
	real(8), allocatable :: err_u(:), err_c(:)
	real(8), allocatable :: Ri(:)
	real(8), allocatable :: root1(:),root2(:)
	real(8), parameter :: kappa=0.4d0, sgm_c=1.3d0

    !open(10,file='u003.d')
    !open(20,file='c003.d')
    !open(30,file='nu_t003.d')
    !open(40,file='X003.d')
    !open(50,file='cb(S=0.01).d')
    open(60,file='Hds(Rp=1,S=0.1).d')
	open(10,file='result.txt')

    Rp=1d0 
    S=0.1d0
	ds=1.0d0/1.0d4
	ks=2.5d0*ds
	Rit=1.0d0/S
	zb=1.0d0*ds/12.0d0

	Fr=S**2d0
	!write(*,*)'Fr',Fr
	n=100000

	allocate (z(1:n), zm(0:n), nu_t(1:n), nu_tm(0:n))
	allocate (nu_tcm(0:n), cof(0:n))
	allocate (u(1:n), uo(1:n), um(0:n), D_u(1:n))
	allocate (dudzm(0:n), dudz(1:n))
	allocate (a_u(1:n), b_u(1:n), c_u(1:n), r_u(1:n), ak(1:n), bk(1:n), ck(1:n))
	allocate (c(1:n), co(1:n), cm(0:n), zeta(0:n), Xm(0:n), X(1:n), fum(0:n))
	allocate (dcdzm(0:n), dcdz(1:n))
	allocate (Ri(1:n))
	allocate (root1(1:n),root2(1:n))
	allocate (err_u(1:n), err_c(1:n))
!   
    

	a1=2.891394478d0
	a2=0.95296d0
	a3=0.056834671d0
	a4=0.0002892046d0
	a5=0.000244647d0

	Rf=dexp(-a1+a2*dlog(Rp)-a3*(dlog(Rp))**2.0d0-a4*(dlog(Rp))**3.0d0+a5*(dlog(Rp))**4.0d0) !Dietrich falling velocity
	!write(*,*)'Rf',Rf
    dt=1.0d-9
	dz=(1.0d0-zb)/dble(n)
	do i=0,n
		zm(i)=zb+dz*dble(i)
	enddo
	z(1:n)=zm(1:n)-0.5d0*dz

	!-------- Initial velocity	-----------------------
do k=0,90
    vs=0.01d0+0.001d0*dble(k)	
	!write(*,*)'vs',vs
	za=zm(100) !value of H=0.001+ds/12
	ips=za-zb !thichness of bottom layer
	h=za/zb 

    u(1:n)=1.0d0/kappa*dlog(30.0d0*z(1:n)/ks)
	c(1:n)=1.0d0-z(1:n)-zb
	X(1:n)=z(1:n)
    
!   
    E=(5.73d0/1.0d3)*((1.0d0/vs)**1.31d0)*(Fr**1.59d0)*(Rp**(-0.86d0))
	write(*,*)'E',E
    alpha=-(vs/kappa)
    ca=1.6d0 
    cb=ca*((h)**(-alpha))
    ub=1.0d0/kappa*dlog(30.0d0*zm(0)/ks)
    Xm(1:n-1)=(X(1:n-1)+X(2:n))/2d0
    Xm(n)=1d0
    Xm(0)=0d0

    maxerr=1.0d0
		do while(maxerr>1d-6)

		!------------  gradient	-------------

		uo(1:n)=u(1:n)
		!write(*,*)'uo', uo(n-3:n)
		um(0)=ub
		um(1:n-1)=(uo(1:n-1)+uo(2:n))/2d0
		um(n)=uo(n)

        co(1:n)=c(1:n)
		cm(0)=cb
		cm(1:n-1)=(co(1:n-1)+co(2:n))/2d0
		cm(n)=0d0


	    dudzm(0)=2d0*(uo(1)-um(0))/dz
		dudzm(1:n-1)=(uo(2:n)-uo(1:n-1))/dz
		dudzm(n)=2d0*(um(n)-uo(n))/dz
		!write(*,*)'dudzm', dudzm(n-3:n)
		dudz(1:n)=(um(1:n)-um(0:n-1))/dz


		dcdzm(0)=2d0*(co(1)-cm(0))/dz
		dcdzm(1:n-1)=(co(2:n)-co(1:n-1))/dz
		dcdzm(n)=2d0*(cm(n)-co(n))/dz
		dcdz(1:n)=(cm(1:n)-cm(0:n-1))/dz


		!-------	eddy viscosity	----------------

    cof(0:n)=(dudzm(0:n)**2d0-1.35d0*Rit*dcdzm(0:n))/(dudzm(0:n)**2d0-14.85d0*Rit*dcdzm(0:n))

  	nu_tm(0:n)=cof(0:n)*(kappa**2d0)*((zm(0:n))**2)*(1d0-zm(0:n))*dabs(dudzm(0:n))
  	  
  	!nu_tm(0:n)=kappa**2*zm(0:n)**2*(1d0-zm(0:n))*dabs(dudzm(0:n))

         
		nu_t(1:n)=(nu_tm(0:n-1)+nu_tm(1:n))/2d0

		!-------    u     ----------------------
        
        !bottom layer equation

        do j=0,100
         um(j)=(1.0d0/kappa)*dlog(zm(j)/zb)+ips*((zm(j)-zb)/(2.0d0*kappa)+zb*(dlog(zm(j)/zb))*(cb-kappa+vs)*(1.0d0/(kappa-vs))&
         &-ips*kappa*cb*zb*((zm(j)/zb)**(1.0d0+alpha)-1.0d0)*(0.5d0/(kappa-vs)**2.0d0))
        enddo
		!write(*,*)'um', um(100)
        !upper layer equation

        !coefficient of dudz equation
        ak(1:n-1)=(kappa**2d0)*((zm(1:n-1))**2)*(1d0-zm(1:n-1))
        bk(1:n-1)=-1.35d0*Rit*(kappa**2d0)*((zm(1:n-1))**2)*(1d0-zm(1:n-1))*dcdzm(1:n-1)-(1d0-Xm(1:n-1))
        ck(1:n-1)=14.85d0*Rit*dcdzm(1:n-1)*(1d0-Xm(1:n-1))

        root1(1:n-1)=(-bk(1:n-1)+(bk(1:n-1)**2.0d0-4d0*ak(1:n-1)*ck(1:n-1))**0.5d0)/(2d0*ak(1:n-1))
 
        root2(1:n-1)=(root1(1:n-1))**0.5d0  !solve Quartic equation
        
        dudz(101:n)=(dudzm(100:n-1)+dudzm(101:n))/2d0
		!write(*,*)'dudz', dudz(n-3:n)
        do i=101,n
        um(i)=um(i-1)+dudz(i)*dz
        enddo
		
        u(1:100)=(um(0:99)+um(1:100))/2.0d0
		!write(*,*)'u', u(n-2:n)
        
		!-------------- C,X --------------------------
       do j=0,100
        cm(j)=cb*(zm(j)/zb)**alpha+ips*vs*cb/(2.0d0*kappa)*(3d0*(zm(j)-zb)+zb*(kappa*cb/(kappa-vs)-1.0d0)*dlog(zm(j)/zb)&
         &-kappa**2d0*cb*zb/(kappa-vs)**2d0*((zm(j)/zb)**(1.0d0+alpha)-1.0d0))*(zm(j)/zb)**alpha
        enddo
      !write(*,*)'cm', cm(0:3)

      Xm(0)=0.0d0
      do j=1,100
       Xm(j)=Xm(j-1)+(cm(j)+cm(j-1))/2.0d0*dz
      enddo

	  !write(*,*)'Xm', Xm(100)

		zeta(100)=0d0
		do i=101,n
			zeta(i)=zeta(i-1)+vs/nu_t(i)*dz
!			write(*,*) i, zeta(i), nu_t(i)
		enddo
		! write(*,*)'zeta', zeta(n-3:n)
		fum(100:n)=dexp(-zeta(100:n))
		! write(*,*)'fum', fum(n-3:n)
		cm(100:n)=cm(100)*fum(100:n)
		!write(*,*)'cm', cm(n-3:n)
		do i=101,n
			Xm(i)=Xm(i-1)+(cm(i-1)+cm(i))*dz/2.0d0
!			write(*,*) i, fum(i)
		enddo



		cb=cb/Xm(n)
		do j=0,100
        cm(j)=cb*(zm(j)/zb)**alpha+ips*vs*cb/(2.0d0*kappa)*(3d0*(zm(j)-zb)+zb*(kappa*cb/(kappa-vs)-1.0d0)*dlog(zm(j)/zb)&
         &-kappa**2d0*cb*zb/(kappa-vs)**2d0*((zm(j)/zb)**(1.0d0+alpha)-1.0d0))*(zm(j)/zb)**alpha
        enddo
		cm(100:n)=cm(100)*fum(100:n)
		c(1:n)=(cm(0:n-1)+cm(1:n))/2.0d0
		Xm(0)=0.0d0
		do i=1,n
		Xm(i)=Xm(i-1)+(cm(i-1)+cm(i))*0.5d0*dz
		enddo
		
		X(1:n)=(Xm(0:n-1)+Xm(1:n))/2d0
		! write(*,*)'X', X(n-2:n)
		!write(*,*)'co', co(0:2),'c',c(0:2)
		!---------  error	-----------------
		err_u(1:n)=dabs(uo(1:n)-u(1:n))
		!write(*,*)'err_u', err_u(0:2)
		err_c(1:n)=dabs(co(1:n)-c(1:n))
		!write(*,*)'err_c', err_c(0:2)
		maxerr=max(maxval(err_u(1:n)),maxval(err_c(1:n)))

		!write(*,*)'maxerr=', maxerr,'k=', k

	enddo
	
	Hds=cm(0)*Rf**2d0/(E*S*vs**2d0)
    !write(*,*)'E', E
    write(10,*)cb, vs
    write(*,*)'H=',Hds,'X_all=',Xm(n)
    !write(50,*)vs, cm(0)
    write(60,*)vs, Hds
enddo


   !do j=0,n
   ! write(10,*)um(j),zm(j)
   ! write(20,*)cm(j),zm(j)
   ! write(30,*)nu_tm(j),zm(j)
   ! write(40,*)Xm(j),zm(j)
   !enddo

    close(10)
    !close(20)
    !close(30)
    !close(40)
    !close(50)
    close(60)
end program main
