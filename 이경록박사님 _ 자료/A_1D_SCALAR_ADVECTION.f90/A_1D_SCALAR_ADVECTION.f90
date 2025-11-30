program RKDG_WENO_1D_SCALAR
  implicit none

  !======================================================================================================================!
  !Counters
  integer :: ord !Index of the degrees of freedom
  integer :: i !Spatial(x-axis) index
  integer :: n !Time index

  !Order of the discontinuous Galerkin method and polynomials
  integer,parameter :: DG_ord=3 !Order of DG solver
  integer,parameter :: Poly_ord=DG_ord-1 !Order of polymomials

  !Useful constants
  real(8),parameter :: pi=4.0d0*datan(1.0d0)
  real(8),parameter :: CFL=0.1d0/dble(2*Poly_ord+1) !CFL number
  real(8) :: TVB_M !TVB constant

  !Spatial variables
  integer,parameter :: imax=50*(2**3) !Number of spatial cells
  integer,parameter :: i_bc=3 !Ghost boundary cells for WENO schemes
  real(8),parameter :: x_left=-1.0d0/1.0d0 !Left boundary on x-axis
  real(8),parameter :: x_right=1.0d0/1.0d0 !Right boundary on x-axis
  real(8),parameter :: x_length=dabs(x_right-x_left) !Length of domain
  real(8),dimension(-i_bc+1:imax+i_bc) :: x !Cell centers
  real(8),dimension(-i_bc+1:imax+i_bc) :: dx !Cell volumes
  real(8):: maxdx !Maximal cell volume
  
  !Time variables
  integer :: RK_ord !Order of Runge-Kutta solver
  integer,parameter :: nmax=34500 !Number of the marching step
  real(8),parameter :: time_out=2.0d0 !Final time of the solution
  real(8) :: dt !Time step size
  real(8) :: flag_0,flag_1 !Time flags

  !Cpu time checkers
  integer :: rate,time_0,time_1

  !Variables for quadratures
  real(8),dimension(1:4,1:2) :: Q 
  
  !Several functions
  real(8),external :: Ini_u,Ext_u

  !Orthogonal polynomials
  real(8),external :: Poly,Dx1_Poly,Dx2_Poly

  !Mass coefficient matrix
  real(8),dimension(0:Poly_ord,0:Poly_ord) :: M

  !Cell interface contributions
  real(8),dimension(0:Poly_ord) :: l_b,r_b

  !Polynomial ideal weights 
  real(8),dimension(0:1,1:4) :: d 

  !Degrees of freedom of the moments by the L2 projection
  real(8),dimension(0:Poly_ord,-i_bc+1:imax+i_bc) :: old_deg 
  real(8),dimension(0:Poly_ord,1:imax) :: new_deg 

  !Calculation variables
  real(8) :: temp
  real(8),dimension(1:imax) :: u,new_u,new_Dx_u,new_Dx2_u
  real(8),dimension(1:imax) :: moment_0,moment_1,moment_2,LIMITER_FLAG
  !======================================================================================================================!
  
  write(*,*) "================================================================================================="
  write(*,*) "                                  NP-RKDG-WENO-LIMITER. (ADVECTION)                              "
  write(*,*) "================================================================================================="
  write(*,*) "                                    Calculations have started.                                   "
  write(*,*) "================================================================================================="
  write(*,*) " "

  !Set several initial settings
  RK_ord=DG_ord-1
  flag_0=0.0d0
  flag_1=flag_0
  call Get_Spatial(i_bc,imax,x_left,x_right,x_length,x,dx,maxdx)
  call Get_Quadrature(Poly_ord,Q)
  call Get_Mass(Poly_ord,maxdx,Q,M)
  call Get_Initial(Poly_ord,Q,M,i_bc,imax,x,dx,maxdx,old_deg,l_b,r_b)
  
  call Get_Weights(Poly_ord,Q,i_bc,imax,x,dx,d)

  !Set the TVB constant
  !if TVB_M is chosen too small(<1.0d0), smooth cells will be declared troubled cells,
  !so unnecessary WENO reconstructions will be appear.

  !if TVB_M is chosen too large(>1.0d0), trouble cells will be declared smooth cells,
  !so suprious oscillations will be appear.
  TVB_M=50.0d0 !1d-2

  !write(*,*) "================================================================================================="
  !write(*,*) "The order of the discontinuous Galerkin solver:",DG_ord
  !write(*,*) "The maximal order of non-polynomials:          ",Poly_ord
  !write(*,*) "The maximal cell size:                         ",maxdx
  !write(*,*) "The selected CFL number:                       ",real(CFL)
  !write(*,*) "The selected TVB constant:                     ",real(TVB_M)
  !write(*,*) "The corresponding Shu's TVB constant:          ",real(TVB_M*maxdx**2)
  !write(*,*) "================================================================================================="
  !write(*,*) " "

  call system_clock(count_rate=rate)
  call system_clock(time_0)

  !============================================================================================================!
  !=======================================Time marching procedure(Start)=======================================!
  !============================================================================================================!

  do n=1,nmax
    flag_1=flag_0
  !Compute the time step size by using the first-derivative of flux
    call Get_Timestep(Poly_ord,i_bc,imax,dx,flag_1,time_out,CFL,dt)
      
  !Update the variable by using TVD RK solver
    call Get_RK(RK_ord,TVB_M,Poly_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,maxdx,dt,old_deg,new_deg,LIMITER_FLAG)

  !Update the real-time
    flag_1=flag_1+dt
    flag_0=flag_1

  !Receive the updated data
    do i=1,imax
      do ord=0,Poly_ord
        old_deg(ord,i)=new_deg(ord,i)
      enddo 
    enddo
    if(flag_0.eq.time_out)then
  !Record the limiter history    
  !    do i=1,imax
  !      write(*,*) i,int(LIMITER_FLAG(i))
  !    enddo
      call Get_Limiter_Flag(i_bc,imax,x,LIMITER_FLAG)

  !Record the numerical result at the finial time
      do i=1,imax
        temp=0.0d0
        do ord=0,Poly_ord
          temp=temp+new_deg(ord,i)*Poly(ord,0.0d0,dx(i))
        enddo
        new_u(i)=temp
      enddo

      do i=1,imax
        temp=0.0d0
        do ord=0,Poly_ord
          temp=temp+new_deg(ord,i)*Dx1_Poly(ord,0.0d0,dx(i))
        enddo
        new_Dx_u(i)=temp
      enddo

      do i=1,imax
        temp=0.0d0
        do ord=0,Poly_ord
          temp=temp+new_deg(ord,i)*Dx2_Poly(ord,0.0d0,dx(i))
        enddo
        new_Dx2_u(i)=temp
      enddo
      goto 10
    endif
  enddo 
  10 continue  

  call system_clock(time_1)

  !============================================================================================================!
  !======================================Time marching procedure(Finish)=======================================!
  !============================================================================================================!

  write(*,*) "=============================Computational and physical time results============================="
  write(*,*) "Number of the time marching step:",n
  write(*,*) "Physical time:                   ",real(flag_0),"//",real(time_out)
  write(*,*) "Computation time:                ",real(time_1-time_0)/real(rate),"sec"
  write(*,*) "================================================================================================="
  write(*,*) " "

  !Set up the exact solution
  open(10,file='Exact_Solution.plt')
  write(10,*)'zone T="file"', 'I=',imax

  do i=1,imax
    u(i)=Ext_u(x(i),time_out)
    write(10,*) x(i),u(i)
  enddo

  close(10)

  !Zero moments
  !call Get_Zero_Moment(i_bc,imax,x,new_deg(0,:))

  !First moments
  !call Get_First_Moment(i_bc,imax,x,new_deg(1,:))
  !call Get_First_Moment(i_bc,imax,x,new_Dx_u)

  !Second moments
  !call Get_Second_Moment(i_bc,imax,x,new_deg(2,:))
  !call Get_Second_Moment(i_bc,imax,x,new_Dx2_u)

  !DG solutions
  call Get_DG_Solution(i_bc,imax,x,new_u)

  !Compute various errors
  write(*,*) "=================================Error results for DG solution==================================="
  call Get_L_2_Error(i_bc,imax,dx,u,new_u)
  call Get_L_Infty_Error(imax,u,new_u)
  write(*,*) "================================================================================================="
  write(*,*) " "

  write(*,*) "================================================================================================="
  write(*,*) "                                   Calculations have finished.                                   "
  write(*,*) "================================================================================================="

  return
endprogram

subroutine Get_Spatial(i_bc,imax,x_left,x_right,x_length,x,dx,maxdx)
  implicit none 

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: x_left,x_right,x_length
  
  real(8),intent(out) :: x(-i_bc+1:imax+i_bc),dx(-i_bc+1:imax+i_bc),maxdx

  !Counter
  integer :: i

  !Calculation variables
  real(8) :: seed
  real(8),dimension(0:imax) :: temp_x
  !-----------------------------------------Calculations have started-----------------------------------------!
  
  !Generate physical cell center locations and volumes of cells
  !Cell interfaces
  do i=0,imax
    temp_x(i)=x_left+dble(i)*dble(x_length/(imax + 1))
  enddo

  !Cell sizes
  maxdx=0.0d0
  do i=1,imax 
    dx(i)=dabs(temp_x(i)-temp_x(i-1))
    maxdx=max(maxdx,dx(i))
  enddo

  !Cell centers
  do i=1,imax
  !TYPE 1. CELL CENTER    
    x(i)=(temp_x(i)+temp_x(i-1))/2.0d0

  !TYPE 2. CELL INTERFACE     
  !  x(i)=temp_x(i)
  enddo
 
  !Set up ghost centers and sizes of cells for Lagrange interpolation
  do i=0,i_bc
	  x(-i)=x(imax-i)-x_length
    x(imax+i)=x(i)+x_length

	  dx(-i)=dx(imax-i)
	  dx(imax+i)=dx(i)
	enddo
  
  return
  !-----------------------------------------Calculations have finished----------------------------------------!
endsubroutine

subroutine Get_Quadrature(Poly_ord,Q)
  implicit none
  
  integer,intent(in) :: Poly_ord
  real(8),intent(out) :: Q(1:4,1:2)

  real(8),parameter :: Jacobian=1.0d0/2.0d0
  !-----------------------------------------Calculations have started------------------------------------------!

  !Quadrature weights are scaled for integrations on [-1.0d0/2.0d0,1.0d0/2.0d0]
  !Jacobian of the affine mapping is producted for each quadrature weights

  !Nodes of 6th-order Gauss-Lobatto quadrature
  Q(1,1)=-1.0d0/1.0d0
  Q(2,1)=-dsqrt(5.0d0)/5.0d0
  Q(3,1)=-Q(2,1)
  Q(4,1)=-Q(1,1)

  !Weights of 6th-order Gauss-Lobatto quadrature
  Q(1,2)=(1.0d0/6.0d0)*Jacobian
  Q(2,2)=(5.0d0/6.0d0)*Jacobian
  Q(3,2)=Q(2,2)
  Q(4,2)=Q(1,2)

  if(Poly_ord.ge.3)then
  !Nodes of 8th-order Gauss-Legendre quadrature
	  Q(1,1)=dsqrt(3.0d0/7.0d0-(2.0d0/7.0d0)*dsqrt(6.0d0/5.0d0)) 
    Q(2,1)=dsqrt(3.0d0/7.0d0+(2.0d0/7.0d0)*dsqrt(6.0d0/5.0d0))
    Q(3,1)=-Q(2,1)
    Q(4,1)=-Q(1,1)

  !Weights of 8th-order Gauss-Legendre quadrature
    Q(1,2)=((18.0d0+dsqrt(30.0d0))/36.0d0)*Jacobian
    Q(2,2)=((18.0d0-dsqrt(30.0d0))/36.0d0)*Jacobian
    Q(3,2)=Q(2,2)
    Q(4,2)=Q(1,2)
  endif
  
  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_Mass(Poly_ord,maxdx,Q,M)
  implicit none
  
  integer,intent(in) :: Poly_ord
  real(8),intent(in) :: maxdx,Q(1:4,1:2) 

  real(8),intent(out) :: M(0:Poly_ord,0:Poly_ord)
 
  !Counters
  integer :: ord_0,ord_1

  !Orthogonal non-polynomials
  real(8),external :: Poly
  !-----------------------------------------Calculations have started------------------------------------------!
  
  !Compute the inverse of the mass matrix
  do ord_0=0,Poly_ord
    do ord_1=0,Poly_ord
      M(ord_0,ord_1)=(Poly(ord_0,Q(1,1)/2.0d0,maxdx)*Poly(ord_1,Q(1,1)/2.0d0,maxdx))*Q(1,2)&
                   &+(Poly(ord_0,Q(2,1)/2.0d0,maxdx)*Poly(ord_1,Q(2,1)/2.0d0,maxdx))*Q(2,2)&
                   &+(Poly(ord_0,Q(3,1)/2.0d0,maxdx)*Poly(ord_1,Q(3,1)/2.0d0,maxdx))*Q(3,2)&
                   &+(Poly(ord_0,Q(4,1)/2.0d0,maxdx)*Poly(ord_1,Q(4,1)/2.0d0,maxdx))*Q(4,2)
      
  !Inversion 
      M(ord_0,ord_1)=1.0d0/(M(ord_0,ord_1)+1d-60)
      if(ord_0.ne.ord_1)then
        M(ord_0,ord_1)=0.0d0
      endif
    enddo
  enddo
  
  write(*,*) "===========================================Mass matrix==========================================="
  do ord_0=0,size(M,1)-1
    write(*,*) real(M(ord_0,:))
  enddo
  write(*,*) "================================================================================================="
  write(*,*) " "

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_Initial(Poly_ord,Q,M,i_bc,imax,x,dx,maxdx,u_0,l_b,r_b)
  implicit none
  
  integer,intent(in) :: Poly_ord
  real(8),intent(in) :: Q(1:4,1:2) 
  real(8),intent(in) :: M(0:Poly_ord,0:Poly_ord)
  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: x(-i_bc+1:imax+i_bc),dx(-i_bc+1:imax+i_bc),maxdx

  real(8),intent(out) :: u_0(0:Poly_ord,-i_bc+1:imax+i_bc)
  real(8),intent(out) :: l_b(0:Poly_ord),r_b(0:Poly_ord)
 
  !Counters
  integer :: ord
  integer :: i

  !Initial data
  real(8),external :: Ini_u

  !Orthogonal polynomials
  real(8),external :: Poly

  !Calculation variables
  real(8) :: temp
  real(8),dimension(1:imax) :: sum_u_0,exact_u_0
  !-----------------------------------------Calculations have started------------------------------------------!
  
  !Compute L2-projection of data into the finite element space
  do i=1,imax
  !Compute degree of freedoms of the moments for the initial data by using quadratures
    do ord=0,Poly_ord
      u_0(ord,i)=Ini_u(x(i)+Q(1,1)*dx(i)/2.0d0)*Poly(ord,Q(1,1)/2.0d0,dx(i))*Q(1,2)&
               &+Ini_u(x(i)+Q(2,1)*dx(i)/2.0d0)*Poly(ord,Q(2,1)/2.0d0,dx(i))*Q(2,2)&
               &+Ini_u(x(i)+Q(3,1)*dx(i)/2.0d0)*Poly(ord,Q(3,1)/2.0d0,dx(i))*Q(3,2)&
               &+Ini_u(x(i)+Q(4,1)*dx(i)/2.0d0)*Poly(ord,Q(4,1)/2.0d0,dx(i))*Q(4,2)
  
  !Product the mass matrix for the orthonormalization                 
      u_0(ord,i)=M(ord,ord)*u_0(ord,i)
    enddo
  enddo

  !Compute the convex summation by using calculated degrees of freedom
  !((x(i)-x(i))/dx(i)=0.0d0: Cell center 
  do i=1,imax
    temp=0.0d0
    do ord=0,Poly_ord
      temp=temp+u_0(ord,i)*Poly(ord,0.0d0,dx(i))
    enddo
    sum_u_0(i)=temp
  enddo
  !call Get_Initial_Condition(i_bc,imax,x,sum_u_0)

  !Set the exact solution of the initial function
  do i=1,imax
    exact_u_0(i)=Ini_u(x(i))
  enddo
  ! call Get_Initial_Condition(i_bc,imax,x,exact_u_0)

  write(*,*) "======================================Initial error result======================================="
  call Get_L_2_Error(i_bc,imax,dx,exact_u_0,sum_u_0)
  call Get_L_Infty_Error(imax,exact_u_0,sum_u_0)
  write(*,*) "================================================================================================="
  write(*,*) " "

  !Compute cell interface contributions
  !((x(i)-dx(i)/2.0d0)-x(i))/dx(i)=-1.0d0/2.0d0: Left cell interface
  !((x(i)+dx(i)/2.0d0)-x(i))/dx(i)=+1.0d0/2.0d0: Right cell interface
  do ord=0,Poly_ord
    l_b(ord)=Poly(ord,-1.0d0/2.0d0,maxdx)
    r_b(ord)=Poly(ord,+1.0d0/2.0d0,maxdx)
  enddo

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_Weights(Poly_ord,Q,i_bc,imax,x,dx,d)
  implicit none
  
  integer,intent(in) :: Poly_ord
  real(8),intent(in) :: Q(1:4,1:2) 
  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: x(-i_bc+1:imax+i_bc)
  real(8),intent(in) :: dx(-i_bc+1:imax+i_bc) 

  real(8),intent(out) :: d(0:1,1:4)
  
  !Counters
  integer :: i
  integer :: l
  integer :: ord
  integer :: quad
 
  !Lagrange coefficient
  real(8),external :: Lagrange_Coefficients

  !Coefficients for global and local polynomials
  real(8),dimension(-1:+0) :: L_K_0
  real(8),dimension(+0:+1) :: L_K_1
  real(8),dimension(-1:+1) :: G_K
  !-----------------------------------------Calculations have started------------------------------------------!
  
  !Ideal weights for the third-order interpolation 
  do quad=1,4
  !Compute global coefficients for each quadrature nodes on [-1,+0,1]
    do l=-1,+1
      G_K(l)=Lagrange_Coefficients(l,-1,+1,i_bc,imax,1,x,dx,Q(quad,1)/2.0d0)
    enddo

  !Compute local coefficients for each quadrature nodes on [-1,+0]
    do l=-1,+0
      L_K_0(l)=Lagrange_Coefficients(l,-1,+0,i_bc,imax,1,x,dx,Q(quad,1)/2.0d0)
    enddo

  !Compute local coefficients for each quadrature nodes on [+0,1]
    do l=+0,+1
      L_K_1(l)=Lagrange_Coefficients(l,+0,+1,i_bc,imax,1,x,dx,Q(quad,1)/2.0d0)
    enddo

  !  write(*,*) quad,real(G_K(-1)),real(G_K(+0)),real(G_K(+1))
  !  write(*,*) quad,real(L_K_0(-1)),real(L_K_0(+0))
  !  write(*,*) quad,real(L_K_1(+0)),real(L_K_1(+1))
    
  !Compute ideal weights
  !  d(0,quad)=G_K(-1)/L_K_0(-1)
  !  d(1,quad)=G_K(+1)/L_K_1(+1)
  !  write(*,*) quad,real(d(0,quad)),real(d(1,quad))
  !  write(*,*) " "
  enddo

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_Timestep(DG_ord,i_bc,imax,dx,flag_1,time_out,CFL,dt)
  implicit none
  
  integer,intent(in) :: DG_ord
  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: dx(-i_bc+1:imax+i_bc)
  real(8),intent(in) :: flag_1,time_out,CFL

  real(8),intent(out) :: dt
  !-----------------------------------------Calculations have started------------------------------------------!
  
  !Determine the time step size
  dt=CFL*(maxval(dx)**(4/3))/1.0d0
  
  !Restriction
  if((flag_1+dt).ge.time_out)then
    dt=(time_out-flag_1)
  endif

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_Boundary(Poly_ord,i_bc,imax,deg)
  implicit none
 
  integer,intent(in) :: Poly_ord
  integer,intent(in) :: i_bc,imax

  real(8),intent(inout) :: deg(0:Poly_ord,-i_bc+1:imax+i_bc)

  !Counters
  integer :: ord
  integer :: i
  !-----------------------------------------Calculations have started------------------------------------------!

  !Apply periodic condition
  do i=1,i_bc
    do ord=0,Poly_ord 
      deg(ord,-i+1)=deg(ord,imax-i+1)
      deg(ord,imax+i)=deg(ord,i)
    enddo
  enddo

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_Original_WENO_Limiter(TVB_M,Poly_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,maxdx,deg,LIMITER_INDEX)
  implicit none
  
  real(8),intent(in) :: TVB_M
  integer,intent(in) :: Poly_ord
  real(8),intent(in) :: Q(1:4,1:2) 
  real(8),intent(in) :: M(0:Poly_ord,0:Poly_ord)
  real(8),intent(in) :: d(0:1,1:4)
  real(8),intent(in) :: l_b(0:Poly_ord),r_b(0:Poly_ord)
  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: x(-i_bc+1:imax+i_bc),dx(-i_bc+1:imax+i_bc),maxdx

  real(8),intent(inout) :: deg(0:Poly_ord,-i_bc+1:imax+i_bc)
  real(8),intent(out) :: LIMITER_INDEX(1:imax)

  !Counters
  integer :: quad
  integer :: ord
  integer :: comp
  integer :: i,i0,i1

  !Boundness for the limiting process
  real(8) :: TVB_C

  !Orthogonal polynomials
  real(8),external :: Poly,Dx1_Poly,Dx2_Poly,Dx3_Poly

  !Lagrange coefficients and polynomials
  real(8),external :: Lagrange_Coefficients 
  real(8),dimension(1:4) :: p_0,p_1

  !Shape parameter
  real(8),dimension(1:4,0:imax+1) :: lambda

  !Several variables for the singularity-cell detectors
  !KXRCF limiter variables
  real(8),dimension(0:imax+1) :: Norm,SI,KXRCF

  !Yang and Yoon's limiter variables
  integer,parameter :: p=1
  real(8),parameter :: xi=1.0d0
  real(8),dimension(-2:2) :: IS
  real(8),dimension(-i_bc+1:imax+i_bc) :: THETA

  !Limiter history
  real(8),dimension(0:imax+1) :: TEMP_LIMITER_INDEX

  !DG polynomial reconstruction variables
  real(8),dimension(0:10) :: temp
  real(8),dimension(1:4,-i_bc+1:imax+i_bc) :: DG_Poly_quad,DG_Dx1_Poly_quad,DG_Dx2_Poly_quad,DG_Dx3_Poly_quad                 
  real(8),dimension(-1:1,1:4,-i_bc+1:imax+i_bc) :: recon_DG_Poly
  real(8),dimension(1:4) :: new_DG_Poly

  !Cell interface values
  real(8),dimension(-i_bc+1:imax+i_bc) :: u_l,u_r,u_c,L_Dx1_u,L_Dx2_u,R_Dx1_u,R_Dx2_u,C_Dx1_u,C_Dx2_u

  !Several free parameters for WENO shceme
  real(8) :: delta
  real(8) :: gamma_0,gamma_1,eta
  real(8) :: exp_0,exp_1,exp_2
  
  !Main WENO reconstruction variables
  integer,parameter :: n=1
  real(8) :: eps
  real(8) :: beta_0,beta_1,tau
  !real(8) :: alpha_0,alpha_1
  real(8) :: omega_0,omega_1
  real(8),dimension(1:4) :: d_0,d_1
  real(8),dimension(1:4) :: alpha_0,alpha_1
  real(8),dimension(1:4) :: WENO
  !-----------------------------------------Calculations have started------------------------------------------!

  !Set up the epsilon parameter
  eps=maxdx**n !maxdx**n !1d-40

  !For given degrees of freedom, compute convex summations at cell interfaces
  do i=0,imax+1
    temp(:)=0.0d0
    do ord=0,Poly_ord
      temp(0)=temp(0)+deg(ord,i)*l_b(ord)
      temp(3)=temp(3)+deg(ord,i)*Dx1_Poly(ord,-1.0d0/2.0d0,dx(i))
      temp(4)=temp(4)+deg(ord,i)*Dx2_Poly(ord,-1.0d0/2.0d0,dx(i))

      temp(2)=temp(2)+deg(ord,i)*Poly(ord,+0.0d0/2.0d0,dx(i))
      temp(7)=temp(7)+deg(ord,i)*Dx1_Poly(ord,+0.0d0/2.0d0,dx(i))
      temp(8)=temp(8)+deg(ord,i)*Dx2_Poly(ord,+0.0d0/2.0d0,dx(i))
      
      temp(1)=temp(1)+deg(ord,i)*r_b(ord)
      temp(5)=temp(5)+deg(ord,i)*Dx1_Poly(ord,+1.0d0/2.0d0,dx(i))
      temp(6)=temp(6)+deg(ord,i)*Dx2_Poly(ord,+1.0d0/2.0d0,dx(i))
    enddo
    u_l(i)=temp(0)
    L_Dx1_u(i)=temp(3)
    L_Dx2_u(i)=temp(4)
    
    u_c(i)=temp(2)
    C_Dx1_u(i)=temp(7)
    C_Dx2_u(i)=temp(8)

    u_r(i)=temp(1)
    R_Dx1_u(i)=temp(5)
    R_Dx2_u(i)=temp(6)
  enddo

  !Set up the quadrature node values 
  do i=-i_bc+1,imax+i_bc
  !Non-derivative values    
    temp(:)=0.0d0
    do ord=0,Poly_ord
      do quad=1,4
        temp(quad)=temp(quad)+deg(ord,i)*Poly(ord,Q(quad,1)/2.0d0,dx(i))
      enddo
    enddo
    do quad=1,4
      DG_Poly_quad(quad,i)=temp(quad)
    enddo

  !First-derivative values
    temp(:)=0.0d0
    do ord=0,Poly_ord
      do quad=1,4
        temp(quad)=temp(quad)+deg(ord,i)*Dx1_Poly(ord,Q(quad,1)/2.0d0,dx(i))
      enddo
    enddo
    do quad=1,4
      DG_Dx1_Poly_quad(quad,i)=temp(quad)
    enddo

  !Second-derivative values
    temp(:)=0.0d0
    do ord=0,Poly_ord
      do quad=1,4
        temp(quad)=temp(quad)+deg(ord,i)*Dx2_Poly(ord,Q(quad,1)/2.0d0,dx(i))
      enddo
    enddo
    do quad=1,4
      DG_Dx2_Poly_quad(quad,i)=temp(quad)
    enddo

  !Third-derivative values
    temp(:)=0.0d0
    do ord=0,Poly_ord
      do quad=1,4
        temp(quad)=temp(quad)+deg(ord,i)*Dx3_Poly(ord,Q(quad,1)/2.0d0,dx(i))
      enddo
    enddo
    do quad=1,4
      DG_Dx3_Poly_quad(quad,i)=temp(quad)
    enddo
  enddo

  !------------------------------------------------------------------------------------------------------------!
  !Activate the vairous Limiters
  !------------------------------------------------------------------------------------------------------------!    
  
  !------------------------------------------------------------------------------------------------------------!    
  !TYPE 1. KXRCF
  !------------------------------------------------------------------------------------------------------------!    
  
  !Prepare the KXRCF shock indicator by using the "density" variable
  do i=0,imax+1
    SI(i)=dabs(u_r(i-1)-u_l(i+0))
    Norm(i)=dabs(maxval(u_c)) !dsqrt((u_c(i)**2)*dx(i))

    KXRCF(i)=SI(i)/(Norm(i)*((dx(i)/2.0d0)**((Poly_ord+1)/2)))
  
  !------------------------------------------------------------------------------------------------------------!      
  !TYPE 2. YANG & YOON'S LIMITER    
  !------------------------------------------------------------------------------------------------------------!    
  
  !local detectors  
  !  IS(-2)=dabs(deg(1,i-1)*dx(i)**1)**p+xi*dabs(deg(2,i-1)*dx(i)**2)**p
  !  IS(+2)=dabs(deg(1,i+1)*dx(i)**1)**p+xi*dabs(deg(2,i+1)*dx(i)**2)**p

    IS(-2)=dabs(((C_Dx1_u(i-1)/dx(i-1)**1)/1.0d0)*dx(i-1)**0)**p&
      &+xi*dabs(((C_Dx2_u(i-1)/dx(i-1)**2)/2.0d0)*dx(i-1)**0)**p

    IS(+2)=dabs(((C_Dx1_u(i+1)/dx(i+1)**1)/1.0d0)*dx(i+1)**0)**p&
      &+xi*dabs(((C_Dx2_u(i+1)/dx(i+1)**2)/2.0d0)*dx(i+1)**0)**p

  !Global detectors      
  !  IS(-1)=dabs((deg(1,i-1)+deg(1,i+0))*dx(i)**1)**p+xi*dabs((deg(2,i-1)+deg(2,i+0))*dx(i)**2)**p
  !  IS(+1)=dabs((deg(1,i+0)+deg(1,i+1))*dx(i)**1)**p+xi*dabs((deg(2,i+0)+deg(2,i+1))*dx(i)**2)**p

    IS(-1)=dabs((((C_Dx1_u(i-1)+C_Dx1_u(i+0))/dx(i-1)**1)/2.0d0)*dx(i-1)**0)**p&
      &+xi*dabs((((C_Dx2_u(i-1)+C_Dx2_u(i+0))/dx(i-1)**2)/4.0d0)*dx(i-1)**0)**p
      
    IS(+1)=dabs((((C_Dx1_u(i+0)+C_Dx1_u(i+1))/dx(i+1)**1)/2.0d0)*dx(i+1)**0)**p&
      &+xi*dabs((((C_Dx2_u(i+0)+C_Dx2_u(i+1))/dx(i+1)**2)/4.0d0)*dx(i+1)**0)**p

  !Set up the singularity detector scheme  
    THETA(i)=(1.0d0/2.0d0)*((min(IS(-1),IS(+1)))/(max(IS(-2),IS(+2))+1d-40))
    
  !------------------------------------------------------------------------------------------------------------!    

  !Conduct the normalization  
    if(real(KXRCF(i)).gt.1.0d0)then  
  !  if(real(THETA(i)).gt.1.0d0)then
  !Singular cell  
      TEMP_LIMITER_INDEX(i)=1.0d0
    else
  !Smooth cell    
      TEMP_LIMITER_INDEX(i)=0.0d0
    endif
    
  !  write(*,*) i,real(x(i)),real(u_c(i)/dx(i)**0),real(C_Dx1_u(i)/dx(i)**1),real(C_Dx2_u(i)/dx(i)**2)
  !  write(*,*) i,real(KXRCF(i)),int(LIMITER_INDEX(i))
  !  write(*,*) i,real(THETA(i)),int(LIMITER_INDEX(i))
  enddo

  do i=1,imax
    LIMITER_INDEX(i)=TEMP_LIMITER_INDEX(i)
  enddo
  
  !call Get_Limiter_Index(i_bc,imax,x,LIMITER_INDEX)
  !stop

  !------------------------------------------------------------------------------------------------------------!
  !Activate the WENO type reconstruction
  !------------------------------------------------------------------------------------------------------------!

  delta=sign(1d-40,deg(0,i+0))
  do i=0,imax+1
  !  if(TEMP_LIMITER_INDEX(i).eq.0.0d0)then  
  !Do not limit
  !    goto 20
  !  else
  !Parameters for first-order generalized undivided differences
      gamma_0=-(deg(0,i+0)-deg(0,i-1))/(deg(0,i+0)+delta)
   
      if(dabs(gamma_0).gt.dx(i)**1) gamma_0=0.0d0
  !    if(gamma_0.ge.0.0d0)then
  !      gamma_0=min(gamma_0,+dx(i)**1)
  !    else
  !      gamma_0=max(gamma_0,-dx(i)**1)
  !    endif

      exp_0=+(gamma_0**0)/1.0d0&
           &+(gamma_0**1)/1.0d0&
           &+(gamma_0**2)/2.0d0&
           &+(gamma_0**3)/6.0d0
      
      gamma_1=-(deg(0,i+1)-deg(0,i+0))/(deg(0,i+0)+delta)
  
      if(dabs(gamma_1).gt.dx(i)**1) gamma_1=0.0d0
  !    if(gamma_1.ge.0.0d0)then
  !      gamma_1=min(gamma_1,+dx(i)**1)
  !    else
  !      gamma_1=max(gamma_1,-dx(i)**1)
  !    endif

      exp_1=+(gamma_1**0)/1.0d0&
           &+(gamma_1**1)/1.0d0&
           &+(gamma_1**2)/2.0d0&
           &+(gamma_1**3)/6.0d0

  !Parameters for the second-order generalized undivided difference
      eta=(deg(0,i-1)-2.0d0*deg(0,i+0)+deg(0,i+1))/(deg(0,i+0)+delta)

      if(dabs(eta).gt.dx(i)**2) eta=0.0d0
  !    if(eta.ge.0.0d0)then
  !      eta=min(eta,+dx(i)**2)
  !    else
  !      eta=max(eta,-dx(i)**2)
  !    endif

      exp_2=+2.0d0*(eta**0)/1.0d0&
           &+2.0d0*(eta**1)/2.0d0&
           &+2.0d0*(eta**2)/24.0d0
  
  !------------------------------------------------------------------------------------------------------------!

  !Local smoothness indicators (WENO-JS,WENO-Z,WENO-P+)
  !    beta_0=(deg(0,i-1)-deg(0,i+0))**2
  !    beta_1=(deg(0,i+0)-deg(0,i+1))**2

  !Local smoothness indicators (WENO-NZ) (Untruncated form)
      beta_0=(deg(0,i-1)-exp_0*deg(0,i+0))**2
      beta_1=(deg(0,i+0)-exp_1*deg(0,i+1))**2
  
  !Local smoothness indicators (WENO-NZ) (Truncated form)
  !    temp(0)=(deg(0,i+0)-deg(0,i-1))**2
  !    temp(0)=(temp(0)/(2.0d0*deg(0,i+0)+dx(i)**2))**2
  !    beta_0=temp(0)*dx(i)**4

  !    temp(0)=(deg(0,i+1)-deg(0,i+0))**2
  !    temp(0)=(temp(0)/(2.0d0*deg(0,i+0)+dx(i)**2))**2
  !    beta_1=temp(0)*dx(i)**4

  !------------------------------------------------------------------------------------------------------------!

  !------------------------------------------------------------------------------------------------------------!

  !Global moothness indicator (WENO-Z)
  !    tau=dabs(beta_1-beta_0)
  
  !Global moothness indicator (WENO-P+)
  !    tau=dabs((beta_0+beta_1)/2.0d0-((deg(0,i-1)-deg(0,i+1))**2)/4.0d0)

  !Global smoothness indicator (WENO-NZ) (Untruncated form)
      tau=(deg(0,i-1)-exp_2*deg(0,i+0)+deg(0,i+1))**2
      
  !Global smoothness indicator (WENO-NZ) (Truncated form)
  !    temp(0)=(deg(0,i-1)-2.0d0*deg(0,i+0)+deg(0,i+1))**2
  !    temp(0)=(temp(0)/(24.0d0*deg(0,i+0)+dx(i)**4))**2
  !    tau=temp(0)*dx(i)**8

  !------------------------------------------------------------------------------------------------------------!
  !Compute nonlinear normalized weights by using WENO schemes
  !    do quad=1,4
  !------------------------------------------------------------------------------------------------------------!
  
  !Set up local kernels and ideal weights on exponential polynomials by using shape paramters on each quadrature nodes
        do quad=1,4
          lambda(quad,i)=DG_Dx3_Poly_quad(quad,i)
          lambda(quad,i)=lambda(quad,i)/(DG_Dx1_Poly_quad(quad,i)+sign(dx(i)**1,DG_Dx1_Poly_quad(quad,i)))
          lambda(quad,i)=(lambda(quad,i)**2)
        enddo

  !      write(*,*) i,real(lambda(1,i)),-real(((4.0d0*datan(1.0d0))**2)*(dx(i)**2)),real(dx(i)**2)&
  !                 &,real(lambda(2,i)),-real(((4.0d0*datan(1.0d0))**2)*(dx(i)**2)),real(dx(i)**2)&
  !                 &,real(lambda(3,i)),-real(((4.0d0*datan(1.0d0))**2)*(dx(i)**2)),real(dx(i)**2)&
  !                 &,real(lambda(4,i)),-real(((4.0d0*datan(1.0d0))**2)*(dx(i)**2)),real(dx(i)**2)        

  !      write(*,*) i,real(DG_Dx1_Poly_quad(1,i)),real(DG_Dx1_Poly_quad(2,i))&
  !                 &,real(DG_Dx1_Poly_quad(3,i)),real(DG_Dx1_Poly_quad(4,i))&
  !                 &,real(dx(i)**1)

  !      write(*,*) i,real(DG_Dx3_Poly_quad(1,i)),real(DG_Dx3_Poly_quad(2,i))&
  !                 &,real(DG_Dx3_Poly_quad(3,i)),real(DG_Dx3_Poly_quad(4,i))&
  !                 &,real(dx(i)**3)      
  !      stop

  !Unnormalized weights (WENO-NZ)
  !1st quadrature node    
        d_0(1)=(+66.0d0/4589.0d0)*(lambda(1,i)**2)-(127.00d0/1068.0d0)*(lambda(1,i)**1)+(907.0d0/1374.0d0)
        d_1(1)=(-66.0d0/4589.0d0)*(lambda(1,i)**2)+(127.00d0/1068.0d0)*(lambda(1,i)**1)+(467.0d0/1374.0d0)

        alpha_0(1)=d_0(1)*(1.0d0+tau/(beta_0+eps)+2.0d0*(beta_0/(tau+eps))**2)
        alpha_1(1)=d_1(1)*(1.0d0+tau/(beta_1+eps)+2.0d0*(beta_1/(tau+eps))**2)

  !2nd quadrature node          
        d_0(2)=(+57.0d0/5741.0d0)*(lambda(2,i)**2)-(123.0d0/1618.0d0)*(lambda(2,i)**1)+(985.0d0/2582.0d0)
        d_1(2)=(-57.0d0/5741.0d0)*(lambda(2,i)**2)+(123.0d0/1618.0d0)*(lambda(2,i)**1)+(441.0d0/713.0d0)

        alpha_0(2)=d_0(2)*(1.0d0+tau/(beta_0+eps)+2.0d0*(beta_0/(tau+eps))**2)
        alpha_1(2)=d_1(2)*(1.0d0+tau/(beta_1+eps)+2.0d0*(beta_1/(tau+eps))**2)

  !3rd quadrature node    
        d_0(3)=(+186.0d0/16139.0d0)*(lambda(3,i)**2)-(824.0d0/8125.0d0)*(lambda(3,i)**1)+(441.0d0/713.0d0)
        d_1(3)=(-186.0d0/16139.0d0)*(lambda(3,i)**2)+(824.0d0/8125.0d0)*(lambda(3,i)**1)+(985.0d0/2582.0d0)

        alpha_0(3)=d_0(3)*(1.0d0+tau/(beta_0+eps)+2.0d0*(beta_0/(tau+eps))**2)
        alpha_1(3)=d_1(3)*(1.0d0+tau/(beta_1+eps)+2.0d0*(beta_1/(tau+eps))**2)

  !4th quadrature node
        d_0(4)=(+49.0d0/4008.0d0)*(lambda(4,i)**2)-(161.0d0/1903.0d0)*(lambda(4,i)**1)+(467.0d0/1374.0d0)
        d_1(4)=(-49.0d0/4008.0d0)*(lambda(4,i)**2)+(161.0d0/1903.0d0)*(lambda(4,i)**1)+(907.0d0/1374.0d0)

        alpha_0(4)=d_0(4)*(1.0d0+tau/(beta_0+eps)+2.0d0*(beta_0/(tau+eps))**2)
        alpha_1(4)=d_1(4)*(1.0d0+tau/(beta_1+eps)+2.0d0*(beta_1/(tau+eps))**2)
                
  !------------------------------------------------------------------------------------------------------------!              
   
  !------------------------------------------------------------------------------------------------------------!

  !Set up adaptive local polynomials
  ![-1,+0]
        p_0(1)=((0.0d0/1.0d0)*(lambda(1,i)**2)+(0.0d0/1.0d0)*(lambda(1,i)**1)-(179.0d0/1053.0d0))*deg(0,i-1)&
             &+((0.0d0/1.0d0)*(lambda(1,i)**2)+(0.0d0/1.0d0)*(lambda(1,i)**1)+(1232.0d0/1053.0d0))*deg(0,i+0)

        p_0(2)=((0.0d0/1.0d0)*(lambda(2,i)**2)+(0.0d0/1.0d0)*(lambda(2,i)**1)-(462.0d0/1073.0d0))*deg(0,i-1)&
             &+((0.0d0/1.0d0)*(lambda(2,i)**2)+(0.0d0/1.0d0)*(lambda(2,i)**1)+(1535.0d0/1073.0d0))*deg(0,i+0)
                
        p_0(3)=((0.0d0/1.0d0)*(lambda(3,i)**2)+(0.0d0/1.0d0)*(lambda(3,i)**1)+(462.0d0/1073.0d0))*deg(0,i-1)&
             &+((0.0d0/1.0d0)*(lambda(3,i)**2)+(0.0d0/1.0d0)*(lambda(3,i)**1)+(611.0d0/1073.0d0))*deg(0,i+0)
                
        p_0(4)=((0.0d0/1.0d0)*(lambda(4,i)**2)+(0.0d0/1.0d0)*(lambda(4,i)**1)+(179.0d0/1053.0d0))*deg(0,i-1)&
             &+((0.0d0/1.0d0)*(lambda(4,i)**2)+(0.0d0/1.0d0)*(lambda(4,i)**1)+(874.0d0/1053.0d0))*deg(0,i+0)

  ![+0,+1]
        p_1(1)=((0.0d0/1.0d0)*(lambda(1,i)**2)+(0.0d0/1.0d0)*(lambda(1,i)**1)+(874.0d0/1053.0d0))*deg(0,i+0)&
             &+((0.0d0/1.0d0)*(lambda(1,i)**2)+(0.0d0/1.0d0)*(lambda(1,i)**1)+(179.0d0/1053.0d0))*deg(0,i+1)
              
        p_1(2)=((0.0d0/1.0d0)*(lambda(2,i)**2)+(0.0d0/1.0d0)*(lambda(2,i)**1)+(611.0d0/1073.0d0))*deg(0,i+0)&
             &+((0.0d0/1.0d0)*(lambda(2,i)**2)+(0.0d0/1.0d0)*(lambda(2,i)**1)+(462.0d0/1073.0d0))*deg(0,i+1)
                
        p_1(3)=((0.0d0/1.0d0)*(lambda(3,i)**2)+(0.0d0/1.0d0)*(lambda(3,i)**1)+(1535.0d0/1073.0d0))*deg(0,i+0)&
             &+((0.0d0/1.0d0)*(lambda(3,i)**2)+(0.0d0/1.0d0)*(lambda(3,i)**1)-(462.0d0/1073.0d0))*deg(0,i+1)
                
        p_1(4)=((0.0d0/1.0d0)*(lambda(4,i)**2)+(0.0d0/1.0d0)*(lambda(4,i)**1)+(1232.0d0/1053.0d0))*deg(0,i+0)&
             &+((0.0d0/1.0d0)*(lambda(4,i)**2)+(0.0d0/1.0d0)*(lambda(4,i)**1)-(179.0d0/1053.0d0))*deg(0,i+1)
  
  ![-1,+0]
  !      p_0(1)=((-45.0d0/9949.0d0)*(lambda(1,i)**2)+(158.0d0/4567.0d0)*(lambda(1,i)**1)-(179.0d0/1053.0d0))*deg(0,i-1)&
  !           &+((+45.0d0/9949.0d0)*(lambda(1,i)**2)-(158.0d0/4567.0d0)*(lambda(1,i)**1)+(1232.0d0/1053.0d0))*deg(0,i+0)

  !      p_0(2)=((-275.0d0/29771.0d0)*(lambda(2,i)**2)+(235.0d0/3076.0d0)*(lambda(2,i)**1)-(462.0d0/1073.0d0))*deg(0,i-1)&
  !           &+((+275.0d0/29771.0d0)*(lambda(2,i)**2)-(235.0d0/3076.0d0)*(lambda(2,i)**1)+(1535.0d0/1073.0d0))*deg(0,i+0)
                
  !      p_0(3)=((+275.0d0/29771.0d0)*(lambda(3,i)**2)-(235.0d0/3076.0d0)*(lambda(3,i)**1)+(462.0d0/1073.0d0))*deg(0,i-1)&
  !           &+((-275.0d0/29771.0d0)*(lambda(3,i)**2)+(235.0d0/3076.0d0)*(lambda(3,i)**1)+(611.0d0/1073.0d0))*deg(0,i+0)
                
  !      p_0(4)=((+45.0d0/9949.0d0)*(lambda(4,i)**2)-(158.0d0/4567.0d0)*(lambda(4,i)**1)+(179.0d0/1053.0d0))*deg(0,i-1)&
  !           &+((-45.0d0/9949.0d0)*(lambda(4,i)**2)+(158.0d0/4567.0d0)*(lambda(4,i)**1)+(874.0d0/1053.0d0))*deg(0,i+0)

  ![+0,+1]
  !      p_1(1)=((+312.0d0/49211.0d0)*(lambda(1,i)**2)-(662.0d0/8529.0d0)*(lambda(1,i)**1)+(874.0d0/1053.0d0))*deg(0,i+0)&
  !           &+((-312.0d0/49211.0d0)*(lambda(1,i)**2)+(662.0d0/8529.0d0)*(lambda(1,i)**1)+(179.0d0/1053.0d0))*deg(0,i+1)
              
  !      p_1(2)=((+149.0d0/15194.0d0)*(lambda(2,i)**2)-(521.0d0/5930.0d0)*(lambda(2,i)**1)+(611.0d0/1073.0d0))*deg(0,i+0)&
  !           &+((-149.0d0/15194.0d0)*(lambda(2,i)**2)+(521.0d0/5930.0d0)*(lambda(2,i)**1)+(462.0d0/1073.0d0))*deg(0,i+1)
                
  !      p_1(3)=((-47.0d0/3841.0d0)*(lambda(3,i)**2)+(949.0d0/4997.0d0)*(lambda(3,i)**1)+(1535.0d0/1073.0d0))*deg(0,i+0)&
  !           &+((+47.0d0/3841.0d0)*(lambda(3,i)**2)-(949.0d0/4997.0d0)*(lambda(3,i)**1)-(462.0d0/1073.0d0))*deg(0,i+1)
                
  !      p_1(4)=((-128.0d0/2537.0d0)*(lambda(4,i)**2)+(317.0d0/13675.0d0)*(lambda(4,i)**1)+(1232.0d0/1053.0d0))*deg(0,i+0)&
  !           &+((+128.0d0/2537.0d0)*(lambda(4,i)**2)-(317.0d0/13675.0d0)*(lambda(4,i)**1)-(179.0d0/1053.0d0))*deg(0,i+1)
  !------------------------------------------------------------------------------------------------------------!           

      do quad=1,4

  !------------------------------------------------------------------------------------------------------------!      
  !Unnormalized weights (WENO-JS)
  !      alpha_0=d(0,quad)/(beta_0+eps)**2
  !      alpha_1=d(1,quad)/(beta_1+eps)**2

  !Unnormalized weights (WENO-Z)
  !      alpha_0=d(0,quad)*(1.0d0+tau/(beta_0+eps))
  !      alpha_1=d(1,quad)*(1.0d0+tau/(beta_1+eps))
    
  !Unnormalized weights (WENO-P+)
  !      alpha_0=d(0,quad)*(1.0d0+tau/(beta_0+eps)+(dx(i)**(1.0d0/6.0d0))*((beta_0+eps)/(tau+eps)))
  !      alpha_1=d(1,quad)*(1.0d0+tau/(beta_0+eps)+(dx(i)**(1.0d0/6.0d0))*((beta_1+eps)/(tau+eps)))

  !      alpha_0=d(0,quad)*(1.0d0+tau/(beta_0+eps)+2.0d0*(beta_0/(tau+eps))**2)
  !      alpha_1=d(1,quad)*(1.0d0+tau/(beta_1+eps)+2.0d0*(beta_1/(tau+eps))**2)      
  !------------------------------------------------------------------------------------------------------------!      

  !Normalized weights
  !      omega_0=alpha_0/(alpha_0+alpha_1)
  !      omega_1=alpha_1/(alpha_0+alpha_1)

        omega_0=alpha_0(quad)/(alpha_0(quad)+alpha_1(quad))
        omega_1=alpha_1(quad)/(alpha_0(quad)+alpha_1(quad))
    
  !Set up adaptive local polynomials
  ![-1,+0]
  !      p_0(quad)=Lagrange_Coefficients(-1,-1,+0,i_bc,imax,1,x,dx,Q(quad,1)/2.0d0)*deg(0,i-1)&
  !              &+Lagrange_Coefficients(+0,-1,+0,i_bc,imax,1,x,dx,Q(quad,1)/2.0d0)*deg(0,i+0)

  ![+0,+1]
  !      p_1(quad)=Lagrange_Coefficients(+0,+0,+1,i_bc,imax,1,x,dx,Q(quad,1)/2.0d0)*deg(0,i+0)&
  !              &+Lagrange_Coefficients(+1,+0,+1,i_bc,imax,1,x,dx,Q(quad,1)/2.0d0)*deg(0,i+1)        

  !Compute convex summations on [-1,+0,+1] by using WENO weights and its local polynomials
        WENO(quad)=omega_0*p_0(quad)&
                 &+omega_1*p_1(quad)
      enddo

  !Evaluate degrees of freedom the moments by using the modified DG polynomials
      do ord=1,Poly_ord
        deg(ord,i)=WENO(1)*Poly(ord,Q(1,1)/2.0d0,dx(i))*Q(1,2)&
                 &+WENO(2)*Poly(ord,Q(2,1)/2.0d0,dx(i))*Q(2,2)&
                 &+WENO(3)*Poly(ord,Q(3,1)/2.0d0,dx(i))*Q(3,2)&
                 &+WENO(4)*Poly(ord,Q(4,1)/2.0d0,dx(i))*Q(4,2)

  !Product the mass matrix for the orthonormalization                 
        deg(ord,i)=M(ord,ord)*deg(ord,i)
      enddo
  !  endif
  !  20 continue
  enddo
  
  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_Simple_WENO_Limiter(TVB_M,Poly_ord,Q,M,l_b,r_b,i_bc,imax,x,dx,maxdx,deg,LIMITER_INDEX)
  implicit none
  
  real(8),intent(in) :: TVB_M
  integer,intent(in) :: Poly_ord
  real(8),intent(in) :: Q(1:4,1:2) 
  real(8),intent(in) :: M(0:Poly_ord,0:Poly_ord)
  real(8),intent(in) :: l_b(0:Poly_ord),r_b(0:Poly_ord)
  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: x(-i_bc+1:imax+i_bc),dx(-i_bc+1:imax+i_bc),maxdx

  real(8),intent(inout) :: deg(0:Poly_ord,-i_bc+1:imax+i_bc)
  real(8),intent(out) :: LIMITER_INDEX(1:imax)

  !Counters
  integer :: quad
  integer :: ord
  integer :: comp
  integer :: i,i0,i1
  integer :: l

  !Boundness for the limiting process
  real(8) :: TVB_C

  !Orthogonal polynomials
  real(8),external :: Poly,Dx1_Poly,Dx2_Poly

  !Linear coefficients for the WENO type reconstruction
  real(8),dimension(-1:1) :: linear_coeff

  !Several variables for the singularity-cell detectors
  !KXRCF limiter variables
  real(8),dimension(1:imax) :: Norm,SI,KXRCF

  !Yang and Yoon's limiter variables
  integer,parameter :: p=1
  real(8),parameter :: xi=1.0d0
  real(8),dimension(-2:2) :: IS
  real(8),dimension(1:imax) :: THETA

  !DG polynomial reconstruction variables
  real(8),dimension(0:10) :: temp
  real(8),dimension(1:4,-i_bc+1:imax+i_bc) :: DG_Poly_quad,DG_Dx1_Poly_quad,DG_Dx2_Poly_quad                 
  real(8),dimension(-1:1,1:4,-i_bc+1:imax+i_bc) :: recon_DG_Poly
  real(8),dimension(1:4) :: new_DG_Poly

  !Cell interface values
  real(8),dimension(-i_bc+1:imax+i_bc) :: u_l,u_r,u_c,L_Dx1_u,L_Dx2_u,R_Dx1_u,R_Dx2_u,C_Dx1_u,C_Dx2_u
  !real(8),dimension(-2:2,-i_bc+1:imax+i_bc) :: Dx1_u,Dx2_u
  
  !Main WENO type reconstruction variables
  integer,parameter :: n=1
  real(8) :: eps
  real(8),dimension(-1:1,1:imax) :: beta
  real(8),dimension(-1:1) :: alpha
  real(8) :: alpha_sum,tau
  real(8) :: max_beta_0,max_beta_1,max_beta_2
  real(8),dimension(-1:1) :: omega
  !-----------------------------------------Calculations have started------------------------------------------!

  !Set up the epsilon parameter
  eps=1d-40 !maxdx**n
  
  !For given degrees of freedom, compute convex summations at cell interfaces
  do i=-i_bc+1,imax+i_bc
    temp(:)=0.0d0
    do ord=0,Poly_ord
      temp(0)=temp(0)+deg(ord,i)*l_b(ord)
      temp(3)=temp(3)+deg(ord,i)*Dx1_Poly(ord,-1.0d0/2.0d0,dx(i))
      temp(4)=temp(4)+deg(ord,i)*Dx2_Poly(ord,-1.0d0/2.0d0,dx(i))

      temp(2)=temp(2)+deg(ord,i)*Poly(ord,+0.0d0/2.0d0,dx(i))
      temp(7)=temp(7)+deg(ord,i)*Dx1_Poly(ord,+0.0d0/2.0d0,dx(i))
      temp(8)=temp(8)+deg(ord,i)*Dx2_Poly(ord,+0.0d0/2.0d0,dx(i))
      
      temp(1)=temp(1)+deg(ord,i)*r_b(ord)
      temp(5)=temp(5)+deg(ord,i)*Dx1_Poly(ord,+1.0d0/2.0d0,dx(i))
      temp(6)=temp(6)+deg(ord,i)*Dx2_Poly(ord,+1.0d0/2.0d0,dx(i))
    enddo
    u_l(i)=temp(0)
    L_Dx1_u(i)=temp(3)
    L_Dx2_u(i)=temp(4)
    
    u_c(i)=temp(2)
    C_Dx1_u(i)=temp(7)
    C_Dx2_u(i)=temp(8)

    u_r(i)=temp(1)
    R_Dx1_u(i)=temp(5)
    R_Dx2_u(i)=temp(6)
  enddo

  !Set up the quadrature node values 
  do i=-i_bc+1,imax+i_bc
  !Non-derivative values    
    temp(:)=0.0d0
    do ord=0,Poly_ord
      do l=1,4
        temp(l)=temp(l)+deg(ord,i)*Poly(ord,Q(l,1)/2.0d0,dx(i))
      enddo
    enddo
    do l=1,4
      DG_Poly_quad(l,i)=temp(l)
    enddo

  !First-derivative values
    temp(:)=0.0d0
    do ord=0,Poly_ord
      do l=1,4
        temp(l)=temp(l)+deg(ord,i)*Dx1_Poly(ord,Q(l,1)/2.0d0,dx(i))
      enddo
    enddo
    do l=1,4
      DG_Dx1_Poly_quad(l,i)=temp(l)
    enddo

  !Second-derivative values
    temp(:)=0.0d0
    do ord=0,Poly_ord
      do l=1,4
        temp(l)=temp(l)+deg(ord,i)*Dx2_Poly(ord,Q(l,1)/2.0d0,dx(i))
      enddo
    enddo
    do l=1,4
      DG_Dx2_Poly_quad(l,i)=temp(l)
    enddo
  enddo

  !------------------------------------------------------------------------------------------------------------!
  !Activate the vairous Limiters
  !------------------------------------------------------------------------------------------------------------!    
  
  !TYPE 1. KXRCF
  !Prepare the KXRCF shock indicator by using the density variable
  do i=1,imax
    SI(i)=dabs(u_r(i-1)-u_l(i+0))
    Norm(i)=dsqrt((u_c(i)**2)*dx(i))

    KXRCF(i)=SI(i)/(Norm(i)*((dx(i)/2.0d0)**((Poly_ord+1)/2))+1d-40)

  !TYPE 2. YANG & YOON'S LIMITER    
  !local detectors  
  !  IS(-2)=dabs(deg(1,i-1)*dx(i)**1)**p+xi*dabs(deg(2,i-1)*dx(i)**2)**p
  !  IS(+2)=dabs(deg(1,i+1)*dx(i)**1)**p+xi*dabs(deg(2,i+1)*dx(i)**2)**p

    IS(-2)=dabs(((C_Dx1_u(i-1)/dx(i-1)**1)/1.0d0)*dx(i-1)**0)**p&
      &+xi*dabs(((C_Dx2_u(i-1)/dx(i-1)**2)/2.0d0)*dx(i-1)**0)**p

    IS(+2)=dabs(((C_Dx1_u(i+1)/dx(i+1)**1)/1.0d0)*dx(i+1)**0)**p&
      &+xi*dabs(((C_Dx2_u(i+1)/dx(i+1)**2)/2.0d0)*dx(i+1)**0)**p

  !Global detectors      
  !  IS(-1)=dabs((deg(1,i-1)+deg(1,i+0))*dx(i)**1)**p+xi*dabs((deg(2,i-1)+deg(2,i+0))*dx(i)**2)**p
  !  IS(+1)=dabs((deg(1,i+0)+deg(1,i+1))*dx(i)**1)**p+xi*dabs((deg(2,i+0)+deg(2,i+1))*dx(i)**2)**p

    IS(-1)=dabs((((C_Dx1_u(i-1)+C_Dx1_u(i+0))/dx(i-1)**1)/2.0d0)*dx(i-1)**0)**p&
      &+xi*dabs((((C_Dx2_u(i-1)+C_Dx2_u(i+0))/dx(i-1)**2)/4.0d0)*dx(i-1)**0)**p
      
    IS(+1)=dabs((((C_Dx1_u(i+0)+C_Dx1_u(i+1))/dx(i+1)**1)/2.0d0)*dx(i+1)**0)**p&
      &+xi*dabs((((C_Dx2_u(i+0)+C_Dx2_u(i+1))/dx(i+1)**2)/4.0d0)*dx(i+1)**0)**p

  !Set up the singularity detector scheme  
    THETA(i)=(1.0d0/2.0d0)*((min(IS(-1),IS(+1)))/(max(IS(-2),IS(+2))+1d-40))

  !Conduct the normalization  
    if(real(KXRCF(i)).gt.1.0d0)then  
  !  if(real(THETA(i)).gt.1.0d0)then
  !Singular cell  
      LIMITER_INDEX(i)=1.0d0
    else
  !Smooth cell    
      LIMITER_INDEX(i)=0.0d0
    endif
    
  !  write(*,*) i,real(x(i)),real(u_c(i)/dx(i)**0),real(C_Dx1_u(i)/dx(i)**1),real(C_Dx2_u(i)/dx(i)**2)
  !  write(*,*) i,real(KXRCF(i)),int(LIMITER_INDEX(i))
  !  write(*,*) i,real(THETA(i)),int(LIMITER_INDEX(i))
  enddo
  
  !call Get_Limiter_Index(i_bc,imax,x,LIMITER_INDEX)
  !stop

  !------------------------------------------------------------------------------------------------------------!
  
  !------------------------------------------------------------------------------------------------------------!
  !Activate the WENO type reconstruction
  !------------------------------------------------------------------------------------------------------------! 
  
  do i=1,imax
    if(LIMITER_INDEX(i).eq.0.0d0)then  
  !Do not limit
      goto 20
    else
  !Set up the linear weights for the WENO type reconstruction
      linear_coeff(-1)=1d-3
      linear_coeff(+1)=1d-3
      linear_coeff(+0)=1.0d0-linear_coeff(-1)-linear_coeff(+1)

  !Enact third-order WENO type limiter
  !Modify the Neumann DG polynomials by changing the cell averages
      do quad=1,4
        do l=-1,1
          recon_DG_Poly(l,quad,i+l)=(DG_Poly_quad(quad,i+l)-deg(0,i+l))+deg(0,i+0)
        enddo
      enddo

  !Compute the smoothness indicators for the high-order WENO reconstruction
  !------------------------------------------------------------------------------------------------------------!    
  
      !The local smoothness indicators
      beta(:,i)=0.0d0
      do l=-1,1
        beta(l,i)=beta(l,i)+&
                &+(dx(i)**0)*(DG_Dx1_Poly_quad(1,i+l)*Q(1,2)&
                            &+DG_Dx1_Poly_quad(2,i+l)*Q(2,2)&
                            &+DG_Dx1_Poly_quad(3,i+l)*Q(3,2)&
                            &+DG_Dx1_Poly_quad(4,i+l)*Q(4,2))**2

        beta(l,i)=beta(l,i)+&
                &+(dx(i)**2)*(DG_Dx2_Poly_quad(1,i+l)*Q(1,2)&
                            &+DG_Dx2_Poly_quad(2,i+l)*Q(2,2)&
                            &+DG_Dx2_Poly_quad(3,i+l)*Q(3,2)&
                            &+DG_Dx2_Poly_quad(4,i+l)*Q(4,2))**2
      
        beta(l,i)=(dx(i)**1)*beta(l,i)
  
      enddo
  !------------------------------------------------------------------------------------------------------------!

  !------------------------------------------------------------------------------------------------------------!    
  
  !The unnormalized nonlinear weights 
      alpha(:)=0.0d0
      do l=-1,1
        alpha(l)=linear_coeff(l)/((beta(l,i)+eps)**2)
      enddo
      alpha_sum=alpha(-1)+alpha(+0)+alpha(+1)
  
  !------------------------------------------------------------------------------------------------------------!
      
  !------------------------------------------------------------------------------------------------------------!    
  
  !The normalized nonlinear weights
      omega(:)=0.0d0
      do l=-1,1
        omega(l)=alpha(l)/alpha_sum
      enddo
  
  !------------------------------------------------------------------------------------------------------------!    

  !Compute the convex summation at each quadrature points to obtain the modified DG polynomials
      new_DG_Poly(:)=0.0d0
      do quad=1,4
        do l=-1,1
          new_DG_Poly(quad)=new_DG_Poly(quad)+omega(l)*recon_DG_Poly(l,quad,i+l)
        enddo
      enddo  
  
  !Evaluate degrees of freedom the moments by using the modified DG polynomials
      do ord=1,Poly_ord
        deg(ord,i)=new_DG_Poly(1)*Poly(ord,Q(1,1)/2.0d0,dx(i))*Q(1,2)&
                 &+new_DG_Poly(2)*Poly(ord,Q(2,1)/2.0d0,dx(i))*Q(2,2)&
                 &+new_DG_Poly(3)*Poly(ord,Q(3,1)/2.0d0,dx(i))*Q(3,2)&
                 &+new_DG_Poly(4)*Poly(ord,Q(4,1)/2.0d0,dx(i))*Q(4,2)

  !Product the mass matrix for the orthonormalization                 
        deg(ord,i)=M(ord,ord)*deg(ord,i)
      enddo
    endif
    20 continue
  enddo
  
  !------------------------------------------------------------------------------------------------------------!

  !Impose the boundary condition to degrees of freedom of the moments
  call Get_Boundary(Poly_ord,i_bc,imax,deg)
  
  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_Residual(TVB_M,Poly_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,maxdx,deg,F_Flux,LIMITER_INDEX)
  implicit none

  real(8),intent(in) :: TVB_M
  integer,intent(in) :: Poly_ord
  real(8),intent(in) :: Q(1:4,1:2) 
  real(8),intent(in) :: M(0:Poly_ord,0:Poly_ord)
  real(8),intent(in) :: d(0:1,1:4)
  real(8),intent(in) :: l_b(0:Poly_ord),r_b(0:Poly_ord)
  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: x(1-i_bc+1:imax+i_bc),dx(-i_bc+1:imax+i_bc),maxdx
  real(8),intent(inout) :: deg(0:Poly_ord,-i_bc+1:imax+i_bc)

  real(8),intent(out) :: F_Flux(0:Poly_ord,1:imax)
  real(8),intent(out) :: LIMITER_INDEX(1:imax)

  !Counters
  integer :: ord
  integer :: i

  !Flux function
  real(8),external :: Flux

  !Orthogonal non-polynomials 
  real(8),external :: Poly,Dx1_Poly

  !Volume integration variables
  real(8),dimension(1:4) :: quad_sum
  real(8),dimension(1:4,1:imax) :: quad_u

  !Volume integration variable
  real(8),dimension(0:Poly_ord,1:imax) :: F_Flux_V

  !Cell interface values
  real(8),dimension(-1:0) :: temp
  real(8),dimension(-1:imax) :: u_p
  real(8),dimension(0:imax+1) :: u_m
  
  !Flux splitting variable
  real(8),dimension(0:imax) :: hat_F

  !Boundary integration variable
  real(8),dimension(0:Poly_ord,1:imax) :: F_Flux_B

  real(8),dimension(1:imax) :: LIMITER_INDEX_1
  !-----------------------------------------Calculations have started------------------------------------------!

  !Calculate the local residual

  !============================================================================================================!
  !Volume integral part
  !============================================================================================================!

  !Compute the volume integral by using quadratures
  do i=1,imax
    quad_sum(:)=0.0d0
    do ord=0,Poly_ord
      quad_sum(1)=quad_sum(1)+deg(ord,i)*Poly(ord,Q(1,1)/2.0d0,dx(i))
      quad_sum(2)=quad_sum(2)+deg(ord,i)*Poly(ord,Q(2,1)/2.0d0,dx(i))
      quad_sum(3)=quad_sum(3)+deg(ord,i)*Poly(ord,Q(3,1)/2.0d0,dx(i))
      quad_sum(4)=quad_sum(4)+deg(ord,i)*Poly(ord,Q(4,1)/2.0d0,dx(i))
    enddo
    quad_u(1,i)=quad_sum(1)
    quad_u(2,i)=quad_sum(2)
    quad_u(3,i)=quad_sum(3)
    quad_u(4,i)=quad_sum(4)
  enddo
   
  do i=1,imax
    do ord=0,Poly_ord
      F_Flux_V(ord,i)=Flux(quad_u(1,i))*Dx1_Poly(ord,Q(1,1)/2.0d0,dx(i))*Q(1,2)&
                    &+Flux(quad_u(2,i))*Dx1_Poly(ord,Q(2,1)/2.0d0,dx(i))*Q(2,2)&
                    &+Flux(quad_u(3,i))*Dx1_Poly(ord,Q(3,1)/2.0d0,dx(i))*Q(3,2)&
                    &+Flux(quad_u(4,i))*Dx1_Poly(ord,Q(4,1)/2.0d0,dx(i))*Q(4,2)
    enddo
  enddo

  !============================================================================================================!
  !Boundary integral part
  !============================================================================================================!
  
  !Impose the boundary condition to degrees of freedom of the moments
  call Get_Boundary(Poly_ord,i_bc,imax,deg)

  !Apply limiting process at the cell interfaces
  call Get_Original_WENO_Limiter(TVB_M,Poly_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,maxdx,deg,LIMITER_INDEX_1)
  !call Get_Simple_WENO_Limiter(TVB_M,Poly_ord,Q,M,l_b,r_b,i_bc,imax,x,dx,maxdx,deg,LIMITER_INDEX_1)
  
  !do i=0,imax+1
  !  write(*,*) i,deg(0,i)
  !enddo
  
  LIMITER_INDEX(:)=LIMITER_INDEX_1(:)

  !Set up the reconstructed boundary contributions
  do i=0,imax+1
    temp(:)=0.0d0
    do ord=0,Poly_ord
      temp(-1)=temp(-1)+deg(ord,i)*l_b(ord)
      temp(+0)=temp(+0)+deg(ord,i)*r_b(ord)
    enddo
    u_p(i-1)=temp(-1)
    u_m(i+0)=temp(+0)
  enddo

  !Take the appropriate numerical trace 
  !As the averaged monotone flux, apply Lax-Friedrich (LF) flux
  do i=0,imax
    hat_F(i)=Flux(u_m(i))+Flux(u_p(i))+1.0d0*(u_m(i)-u_p(i))

    hat_F(i)=hat_F(i)/2.0d0
  enddo

  !Compute the boundary integral by using the fundamental theorem of calculus
  do i=1,imax
    do ord=0,Poly_ord  
      F_Flux_B(ord,i)=hat_F(i-1)*l_b(ord)-hat_F(i+0)*r_b(ord)
    enddo
  enddo

  !============================================================================================================!
  !============================================================================================================!

  !!Compute local residuals by using calculated integral values
  do i=1,imax
    do ord=0,Poly_ord
      F_Flux(ord,i)=M(ord,ord)*(F_Flux_B(ord,i)+F_Flux_V(ord,i))/dx(i)
    enddo
  enddo
  
  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_RK(RK_ord,TVB_M,Poly_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,maxdx,dt,old_deg,new_deg,LIMITER_FLAG)
  implicit none
  
  integer,intent(in) :: RK_ord
  real(8),intent(in) :: TVB_M
  integer,intent(in) :: Poly_ord
  real(8),intent(in) :: Q(1:4,1:2) 
  real(8),intent(in) :: M(0:Poly_ord,0:Poly_ord)
  real(8),intent(in) :: d(0:1,1:4)
  real(8),intent(in) :: l_b(0:Poly_ord),r_b(0:Poly_ord)
  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: x(-i_bc+1:imax+i_bc),dx(-i_bc+1:imax+i_bc),maxdx
  real(8),intent(in) :: dt
  real(8),intent(inout) :: old_deg(0:Poly_ord,-i_bc+1:imax+i_bc)

  real(8),intent(out) :: new_deg(0:Poly_ord,1:imax)
  real(8),intent(out) :: LIMITER_FLAG(1:imax)

  !Counters
  integer :: ord
  integer :: i

  !Calculation variables
  real(8),dimension(0:Poly_ord,-i_bc+1:imax+i_bc,1:3) :: temp_deg
  real(8),dimension(0:Poly_ord,1:imax) :: F_Flux
  real(8),dimension(1:4,1:imax) :: LIMITER_INDEX
  !-----------------------------------------Calculations have started------------------------------------------!
   
  if(RK_ord.eq.1)then

  !1st-order TVD RK solver
  !============================================================================================================!
  !Stage 1
  !============================================================================================================!

  !Compute the residual
    call Get_Residual(TVB_M,Poly_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,maxdx,old_deg,F_Flux,LIMITER_INDEX(1,:))
    do i=1,imax
      do ord=0,Poly_ord
        new_deg(ord,i)=+1.0d0*old_deg(ord,i)&
                      &+1.0d0*dt*F_Flux(ord,i)
      enddo
    enddo

    LIMITER_FLAG(:)=LIMITER_INDEX(1,:)

  elseif(RK_ord.eq.3)then

  !3rd-order TVD RK solver
  !============================================================================================================!
  !Stage 1
  !============================================================================================================!

  !Compute the residual
    call Get_Residual(TVB_M,Poly_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,maxdx,old_deg,F_Flux,LIMITER_INDEX(1,:))
    do i=1,imax
      do ord=0,Poly_ord
        temp_deg(ord,i,1)=+1.0d0*old_deg(ord,i)&
                         &+1.0d0*dt*F_Flux(ord,i)
      enddo
    enddo

  !============================================================================================================!
  !Stage 2
  !============================================================================================================!

  !Compute the residual
    call Get_Residual(TVB_M,Poly_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,maxdx,temp_deg(:,:,1),F_Flux,LIMITER_INDEX(2,:))
    do i=1,imax
      do ord=0,Poly_ord
        temp_deg(ord,i,2)=+3.0d0*old_deg(ord,i)&
                         &+1.0d0*temp_deg(ord,i,1)&
                         &+1.0d0*dt*F_Flux(ord,i)

        temp_deg(ord,i,2)=temp_deg(ord,i,2)/4.0d0
      enddo
    enddo

  !============================================================================================================!
  !Stage 3
  !============================================================================================================!

  !Compute the residual
    call Get_Residual(TVB_M,Poly_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,maxdx,temp_deg(:,:,2),F_Flux,LIMITER_INDEX(3,:))
    do i=1,imax
      do ord=0,Poly_ord
        new_deg(ord,i)=+1.0d0*old_deg(ord,i)&
                      &+2.0d0*temp_deg(ord,i,2)&
                      &+2.0d0*dt*F_Flux(ord,i)

        new_deg(ord,i)=new_deg(ord,i)/3.0d0
      enddo
    enddo

    LIMITER_FLAG(:)=LIMITER_INDEX(3,:)

  elseif(RK_ord.eq.4)then

  !4th-order non-TVD RK solver
  !============================================================================================================!
  !Stage 1
  !============================================================================================================!
  
  !Compute the residual
    call Get_Residual(TVB_M,Poly_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,maxdx,old_deg,F_Flux,LIMITER_INDEX(1,:))
    do i=1,imax
      do ord=0,Poly_ord
        temp_deg(ord,i,1)=+1.0d0*old_deg(ord,i)&
                         &+0.5d0*dt*F_Flux(ord,i)
      enddo
    enddo

  !============================================================================================================!
  !Stage 2
  !============================================================================================================!
  
  !Compute the residual
    call Get_Residual(TVB_M,Poly_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,maxdx,temp_deg(:,:,1),F_Flux,LIMITER_INDEX(2,:))
    do i=1,imax
      do ord=0,Poly_ord
        temp_deg(ord,i,2)=+1.0d0*old_deg(ord,i)&
                         &+(dt/2.0d0)*F_Flux(ord,i)
      enddo
    enddo

  !============================================================================================================!
  !Stage 3
  !============================================================================================================!

  !Compute the residual
    call Get_Residual(TVB_M,Poly_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,maxdx,temp_deg(:,:,2),F_Flux,LIMITER_INDEX(3,:)) 
    do i=1,imax
      do ord=0,Poly_ord
        temp_deg(ord,i,3)=+1.0d0*old_deg(ord,i)&
                         &+1.0d0*dt*F_Flux(ord,i)
      enddo 
    enddo
  
  !============================================================================================================!
  !Stage 4
  !============================================================================================================!

  !Compute the residual
    call Get_Residual(TVB_M,Poly_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,maxdx,temp_deg(:,:,3),F_Flux,LIMITER_INDEX(4,:))
    do i=1,imax
      do ord=0,Poly_ord
        new_deg(ord,i)=-1.0d0*old_deg(ord,i)&
                      &+1.0d0*temp_deg(ord,i,1)&
                      &+2.0d0*temp_deg(ord,i,2)&
                      &+1.0d0*temp_deg(ord,i,3)&
                      &+(dt/2.0d0)*F_Flux(ord,i)

        new_deg(ord,i)=new_deg(ord,i)/3.0d0
      enddo
    enddo

    LIMITER_FLAG(:)=LIMITER_INDEX(4,:)

  endif
  
  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_L_2_Error(i_bc,imax,dx,exact,numerical)
  implicit none

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: dx(-i_bc+1:imax+i_bc)
  real(8),intent(in) :: exact(1:imax),numerical(1:imax)

  !Counter
  integer :: i

  !Error calculation variables
  real(8),dimension(1:imax) :: e
  !-----------------------------------------Calculations have started------------------------------------------!

  !Compute square differences
  do i=1,imax
    e(i)=(exact(i)-numerical(i))**2
    e(i)=e(i)*dx(i)
  enddo

  !Compute the square root average sum
  write(*,*) "L_2 error:",real(dsqrt(sum(e)))

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_L_Infty_Error(imax,exact,numerical)
  implicit none

  integer,intent(in) :: imax
  real(8),intent(in) :: exact(1:imax),numerical(1:imax)

  !Counters
  integer :: i

  !Error calculation variables
  real(8),dimension(1:imax) :: e
  !-----------------------------------------Calculations have started------------------------------------------!

  !Compute absolute differences
  do i=1,imax
    e(i)=dabs(exact(i)-numerical(i))
  enddo

  !Take the maximum value
  write(*,*) "L_infty error:",real(maxval(e))
  write(*,*) "Maximum error location:",maxloc(e)

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_Initial_Condition(i_bc,imax,x,u)
  implicit none

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: x(-i_bc+1:imax+i_bc)
  real(8),intent(in) :: u(1:imax)

  !Counter
  integer :: i
  !-----------------------------------------Calculations have started------------------------------------------!

  open(10,file='Initial_Condition.plt')
  write(10,*)'zone T="file"', 'I=',imax
  
  do i=1,imax
    write(10,*) real(x(i)),real(u(i))
  enddo

  close(10)

  return
  !-----------------------------------------Calculations have finished----------------------------------------!
endsubroutine

subroutine Get_Limiter_Index(i_bc,imax,x,u)
  implicit none

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: x(-i_bc+1:imax+i_bc)
  real(8),intent(in) :: u(1:imax)

  !Counter
  integer :: i
  !-----------------------------------------Calculations have started------------------------------------------!
  
  open(20,file='Limiter_Index.plt')
  write(20,*)'zone T="file"', 'I=',imax

  do i=1,imax
    write(20,*) real(x(i)),real(u(i))
  enddo

  close(20)

  return
  !-----------------------------------------Calculations have finished----------------------------------------!
endsubroutine

subroutine Get_Limiter_Flag(i_bc,imax,x,u)
  implicit none

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: x(-i_bc+1:imax+i_bc)
  real(8),intent(in) :: u(1:imax)

  !Counter
  integer :: i
  !-----------------------------------------Calculations have started------------------------------------------!

  open(11,file='Limiter_Flag.plt')
  write(11,*)'zone T="file"', 'I=',imax
  
  do i=1,imax
    write(11,*) real(x(i)),real(u(i))
  enddo

  close(11)

  return
  !-----------------------------------------Calculations have finished----------------------------------------!
endsubroutine

subroutine Get_Zero_Moment(i_bc,imax,x,u)
  implicit none

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: x(-i_bc+1:imax+i_bc)
  real(8),intent(in) :: u(1:imax)

  !Counter
  integer :: i
  !-----------------------------------------Calculations have started------------------------------------------!

  open(0,file='Zero_Moment.plt')
  write(0,*)'zone T="file"', 'I=',imax
  
  do i=1,imax
    write(0,*) real(x(i)),real(u(i))
  enddo

  close(0)

  return
  !-----------------------------------------Calculations have finished----------------------------------------!
endsubroutine

subroutine Get_First_Moment(i_bc,imax,x,u)
  implicit none

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: x(-i_bc+1:imax+i_bc)
  real(8),intent(in) :: u(1:imax)

  !Counter
  integer :: i
  !-----------------------------------------Calculations have started------------------------------------------!
  
  open(1,file='First_Moment.plt')
  write(1,*)'zone T="file"', 'I=',imax
  
  do i=1,imax
    write(1,*) real(x(i)),real(u(i))
  enddo

  close(1)

  return
  !-----------------------------------------Calculations have finished----------------------------------------!
endsubroutine

subroutine Get_Second_Moment(i_bc,imax,x,u)
  implicit none

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: x(-i_bc+1:imax+i_bc)
  real(8),intent(in) :: u(1:imax)

  !Counter
  integer :: i
  !-----------------------------------------Calculations have started------------------------------------------!
  
  open(2,file='Second_Moment.plt')
  write(2,*)'zone T="file"', 'I=',imax

  do i=1,imax
    write(2,*) real(x(i)),real(u(i))
  enddo

  close(2)

  return
  !-----------------------------------------Calculations have finished----------------------------------------!
endsubroutine

subroutine Get_DG_Solution(i_bc,imax,x,u)
  implicit none

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: x(-i_bc+1:imax+i_bc)
  real(8),intent(in) :: u(1:imax)

  !Counter
  integer :: i
  !-----------------------------------------Calculations have started------------------------------------------!

  open(3,file='DG_Soluttion.plt')
  write(3,*)'zone T="file"', 'I=',imax
  
  do i=1,imax
    write(3,*) real(x(i)),real(u(i))
  enddo

  close(3)

  return
  !-----------------------------------------Calculations have finished----------------------------------------!
endsubroutine

real(8) function Ini_u(x)
  implicit none

  real(8),intent(in) :: x

  !Useful constant
  real(8),parameter :: pi=4.0d0*datan(1.0d0)
  !-----------------------------------------Calculations have started------------------------------------------!
  
  !Set up the initial function
  !Smooth function. DOMAIN: [-1.0d0/1.0d0,+1.0d0/1.0d0]
  !Ini_u=dsin(pi*x)**1
  Ini_u=dcos(pi*x)**4
  !Ini_u=dsin(pi*x-dsin(pi*x)/pi)
  
  !Nonsmooth function
  !Square wave. DOMAIN: [-1.0d0/1.0d0,+1.0d0/1.0d0]
  !if((x.gt.-0.5d0).and.(x.lt.0.5d0))then
  !  Ini_u=1.0d0
  !else
  !  Ini_u=0.0d0
  !endif

  !Piecewise smooth function
  !TYPE 1. DOMAIN: [-1.0d0/1.0d0,+1.0d0/1.0d0]
  !Ini_u=max(0.5d0-1.0d0*dabs(x)**3,0.0d0)
  
  !TYPE 2. DOMAIN: [-1.0d0/2.0d0,+1.0d0/2.0d0]
  !if((x.ge.-1.0d0).and.(x.le.-1.0d0/3.0d0))then
  !  Ini_u=-x*dsin(((3.0d0*pi)/2.0d0)*x**2)
  !elseif((x.gt.-1.0d0/3.0d0).and.(x.le.1.0d0/3.0d0))then
  !  Ini_u=dabs(dsin(2.0d0*pi*x))
  !else
  !  Ini_u=2.0d0*x-1.0d0-(1.0d0/6.0d0)*dsin(3.0d0*pi*x)
  !endif

  !TYPE 3. DOMAIN: [-1.0d0/2.0d0,+1.0d0/2.0d0]
  !if((x.ge.-1.0d0).and.(x.le.-1.0d0/3.0d0))then
  !  Ini_u=-x*dsin(((3.0d0*pi)/2.0d0)*x**2)-1.0d0/6.0d0+dsqrt(3.0d0)/2.0d0
  !elseif((x.gt.-1.0d0/3.0d0).and.(x.le.1.0d0/3.0d0))then
  !  Ini_u=dabs(dsin(2.0d0*pi*x))
  !else
  !  Ini_u=2.0d0*x-2.0d0/3.0d0+dsqrt(3.0d0)/2.0d0-(1.0d0/6.0d0)*dsin(3.0d0*pi*x)
  !endif

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endfunction

real(8) function Ext_u(x,time_out)
  implicit none

  real(8),intent(in) :: x
  real(8),intent(in) :: time_out

  !Useful constant
  real(8),parameter :: pi=4.0d0*datan(1.0d0)
  !-----------------------------------------Calculations have started------------------------------------------!
  
  !Set up the exact solution
  !Smooth functions
  !Ext_u=dsin(pi*(x-time_out))**1
  Ext_u=dcos(pi*(x-time_out))**4
  !Ext_u=dsin(pi*(x-time_out)-dsin(pi*(x-time_out))/pi)

  !Nonsmooth function
  !Square wave
  !if((x.gt.-0.5d0).and.(x.lt.0.5d0))then
  !  Ext_u=1.0d0
  !else
  !  Ext_u=0.0d0
  !endif

  !Piecewise smooth function
  !TYPE 1.
  !Ext_u=max(0.5d0-1.0d0*dabs(x)**1,0.0d0)

  !TYPE 2.
  !if((x.ge.-1.0d0).and.(x.le.-1.0d0/3.0d0))then
  !  Ext_u=-x*dsin(((3.0d0*pi)/2.0d0)*x**2)
  !elseif((x.gt.-1.0d0/3.0d0).and.(x.le.1.0d0/3.0d0))then
  !  Ext_u=dabs(dsin(2.0d0*pi*x))
  !else
  !  Ext_u=2.0d0*x-1.0d0-(1.0d0/6.0d0)*dsin(3.0d0*pi*x)
  !endif

  !TYPE 3.
  !if((x.ge.-1.0d0).and.(x.le.-1.0d0/3.0d0))then
  !  Ini_u=-x*dsin(((3.0d0*pi)/2.0d0)*x**2)-1.0d0/6.0d0+dsqrt(3.0d0)/2.0d0
  !elseif((x.gt.-1.0d0/3.0d0).and.(x.le.1.0d0/3.0d0))then
  !  Ini_u=dabs(dsin(2.0d0*pi*x))
  !else
  !  Ini_u=2.0d0*x-2.0d0/3.0d0+dsqrt(3.0d0)/2.0d0-(1.0d0/6.0d0)*dsin(3.0d0*pi*x)
  !endif

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endfunction

real(8) function Flux(u)
  implicit none

  real(8),intent(in) :: u
  !-----------------------------------------Calculations have started------------------------------------------!
  
  !Set up the flux
  Flux=u

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endfunction

real(8) function Poly(ord,x,dx)
  implicit none

  integer,intent(in) :: ord
  real(8),intent(in) :: x,dx

  real(8) :: lambda
  !-----------------------------------------Calculations have started------------------------------------------!
   
  !Set up scaled-orthoronal bases
  !Legendre polynomials
  if(ord.eq.0)then
	  Poly=1.0d0
  elseif(ord.eq.1)then
    Poly=x
  elseif(ord.eq.2)then
    Poly=x**2-(1.0d0/12.0d0)
  elseif(ord.eq.3)then
    Poly=x**3-(3.0d0/20.0d0)*x
  endif

  !Determine shpae parameter
  !lambda=dx**2
  
  !Exponential polynomials
  !if(ord.eq.0)then
	!  Poly=1.0d0
  !elseif(ord.eq.1)then
  !  Poly=(lambda*x)**5+(5.0d0*(lambda*x)**3)+(15.0d0*(lambda*x))/2.0d0
  !elseif(ord.eq.2)then
  !  Poly=(lambda*x)**4+(3.0d0*(lambda*x)**2)-(3.0d0/2.0d0)*((lambda**4)/120.0d0+(lambda**2)/6.0d0)
  !endif

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endfunction

real(8) function Dx1_Poly(ord,x,dx)
  implicit none
  
  integer,intent(in) :: ord
  real(8),intent(in) :: x,dx

  real(8) :: lambda
  !-----------------------------------------Calculations have started------------------------------------------!
  
  !Set up derivatives of scaled-orthoronal bases
  !Legendre polynomials
  if(ord.eq.0)then
	  Dx1_Poly=0.0d0
  elseif(ord.eq.1)then
    Dx1_Poly=1.0d0
  elseif(ord.eq.2)then
    Dx1_Poly=2.0d0*x
  elseif(ord.eq.3)then
    Dx1_Poly=3.0d0*x**2-(3.0d0/20.0d0)
  endif

  !Determine shpae parameter
  !lambda=dx**2

  !Exponential polynomials
  !if(ord.eq.0)then
	!  Dx1_Poly=0.0d0
  !elseif(ord.eq.1)then
  !  Dx1_Poly=5.0d0*(lambda**5)*(x**4)+15.0d0*(lambda**3)*(x**2)+(15.0d0*lambda)/2.0d0
  !elseif(ord.eq.2)then
  !  Dx1_Poly=4.0d0*(lambda**4)*(x**3)+6.0d0*(lambda**2)*x
  !endif

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endfunction

real(8) function Dx2_Poly(ord,x,dx)
  implicit none
  
  integer,intent(in) :: ord
  real(8),intent(in) :: x,dx

  real(8) :: lambda
  !-----------------------------------------Calculations have started------------------------------------------!
  
  !Set up derivatives of scaled-orthoronal bases
  !Legendre polynomials
  if(ord.eq.0)then
	  Dx2_Poly=0.0d0
  elseif(ord.eq.1)then
    Dx2_Poly=0.0d0
  elseif(ord.eq.2)then
    Dx2_Poly=2.0d0
  elseif(ord.eq.3)then
    Dx2_Poly=6.0d0*x
  endif

  !Determine shpae parameter
  !lambda=dx**2

  !Exponential polynomials
  !if(ord.eq.0)then
	!  Dx2_Poly=0.0d0
  !elseif(ord.eq.1)then
  !  Dx2_Poly=20.0d0*(lambda**5)*(x**3)+30.0d0*(lambda**3)*x
  !elseif(ord.eq.2)then
  !  Dx2_Poly=12.0d0*(lambda**4)*(x**2)+6.0d0*(lambda**2)
  !endif

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endfunction

real(8) function Dx3_Poly(ord,x,dx)
  implicit none
  
  integer,intent(in) :: ord
  real(8),intent(in) :: x,dx

  real(8) :: lambda
  !-----------------------------------------Calculations have started------------------------------------------!
  
  !Set up derivatives of scaled-orthoronal bases
  !Legendre polynomials
  if(ord.eq.0)then
	  Dx3_Poly=0.0d0
  elseif(ord.eq.1)then
    Dx3_Poly=0.0d0
  elseif(ord.eq.2)then
    Dx3_Poly=0.0d0
  elseif(ord.eq.3)then
    Dx3_Poly=6.0d0
  endif

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endfunction

real(8) function Lagrange_Coefficients(l,k_l,k_r,i_bc,imax,i,x,dx,x_p)
  implicit none

  integer,intent(in) :: l !Cell counter
  integer,intent(in) :: k_l !Minimal index of the left cell
  integer,intent(in) :: k_r !Maximal index of the right cell
  integer,intent(in) :: i_bc,imax,i
  real(8),intent(in) :: x(-i_bc+1:imax+i_bc)
  real(8),intent(in) :: dx(-i_bc+1:imax+i_bc) 
  real(8),intent(in) :: x_p !Present location
  
  !Counters
  integer :: j,k,s

  !Calculation variables
  real(8) :: x_n !Neighbor cell
  real(8) :: c1,c2,c3
  !-----------------------------------------Calculations have started------------------------------------------!
  
  !Initialization
  Lagrange_Coefficients=0.0d0

  !Set up the neighbor cell
  x_n=x(i)+x_p*dx(i)

  !Compute coefficients for Lagrange interpolation
	do k=l,k_r
	  c1=1.0d0
	  do j=k_l-1,k_r
	    if(j.ne.k)then
        c1=c1*(x(k+i)+(dx(k+i)-dx(j+i))/2.0d0-x(j+i))
      endif
    enddo
	  c2=0.0d0
	  do j=k_l-1,k_r
	    if(j.ne.k) then
	      c3=1.0d0
	      do s=k_l-1,k_r
	        if((s.ne.j).and.(s.ne.k))then
            c3=(x_n-x(s+i)-dx(s+i)/2.0d0)*c3
          endif 
        enddo
	      c2=c2+c3
      endif
	  enddo
	  Lagrange_Coefficients=Lagrange_Coefficients+c2/c1
	enddo  
	Lagrange_Coefficients=Lagrange_Coefficients*dx(i+l)
  
  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endfunction