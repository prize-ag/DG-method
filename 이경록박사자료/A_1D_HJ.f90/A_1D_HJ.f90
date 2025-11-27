program HJ_WENO_1D
  implicit none

  !======================================================================================================================!
  !Counters
  integer :: i !Spatial(x-axis) index
  integer :: n !Time index

  !Useful constants
  real(8),parameter :: pi=4.0d0*datan(1.0d0)
  real(8),parameter :: CFL=0.4d0 !Courant–Friedrichs–Lewy number

  !Spatial variables
  integer,parameter :: p=8
  integer,parameter :: imax=10*(2**p) !Number of spatial nodes

  integer,parameter :: imax_REF=100*(2**8) !Number of spatial nodes
  
  integer,parameter :: i_bc=5 !Ghost boundary nodes for WENO schemes

  real(8),parameter :: x_left=-1.0d0 !Left boundary point on x-axis
  real(8),parameter :: x_right=1.0d0 !Right boundary point on x-axis

  real(8),parameter :: x_length=dabs(x_right-x_left) !Length of domain

  real(8),parameter :: dx=x_length/dble(imax) !Cell size of domain

  real(8),parameter :: dx_REF=x_length/dble(imax_REF) !Cell size of domain
  
  integer,parameter :: ini=1 !Initial condition controller

  !Time variables
  integer,parameter :: nmax=543000 !Number of the marching step
  real(8),parameter :: time_out=0.5d0/(pi**2) !Final time of the solution
  real(8) :: dt !Time step size
  real(8) :: flag_0,flag_1

  !Calculation variables
  real(8),dimension(-i_bc+0:imax+i_bc) :: old_phi
  real(8),dimension(0:imax) :: new_phi
  real(8),dimension(-i_bc+0:imax_REF+i_bc) :: old_phi_REF
  real(8),dimension(0:imax_REF) :: new_phi_REF
  !======================================================================================================================!
  
  write(*,*) "================================================================================================="
  write(*,*) "                               HAMILTON-JACOBI EQUATION. (NONLINEAR)                             "
  write(*,*) "================================================================================================="
  write(*,*) "                                    Calculations have started.                                   "
  write(*,*) "================================================================================================="
  
  !Set up the initial conditions
  flag_0=0.0d0
  flag_1=flag_0
  call Get_Initial(ini,i_bc,imax,dx,x_left,old_phi(0:imax))

  !!Import the reference solution
  open(0,file='REF_SOL.txt')

  do i=0,imax_REF
    read(0,*) new_phi_REF(i)
  enddo
  
  close(0)

  !============================================================================================================!
  !=======================================Time marching procedure(Start)=======================================!
  !============================================================================================================!

  do n=1,nmax
  !Compute the time step size
    call Get_Timestep(i_bc,imax,dx,flag_1,time_out,CFL,old_phi,dt)
     
  !Update the variable by using RK solver
    call Get_RK(i_bc,imax,dx,dt,old_phi,new_phi)
    
  !Update the real-time
    flag_1=flag_1+dt
    flag_0=flag_1

  !Receive the updated data
    do i=0,imax
      old_phi(i)=new_phi(i)
    enddo

    if(flag_0.eq.time_out)then
  !    call Get_Results(imax,dx,x_left,new_phi)
      goto 10
    endif
  enddo
  10 continue 

  !============================================================================================================!
  !======================================Time marching procedure(Finish)=======================================!
  !============================================================================================================!

  
  !write(*,*) "================================================================================================="
  !write(*,*) "                               Two phase calculations have started.                              "
  !write(*,*) "================================================================================================="

  !Set up the initial conditions
  !flag_0=0.0d0
  !flag_1=flag_0
  !call Get_Initial(ini,i_bc,imax_REF,dx_REF,x_left,old_phi_REF(0:imax_REF))

  !============================================================================================================!
  !=======================================Time marching procedure(Start)=======================================!
  !============================================================================================================!

  !do n=1,nmax
  !Compute the time step size
  !  call Get_Timestep_REF(i_bc,imax_REF,dx_REF,flag_1,time_out,CFL,old_phi_REF,dt)

  !Update the variable by using RK solver
  !  call Get_RK_REF(i_bc,imax_REF,dx_REF,dt,old_phi_REF,new_phi_REF)  
  
  !Update the real-time
  !  flag_1=flag_1+dt
  !  flag_0=flag_1

  !Receive the updated data
  !  do i=0,imax_REF
  !    old_phi_REF(i)=new_phi_REF(i)
  !  enddo  
  
  !  if(flag_0.eq.time_out)then
  !    call Get_Results_REF(imax_REF,dx_REF,x_left,new_phi_REF)
  !    goto 20
  !  endif
  !enddo
  !20 continue

  !============================================================================================================!
  !======================================Time marching procedure(Finish)=======================================!
  !============================================================================================================!

  !Compute various errors
  write(*,*) "=========================================Error results==========================================="
  call Get_L_1_Error(imax,imax_REF,p,dx,new_phi_REF,new_phi)
  call Get_L_Infty_Error(imax,imax_REF,p,dx,new_phi_REF,new_phi)
  write(*,*) "================================================================================================="
  write(*,*) " "

  write(*,*) "================================================================================================="
  write(*,*) "                                   Calculations have finished.                                   "
  write(*,*) "================================================================================================="

  return
endprogram

!============================================================================================================!

!============================================================================================================!

subroutine Get_Initial(ini,i_bc,imax,dx,x_left,phi)
  implicit none

  integer,intent(in) :: ini
  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: dx,x_left

  real(8),intent(out) :: phi(0:imax)

  !Counter
  integer :: i

  !Useful constant
  real(8),parameter :: pi=4.0d0*datan(1.0d0)
  !-----------------------------------------Calculations have started-----------------------------------------!

  if(ini.eq.1)then
    do i=0,imax
      phi(i)=-dcos(pi*(x_left+dble(i)*dx))
    enddo
  endif

  return
  !-----------------------------------------Calculations have finished----------------------------------------!
endsubroutine

subroutine Get_Boundary(i_bc,imax,phi)
  implicit none

  integer,intent(in) :: i_bc,imax

  real(8),intent(inout) :: phi(-i_bc+0:imax+i_bc)

  !Counter
  integer :: i
  !-----------------------------------------Calculations have started------------------------------------------!

  !Apply periodic condition
  do i=1,i_bc
    phi(-i)=phi(imax-i)
    phi(imax+i)=phi(i)
  enddo

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

!============================================================================================================!

!============================================================================================================!

subroutine Get_Timestep(i_bc,imax,dx,flag_1,time_out,CFL,phi,dt)
  implicit none

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: dx
  real(8),intent(in) :: flag_1,time_out,CFL
  real(8),intent(inout) :: phi(-i_bc+0:imax+i_bc)

  real(8),intent(out) :: dt

  !Counter
  integer :: i

  !Calcuation variables
  real(8),dimension(0:imax) :: F_phi_dx,B_phi_dx
  real(8) :: max_DH !Artificial Hamiltonian speed
  real(8) :: S_max !Maximal wave speed
  !-----------------------------------------Calculations have started------------------------------------------!
  
  !Impose the boundary condition
  call Get_Boundary(i_bc,imax,phi)
 
  !Compute the artificial Hamiltonian speed
  S_max=0.0d0
  do i=0,imax
    F_phi_dx(i)=(phi(i+1)-phi(i+0))/dx
    B_phi_dx(i)=(phi(i+0)-phi(i-1))/dx
    max_DH=max(dabs(dsin(F_phi_dx(i)+1.0d0)),dabs(dsin((B_phi_dx(i))+1.0d0)))
    S_max=max(S_max,max_DH)
  enddo

  !Determine the time step size
  dt=(CFL*dx)/S_max

  !Restriction
  if((flag_1+dt).ge.time_out)then
    dt=(time_out-flag_1)
  endif

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_Hat_Hx(i_bc,imax,dx,phi,hat_H)
  implicit none

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: dx
  real(8),intent(inout) :: phi(-i_bc+0:imax+i_bc)

  real(8),intent(out) :: hat_H(0:imax)

  !Counter
  integer :: i

  !Reconstructed variable
  real(8),dimension(0:imax,1:2) :: WENO

  !Canculation variables
  real(8),dimension(0:imax) :: P_phi_dx,M_phi_dx
  real(8),dimension(0:imax) :: L_max_alpha
  real(8) :: G_max_alpha
  !-----------------------------------------Calculations have started------------------------------------------!

  !Impose the boundary condition
  call Get_Boundary(i_bc,imax,phi)

  !Apply the WENO scheme
  call Get_WENO(i_bc,imax,dx,phi,WENO)

  !Seperate reconstructed arraies into two parts
  do i=0,imax
    M_phi_dx(i)=WENO(i,1)
    P_phi_dx(i)=WENO(i,2)
  enddo

  !Compute alpha in Lax-Friedrich flux
  !In this case, |H'(phi_dx)|=|dsin(phi_dx+1.0d0)|
  do i=0,imax
    L_max_alpha(i)=max(dabs(dsin(M_phi_dx(i)+1.0d0)),dabs(dsin((P_phi_dx(i))+1.0d0)))
  enddo
  G_max_alpha=maxval(L_max_alpha)

  !Apply Lax-Friedrich flux to the reconstructed Hamiltonian
  do i=0,imax
    hat_H(i)=-dcos((P_phi_dx(i)+M_phi_dx(i))/2.0d0+1.0d0)-G_max_alpha*(P_phi_dx(i)-M_phi_dx(i))

    hat_H(i)=-hat_H(i)/1.0d0
  enddo

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_WENO(i_bc,imax,dx,phi,WENO)
  implicit none

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: dx
  real(8),intent(in) :: phi(-i_bc+0:imax+i_bc)

  real(8),intent(out) :: WENO(0:imax,1:2)

  !Counters
  integer :: i,i0,i1
  integer :: l

  !Derivatives of the flux
  real(8) :: D_phi_1,D_phi_2,D_phi_3

  !Shape parameter
  real(8) :: lambda_dx,lambda2_dx2

  !Local coefficients
  real(8),dimension(0:imax,1:2) :: a_0,a_1,aL_0,aL_1 
  real(8),dimension(0:imax,1:2) :: b_1,b_2,bR_1,bR_2

  !Ideal weights
  real(8),dimension(0:imax,1:2) :: d_0,d_1

  !Calculation variables
  real(8),dimension(-i_bc+1:imax+i_bc) :: phi_dx
  real(8),dimension(0:imax,1:2,1:5) :: trans_phi_dx

  !Several parameters for WENO shceme
  real(8),external :: B3_Spline
  real(8) :: temp,delta
  real(8) :: gamma_0,gamma_1,eta
  real(8) :: exp_0,exp_1,exp_2
  
  !Main WENO reconstruction variables
  integer,parameter :: n=2 !Degenerate number: 10
  real(8) :: eps
  real(8) :: xi
  real(8) :: tau
  real(8) :: beta_0,beta_1
  real(8) :: mu
  real(8) :: alpha_0,alpha_1
  real(8) :: omega_0,omega_1
  !---------------------------------------------------Calculations have started----------------------------------------------------!

  !Set up some free parameters
  eps=dx**n !dx**n !1d-40
  xi=eps

  !Compute the first derivative of phi using the first-order difference
  do i=-i_bc+1,imax+i_bc
    phi_dx(i)=(phi(i+0)-phi(i-1))/dx
  enddo

  !Transform computed first derivative of phi in each directions
  do l=1,5
    i0=l-3
    i1=3-l
    do i=0,imax
  !Upwind part [j-1,j+0,j+1]    
      trans_phi_dx(i,1,l)=phi_dx(i+i0+0) 
  !Downwind part [j+2,j+1,j+0]
      trans_phi_dx(i,2,l)=phi_dx(i+i1+1) 
    enddo
  enddo
  
  do l=1,2
    do i=0,imax
  
  !============================================================================================================!

  !First derivative value of phi
      D_phi_1=-(12.0d0/30.0d0)*trans_phi_dx(i,l,2)&
             &-(12.0d0/60.0d0)*trans_phi_dx(i,l,3)&
             &+(12.0d0/20.0d0)*trans_phi_dx(i,l,4)

  !Third derivative value of phi
      D_phi_3=-1.0d0*trans_phi_dx(i,l,2)&
             &+3.0d0*trans_phi_dx(i,l,3)&
             &-3.0d0*trans_phi_dx(i,l,4)&
             &+1.0d0*trans_phi_dx(i,l,5)
   
  !Shape parameter [j-1,j+0,j+1]
      lambda2_dx2=+D_phi_3/(D_phi_1+sign(dx**4,D_phi_1))    

  !Local coefficients
      a_0(i,l)=-1.0d0/2.0d0+lambda2_dx2/12.0d0
      a_1(i,l)=+3.0d0/2.0d0-lambda2_dx2/12.0d0

      b_1(i,l)=+1.0d0/2.0d0+lambda2_dx2/12.0d0
      b_2(i,l)=+1.0d0/2.0d0-lambda2_dx2/12.0d0
              
  !Ideal weights
      d_0(i,l)=+1.0d0/3.0d0-lambda2_dx2/90.0d0
      d_1(i,l)=+2.0d0/3.0d0+lambda2_dx2/90.0d0
  
  !============================================================================================================!

  !First derivative value of phi
      D_phi_1=-(1.0d0)*trans_phi_dx(i,l,2)&
             &+(1.0d0)*trans_phi_dx(i,l,3)

  !Second derivative value of phi
      D_phi_2=+(1.0d0)*trans_phi_dx(i,l,2)&
             &-(2.0d0)*trans_phi_dx(i,l,3)&
             &+(1.0d0)*trans_phi_dx(i,l,4)  

  !Shape parameter on [j-1,j+0]
      lambda_dx=-D_phi_2/(D_phi_1+sign(dx**4,D_phi_1))      

  !Local coefficients
      aL_0(i,l)=-1.0d0/2.0d0+lambda_dx/3.0d0
      aL_1(i,l)=+3.0d0/2.0d0-lambda_dx/3.0d0

  !============================================================================================================!
  
  !First derivative value of phi
      D_phi_1=(1.0d0/12.0d0)*trans_phi_dx(i,l,2)&
            &-(5.0d0/4.00d0)*trans_phi_dx(i,l,3)&
            &+(5.0d0/4.00d0)*trans_phi_dx(i,l,4)&
            &-(1.0d0/12.0d0)*trans_phi_dx(i,l,5)

  !Second derivative value of phi
      D_phi_2=(1.0d0/2.0d0)*trans_phi_dx(i,l,2)&
            &-(1.0d0/2.0d0)*trans_phi_dx(i,l,3)&
            &-(1.0d0/2.0d0)*trans_phi_dx(i,l,4)&
            &+(1.0d0/2.0d0)*trans_phi_dx(i,l,5)

  !Shape parameter on [j+0,j+1]
      lambda_dx=-D_phi_2/(D_phi_1+sign(dx**4,D_phi_1))    
  
      bR_1(i,l)=+1.0d0/2.0d0-lambda_dx/6.0d0
      bR_2(i,l)=+1.0d0/2.0d0+lambda_dx/6.0d0

  !============================================================================================================! 
  
  !------------------------------------------------------------------------------------------------------------!
  !Activate the WENO reconstruction
  !------------------------------------------------------------------------------------------------------------!

  !Parameters for first-order generalized undivided differences
      delta=sign(dx**2,trans_phi_dx(i,l,3))  

      gamma_0=-(trans_phi_dx(i,l,3)-trans_phi_dx(i,l,2))/(trans_phi_dx(i,l,3)+delta)

  !    if(dabs(gamma_0).gt.dx**1) gamma_0=0.0d0

      exp_0=+(gamma_0**0)/1.0d0&
           &+(gamma_0**1)/1.0d0&
           &+(gamma_0**2)/2.0d0&
           &+(gamma_0**3)/6.0d0

      gamma_1=-(trans_phi_dx(i,l,4)-trans_phi_dx(i,l,3))/(trans_phi_dx(i,l,3)+delta)

  !    if(dabs(gamma_1).gt.dx**1) gamma_1=0.0d0

      exp_1=+(gamma_1**0)/1.0d0&
           &+(gamma_1**1)/1.0d0&
           &+(gamma_1**2)/2.0d0&
           &+(gamma_1**3)/6.0d0

  !Parameters for the second-order generalized undivided difference
      eta=(trans_phi_dx(i,l,2)-2.0d0*trans_phi_dx(i,l,3)+trans_phi_dx(i,l,4))/(trans_phi_dx(i,l,3)+delta)

  !    if(dabs(eta).gt.dx**1) eta=0.0d0

      exp_2=+2.0d0*(eta**0)/1.0d0&
           &+2.0d0*(eta**1)/2.0d0&
           &+2.0d0*(eta**2)/24.0d0
  
  !------------------------------------------------------------------------------------------------------------!
  
  !Local smoothness indicators (WENO-JS,WENO-M,WENO-Z,WENO-P+)
  !    beta_0=(trans_phi_dx(i,l,3)-trans_phi_dx(i,l,2))**2
  !    beta_1=(trans_phi_dx(i,l,4)-trans_phi_dx(i,l,3))**2

  !Local smoothness indicators (WENO-NZ)
  !    beta_0=(trans_phi_dx(i,l,2)-exp_0*trans_phi_dx(i,l,3))**2
  !    beta_1=(trans_phi_dx(i,l,3)-exp_1*trans_phi_dx(i,l,4))**2

  !Local smoothness indicators (Truncated form)
      temp=(trans_phi_dx(i,l,3)-trans_phi_dx(i,l,2))**2
      temp=(temp/(2.0d0*trans_phi_dx(i,l,3)+dx**2))**2
      beta_0=temp*dx**4

      temp=(trans_phi_dx(i,l,4)-trans_phi_dx(i,l,3))**2
      temp=(temp/(2.0d0*trans_phi_dx(i,l,3)+dx**2))**2
      beta_1=temp*dx**4
  
  !------------------------------------------------------------------------------------------------------------!

  !------------------------------------------------------------------------------------------------------------! 
  
  !Global smoothness indicator (WENO-Z)
  !    tau=dabs(beta_1-beta_0)
 
  !Global smoothness indicator (WENO-P,WENO-P+)
  !    tau=dabs((beta_0+beta_1)/2.0d0-(trans_phi_dx(i,l,2)-trans_phi_dx(i,l,4))**2/4.0d0)

  !Global smoothness indicator (WENO-NZ)
  !    tau=(trans_phi_dx(i,l,2)-exp_2*trans_phi_dx(i,l,3)+trans_phi_dx(i,l,4))**2

  !Global smoothness indicator (Truncated form)
      temp=(trans_phi_dx(i,l,2)-2.0d0*trans_phi_dx(i,l,3)+trans_phi_dx(i,l,2))**2
      temp=(temp/(24.0d0*trans_phi_dx(i,l,3)+dx**2))**2
      tau=temp*dx**8      

  !------------------------------------------------------------------------------------------------------------!  
      
  !------------------------------------------------------------------------------------------------------------!  

  !Unnormalized weights (WENO-JS,M)
  !    alpha_0=d_0(i,l)/(beta_0+eps)**2
  !    alpha_1=d_1(i,l)/(beta_1+eps)**2

  !Unnormalized weights (WENO-Z,WENO-P)
  !    alpha_0=d_0(i,l)*(1.0d0+tau/(beta_0+eps))
  !    alpha_1=d_1(i,l)*(1.0d0+tau/(beta_1+eps))
    
  !Unnormalized weights (WENO-P+)
  !    alpha_0=d_0(i,l)*(1.0d0+tau/(beta_0+eps)+((dx**(1.0d0/6.0d0)))*((beta_0+eps)/(tau+eps)))
  !    alpha_1=d_1(i,l)*(1.0d0+tau/(beta_1+eps)+((dx**(1.0d0/6.0d0)))*((beta_1+eps)/(tau+eps)))      

  !Unnormalized weights (WENO-NZ)
      alpha_0=d_0(i,l)*(1.0d0+(tau/(beta_0+eps))+2.0d0*((beta_0/(tau+eps))**2))
      alpha_1=d_1(i,l)*(1.0d0+(tau/(beta_1+eps))+2.0d0*((beta_1/(tau+eps))**2))

  !------------------------------------------------------------------------------------------------------------!
     
  !------------------------------------------------------------------------------------------------------------!

  !Normalized weights
      omega_0=alpha_0/(alpha_0+alpha_1)
      omega_1=alpha_1/(alpha_0+alpha_1)

  !------------------------------------------------------------------------------------------------------------!
  
  !The below procedures are valid only for WENO-M   
  !Unnormalized nonlinear weights (WENO-M)
  !    alpha_0=omega_0*(d_0(i,l)+(d_0(i,l)**2)-3.0d0*d_0(i,l)*omega_0+omega_0**2)
  !    alpha_0=alpha_0/((d_0(i,l)**2)+omega_0*(1.0d0-2.0d0*d_0(i,l)))

  !    alpha_1=omega_1*(d_1(i,l)+(d_1(i,l)**2)-3.0d0*d_1(i,l)*omega_1+omega_1**2)
  !    alpha_1=alpha_1/((d_1(i,l)**2)+omega_1*(1.0d0-2.0d0*d_1(i,l)))
  
  !Normalized nonlinear weights (WENO-M)
  !    omega_0=alpha_0/(alpha_0+alpha_1)
  !    omega_1=alpha_1/(alpha_0+alpha_1)      

  !------------------------------------------------------------------------------------------------------------!
  
  !------------------------------------------------------------------------------------------------------------!    
  
  !The below procedures are valid only for WENO-NZM  
  !Set up the mu function
      if(dabs(tau).lt.xi)then
        mu=1.0d0
      elseif((xi.ge.tau).and.(tau.lt.(eps+xi)))then
        mu=(3.0d0/2.0d0)*(1.0d0-eps)*B3_Spline((2.0d0/eps)*(tau-xi))
      else
        mu=eps
      endif
  
  !Filtered Normalized weights    
      omega_0=mu*d_0(i,l)+(1.0d0-mu)*omega_0
      omega_1=mu*d_1(i,l)+(1.0d0-mu)*omega_1

  !------------------------------------------------------------------------------------------------------------!

  !------------------------------------------------------------------------------------------------------------!    
  
  !Set up the local coefficients
  !Singularity location: [j-1,j+0]    
      if(dabs(beta_0).ge.1.0d0)then
        b_1(i,l)=bR_1(i,l)
        b_2(i,l)=bR_2(i,l)       
      endif
  
  !Singularity location: [j+0,j+1]        
      if(dabs(beta_1).ge.1.0d0)then
        a_0(i,l)=aL_0(i,l)
        a_1(i,l)=aL_1(i,l)
      endif

  !Convex summation  
      WENO(i,l)=omega_0*(a_0(i,l)*trans_phi_dx(i,l,2)+a_1(i,l)*trans_phi_dx(i,l,3))&
              &+omega_1*(b_1(i,l)*trans_phi_dx(i,l,3)+b_2(i,l)*trans_phi_dx(i,l,4)) 
            
  !Modification at the critical points   
      if(dabs(D_phi_1).le.dx)then  
         WENO(i,l)=-(1.0d0/6.0d0)*trans_phi_dx(i,l,2)&
                  &+(5.0d0/6.0d0)*trans_phi_dx(i,l,3)&
                  &+(2.0d0/6.0d0)*trans_phi_dx(i,l,4)&
                  &-(1.0d0/12.0d0)*D_phi_3
      endif               
  
  !Singularity location: [j-1,j+0]    
  !    if(dabs(beta_0).ge.1.0d0)then
        WENO(i,l)=bR_1(i,l)*trans_phi_dx(i,l,3)&
                &+bR_2(i,l)*trans_phi_dx(i,l,4)
  
  !Modification at the critical points        
        if(dabs(D_phi_1).le.dx**2)then
          WENO(i,l)=+(1.0d0/2.0d0)*trans_phi_dx(i,l,3)&
                   &+(1.0d0/2.0d0)*trans_phi_dx(i,l,4)&
                   &-(1.0d0/6.0d0)*D_phi_2
        endif           
  !    endif
  
  !Singularity location: [j+0,j+1]        
  !    if(dabs(beta_1).ge.1.0d0)then
  !      WENO(i,l)=aL_0(i,l)*trans_phi_dx(i,l,2)&
  !               +aL_1(i,l)*trans_phi_dx(i,l,3)
  
  !Modification at the critical points        
  !      if(dabs(D_phi_1).le.dx**2)then
  !        WENO(i,l)=-(1.0d0/2.0d0)*trans_phi_dx(i,l,2)&
  !                 &+(3.0d0/2.0d0)*trans_phi_dx(i,l,3)&
  !                 &+(1.0d0/3.0d0)*D_phi_2
  !      endif
  !    endif
  
    enddo
  enddo

  return
  !---------------------------------------------------Calculations have finished---------------------------------------------------!
endsubroutine

subroutine Get_RK(i_bc,imax,dx,dt,old_phi,new_phi)
  implicit none

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: dx,dt
  real(8),intent(inout) :: old_phi(-i_bc+0:imax+i_bc)

  real(8),intent(out) :: new_phi(-i_bc+0:imax+i_bc)

  !Counter
  integer :: i

  !Claculation variables
  real(8),dimension(0:imax,1:4) :: hat_H
  real(8),dimension(-i_bc+0:imax+i_bc,1:3) :: phi
  !-----------------------------------------Calculations have started------------------------------------------!

  !============================================================================================================!
  !Stage 1
  !============================================================================================================!
  
  !Compute the Hamiltonian
  call Get_Hat_Hx(i_bc,imax,dx,old_phi,hat_H(:,1))

  !Update variables
  do i=0,imax
    phi(i,1)=old_phi(i)+dt*hat_H(i,1)/2.0d0
  enddo

  !============================================================================================================!
  !Stage 2
  !============================================================================================================!

  !Compute the Hamiltonian
  call Get_Hat_Hx(i_bc,imax,dx,phi(:,1),hat_H(:,2))

  !Update variables
  do i=0,imax
    phi(i,2)=phi(i,1)+dt*(-1.0d0*hat_H(i,1)+1.0d0*hat_H(i,2))/2.0d0
  enddo

  !============================================================================================================!
  !Stage 3
  !============================================================================================================!
  
  !Compute the Hamiltonian
  call Get_Hat_Hx(i_bc,imax,dx,phi(:,2),hat_H(:,3))

  !Update variables
  do i=0,imax
    phi(i,3)=phi(i,2)+dt*(-1.0d0*hat_H(i,2)+2.0d0*hat_H(i,3))/2.0d0
  enddo

  !============================================================================================================!
  !Stage 4
  !============================================================================================================!

  !Compute the Hamiltonian
  call Get_Hat_Hx(i_bc,imax,dx,phi(:,3),hat_H(:,4))

  !Update variables
  do i=0,imax
    new_phi(i)=phi(i,3)+dt*(1.0d0*hat_H(i,1)+2.0d0*hat_H(i,2)-4.0d0*hat_H(i,3)+1.0d0*hat_H(i,4))/6.0d0
  enddo

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_Results(imax,dx,x_left,phi)
  implicit none

  integer,intent(in) :: imax
  real(8),intent(in) :: dx,x_left
  real(8),intent(in) :: phi(0:imax)

  !Counter
  integer :: i
  !-----------------------------------------Calculations have started------------------------------------------!

  open(1,file='Numerical.plt')
  write(1,*)'zone T="file"', 'I=',imax+1

  do i=0,imax
    write(1,*) x_left+dble(i)*dx,phi(i)
  enddo

  close(1)

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

!============================================================================================================!

!============================================================================================================!

subroutine Get_Timestep_REF(i_bc,imax,dx,flag_1,time_out,CFL,phi,dt)
  implicit none

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: dx
  real(8),intent(in) :: flag_1,time_out,CFL
  real(8),intent(inout) :: phi(-i_bc+0:imax+i_bc)

  real(8),intent(out) :: dt

  !Counter
  integer :: i

  !Calcuation variables
  real(8),dimension(0:imax) :: F_phi_dx,B_phi_dx
  real(8) :: max_DH !Artificial Hamiltonian speed
  real(8) :: S_max !Maximal wave speed
  !-----------------------------------------Calculations have started------------------------------------------!
  
  !Impose the boundary condition
  call Get_Boundary(i_bc,imax,phi)

  !Compute the artificial Hamiltonian speed
  S_max=0.0d0
  do i=0,imax
    F_phi_dx(i)=(phi(i+1)-phi(i+0))/dx
    B_phi_dx(i)=(phi(i+0)-phi(i-1))/dx
    max_DH=max(dabs(dsin(F_phi_dx(i)+1.0d0)),dabs(dsin((B_phi_dx(i))+1.0d0)))
    S_max=max(S_max,max_DH)
  enddo

  !Determine the time step size
  dt=(CFL*(dx**(5/4)))/S_max

  !Restriction
  if((flag_1+dt).ge.time_out)then
    dt=(time_out-flag_1)
  endif

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_Hat_Hx_REF(i_bc,imax,dx,phi,hat_H)
  implicit none

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: dx
  real(8),intent(inout) :: phi(-i_bc+0:imax+i_bc)

  real(8),intent(out) :: hat_H(0:imax)

  !Counter
  integer :: i

  !Reconstructed variable
  real(8),dimension(0:imax,1:2) :: WENO

  !Canculation variables
  real(8),dimension(0:imax) :: P_phi_dx,M_phi_dx
  real(8),dimension(0:imax) :: L_max_alpha
  real(8) :: G_max_alpha
  !-----------------------------------------Calculations have started------------------------------------------!

  !Impose the boundary condition
  call Get_Boundary(i_bc,imax,phi)

  !Apply the WENO scheme
  call Get_WENO_REF(i_bc,imax,dx,phi,WENO)

  !Seperate reconstructed arraies into two parts
  do i=0,imax
    M_phi_dx(i)=WENO(i,1)
    P_phi_dx(i)=WENO(i,2)
  enddo

  !Compute alpha in Lax-Friedrich flux
  !In this case, |H'(phi_dx)|=dsin(|phi_dx+1.0d0|)
  do i=0,imax
    L_max_alpha(i)=max(dabs(dsin(M_phi_dx(i)+1.0d0)),dabs(dsin((P_phi_dx(i))+1.0d0)))
  enddo
  G_max_alpha=maxval(L_max_alpha)

  !Apply Lax-Friedrich flux to the reconstructed Hamiltonian
  do i=0,imax
    hat_H(i)=-dcos((P_phi_dx(i)+M_phi_dx(i))/2.0d0+1.0d0)-G_max_alpha*(P_phi_dx(i)-M_phi_dx(i))

    hat_H(i)=-hat_H(i)/1.0d0
  enddo

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_WENO_REF(i_bc,imax,dx,phi,WENO)
  implicit none

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: dx
  real(8),intent(in) :: phi(-i_bc+0:imax+i_bc)

  real(8),intent(out) :: WENO(0:imax,1:2)

  !Counters
  integer :: i,i0,i1
  integer :: l

  !Derivatives of the flux
  real(8) :: D_phi_1,D_phi_2,D_phi_3

  !Shape parameter
  real(8) :: lambda_dx,lambda2_dx2

  !Local coefficients
  real(8),dimension(0:imax,1:2) :: a_0,a_1,a_2 
  real(8),dimension(0:imax,1:2) :: b_1,b_2,b_3 
  real(8),dimension(0:imax,1:2) :: c_2,c_3,c_4 

  !Ideal weights
  real(8),dimension(0:imax,1:2) :: d_0,d_1,d_2

  !Calculation variables
  real(8),dimension(-i_bc+1:imax+i_bc) :: phi_dx
  real(8),dimension(0:imax,1:2,1:5) :: trans_phi_dx
  
  !Main WENO reconstruction variables
  real(8) :: eps
  real(8) :: beta_0,beta_1,beta_2
  real(8) :: alpha_0,alpha_1,alpha_2
  real(8) :: omega_0,omega_1,omega_2
  !---------------------------------------------------Calculations have started----------------------------------------------------!

  !Set up the free parameter
  eps=1d-6

  !Compute the first derivative of phi using the first-order difference
  do i=-i_bc+1,imax+i_bc
    phi_dx(i)=(phi(i+0)-phi(i-1))/dx
  enddo

  !Transform computed first derivative of phi in each directions
  do l=1,5
    i0=l-3
    i1=3-l
    do i=0,imax
  !Upwind part [j-1,j+0,j+1]    
      trans_phi_dx(i,1,l)=phi_dx(i+i0+0) 
  !Downwind part [j+2,j+1,j+0]
      trans_phi_dx(i,2,l)=phi_dx(i+i1+1) 
    enddo
  enddo
  
  do l=1,2
    do i=0,imax
  
  !============================================================================================================!

  !Local coefficients
      a_0(i,l)=+2.0d0/6.0d0
      a_1(i,l)=-7.0d0/6.0d0
      a_2(i,l)=+11.0d0/6.0d0

      b_1(i,l)=-1.0d0/6.0d0
      b_2(i,l)=+5.0d0/6.0d0
      b_3(i,l)=+2.0d0/6.0d0

      c_2(i,l)=+2.0d0/6.0d0
      c_3(i,l)=+5.0d0/6.0d0
      c_4(i,l)=-1.0d0/6.0d0
              
  !Ideal weights
      d_0(i,l)=+1.0d0/10.0d0
      d_1(i,l)=+6.0d0/10.0d0
      d_2(i,l)=+3.0d0/10.0d0
  
  !============================================================================================================! 
  
  !------------------------------------------------------------------------------------------------------------!
  !Activate the WENO reconstruction
  !------------------------------------------------------------------------------------------------------------!
  
  !------------------------------------------------------------------------------------------------------------!
  
  !Local smoothness indicators 
      beta_0=+(1.0d0/4.0d0)*(1.0d0*trans_phi_dx(i,l,1)-4.0d0*trans_phi_dx(i,l,2)+3.0d0*trans_phi_dx(i,l,3))**2&
          &+(13.0d0/12.0d0)*(1.0d0*trans_phi_dx(i,l,1)-2.0d0*trans_phi_dx(i,l,2)+1.0d0*trans_phi_dx(i,l,3))**2
      
      beta_1=+(1.0d0/4.0d0)*(1.0d0*trans_phi_dx(i,l,2)+0.0d0*trans_phi_dx(i,l,3)-1.0d0*trans_phi_dx(i,l,4))**2&
          &+(13.0d0/12.0d0)*(1.0d0*trans_phi_dx(i,l,2)-2.0d0*trans_phi_dx(i,l,3)+1.0d0*trans_phi_dx(i,l,4))**2    
   
      beta_2=+(1.0d0/4.0d0)*(3.0d0*trans_phi_dx(i,l,3)-4.0d0*trans_phi_dx(i,l,4)+1.0d0*trans_phi_dx(i,l,5))**2&
          &+(13.0d0/12.0d0)*(1.0d0*trans_phi_dx(i,l,3)-2.0d0*trans_phi_dx(i,l,4)+1.0d0*trans_phi_dx(i,l,5))**2    
  
  !------------------------------------------------------------------------------------------------------------!
      
  !------------------------------------------------------------------------------------------------------------!  

  !Unnormalized weights 
      alpha_0=d_0(i,l)/(beta_0+eps)**2
      alpha_1=d_1(i,l)/(beta_1+eps)**2
      alpha_2=d_2(i,l)/(beta_2+eps)**2

  !------------------------------------------------------------------------------------------------------------!
     
  !------------------------------------------------------------------------------------------------------------!

  !Normalized weights
      omega_0=alpha_0/(alpha_0+alpha_1+alpha_2)
      omega_1=alpha_1/(alpha_0+alpha_1+alpha_2)
      omega_2=alpha_2/(alpha_0+alpha_1+alpha_2)

  !------------------------------------------------------------------------------------------------------------!

  !------------------------------------------------------------------------------------------------------------!
  
  !Convex summation  
      WENO(i,l)=omega_0*(a_0(i,l)*trans_phi_dx(i,l,1)+a_1(i,l)*trans_phi_dx(i,l,2)+a_2(i,l)*trans_phi_dx(i,l,3))&
              &+omega_1*(b_1(i,l)*trans_phi_dx(i,l,2)+b_2(i,l)*trans_phi_dx(i,l,3)+b_3(i,l)*trans_phi_dx(i,l,4))&
              &+omega_2*(c_2(i,l)*trans_phi_dx(i,l,3)+c_3(i,l)*trans_phi_dx(i,l,4)+c_4(i,l)*trans_phi_dx(i,l,5))
            
  !------------------------------------------------------------------------------------------------------------!    

    enddo
  enddo

  return
  !---------------------------------------------------Calculations have finished---------------------------------------------------!
endsubroutine

subroutine Get_RK_REF(i_bc,imax,dx,dt,old_phi,new_phi)
  implicit none

  integer,intent(in) :: i_bc,imax
  real(8),intent(in) :: dx,dt
  real(8),intent(inout) :: old_phi(-i_bc+0:imax+i_bc)

  real(8),intent(out) :: new_phi(0:imax)

  !Counter
  integer :: i

  !Claculation variables
  real(8),dimension(0:imax,1:4) :: hat_H
  real(8),dimension(-i_bc+0:imax+i_bc,1:3) :: phi
  !-----------------------------------------Calculations have started------------------------------------------!

  !============================================================================================================!
  !Stage 1
  !============================================================================================================!

  !Compute the Hamiltonian
  call Get_Hat_Hx_REF(i_bc,imax,dx,old_phi,hat_H(:,1))

  !Update variables
  do i=0,imax
    phi(i,1)=old_phi(i)+dt*hat_H(i,1)/2.0d0
  enddo

  !============================================================================================================!
  !Stage 2
  !============================================================================================================!

  !Compute the Hamiltonian
  call Get_Hat_Hx_REF(i_bc,imax,dx,phi(:,1),hat_H(:,2))

  !Update variables
  do i=0,imax
    phi(i,2)=phi(i,1)+dt*(-1.0d0*hat_H(i,1)+1.0d0*hat_H(i,2))/2.0d0
  enddo

  !============================================================================================================!
  !Stage 3
  !============================================================================================================!
  
  !Compute the Hamiltonian
  call Get_Hat_Hx_REF(i_bc,imax,dx,phi(:,2),hat_H(:,3))

  !Update variables
  do i=0,imax
    phi(i,3)=phi(i,2)+dt*(-1.0d0*hat_H(i,2)+2.0d0*hat_H(i,3))/2.0d0
  enddo

  !============================================================================================================!
  !Stage 4
  !============================================================================================================!

  !Compute the Hamiltonian
  call Get_Hat_Hx_REF(i_bc,imax,dx,phi(:,3),hat_H(:,4))

  !Update variables
  do i=0,imax
    new_phi(i)=phi(i,3)+dt*(1.0d0*hat_H(i,1)+2.0d0*hat_H(i,2)-4.0d0*hat_H(i,3)+1.0d0*hat_H(i,4))/6.0d0
  enddo

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_Results_REF(imax,dx,x_left,phi)
  implicit none

  integer,intent(in) :: imax
  real(8),intent(in) :: dx,x_left
  real(8),intent(in) :: phi(0:imax)

  !Counter
  integer :: i
  !-----------------------------------------Calculations have started------------------------------------------!
  
  open(0,file='REF_SOL.txt')

  do i=0,imax
    write(0,*) phi(i)
  enddo

  close(0)

  open(2,file='Numerical_REF.plt')
  write(2,*)'zone T="file"', 'I=',imax+1

  do i=0,imax
    write(2,*) x_left+dble(i)*dx,phi(i)
  enddo

  close(2)

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

!============================================================================================================!

!============================================================================================================!

subroutine Get_L_1_Error(imax,imax_REF,p,dx,exact,numer)
  implicit none

  integer,intent(in) :: imax,imax_REF,p
  real(8),intent(in) :: dx,exact(0:imax_REF),numer(0:imax)

  !Counter
  integer :: i

  !Calculation variable
  real(8),dimension(0:imax) :: e
  !-----------------------------------------Calculations have started------------------------------------------!

  !Compute absolute differences
  do i=0,imax
    e(i)=dabs(exact(i*(1280/(2**(p-1))))-numer(i))
  enddo
 

  !Compute the average sum
  write(*,*) "L_1 error:",real(sum(e)/dble(imax))

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

subroutine Get_L_Infty_Error(imax,imax_REF,p,dx,exact,numer)
  implicit none

  integer,intent(in) :: imax,imax_REF,p
  real(8),intent(in) :: dx,exact(0:imax_REF),numer(0:imax)

  !Counters
  integer :: i

  !Calculation variable
  real(8),dimension(0:imax) :: e
  !-----------------------------------------Calculations have started------------------------------------------!

  !Compute absolute differences
  do i=0,imax
    e(i)=dabs(exact(i*(1280/(2**(p-1))))-numer(i))
  enddo

  !Take the maximum value
  !write(*,*) "Maximum error location:",maxloc(e)
  write(*,*) "L_infty error:",real(maxval(e))

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endsubroutine

real(8) function B3_Spline(x)
  implicit none

  real(8),intent(in) :: x
  !-----------------------------------------Calculations have started------------------------------------------!

  if((0.0d0.le.x).and.(x.lt.1.0d0))then
    B3_Spline=(2.0d0/3.0d0)-(x**2)*(1.0d0-x/2.0d0)
  elseif((1.0d0.le.x).and.(x.le.2.0d0))then
    B3_Spline=(1.0d0/6.0d0)*(2.0d0-x)**3
  endif

  return
  !-----------------------------------------Calculations have finished-----------------------------------------!
endfunction
