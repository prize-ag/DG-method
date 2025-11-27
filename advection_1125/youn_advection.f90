program RKDG_WENO_SCALAR_LEGENDRE
    implicit none

    !Counters
    integer :: ord !Index of the degree for Legendre polymomials 
    ! <<>> 이건 언제 사용하는 걸까나? 
    integer :: i !Spatial(x-axis) index
    integer :: n !Time index

    !Order of the discontinuous Galerkin method and polynomials
    integer,parameter :: DG_ord = 3 !Order of DG solver (2,3,4,5)
    integer,parameter :: Legen_ord = DG_ord - 1 !Order of polymomials (1,2,3,4)

    !Useful constants
    double precision,parameter :: pi = 4.0d0 * datan(1.0d0)
    double precision,parameter :: CFL = 0.5d0 / dble(2 * Legen_ord + 1) !CFL number
    double precision :: TVB_M !TVB constant

    !Spatial variables
    integer,parameter :: imax= 200 * (2 ** 0) !Number of spatial cells
    integer,parameter :: i_bc = 2 * Legen_ord + 1 !Ghost boundary cells for WENO schemes
    double precision,parameter :: x_left  = -1.0d0 !Left boundary on x-axis
    double precision,parameter :: x_right =  1.0d0 !Right boundary on x-axis
    double precision,parameter :: x_length= dabs(x_right - x_left) !Length of domain
    double precision,dimension(-i_bc + 1: imax + i_bc) :: x !Cell centers
    double precision,dimension(-i_bc + 1: imax + i_bc) :: dx !Cell volumes
    
    !Time variables
    integer :: RK_ord !Order of Runge-Kutta solver
    integer,parameter :: nmax = 999999999 !Number of the marching step
    double precision,parameter :: time_out = 2.0d0 !Final time of the solution
    double precision :: dt !Time step size
    double precision :: flag_0,flag_1 !Time flags

    !Cpu time checkers
    integer :: rate,time_0,time_1

    !Variables for quadratures
    double precision,dimension(1:6,1:2) :: Q 
    
    !Several functions
    double precision,external :: Ini_u , Ext_u

    !Legendre polynomials
    double precision,external :: Le_poly

    !Mass coefficient matrices
    double precision,dimension(0 : Legen_ord , 0 : Legen_ord) :: M

    !Cell interface contributions
    double precision,dimension(0 : Legen_ord) :: l_b , r_b
    
    !Polynomial ideal weights 
    double precision,dimension(0 : Legen_ord , 1 : Legen_ord + 2) :: d 

    !Degrees of freedom in the moments by the L2 projection
    double precision,dimension(0 : Legen_ord,-i_bc + 1 : imax + i_bc) :: old_u 
    double precision,dimension(0 : Legen_ord, 1 : imax) :: new_u 

    !Convex summation of degrees of freedom and Legendre polynommials 
    double precision :: temp_sum
    double precision,dimension(1 : imax) :: sum_ext_u , sum_new_u

    write(*,*) "================================================================="
    write(*,*) "        Scalar conservation laws solver. (Linear advection)      "
    write(*,*) "================================================================="
    write(*,*) "                    Calculations have started.                   "
    write(*,*) "================================================================="
    write(*,*) " "

    !Set several initial settings
    n=0
    flag_0=0.0d0
    flag_1=0.0d0
    call Get_Quadrature(7 , Q)
    call Get_Mass(Legen_ord,Q,M)
    call Get_Spatial(i_bc,imax,x_left,x_length,x,dx)
    call Get_Initial(Legen_ord,Q,M,i_bc,imax,x,dx,old_u,l_b,r_b)
    !call Get_Weights(Legen_ord,Q,i_bc,imax,x,dx,d)

    d = 0.0d0 

    !Set the TVB constant
    !if TVB_M is chosen too small(<1.0d0), smooth cells will be declared troubled cells,
    !so unnecessary WENO reconstructions will be appear.

    !if TVB_M is chosen too large(>1.0d0), trouble cells will be declared smooth cells,
    !so suprious oscillations will be appear.
    TVB_M=1.0d0
    
    !Smooth case
    RK_ord=DG_ord

    !if(RK_ord.eq.5)then
    !    RK_ord=RK_ord-1
    !endif

    !Nonsmooth case
    !RK_ord=3

    !write(*,*) "================================================================="
    !write(*,*) "The order of the discontinuous Galerkin solver:", DG_ord
    !write(*,*) "The order of the stage of Runge-Kutta solver:  ", RK_ord
    !write(*,*) "The maximal order of Legendre polynomials:     ", Legen_ord
    !write(*,*) "The selected CFL number:                       ", real(CFL)
    !write(*,*) "The selected TVB constant:                     ", real(TVB_M)
    !write(*,*) "================================================================="
    !write(*,*) " "

    open(0,file='Density(Exact).plt')
    open(1,file='Density(Numerical).plt')
    
    !write(0,*)'zone T="file"', 'I=',imax
    !write(1,*)'zone T="file"', 'I=',imax

    
    call system_clock(count_rate = rate)
    call system_clock(time_0)

    !============================================================================!
    !=======================Time marching procedure(Start)=======================!
    !============================================================================!
      !flag_0,1 ㅣ 둘다 필요할까..???
    ! do while (flag_0 <= time_out)

    !   call Get_Timestep(DG_ord,i_bc,imax,dx,flag_1,time_out,CFL,dt)
    !   call Get_RK(RK_ord,TVB_M,Legen_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,dt&
    !                                                         &,old_u,new_u)

    !   flag_0 = flag_0 + dt
    !                         !Receive the updated data
    !   do i=1,imax
    !     do ord=0,Legen_ord
    !       old_u(ord,i)=new_u(ord,i)
    !     enddo 
    !   enddo  
    !   flag_1 = flag_0 

    ! end do

    do n=1,nmax
        if(flag_0 < time_out)then
          flag_1=flag_0

          !Compute the time step size by using the first-derivative of flux
          call Get_Timestep(DG_ord,i_bc,imax,dx,flag_1,time_out,CFL,dt)
        
          !Update the variable by using TVD RK solver
          call Get_RK(RK_ord,TVB_M,Legen_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,dt&
                                                            &,old_u,new_u)
          !Update the real-time
          flag_1=flag_1+dt
          flag_0=flag_1

          !Receive the updated data
          do i=1,imax
              do ord=0,Legen_ord
                old_u(ord,i)=new_u(ord,i)
              enddo 
          enddo
          
        elseif(flag_0 == time_out)then
         !Record the numerical result at the finial time
          do i=1,imax
              temp_sum=0.0d0
              do ord=0,Legen_ord
                temp_sum=temp_sum+new_u(ord,i)*Le_poly(ord,0.0d0)
                sum_new_u(i)=temp_sum
              enddo
          enddo
    
          !Record the final numerical result
          goto 10
        endif
    enddo 
    10 continue 
    call Get_Result(i_bc,imax,x,sum_new_u(1:imax))

    !============================================================================!
    !======================Time marching procedure(Finish)=======================!
    !============================================================================!

    call system_clock(time_1)
    write(*,*) "=============Computational and physical time results============="
    write(*,*) "Physical time:", real(flag_0),"//   ",real(time_out)
    write(*,*) "Computation time:",real(time_1-time_0)/real(rate),"sec"
    write(*,*) "================================================================="
    write(*,*) " "

    !Set up the exact solution
    do i=1,imax
        sum_ext_u(i)=Ext_u(x(i),time_out)
        write(0,*) x(i),sum_ext_u(i)
    enddo

    close(0)
    close(1)

    !Compute various errors
    write(*,*) "==========================Error results=========================="
    call Get_L_2_Error(i_bc,imax,dx,sum_ext_u,sum_new_u)
    call Get_L_Infty_Error(imax,sum_ext_u,sum_new_u)
    write(*,*) "================================================================="
    write(*,*) " "

    write(*,*) "================================================================="
    write(*,*) "                   Calculations have finished.                   "
    write(*,*) "================================================================="
    return
endprogram

subroutine Get_Quadrature(ord,Q) 

    implicit none

    integer,intent(in) :: ord 

    double precision,intent(out) ::Q(1:6,1:2)

    Q = 0.0d0
    
    select case(ord)

        ! case 2 3 4 는 가우시안 쿼드락쳐 를 사용한 것 
        ! 2 3 4 는 각각 점의 개수 
        case (2) 
            Q(1,1) = -dsqrt(3.0d0) / 3.0d0 / 2.0d0
            Q(2,1) = +dsqrt(3.0d0) / 3.0d0 / 2.0d0
            
            Q(1,2) = 1.0d0 / 2.0d0
            Q(2,2) = 1.0d0 / 2.0d0

        case (3) 
            Q(1,1) = -dsqrt(15.0d0) / 5.0d0 / 2.0d0
            Q(2,1) = 0.0d0
            Q(3,1) = +dsqrt(15.0d0) / 5.0d0 / 2.0d0

            Q(1,2) = 5.0d0 / 9.0d0 / 2.0d0
            Q(2,2) = 8.0d0 / 9.0d0 / 2.0d0
            Q(3,2) = 5.0d0 / 9.0d0 / 2.0d0

        case(4) 
            Q(1,1) = -dsqrt(3.0d0 / 7.0d0 + 2.0d0 / 7.0d0 * dsqrt(6.0d0 / 5.0d0)) / 2.0d0
            Q(2,1) = -dsqrt(3.0d0 / 7.0d0 - 2.0d0 / 7.0d0 * dsqrt(6.0d0 / 5.0d0)) / 2.0d0
            Q(3,1) = +dsqrt(3.0d0 / 7.0d0 - 2.0d0 / 7.0d0 * dsqrt(6.0d0 / 5.0d0)) / 2.0d0
            Q(4,1) = +dsqrt(3.0d0 / 7.0d0 + 2.0d0 / 7.0d0 * dsqrt(6.0d0 / 5.0d0)) / 2.0d0

            Q(1,2) = (18.0d0 - dsqrt(30.0d0))/36.0d0 / 2.0d0
            Q(2,2) = (18.0d0 + dsqrt(30.0d0))/36.0d0 / 2.0d0
            Q(3,2) = (18.0d0 + dsqrt(30.0d0))/36.0d0 / 2.0d0
            Q(4,2) = (18.0d0 - dsqrt(30.0d0))/36.0d0 / 2.0d0

        !점 4개를 사용하는 가우시안 로바또 쿼드락쳐 노드와 가중치 
        case(5) 
            Q(1,1) = -1.0d0/2.0d0
            Q(2,1) = -dsqrt(5.0d0)/5.0d0 / 2
            Q(3,1) = +dsqrt(5.0d0)/5.0d0 / 2
            Q(4,1) = +1.0d0/2.0d0

            Q(1,2) = 1.0d0/12.0d0
            Q(2,2) = 5.0d0/12.0d0
            Q(3,2) = 5.0d0/12.0d0
            Q(4,2) = 1.0d0/12.0d0
        ! 점 5개를 사용하는 가우시안 로바또 쿼드락쳐 노드와 가중치 
        case(6)  
            Q(1,1) = -1.0d0/2.0d0
            Q(2,1) = -sqrt(5.0d0/7.0d0) / 2.0d0
            Q(3,1) =  0.0d0
            Q(4,1) = +sqrt(5.0d0/7.0d0) / 2.0d0
            Q(5,1) = +1.0d0/2.0d0

            Q(1,2) = 1.0d0/20.0d0
            Q(2,2) = 49.0d0/90.0d0
            Q(3,2) = 64.0d0/90.0d0
            Q(4,2) = 49.0d0/90.0d0
            Q(5,2) = 1.0d0/20.0d0

        case(7) 
            Q(1,1) = -1.0d0/2.0d0
            Q(2,1) = -dsqrt(147.0d0 + 42.0d0*dsqrt(7.0d0)) / 42.0d0
            Q(3,1) = -dsqrt(147.0d0 - 42.0d0*dsqrt(7.0d0)) / 42.0d0
            Q(4,1) =  dsqrt(147.0d0 - 42.0d0*dsqrt(7.0d0)) / 42.0d0
            Q(5,1) =  dsqrt(147.0d0 + 42.0d0*dsqrt(7.0d0)) / 42.0d0
            Q(6,1) =  1.0d0/2.0d0


            Q(1,2) = 1.0d0/30.0d0
            Q(2,2) = dsqrt(7.0d0)*(7.0d0 + dsqrt(7.0d0)) * (-7.0d0 + 5.0d0*dsqrt(7.0d0)) / 840.0d0
            Q(3,2) = dsqrt(7.0d0) * (7.0d0 - dsqrt(7.0d0)) * (7.0d0 + 5.0d0*dsqrt(7.0d0)) / 840.0d0
            Q(4,2) = Q(3,2)
            Q(5,2) = Q(2,2)
            Q(6,2) = Q(1,2)

        case default 
            print * , "order는 숫자여야하고, 7이하여야합니다. "
        
    end select

    return 
endsubroutine

subroutine Get_Mass(Legen_ord,Q,M)

    implicit none 
    integer, intent(in) :: Legen_ord
    double precision, intent(in) :: Q(1:6,1:2)
    double precision, intent(out) :: M(0:Legen_ord,0:Legen_ord)

    integer :: i , j ,k 

    double precision,external :: Le_poly

    M = 0.0d0 

    ! -1/2에서 1/2 구간에서의 질량 행열 

    do i = 0,Legen_ord
        do j = 0,Legen_ord
            do k = 1,6
                M(i,j) = M(i,j) + Le_poly(i,Q(k,1)) * Le_poly(j,Q(k,1)) * Q(k,2)
            enddo

            if (M(i,j) > 1d-10) then
                M(i,j) = 1.0d0 / M(i,j)
            else 
                M(i,j) = 0.0d0 
            endif 
        enddo   
    enddo

    write(*,*) "===========================Mass matrix==========================="
    do i=0,Legen_ord
        write(*,'(4G16.6)') M(i,:)
    enddo
    write(*,*) "================================================================="
    write(*,*) " "
    
    return
endsubroutine

subroutine Get_Spatial(i_bc,imax,x_left,x_length,x,dx)
    implicit none 

    integer,intent(in) :: i_bc , imax
    double precision,intent(in) :: x_left , x_length
    
    double precision,intent(out) :: x(1 - i_bc : imax + i_bc) , dx(1 - i_bc : imax + i_bc)

    !Counter
    integer :: i

    !Calculation variable
    double precision,dimension(0:imax) :: temp_x 
    !-------------------------Calculations have started--------------------------!
    
    !Generate physical cell center locations and volumes of cells
    do i=0,imax
        temp_x(i) = x_left + dble(i) * dble(x_length)/dble(imax)
    enddo
    
    do i=1,imax
        x(i)= (temp_x(i-1) + temp_x(i))/2.0d0
        dx(i) = dabs(temp_x(i) - temp_x(i-1))
    enddo

    !Generate ghost cell center locations and volumes of cells 
    !for Lagrange interpolation at boundaries
    do i = 1 - i_bc , 0  
        x(i) =   x(imax - i) - x_length 
        dx(i) =  dx(imax -i)
    enddo 

    do i = 1 , i_bc 
        x (imax + i) = x(i) + x_length 
        dx(imax + i) = dx(i) 
    enddo 
    
    return
  !-------------------------Calculations have finished-------------------------!
endsubroutine

! u_0 : u의 초기값을 이렇게 설정. 이 서브루틴은 (초기값 + 양쪽 끝 값) 을 가져오는 루틴 
subroutine Get_Initial(Legen_ord,Q,M,i_bc,imax,x,dx,u_0,l_b,r_b)
    implicit none
    
    integer,intent(in) :: Legen_ord
    double precision,intent(in) :: Q(1:6,1:2) 
    double precision,intent(in) :: M(0:Legen_ord,0:Legen_ord)
    integer,intent(in) :: i_bc,imax
    double precision,intent(in) :: x(-i_bc+1:imax+i_bc),dx(-i_bc+1:imax+i_bc)

    double precision,intent(out) :: u_0(0:Legen_ord,-i_bc+1:imax+i_bc)
    double precision,intent(out) :: l_b(0:Legen_ord),r_b(0:Legen_ord)
    
    !Counters
    integer :: ord
    integer :: i , j , k , r

    !Initial data
    double precision,external :: Ini_u

    !Legendre polynomials
    double precision,external :: Le_poly

    !Calculation variables
    double precision,dimension(1:6) :: xq 
    double precision :: temp_sum
    double precision,dimension(1:imax) :: sum_u_0,exact_u_0
    !-------------------------Calculations have started--------------------------!
    
    u_0 = 0.0d0

    !Compute L2-projection of initial data into the finite element space
    do i = 1 , imax
    !Transform nodes by using the liner transform
        do k = 1,6
            xq(k) =  dx(i) * Q(k,1) + x(i) 
        enddo

    !Compute degree of freedoms of initial data by using quadratures
        do ord = 0 , Legen_ord           
            do r = 1,6 
                u_0(ord,i) = u_0(ord,i) + Ini_u(xq(r)) * Le_poly(ord,Q(r,1)) * Q(r,2) 
                ! weight를 바꾸야하지 않나..? 자코비안 곱해야할 것 같은데.... & 노드를 그냥 원래 노드를 넣어도 될까? 
                ! --> 무척 신기하네.. 다 포함되어 있는 것임... Q(1,2)를 넣으므로 써 자코비안이 이미 곱해진 형태를 가지고 있다. 
            enddo

                !Product the mass matrix for the orthonormalization                 
            u_0(ord,i)=M(ord,ord)*u_0(ord,i)
        enddo
    enddo              
                                    
    !Compute the convex summation by using calculated degrees of freedom
    !((x(i)-x(i))/dx(i)=0.0d0: Cell center 

    do i=1,imax

        temp_sum=0.0d0

        do ord=0,Legen_ord
            temp_sum = temp_sum + u_0(ord,i) * Le_poly(ord,0.0d0)
        enddo   

        sum_u_0(i) = temp_sum
    enddo

    !call Get_Result(i_bc,imax,x,sum_u_0)

    !Set the exact solution of the initial function
    do i=1 , imax
        exact_u_0(i)=Ini_u(x(i))
    enddo
    
    write(*,*) "======================Initial error results======================"
    call Get_L_2_Error(i_bc,imax,dx,exact_u_0,sum_u_0)
    call Get_L_Infty_Error(imax,exact_u_0,sum_u_0)
    write(*,*) "================================================================="
    write(*,*) " "

    !Compute cell interface contributions
    !((x(i)-dx(i)/2.0d0)-x(i))/dx(i)=-1.0d0/2.0d0: Left cell interface
    !((x(i)+dx(i)/2.0d0)-x(i))/dx(i)=+1.0d0/2.0d0: Right cell interface
    ! <<>> 이거 왜 필요하지? 무슨 의미가 있는 거지? 
    ! 나중에 값을 넣을건데,, 경계에서의 u_h값이 필요하다 그러면 다시 프로젝션을 시켜야하는데, 그때 계산을 위한..
    do ord=0 , Legen_ord
        l_b(ord)=Le_poly(ord , -1.0d0/2.0d0)
        r_b(ord)=Le_poly(ord , +1.0d0/2.0d0)
    enddo

    return
  !-------------------------Calculations have finished-------------------------!
endsubroutine

subroutine Get_Timestep(DG_ord,i_bc,imax,dx,flag_1,time_out,CFL,dt)
  implicit none
  
  integer,intent(in) :: DG_ord
  integer,intent(in) :: i_bc,imax
  double precision,intent(in) :: dx(-i_bc+1:imax+i_bc)
  double precision,intent(in) :: flag_1,time_out,CFL

  double precision,intent(out) :: dt
  !-------------------------Calculations have started--------------------------!

  ! 현재 CFL 은 0.1로 설정되어 있음 
  dt=CFL*maxval(dx)/1.0d0

  !Smooth case
  if(DG_ord > 4)then
    dt=CFL*(maxval(dx)**(DG_ord/4))/1.0d0
  endif

  !Nonsmooth case
  !if(DG_ord.ge.3)then
  !  dt=CFL*(maxval(dx)**(DG_ord/3))/1.0d0
  !endif

  if( flag_1 + dt >= time_out)then
    dt= time_out - flag_1
  endif

  return
  !-------------------------Calculations have finished-------------------------!
endsubroutine

! <<>> 이거는 그냥 고스트셀 채우는 용도. periodic 하게 채운다.
subroutine Get_Boundary(Legen_ord,i_bc,imax,u)
    implicit none
    
    integer,intent(in) :: Legen_ord
    integer,intent(in) :: i_bc,imax

    double precision,intent(inout) :: u(0:Legen_ord,-i_bc+1:imax+i_bc)

    !Counters
    integer :: ord
    integer :: i
    !-------------------------Calculations have started--------------------------!

    !Apply periodic conditions to ghost cells
    do i=1,i_bc
        do ord=0,Legen_ord 
        u(ord,-i+1)=u(ord,imax-i+1)
        u(ord,imax+i)=u(ord,i)
        enddo
    enddo

    return
    !-------------------------Calculations have finished-------------------------!
endsubroutine

! 아직 안봄 
subroutine Get_Minmod(TVB_M,Legen_ord,l_b,r_b,i_bc,imax,dx,u)
  implicit none
  
  double precision,intent(in) :: TVB_M
  integer,intent(in) :: Legen_ord  
  double precision,intent(in) :: l_b(0:Legen_ord),r_b(0:Legen_ord)
  integer,intent(in) :: i_bc,imax
  double precision,intent(in) :: dx(1:imax)
  double precision,intent(inout) :: u(0:Legen_ord,-i_bc+1:imax+i_bc)

  !Counters
  integer :: ord
  integer :: i

  !First-order undivived difference of the mean
  double precision,dimension(-i_bc+1:imax+i_bc) :: Forward_u_dx

  !Variables for the limiting process
  double precision :: temp
  double precision,dimension(0:2) :: s !Sign values of contributions
  double precision :: sgn_s
  double precision :: u_left,u_right
  double precision :: min_candi_u !A minmum of contributions
  double precision,dimension(0:2) :: candi_u
  !-------------------------Calculations have started--------------------------!

  !Compute first-order undivived difference of the mean
  do i=-1,imax+1
    Forward_u_dx(i)=u(0,i+1)-u(0,i)
  enddo
  
  !Apply the minmod limiting procedure to detect trouble cells
  do i=0,imax+1

  !============================================================================!
  !Limiting procedure at the right boundary (Left flux direction)
  !============================================================================!
  
  !Set up candidates and its signs for the minmod function
    candi_u(1)=Forward_u_dx(i)
    candi_u(2)=Forward_u_dx(i-1)
    s(1)=sign(1.0d0,candi_u(1))
    s(2)=sign(1.0d0,candi_u(2))

    u_left=0.0d0
    u_right=0.0d0
    do ord=0,Legen_ord
      u_left=u_left+u(ord,i)*r_b(ord)
      u_right=u_right+u(ord,i)*l_b(ord)
    enddo

  !Set up the main candidate for the minmod function
    candi_u(0)=u_left-u(0,i)

  !Enact the Shu's minmod function  
    if(dabs(candi_u(0)).le.TVB_M)then
  !Do nothing
      goto 10
    else
  !Enact the classical minmod function
      s(0)=sign(1.0d0,candi_u(0))
      sgn_s=(s(0)+s(1)+s(2))/3.0d0
    
      min_candi_u=min(dabs(candi_u(0)),dabs(candi_u(1)))
      min_candi_u=min(min_candi_u,dabs(candi_u(2)))

      do ord=0,Legen_ord
        if(dabs(sgn_s).eq.1.0d0)then
          u(ord,i)=sgn_s*min_candi_u
        else
          u(ord,i)=0.0d0
        endif
      enddo
    endif
    10 continue

  !============================================================================!
  !Limiting procedure at the left boundary (Right flux direction)
  !============================================================================!
   
  !Set up the main candidate for the minmod function
    candi_u(0)=u(0,i)-u_right

  !Enact the Shu's minmod function  
    if(dabs(candi_u(0)).le.TVB_M)then
  !Do nothing
      goto 20
    else
  !Enact the classical minmod function
      s(0)=sign(1.0d0,candi_u(0))
      sgn_s=(s(0)+s(1)+s(2))/3.0d0
   
      min_candi_u=min(dabs(candi_u(0)),dabs(candi_u(1)))
      min_candi_u=min(min_candi_u,dabs(candi_u(2)))

      do ord=0,Legen_ord
        if(dabs(sgn_s).eq.1.0d0)then
          u(ord,i)=sgn_s*min_candi_u
        else
          u(ord,i)=0.0d0
        endif
      enddo
    endif
    20 continue
  enddo

  return
  !-------------------------Calculations have finished-------------------------!
endsubroutine

! 아직 안봄 
subroutine Get_Minmod_WENO_JS(TVB_M,Legen_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,u)
  implicit none
  
  double precision,intent(in) :: TVB_M
  integer,intent(in) :: Legen_ord
  double precision,intent(in) :: Q(1:6,1:2) 
  double precision,intent(in) :: M(0:Legen_ord,0:Legen_ord)
  double precision,intent(in) :: d(0:Legen_ord,1:Legen_ord+2)
  double precision,intent(in) :: l_b(0:Legen_ord),r_b(0:Legen_ord)
  integer,intent(in) :: i_bc,imax
  double precision,intent(in) :: x(-i_bc+1:imax+i_bc),dx(-i_bc+1:imax+i_bc)
  double precision,intent(inout) :: u(0:Legen_ord,-i_bc+1:imax+i_bc)

  !Counters
  integer :: quad
  integer :: ord
  integer :: i

  !Legendre polynomials
  double precision,external :: Le_poly

  !Cell interface values
  double precision,dimension(0:1,0:imax+1) :: til_u

  !Lagrange coefficients and polynomials
  double precision,external :: Lag_C,Lag_P 
  double precision,dimension(1:Legen_ord+2) :: p_0,p_1,p_2
  
  !Main WENO reconstruction variables
  integer,parameter :: n=0
  double precision :: eps
  double precision :: beta_0,beta_1,beta_2
  double precision :: alpha_0,alpha_1,alpha_2,alpha_sum
  double precision :: omega_0,omega_1,omega_2
  double precision,dimension(1:Legen_ord+2) :: WENO
  !-------------------------Calculations have started--------------------------!

  !Set up the epsilon parameter
  eps=maxval(dx)**n
  !eps=1d-6

  !Set up the sum of degrees of freedom except the mean
  do i=0,imax+1
    do ord=0,Legen_ord
      til_u(0,i)=til_u(0,i)+u(ord,i)*l_b(ord)
      til_u(1,i)=til_u(1,i)+u(ord,i)*r_b(ord)
    enddo
    til_u(0,i)=dabs(u(0,i)-til_u(0,i))
    til_u(1,i)=dabs(til_u(1,i)-u(0,i))
  enddo
  
  !Apply the WENO limiting procedure by using Shu's minmod trouble cell detector
  do i=0,imax+1
    if((til_u(0,i).le.TVB_M).and.(til_u(1,i).le.TVB_M))then
   !Do nothing
      goto 20
    else
   !Enact WENO limiter
   !Fifth-order WENO scheme at four point Gauss-Lobatto quadrature
      do quad=1,Legen_ord+2
   !Set up adaptive local polynomials
   ![-2,-1,0]
        p_0(quad)=Lag_P(-2,0,Legen_ord,i_bc,imax,i,x,dx,Q(quad,1),u)
   ![-1,0,1]
        p_1(quad)=Lag_P(-1,1,Legen_ord,i_bc,imax,i,x,dx,Q(quad,1),u)
   ![0,1,2]
        p_2(quad)=Lag_P(0,2,Legen_ord,i_bc,imax,i,x,dx,Q(quad,1),u) 
  
   !Compute nonlinear normalized weights by using WENO scheme
   !Local smoothness indicator
        beta_0=(13.0d0/12.0d0)*(1.0d0*u(0,i-2)&
                              &-2.0d0*u(0,i-1)&
                              &+1.0d0*u(0,i+0))**2&
               &+(1.0d0/4.0d0)*(1.0d0*u(0,i-2)&
                              &-4.0d0*u(0,i-1)&
                              &+3.0d0*u(0,i+0))**2
          
        beta_1=(13.0d0/12.0d0)*(1.0d0*u(0,i-1)&
                              &-2.0d0*u(0,i+0)&
                              &+1.0d0*u(0,i+1))**2&
               &+(1.0d0/4.0d0)*(1.0d0*u(0,i-1)&
                              &+0.0d0*u(0,i+0)&
                              &-1.0d0*u(0,i+1))**2

        beta_2=(13.0d0/12.0d0)*(1.0d0*u(0,i+0)&
                              &-2.0d0*u(0,i+1)&
                              &+1.0d0*u(0,i+2))**2&
               &+(1.0d0/4.0d0)*(3.0d0*u(0,i+0)&
                              &-4.0d0*u(0,i+1)&
                              &+1.0d0*u(0,i+2))**2
   !Unnormalized weights
        alpha_0=d(0,quad)/(eps+beta_0)**2
        alpha_1=d(1,quad)/(eps+beta_1)**2
        alpha_2=d(2,quad)/(eps+beta_2)**2
        alpha_sum=alpha_0+alpha_1+alpha_2

   !Normalized weights
        omega_0=alpha_0/alpha_sum
        omega_1=alpha_1/alpha_sum
        omega_2=alpha_2/alpha_sum

   !Convex summation
        WENO(quad)=omega_0*p_0(quad)&
                 &+omega_1*p_1(quad)&
                 &+omega_2*p_2(quad)
      enddo
   !Reconstruct degrees of freedom by comparing degrees of polynomials
      u(1,i)=M(1,1)*(Q(1,2)*WENO(1)*Le_Poly(1,Q(1,1))&
                   &+Q(2,2)*WENO(2)*Le_Poly(1,Q(2,1))&
                   &+Q(3,2)*WENO(3)*Le_Poly(1,Q(3,1))&
                   &+Q(4,2)*WENO(4)*Le_Poly(1,Q(4,1)))

      u(2,i)=M(2,2)*(Q(1,2)*WENO(1)*Le_Poly(2,Q(1,1))&
                   &+Q(2,2)*WENO(2)*Le_Poly(2,Q(2,1))&
                   &+Q(3,2)*WENO(3)*Le_Poly(2,Q(3,1))&
                   &+Q(4,2)*WENO(4)*Le_Poly(2,Q(4,1)))
    endif
    20 continue
  enddo
  
  return
  !-------------------------Calculations have finished-------------------------!
endsubroutine

subroutine Get_Residual(TVB_M,Legen_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx&
                                                                    &,u,F_Flux)
    implicit none

    double precision,intent(in) :: TVB_M
    integer,intent(in) :: Legen_ord
    double precision,intent(in) :: Q(1:6,1:2) 
    double precision,intent(in) :: M(0:Legen_ord,0:Legen_ord)
    double precision,intent(in) :: d(0:Legen_ord,1:Legen_ord+3)
    double precision,intent(in) :: l_b(0:Legen_ord),r_b(0:Legen_ord)
    integer,intent(in) :: i_bc,imax
    double precision,intent(in) :: x(1-i_bc+1:imax+i_bc),dx(-i_bc+1:imax+i_bc)
    double precision,intent(inout) :: u(0:Legen_ord,-i_bc+1:imax+i_bc)

    double precision,intent(out) :: F_Flux(0:Legen_ord,1:imax)

    !Counters
    integer :: ord
    integer :: i

    !Flux function
    double precision,external :: Flux

    !Legendre polynomials
    double precision,external :: Le_poly,D_Le_poly

    !Cell interface values
    double precision :: temp
    double precision,dimension(0:imax+1) :: u_m
    double precision,dimension(-1:imax) :: u_p
    double precision,dimension(0:1,0:imax+1) :: til_m
    
    !Flux splitting variables
    double precision :: F_max
    double precision,dimension(0:imax) :: hat_F

    !Boundary integration variable
    double precision,dimension(0:Legen_ord,1:imax) :: F_Flux_B

    !Volume integration variables
    double precision :: uh_1,uh_2,uh_3,uh_4,uh_5,uh_6
    double precision,dimension(0:Legen_ord,1:imax) :: F_Flux_V

    !Lagrange polynomials
    double precision,external :: Lag_P
    !-------------------------Calculations have started--------------------------!

    !Calculate the local residual

    !============================================================================!
    !Volume integral part
    !============================================================================!

    !Compute the volume integral by using quadratures
    do i=1,imax
        uh_1=0.0d0
        uh_2=0.0d0
        uh_3=0.0d0
        uh_4=0.0d0
        uh_5=0.0d0
        uh_6=0.0d0
        ! 각 셀에서 u_h 찾기 
        do ord=0,Legen_ord
        uh_1=uh_1+u(ord,i)*Le_poly(ord,Q(1,1))
        uh_2=uh_2+u(ord,i)*Le_poly(ord,Q(2,1))
        uh_3=uh_3+u(ord,i)*Le_poly(ord,Q(3,1))
        uh_4=uh_4+u(ord,i)*Le_poly(ord,Q(4,1))
        uh_5=uh_5+u(ord,i)*Le_poly(ord,Q(5,1))
        uh_6=uh_6+u(ord,i)*Le_poly(ord,Q(6,1))
        enddo

        do ord=0,Legen_ord
          F_Flux_V(ord,i)=Flux(uh_1)*D_Le_poly(ord,Q(1,1))*Q(1,2)&
                          &+Flux(uh_2)*D_Le_poly(ord,Q(2,1))*Q(2,2)&
                          &+Flux(uh_3)*D_Le_poly(ord,Q(3,1))*Q(3,2)&
                          &+Flux(uh_4)*D_Le_poly(ord,Q(4,1))*Q(4,2)&
                          &+Flux(uh_5)*D_Le_poly(ord,Q(5,1))*Q(5,2)&
                          &+Flux(uh_6)*D_Le_poly(ord,Q(6,1))*Q(6,2)
        enddo
    enddo
    
    !============================================================================!
    !============================================================================!

    !============================================================================!
    !Boundary integral part
    !============================================================================!

    !Impose the periodic boundary condition to degrees of freedom
    call Get_Boundary(Legen_ord,i_bc,imax,u)
    
    !Apply limiting process   
    !<<>> 일단 밑에 서브루틴들은 u를 limiting 시키는데, 나는 아직 모르겠음 
    !call Get_Minmod(TVB_M,Legen_ord,l_b,r_b,i_bc,imax,dx,u)
    !call Get_Minmod_WENO_JS(TVB_M,Legen_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,u)
    
    !For given degrees of freedom, compute convex summations at cell interfaces
    do i=0 , imax+1
    !Left limit (-) at the right boundary 
      temp=0.0d0
      ! 여기서 미리 구한 의미가 있는거죠~ 굿 
      do ord=0,Legen_ord
        temp=temp+u(ord,i)*r_b(ord)
      enddo
      ! 이거 왜 이렇게 해놨지? u_m??? why? 
      u_m(i+0)=temp

    !Right limit (+) at the left boundary
      temp=0.0d0
      do ord=0,Legen_ord
        temp=temp+u(ord,i)*l_b(ord)
      enddo
      u_p(i-1)=temp
    enddo

    !Take the apropriate numerical trace 
    !As the averaged monotone flux, apply the Lax-Friedrich flux splitting method
    F_max=1.0d0
    do i=0,imax
        hat_F(i)=Flux(u_m(i))+Flux(u_p(i))
        hat_F(i)=(hat_F(i)+F_max*(u_m(i)-u_p(i)))/2.0d0
    enddo

    !Compute the boundary integral by using the fundamental theorem of calculus
    do i=1,imax
      do ord=0,Legen_ord  
          F_Flux_B(ord,i)=hat_F(i-1)*l_b(ord)-hat_F(i+0)*r_b(ord)
      enddo
    enddo

    !============================================================================!
    !============================================================================!
    
    !Compute local residuals by using calculated integral values
    do i=1,imax
      do ord=0,Legen_ord
          F_Flux(ord,i)=M(ord,ord)*(F_Flux_B(ord,i)+F_Flux_V(ord,i))/dx(i)
      enddo
    enddo
    
    return
    !-------------------------Calculations have finished-------------------------!
endsubroutine

subroutine Get_RK(RK_ord,TVB_M,Legen_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx,dt&
                                                           &,old_u,new_u)
  implicit none
  
  integer,intent(in) :: RK_ord
  double precision,intent(in) :: TVB_M
  integer,intent(in) :: Legen_ord
  double precision,intent(in) :: Q(1:6,1:2) 
  double precision,intent(in) :: M(0:Legen_ord,0:Legen_ord)
  double precision,intent(in) :: d(0:Legen_ord,1:Legen_ord+3)
  double precision,intent(in) :: l_b(0:Legen_ord),r_b(0:Legen_ord)
  integer,intent(in) :: i_bc,imax
  double precision,intent(in) :: x(-i_bc+1:imax+i_bc),dx(-i_bc+1:imax+i_bc)
  double precision,intent(in) :: dt
  double precision,intent(inout) :: old_u(0:Legen_ord,-i_bc+1:imax+i_bc)

  double precision,intent(out) :: new_u(0:Legen_ord,1:imax)

  !Counters
  integer :: ord
  integer :: i

  !Calculation variables
  double precision,dimension(0:Legen_ord,-i_bc+1:imax+i_bc,1:3) :: temp_u
  double precision,dimension(0:Legen_ord,1:imax) :: F_Flux
  !-------------------------Calculations have started--------------------------!

  if(RK_ord.eq.1)then
  !1st-order TVD RK solver
  !============================================================================!
  !Stage 1
  !============================================================================!

  !Compute the residual
    call Get_Residual(TVB_M,Legen_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx&
                                                    &,old_u,F_Flux)

    do i=1,imax
      do ord=0,Legen_ord
        new_u(ord,i)=+1.0d0*old_u(ord,i)&
                    &+1.0d0*dt*F_Flux(ord,i)
      enddo
    enddo

  !============================================================================!
  !============================================================================!
  elseif(RK_ord.eq.2)then
  !2nd-order TVD RK solver
  !============================================================================!
  !Stage 1
  !============================================================================!

  !Compute the residual
    call Get_Residual(TVB_M,Legen_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx&
                                                    &,old_u,F_Flux)

    do i=1,imax
      do ord=0,Legen_ord
        temp_u(ord,i,1)=+1.0d0*old_u(ord,i)&
                       &+1.0d0*dt*F_Flux(ord,i)
      enddo
    enddo

  !============================================================================!
  !============================================================================!

  !============================================================================!
  !Stage 2
  !============================================================================!

  !Compute the residual
    call Get_Residual(TVB_M,Legen_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx&
                                            &,temp_u(:,:,1),F_Flux)
    do i=1,imax
      do ord=0,Legen_ord
        new_u(ord,i)=+1.0d0*old_u(ord,i)&
                    &+1.0d0*temp_u(ord,i,1)&
                    &+1.0d0*dt*F_Flux(ord,i)
        new_u(ord,i)=new_u(ord,i)/2.0d0
      enddo
    enddo

  !============================================================================!
  !============================================================================!
  !3rd-order TVD RK solver
  elseif(RK_ord.eq.3)then
  !============================================================================!
  !Stage 1
  !============================================================================!

  !Compute the residual
    call Get_Residual(TVB_M,Legen_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx&
                                                    &,old_u,F_Flux)
    do i=1,imax
      do ord=0,Legen_ord
        temp_u(ord,i,1)=+1.0d0*old_u(ord,i)&
                       &+1.0d0*dt*F_Flux(ord,i)
      enddo
    enddo

  !============================================================================!
  !============================================================================!

  !============================================================================!
  !Stage 2
  !============================================================================!

  !Compute the residual
    call Get_Residual(TVB_M,Legen_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx&
                                            &,temp_u(:,:,1),F_Flux)
    do i=1,imax
      do ord=0,Legen_ord
        temp_u(ord,i,2)=+3.0d0*old_u(ord,i)&
                        &+1.0d0*temp_u(ord,i,1)&
                        &+1.0d0*dt*F_Flux(ord,i)
        temp_u(ord,i,2)=temp_u(ord,i,2)/4.0d0
      enddo
    enddo

  !============================================================================!
  !============================================================================!

  !============================================================================!
  !Stage 3
  !============================================================================!

  !Compute the residual
    call Get_Residual(TVB_M,Legen_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx&
                                            &,temp_u(:,:,2),F_Flux)
    do i=1,imax
      do ord=0,Legen_ord
        new_u(ord,i)=+1.0d0*old_u(ord,i)&
                    &+2.0d0*temp_u(ord,i,2)&
                    &+2.0d0*dt*F_Flux(ord,i)
        new_u(ord,i)=new_u(ord,i)/3.0d0
      enddo
    enddo
  
  !============================================================================!
  !============================================================================!
  elseif(RK_ord.eq.4)then
  !4th-order non-TVD RK solver
  !============================================================================!
  !Stage 1
  !============================================================================!
  
  !Compute the residual
    call Get_Residual(TVB_M,Legen_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx&
                                                    &,old_u,F_Flux)
    do i=1,imax
      do ord=0,Legen_ord
        temp_u(ord,i,1)=+1.0d0*old_u(ord,i)&
                       &+0.5d0*dt*F_Flux(ord,i)
      enddo
    enddo

  !============================================================================!
  !============================================================================!

  !============================================================================!
  !Stage 2
  !============================================================================!
  
  !Compute the residual
    call Get_Residual(TVB_M,Legen_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx&
                                            &,temp_u(:,:,1),F_Flux)
    do i=1,imax
      do ord=0,Legen_ord
        temp_u(ord,i,2)=+1.0d0*old_u(ord,i)&
                       &+0.5d0*dt*F_Flux(ord,i)
      enddo
    enddo

  !============================================================================!
  !============================================================================!

  !============================================================================!
  !Stage 3
  !============================================================================!

  !Compute the residual
    call Get_Residual(TVB_M,Legen_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx&
                                            &,temp_u(:,:,2),F_Flux) 
    do i=1,imax
      do ord=0,Legen_ord
        temp_u(ord,i,3)=+1.0d0*old_u(ord,i)&
                       &+1.0d0*dt*F_Flux(ord,i)
      enddo
    enddo
  
  !============================================================================!
  !============================================================================!

  !============================================================================!
  !Stage 4
  !============================================================================!

  !Compute the residual
    call Get_Residual(TVB_M,Legen_ord,Q,M,d,l_b,r_b,i_bc,imax,x,dx&
                                            &,temp_u(:,:,3),F_Flux)
    do i=1,imax
      do ord=0,Legen_ord
        new_u(ord,i)=-1.0d0*old_u(ord,i)&
                    &+1.0d0*temp_u(ord,i,1)&
                    &+2.0d0*temp_u(ord,i,2)&
                    &+1.0d0*temp_u(ord,i,3)&
                    &+0.5d0*dt*F_Flux(ord,i)
        new_u(ord,i)=new_u(ord,i)/3.0d0
      enddo
    enddo
  
  !============================================================================!
  !============================================================================!
  endif
  
  return
  !-------------------------Calculations have finished-------------------------!
endsubroutine

subroutine Get_L_1_Error(i_bc,imax,dx,exact,numerical)
  implicit none

  integer,intent(in) :: i_bc,imax
  double precision,intent(in) :: dx(1-i_bc+1:imax+i_bc)
  double precision,intent(in) :: exact(1:imax),numerical(1:imax)

  !Counter
  integer :: i

  !Error calculation variables
  double precision,dimension(1:imax) :: e
  !-------------------------Calculations have started--------------------------!

  !Compute absolute differences
  do i=1,imax
    e(i)=dabs(exact(i)-numerical(i))
    e(i)=e(i)*dx(i)
  enddo

  !Compute the average sum
  write(*,*) "L_1 error:",sum(e)

  return
  !-------------------------Calculations have finished-------------------------!
endsubroutine

subroutine Get_L_2_Error(i_bc,imax,dx,exact,numerical)
  implicit none

  integer,intent(in) :: i_bc,imax
  double precision,intent(in) :: dx(-i_bc+1:imax+i_bc)
  double precision,intent(in) :: exact(1:imax),numerical(1:imax)

  !Counter
  integer :: i

  !Error calculation variables
  double precision,dimension(1:imax) :: e
  !-------------------------Calculations have started--------------------------!

  !Compute square differences
  do i=1,imax
    e(i)=(exact(i)-numerical(i))**2
    e(i)=e(i)*dx(i)
  enddo

  !Compute the square root average sum
  write(*,*) "L_2 error:",dsqrt(sum(e))

  return
  !-------------------------Calculations have finished-------------------------!
endsubroutine

subroutine Get_L_Infty_Error(imax,exact,numerical)
  implicit none

  integer,intent(in) :: imax
  double precision,intent(in) :: exact(1:imax),numerical(1:imax)

  !Counters
  integer :: i

  !Error calculation variables
  double precision,dimension(1:imax) :: e
  !-------------------------Calculations have started--------------------------!

  !Compute absolute differences
  do i=1,imax
    e(i)=dabs(exact(i)-numerical(i))
  enddo

  !Take the maximum value
  write(*,*) "L_infty error:",maxval(e)
  write(*,*) "Maximum error location:",maxloc(e)

  return
  !-------------------------Calculations have finished-------------------------!
endsubroutine

subroutine Get_Result(i_bc,imax,x,u)
  implicit none

  integer,intent(in) :: i_bc,imax
  double precision,intent(in) :: x(-i_bc+1:imax+i_bc)
  double precision,intent(in) :: u(1:imax)

  !Counter
  integer :: i
  !-------------------------Calculations have started--------------------------!
  
  do i=1,imax
    write(1,*) x(i),u(i)
  enddo

  return
  !-------------------------Calculations have finished-------------------------!
endsubroutine

double precision function Ini_u(x)
  implicit none

  double precision,intent(in) :: x

  !Useful constant
  double precision,parameter :: pi=4.0d0*datan(1.0d0)
  !-------------------------Calculations have started--------------------------!
  
  !Set up the initial function
  !Smooth functions
  !Ini_u=dcos(pi*x)**1
  Ini_u=dcos(pi*x)**4
  !Ini_u=dsin(pi*x)**1
  !Ini_u=dsin(pi*x-dsin(pi*x)/pi)
  !Ini_u = 1.0d0/4.0d0 + 1.0d0/2.0d0 * dsin(pi*x)

  !Nonsmooth function
  !Square wave
  !if((x.gt.-0.5d0).and.(x.lt.0.5d0))then
  !  Ini_u=1.0d0
  !else
  !  Ini_u=0.0d0
  !endif

  return
  !-------------------------Calculations have finished-------------------------!
endfunction

double precision function Ext_u(x,time_out)
  implicit none

  double precision,intent(in) :: x
  double precision,intent(in) :: time_out

  !Useful constant
  double precision,parameter :: pi=4.0d0*datan(1.0d0)
  !-------------------------Calculations have started--------------------------!
  
  !Set up the exact solution
  !Smooth functions
  !Ext_u=dcos(pi*(x-time_out))**1
  Ext_u=dcos(pi*(x-time_out))**4
  !Ext_u=dsin(pi*(x-time_out))**1
  !Ext_u=dsin(pi*(x-time_out)-dsin(pi*(x-time_out))/pi)

  !Nonsmooth function
  !Square wave
  !if((x.gt.-0.5d0).and.(x.lt.0.5d0))then
  !  Ext_u=1.0d0
  !else
  !  Ext_u=0.0d0
  !endif

  return
  !-------------------------Calculations have finished-------------------------!
endfunction

double precision function Flux(u)
  implicit none

  double precision,intent(in) :: u
  !-------------------------Calculations have started--------------------------!
  
  !Set up the flux
  Flux=  u


  return
  !-------------------------Calculations have finished-------------------------!
endfunction

double precision function Le_poly(ord, x) 
    implicit none

    integer,intent(in) :: ord 
    double precision,intent(in) :: x 

    select case(ord) 
    case(0) 
        Le_poly = 1 
    case(1) 
        Le_poly = x**1
    case(2)
        Le_poly = x**2 - (1.0d0/12.0d0) 
    case(3) 
        Le_poly = x**3 - (3.0d0 /20.0d0) * x 
    case(4) 
        Le_poly = (x**4) - ((3.0d0 / 14.0d0) * x**2) + (3.0d0 / 560.0d0)
    case default 
        Le_poly = 0.0d0
    end select

    return

endfunction

double precision function D_Le_poly(ord, x)
    implicit none

    integer, intent(in) :: ord
    double precision, intent(in) :: x

    !-------------------------Calculations have started--------------------------!

    select case (ord)
    case (0)
        D_Le_poly = 0.0d0
    case (1)
        D_Le_poly = 1.0d0
    case (2)
        D_Le_poly = 2.0d0 * x
    case (3)
        D_Le_poly = 3.0d0 * x**2 - (3.0d0 / 20.0d0)
    case (4)
        D_Le_poly = 4.0d0 * x**3 - (3.0d0 / 7.0d0) * x
    case default
        D_Le_poly = 0.0d0
    end select

    return
    !-------------------------Calculations have finished--------------------------!
endfunction 

double precision function Lag_C(l,k_l,k_r,i_bc,imax,i,x,dx,x_p)
  implicit none

  integer,intent(in) :: l !Cell counter
  integer,intent(in) :: k_l !Minimal index of the left cell
  integer,intent(in) :: k_r !Maximal index of the right cell
  integer,intent(in) :: i_bc,imax,i
  double precision,intent(in) :: x(-i_bc+1:imax+i_bc)
  double precision,intent(in) :: dx(-i_bc+1:imax+i_bc) 
  double precision,intent(in) :: x_p !Present location
  
  !Counters
  integer :: j,k,s

  !Calculation variables
  double precision :: x_n !Neighbor cell
  double precision :: c1,c2,c3
  !-------------------------Calculations have started--------------------------!
  
  !Initialization
  Lag_C=0.0d0

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
	  Lag_C=Lag_C+c2/c1
	enddo  
	Lag_C=Lag_C*dx(i+l)
  
  return
  !-------------------------Calculations have finished-------------------------!
endfunction

double precision function Lag_P(k_l,k_r,Legen_ord,i_bc,imax,i,x,dx,x_p,u)
  implicit none

  integer,intent(in) :: k_l !Minimal index of the left cell
  integer,intent(in) :: k_r !Maximal index of the right cell
  integer,intent(in) :: Legen_ord
  integer,intent(in) :: i_bc,imax
  integer,intent(in) :: i !Present index of the cell
  double precision,intent(in) :: x(1:imax)
  double precision,intent(in) :: dx(1:imax) 
  double precision,intent(in) :: x_p !Present location
  double precision,intent(inout) :: u(0:Legen_ord,-i_bc+1:imax+i_bc) 

  !Counter
  integer :: l

  !Lagrange coefficients
  double precision,external :: Lag_C
  !-------------------------Calculations have started--------------------------!
  
  !Initialization
  Lag_P=0.0d0

  !Compute the convex summation of Lagrange basis by using mean values
  do l=k_l,k_r
    Lag_P=Lag_P+Lag_C(l,k_l,k_r,i_bc,imax,i,x,dx,x_p)*u(0,l+i)
 	enddo

 	return
  !-------------------------Calculations have finished-------------------------!
endfunction