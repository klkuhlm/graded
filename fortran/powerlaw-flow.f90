!
! Copyright (c) 2021-2022 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
!

program powerlaw

  use constants, only : DP, PISQ, ULP !, PIOV2 ! PI
  use utility, only : logspace !, operator(.SUM.)
  use invlap
  use complex_bessel, only : cbesk

  implicit none

  ! SPECQ is [user-specified] specified flowrate (Q)?
  ! CALCT is [user-specified] calculate time vector?
  ! HEADER is [user-specified] output header (echo of input) in output?
  ! DERIV is [user-specified] compute derivative of solution wrt time?
  ! WBSTORAGE is [computed] whether wellbore storage is used?
  ! ATSOURCE is [computed] whether solution is inside source boreholes?
  ! COMPUTEP is [user-specified] whether computing pressure solution for specified
  !  pressure (only makes sense outside source borehole)?
  logical :: SPECQ, CALCT, HEADER, DERIV, WBSTORAGE, ATSOURCE, COMPUTEP

  ! kappa is permeability exponent
  ! eta is porosity exponent
  ! mdim is the dimension factor (0=cartesian, 1=cylindrical, 2=spherical)
  ! alpha = (1 + kappa - m) / 2
  ! gamma = (2 + kappa - eta) / 2
  ! nu = alpha / gamma
  ! lambda is WR interporosity transfer coefficient
  real(DP) :: eta, alpha, gamma, kappa, mdim, nu, lambda
  
  ! omega is WR fracture domain storage ratio (omega_f in Kuhlman et al. 2015)
  ! omomega = 1-omega is matrix domain storage ratio (omega_i in Kuhlman et al. 2015)
  ! sigma is wellbore storage factor
  ! rD is radius of calcuation point
  real(DP) :: sigma, omega, omomega, dummy, rD

  ! min and maximum log t for computing solution
  ! Amos bessl function solutions are computed scaled, unscale is to convert back
  real(DP) :: minlogt, maxlogt, tee, r_term

  ! NT num times, NP num laplace parameters, NA number Kazemi series terms, PORO is model flag
  integer :: NT, NP, i, j, n, PORO, NA, ii
  integer :: nz, ierr, timeFlag
  real(DP), allocatable :: timePar(:)
  character(256) :: out_file, time_file
  real(DP), dimension(2) :: order(-1:0)
  complex(DP), dimension(2) :: K(-1:0)
  complex(DP), dimension(1) :: Ksigma
  real(DP) :: ft, dft
  real(DP), allocatable :: t(:), ivsq(:)
  complex(DP), allocatable :: p(:), fp(:), barft(:), beta(:), arg(:,:), betgam(:), unscale(:)
  type(invLaplace) :: lap

  open(unit=20,file='powerlaw.in',action='read',status='old')

  ! SPECQ = specified Q? (False => specified head)
  ! PORO = 1=single, 2=warren-root, 3=kazemi (multiporosity), 4=analytical sum Kazemi
  ! HEADER = print header? (easier to read in Matlab w/o header)
  ! DERIV = compute log-time derivative?
  read(20,*,iostat=ierr) SPECQ, PORO, HEADER, DERIV, COMPUTEP

  if (ierr /= 0) then
     write(*,*) 'error reading line 1:',SPECQ,PORO,HEADER,DERIV,COMPUTEP
     stop 100
  end if

  read(20,*,iostat=ierr) eta, kappa, mdim, NA
  if (ierr /= 0) then
     write(*,*) 'error reading line 2: ',eta,kappa,mdim,NA
     stop 102
  end if

  ! sigma, wellbore storage coefficient (0= no WB storage)
  ! lambda, omega: inter-porosity exchange coeff, fracture porosity ratio
  ! rD: radius for calculation point (anything <= 1 is in source well)
  read(20,*,iostat=ierr) sigma, lambda, omega, dummy, rD
  if (ierr /= 0) then
     write(*,*) 'error reading line 3: ',sigma,lambda,omega,"not used",rD
     stop 103
  end if

  ! calcT = compute time vector (read in vector if False)
  ! min/maxlogt = like used in logspace() numpy function
  ! nt = number of times in vector
  ! time_file = filename to read times from (if calcT == False)
  read(20,*,iostat=ierr) CALCT, minlogt, maxlogt, NT, time_file
  if (ierr /= 0) then
     write(*,*) 'error reading line 4: ',CALCT,minlogt,maxlogt,NT,trim(time_file)
     stop 104
  end if

  ! deHoog Laplace transform parameters
  read(20,*,iostat=ierr) lap%alpha, lap%tol, lap%M
  if (ierr /= 0) then
     write(*,*) 'error reading line 5: ',lap%alpha,lap%tol,lap%M
     stop 105
  end if

  ! filename for output
  read(20,*,iostat=ierr) out_file
  if (ierr /= 0) then
     write(*,*) 'error reading line 6: ',trim(out_file)
     stop 106
  end if

  ! integer flag indicating time behavior (see below)
  read(20,*,iostat=ierr) timeFlag
  if (ierr /= 0) then
     write(*,*) 'error reading line 7: ',timeFlag
     stop 107
  end if

  allocate(timePar(2*abs(timeFlag)+1))
  timePar = -999.

  ! timeFlag == 0 simple step on at t=0 behavior (timePar not used)
  ! timeFlag < 0 piecewise constant arbitrary time behavior
  ! timeFlag > 0 piecewise linear arbitrary time behevior

  ! for piecewise constant|linear behavior, there should be 2n+1 values in
  ! timePar for n discrete sections.

  ! first n values are initial times for a section, n+1 is final time
  ! n+2:2n+1 are the strengths associated with each initial time.
  ! piecewise constant means horizontal lines until next time, while
  ! piecwise linear means a straight line connecting points.

  ! EXAMPLE
  ! to specify a unit-height tent function on t={1,3}, centered on t=2:
  ! timeFlag = +2
  ! timePar = 1.0 2.0 3.0  0.0 1.0
  ! the first 3 values are times, the last 2 numbers are height, with
  ! the 6th 0.0 assumed to occur at t=3.0

  if (timeFlag /= 0) then
     backspace(20)
     read(20,*,iostat=ierr) timeFlag,timePar(:)
     if (ierr /= 0) then
        print *, 'error reading line 8: ', timePar(:)
        stop 108
     end if
  end if
  close(20)

  allocate(t(NT))
  if (CALCT) then
     if (minlogt > maxlogt .or. abs(minlogt - maxlogt) <= 0.0_DP) then
        print *, 'invalid (minlogt,maxlogt): ',minlogt,maxlogt
        stop 109
     end if
     t = logspace(minlogt,maxlogt,NT)
  else
     open(unit=50,file=trim(time_file),action='read',status='old')
     do i=1,NT
        read(50,*,iostat=ierr) t(i)
        if (ierr /= 0) then
           write(*,*) 'error reading time in row=',i,' from: ',trim(time_file)
           write(*,*) t(:)
           stop 110
        end if
     end do
     if (any(t <= 0.0)) then
        print *, 'ERROR: times read from ',trim(time_file),' must all be >0:'
        print *, 't vector: ',t
        stop 111
     end if
     if (any((t(2:NT)-t(1:NT-1)) <= 0.0_DP)) then
        print *, 'ERROR: times read from ',trim(time_file),&
             &' must be monotonically increasing:',t(2:NT)-t(1:NT-1)
        stop 112
     end if
  end if

  !! input error checking
  !! ====================================================

  if (PORO < 1 .or. PORO > 4) then
     write(*,*) 'PORO must be in {1,2,3,4}, not ',PORO
     stop 201
  else if (PORO == 3) then
     if (NA < 1) then
        write(*,*) 'NA (# terms in Kazemi series) must be >= 1, not: ', NA
        stop 202
     end if
     ! finite series multiporosity approximation to Kazemi
     allocate(ivsq(NA))
  end if

  if (rD < 1.0_DP) then
     print *, "*******************************"
     print *, 'WARNING: dimensionless radius < 1.0 (i.e., inside'//&
          &' well), resetting to 1.0 (borehole wall): ', rD
     print *, "*******************************"
     rD = 1.0_DP
  end if

  if (any([eta,kappa,mdim,sigma,omega,lambda,lap%alpha,lap%tol] < 0.0)) then
     ! non-positive input parameters
     print *, 'ERROR: parameters [eta,kappa,mdim,sigma,omega,lambda,'//&
          &'lap%alpha,lap%tol] must be >= 0.0: ',&
          &[eta,kappa,mdim,sigma,omega,lambda,lap%alpha,lap%tol]
     stop 203
  else if (omega > 1.0) then
     print *, 'ERROR: only 0<=omega<=1, not: ',omega
     stop 204
  end if

  if (lap%M < 1) then
     print *, 'ERROR: invalid number of terms in deHoog et al. approxation: ',lap%M
     stop 205
  end if

  ! matrix capacity for double-porosity case
  ! omega = omega_f
  ! omomega = omega_m
  omomega = 1.0_DP - omega

  ! compute logical flags
  if(sigma*ULP < spacing(0.0_DP)) then
     WBSTORAGE = .false. ! sigma = 0
  else
     WBSTORAGE = .true. ! sigma > 0
  end if

  if (abs(rD - 1.0_DP)*ULP < spacing(1.0_DP)) then
     ATSOURCE = .true. ! rD = 1
  else
     ATSOURCE = .false. ! rD > 1
  end if

  NP = 2*lap%M+1  ! number of laplace parameters per time
  allocate(p(NP),fp(NP),barft(NP),beta(NP),arg(NP,-1:0),betgam(NP),unscale(NP))

  alpha = (1.0_DP + kappa - mdim) / 2.0_DP
  gamma = (2.0_DP + kappa - eta) / 2.0_DP

  ! not nu = sqrt(alpha**2/gamma**2), which is always positive
  ! eta=3, kappa=0 is wrong for m={0,1,2}, seems like related to sign of gamma
  nu = alpha/gamma

  open(unit=30,file=trim(out_file),action='write',status='replace')

  if (HEADER) then
     ! echo key input to header for documentation (could be used as input file, save last line)
     write(30,'(A,L1,1X,I1,3(1X,L1))') '# ',SPECQ,PORO,HEADER,DERIV,COMPUTEP
     write(30,'(A,3(ES14.7,1X),I0)') '# ',eta,kappa,mdim,NA
     write(30,'(A,5(ES14.7,1X))') '# ',sigma,lambda,omega,dummy,rD
     write(30,'(A,L1,1X,2(ES14.7,1X),I0,1X,A)') '# ',CALCT,minlogt,maxlogt,NT,trim(time_file)
     write(30,'(A,2(ES14.7,1X),I0)') '# ',lap%alpha,lap%tol,lap%M
     write(30,'(2A)') '# ',trim(out_file)
     if (timeFlag == 0) then
        write(30,'(A,I0)') '# ',timeFlag
     else
        write(30,'(A,I0,1X,999(ES14.7,1X))') '# L7: ',timeFlag,timepar(:)
     end if
     write(30,'(A,3(1X,ES14.7),2(1X,L1))') '## alpha,gamma,nu,WBSTOR,ATSRC', &
          & alpha,gamma,nu,WBSTORAGE,ATSOURCE
     if (SPECQ) then
        write(30,'(A)') '#      t_D                 p_D(r_D)'//&
             &'             deriv wrt ln(t)'
     else
        write(30,'(A)') '#      t_D                dp_D(r_D) '//&
             &'            deriv wrt ln(t)'
     end if
     write(30,'(A)') '#======================================='//&
          &'=========================='
  end if

  if (PORO == 3) then
     ! things from Kazemi (in)finite sum that don't depend on p
     ivsq = PISQ / 4.0_DP * real([((2*ii - 1)**2, ii = 1, NA )],DP)
  end if

  do n = 1, NT
     tee = t(n)*2.0_DP
     p(1:NP) = deHoog_pvalues(tee,lap)
     barft(1:NP) = time_pvect(p(:),timePar(:),timeFlag)

     select case (PORO)
     case (1)
       ! single porosity, barker-like
       beta(1:NP) = 1.0_DP
     case (2)
       ! warren-root double porosity
       beta(1:NP) = omega + lambda / (lambda / omomega + p(:))
     case (3)
       ! approximation to Kazemi via Kuhlman et al. (2015)
       ! sum() is \bar{g}(p) (i.e., matrix memory kernel)
       beta(1:NP) = omega + sum(2.0_DP * lambda / &
            (spread(ivsq(1:NA) * lambda / omomega, dim=2, ncopies=NP) + &
            spread(p(1:NP), dim=1, ncopies=NA)), dim=1)
     case (4)
       ! infinite sum evaluated using Mathematica
       beta(1:NP) = omega +  &
            & sqrt(lambda * omomega / p(:)) * &
            & tanh(sqrt(p(:) * omomega / lambda))
     end select
     ! apply p and take square-root
     betgam(1:NP) = sqrt(p(:) * beta(:))  ! beta*gamma
     beta(1:NP) = betgam/gamma ! just beta (1/gamma^2 factored out of square root)

     ! K_nu(z) = K_-nu(z) for complex nu (DLMF 10.27.3)
     ! amos library only works for nu >= 0
     order(-1) = abs(nu - 1.0_DP)
     ! when nu = 1/2, both are equal (-> 1/sqrt(pi*t))
     order(0) = abs(nu)

     arg(1:NP,-1:0) = spread(beta(:),2,2)
     if (ATSOURCE) then
        r_term = 1.0_DP
        unscale = cmplx(1.0,0.0,DP)
     else ! for rD > 1
        if (SPECQ) then
           ! type II/III BC, beta*r**gamma in K_nu term
           r_term = rD ** alpha
           arg(1:NP,0) = beta(:) * rD ** gamma ! nu term in numerator
        else
           if (.not. COMPUTEP) then
              ! type I BC, predict flux with beta*r**gamma in K_nu-1 term
              r_term = rD ** (alpha + gamma - 1.0_DP)
              arg(1:NP,-1) = beta(:) * rD ** gamma ! nu-1 term in numerator
           else
              ! type I BC, predict pressure for rD>1
              r_term = rD ** alpha
              arg(1:NP,-1) = beta(:) * rD ** gamma ! nu term in numerator with r-dependnce
           end if
        end if
        unscale = exp(-beta(:)*(rD ** gamma - 1.0_DP)) ! account for different scaling when rD > 1
     end if

     if (SPECQ) then
        ! solutions for pressure at wellbore given specified flowrate
        do j = 1,NP
           do i = -1,0
              call cbesk(Z=arg(j,i), FNU=order(i), KODE=2, N=1, CY=K(i), NZ=nz, IERR=ierr)
              if (ierr /= 0 .and. ierr /= 3) then
                 print '(A,I0,A,"(",ES10.3,",",ES10.3,")",A,ES10.3)', &
                      & 'ERROR: cbesk K_nu(z) (Q) ierr:',&
                      & ierr,' z:',arg(j,i),' nu:',order(i)
                 stop 901
              end if
           end do

           ! wellbore storage _and_ computing solution away from source well
           if (WBSTORAGE) then ! non-zero wellbore storage
              if (.not. ATSOURCE) then ! solution away from source well
                 call cbesk(Z=beta(j), FNU=order(0), KODE=2, N=1, CY=Ksigma(1), NZ=nz, IERR=ierr)
                 if (ierr /= 0 .and. ierr /= 3) then
                    print '(A,I0,A,"(",ES10.3,",",ES10.3,")",A,ES10.3)', &
                         & 'ERROR: cbesk K_nu(beta) (Q) ierr:',&
                         & ierr,' z:',beta(j),' nu:',order(0)
                    stop 902
                 end if
              else ! ATSOURCE
                 ! wellbore storage and at source well
                 Ksigma(1) = K(0)
              end if
           else ! .not. WBSTORAGE
              Ksigma(1) = cmplx(0.0,0.0,DP)
           end if
           fp(j) = (barft(j) * K(0))/(betgam(j) * K(-1) + sigma*p(j) * Ksigma(1))
        end do

     else ! .not. SPECQ
        if (.not. COMPUTEP) then
           ! solutions for gradient given specified pressure
           do j = 1,NP
              do i = -1,0
                 call cbesk(Z=arg(j,i), FNU=order(i), KODE=2, N=1, CY=K(i), NZ=nz, IERR=ierr)
                 if (ierr /= 0 .and. ierr /= 3) then
                    print '(A,I0,A,"(",ES10.3,",",ES10.3,")",A,ES10.3)', &
                         & 'ERROR: cbesk K_nu(z) (P) ierr:',&
                         & ierr,' z:',arg(j,i),' nu:',order(i)
                    stop 903
                 end if
              end do
              fp(j) = (barft(j) * betgam(j) * K(-1))/K(0) 
           end do
        else ! COMPUTEP
           ! solution for pressure given specified pressure (only interesting for rD>1)
           do j = 1,NP
              do i = -1,0
                 ! both numerator and denominator terms are order nu (not nu-1)
                 call cbesk(Z=arg(j,i), FNU=order(0), KODE=2, N=1, CY=K(i), NZ=nz, IERR=ierr)
                 if (ierr /= 0 .and. ierr /= 3) then
                    print '(A,I0,A,"(",ES10.3,",",ES10.3,")",A,ES10.3)', &
                         & 'ERROR: cbesk K_nu(z) (P) ierr:',&
                         & ierr,' z:',arg(j,i),' nu:',order(0)
                    stop 904
                 end if
              end do
              fp(j) = (barft(j) * K(-1))/K(0)
           end do
        end if
     end if

     fp(1:np) = (r_term * unscale(:)) * fp(:) ! apply numerator r_term and unscaling factor last
     
     ft = deHoog_invlap(t(n),tee,fp(:),lap)
     if (DERIV) then
        fp = fp*p
        dft = deHoog_invlap(t(n),tee,fp(:),lap)*t(n)
        write(30,'(3(ES20.12E3,2X))') t(n),ft,dft
     else
        write(30,'(2(ES20.12E3,2X))') t(n),ft
     end if

  end do

  deallocate(t,p,fp,beta,barft,timePar,unscale)
  if (allocated(ivsq)) then
     deallocate(ivsq)
  end if

contains

  function time_pvect(p,par,flag) result(mult)
    use constants, only : DP
    use utility, only : operator(.X.)
    implicit none

    complex(DP), dimension(:), intent(in) :: p
    integer, intent(in) :: flag
    real(DP), dimension(:), intent(in) :: par
    complex(DP), dimension(size(p,1)) :: mult

    real(DP), allocatable :: ti(:), y(:), dy(:), W(:), b(:), deltaW(:)
    real(DP), allocatable :: denom(:), numer(:)

    real(DP) :: tf, yf, deltaWf, bf
    integer :: n, NP
    NP = size(p,1)

    select case (flag)
    case (0)
       ! step on at t=0 (default)
       mult(1:NP) = 1.0_DP/p

    case (:-1)

       !! arbitrary piecewise constant time behavior
       !! with n steps, from ti(1) to tf
       n = -flag
       allocate(ti(n),y(0:n),dy(n))

       ! unpack initial times, pumping rates and final time
       ti(1:n) = par(1:n)
       tf = par(n+1)
       y(0) = 0.0_DP
       y(1:n) = par(n+2:2*n+1)
       dy = y(1:n) - y(0:n-1)
       yf = sum(dy(:))

       mult(1:NP) = (sum(spread(dy,2,NP)*exp(-ti .X. p),1)-yf*exp(-tf*p))/p

       deallocate(ti,y,dy)

    case (1:)
       !! piecewise linear pumping rate with n steps, from ti(1) to tf
       !! no jumps in value (no vertical slopes)
       n = flag
       allocate(ti(n),W(0:n+1),y(1:n+1),denom(n),numer(n),&
            & b(n),deltaW(n))

       ! unpack initial times, pumping rates and final time
       ti(1:n) = par(1:n)
       tf = par(n+1)
       y(1:n) = par(n+2:2*n+1)

       ! compute slope between each pair of points
       W(0) = 0.0_DP  ! 0 and n+1 are ghost slopes for FD calc
       W(n+1) = 0.0_DP
       y(n+1) = 0.0_DP  ! assumed specified BC at tf
       denom = [ti(2:n),tf] - ti(1:n) ! run
       numer = y(2:n+1) - y(1:n) ! rise
       where (abs(denom) < epsilon(abs(numer)))
          ! ensure no divide by zero errors
          denom = denom + epsilon(abs(numer))
       end where
       W(1:n) = numer/denom ! slope = rise/run

       deltaW(1:n) = W(1:n) - W(0:n-1)
       deltaWf = sum(deltaW)

       ! intercept of lines shifted down to start on x-axis
       b(1:n) = -deltaW(1:n)*ti(1:n)
       bf = -deltaWf*tf

       mult(1:NP) = sum((spread(deltaW,2,NP)/spread(p**2,1,n) + &
            & spread(b,2,NP))*exp(-ti .X. p),dim=1) - &
            & (deltaWf/p**2 + bf)*exp(-tf*p)

       deallocate(ti,W,y,denom,deltaW,b)

    end select
  end function time_pvect

end program powerlaw
