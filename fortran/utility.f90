!
! Copyright (c) 2014-2022 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
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

module utility
  implicit none
  private
  public :: logspace, linspace, dbessel_j0, ddbessel_j0, ctanh, &
       & dbessel_y0, ddbessel_y0, solve_tridiag, operator(.X.), operator(.SUM.)
  
  interface linspace
     module procedure linspacedp,linspaceep
  end interface linspace

  interface logspace
     module procedure logspace_int, logspace_real
  end interface logspace

  interface operator(.X.)
     module procedure outerprod_zd, outerprod_dz, outerprod_dd, outerprod_zz, &
          & outerprod_dpzd, outerprod_dpdz, outerprod_dpdd, outerprod_dpzz
  end interface operator(.X.)
  
  interface operator(.SUM.)
     module procedure outersum_zz, outersum_zd, outersum_dz, &
          & outersum_dpzz, outersum_dpzd, outersum_dpdz
  end interface operator(.SUM.)
  
contains
  
  function linspacedp(lo,hi,num) result(v)
    use constants, only : DP
    real(DP), intent(in) :: lo,hi
    integer, intent(in) :: num
    real(DP), dimension(num) :: v
    integer :: i
    real(DP) :: dx

    if (num == 1) then
       ! avoid div-by-zero in dx
       v = [lo]
    else
       dx = (hi-lo)/(num-1)
       do i=1,num
          v(i) = lo + (i-1)*dx
       end do
    end if
  end function linspacedp

  function linspaceep(lo,hi,num) result(v)
    use constants, only : EP
    real(EP), intent(in) :: lo,hi
    integer, intent(in) :: num
    real(EP), dimension(num) :: v
    integer :: i
    real(EP) :: dx
    
    if (num == 1) then
       v = [lo]
    else
       dx = (hi-lo)/(num-1)
       do i=1,num 
          v(i) = lo + (i-1)*dx
       end do
    end if
  end function linspaceep

  function logspace_int(lo,hi,num) result(v)
    use constants, only : DP
    integer, intent(in) :: lo,hi,num
    real(DP), dimension(num) :: v
    v = 10.0_DP**linspace(real(lo,DP),real(hi,DP),num)
  end function logspace_int

  function logspace_real(lo,hi,num) result(v)
    use constants, only : DP
    real(DP), intent(in) :: lo,hi
    integer, intent(in) :: num
    real(DP), dimension(num) :: v
    v = 10.0_DP**linspace(lo,hi,num)
  end function logspace_real

  elemental function dbessel_j0(x,a) result(df)
    use constants, only : EP
    real(EP), intent(in) :: x,a
    real(EP) :: df
    intrinsic :: bessel_j1

    df = -a*bessel_j1(x*a)

  end function dbessel_j0

  elemental function ddbessel_j0(x,a) result(ddf)
    use constants, only : EP
    real(EP), intent(in) :: x,a
    real(EP) :: ddf,ax
    intrinsic :: bessel_j0, bessel_j1

    if (abs(x) > epsilon(0.0)) then
       ax = a*x
       ! http://dlmf.nist.gov/10.6#E2
       ddf = -a**2*(bessel_j0(ax) - bessel_j1(ax)/ax)
    else
       ddf = -0.5_EP ! j1 has 1/2 slope at origin
    end if
    
  end function ddbessel_j0

  elemental function dbessel_y0(x,a) result(df)
    use constants, only : EP
    real(EP), intent(in) :: x,a
    real(EP) :: df
    intrinsic :: bessel_y1

    df = -a*bessel_y1(x*a)

  end function dbessel_y0

  elemental function ddbessel_y0(x,a) result(ddf)
    use constants, only : EP
    real(EP), intent(in) :: x,a
    real(EP) :: ddf,ax
    intrinsic :: bessel_y0, bessel_y1

    if (abs(x) > epsilon(0.0)) then
       ax = a*x
       ddf = -a**2*(bessel_y0(ax) - bessel_y1(ax)/ax)
    else
       ddf = -(huge(1.0_EP) + 1.0) ! y1 ->  negative infinity
    end if
    
  end function ddbessel_y0

  pure subroutine solve_tridiag(a,b,c,v,x)
    use constants, only : EP
    implicit none
    ! modified from en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    !      a - sub-diagonal (below main diagonal)
    !      b - the main diagonal
    !      c - sup-diagonal (above main diagonal)
    !      v - right part
    !      x - the answer

    complex(EP),dimension(:,:,:),intent(in) :: a,b,c,v
    complex(EP),dimension(size(a,1),size(a,2),size(a,3)),intent(out) :: x
    complex(EP),dimension(size(a,1),size(a,2),size(a,3)) :: cp,vp
    complex(EP),dimension(size(a,2),size(a,3)) :: m
    integer :: i,np,nj,n

    n = size(a,1)
    nj = size(a,2)
    np = size(a,3)

    ! initialize c-prime and v-prime
    cp(1,:,:) = c(1,:,:)/b(1,:,:)
    vp(1,:,:) = v(1,:,:)/b(1,:,:)

    !solve for vectors c-prime and v-prime
    do i = 2,n
       m = b(i,:,:) - cp(i-1,:,:)*a(i,:,:) 
       cp(i,:,:) = c(i,:,:)/m
       vp(i,:,:) = (v(i,:,:) - vp(i-1,:,:)*a(i,:,:))/m
    end do

    !initialize x
    x(n,:,:) = vp(n,:,:)

    !solve for x from vectors c-prime and d-prime
    do i = n-1, 1, -1
       x(i,:,:) = vp(i,:,:) - c(i,:,:)*x(i+1,:,:)
    end do

  end subroutine solve_tridiag

  pure function outerprod_dpdd(da,db) result(c)
    use constants, only : DP
    real(DP), intent(in), dimension(:) :: da,db
    real(DP), dimension(size(da),size(db)) :: c
    c = spread(da,dim=2,ncopies=size(db))*spread(db,dim=1,ncopies=size(da))
  end function outerprod_dpdd

  pure function outerprod_dpzd(za,db) result(c)
    use constants, only : DP
    complex(DP), intent(in), dimension(:) :: za
    real(DP), intent(in), dimension(:) :: db
    complex(DP), dimension(size(za),size(db)) :: c
    c = spread(za,dim=2,ncopies=size(db))*spread(db,dim=1,ncopies=size(za))
  end function outerprod_dpzd

  pure function outerprod_dpdz(da,zb) result(c)
    use constants, only : DP
    real(DP), intent(in), dimension(:) :: da
    complex(DP), intent(in), dimension(:) :: zb
    complex(DP), dimension(size(da),size(zb)) :: c
    c = spread(da,dim=2,ncopies=size(zb))*spread(zb,dim=1,ncopies=size(da))
  end function outerprod_dpdz

  pure function outerprod_dpzz(za,zb) result(c)
    use constants, only : DP
    complex(DP), intent(in), dimension(:) :: za,zb
    complex(DP), dimension(size(za),size(zb)) :: c
    c = spread(za,dim=2,ncopies=size(zb))*spread(zb,dim=1,ncopies=size(za))
  end function outerprod_dpzz

  !===

  pure function outerprod_dd(da,db) result(c)
    use constants, only : EP
    real(EP), intent(in), dimension(:) :: da,db
    real(EP), dimension(size(da),size(db)) :: c
    c = spread(da,dim=2,ncopies=size(db))*spread(db,dim=1,ncopies=size(da))
  end function outerprod_dd

  pure function outerprod_zd(za,db) result(c)
    use constants, only : EP
    complex(EP), intent(in), dimension(:) :: za
    real(EP), intent(in), dimension(:) :: db
    complex(EP), dimension(size(za),size(db)) :: c
    c = spread(za,dim=2,ncopies=size(db))*spread(db,dim=1,ncopies=size(za))
  end function outerprod_zd

  pure function outerprod_dz(da,zb) result(c)
    use constants, only : EP
    real(EP), intent(in), dimension(:) :: da
    complex(EP), intent(in), dimension(:) :: zb
    complex(EP), dimension(size(da),size(zb)) :: c
    c = spread(da,dim=2,ncopies=size(zb))*spread(zb,dim=1,ncopies=size(da))
  end function outerprod_dz

  pure function outerprod_zz(za,zb) result(c)
    use constants, only : EP
    complex(EP), intent(in), dimension(:) :: za,zb
    complex(EP), dimension(size(za),size(zb)) :: c
    c = spread(za,dim=2,ncopies=size(zb))*spread(zb,dim=1,ncopies=size(za))
  end function outerprod_zz

  !===

  pure function outersum_zz(za,zb) result(c)
    use constants, only : EP
    complex(EP), intent(in), dimension(:) :: za,zb
    complex(EP), dimension(size(za),size(zb)) :: c
    c = spread(za,dim=2,ncopies=size(zb))+spread(zb,dim=1,ncopies=size(za))
  end function outersum_zz

  pure function outersum_zd(za,db) result(c)
    use constants, only : EP
    complex(EP), intent(in), dimension(:) :: za
    real(EP), intent(in), dimension(:) :: db
    complex(EP), dimension(size(za),size(db)) :: c
    c = spread(za,dim=2,ncopies=size(db))+spread(db,dim=1,ncopies=size(za))
  end function outersum_zd

  pure function outersum_dz(da,zb) result(c)
    use constants, only : EP
    real(EP), intent(in), dimension(:) :: da
    complex(EP), intent(in), dimension(:) :: zb
    complex(EP), dimension(size(da),size(zb)) :: c
    c = spread(da,dim=2,ncopies=size(zb))+spread(zb,dim=1,ncopies=size(da))
  end function outersum_dz

  !===

  pure function outersum_dpzz(za,zb) result(c)
    use constants, only : DP
    complex(DP), intent(in), dimension(:) :: za,zb
    complex(DP), dimension(size(za),size(zb)) :: c
    c = spread(za,dim=2,ncopies=size(zb))+spread(zb,dim=1,ncopies=size(za))
  end function outersum_dpzz

  pure function outersum_dpzd(za,db) result(c)
    use constants, only : DP
    complex(DP), intent(in), dimension(:) :: za
    real(DP), intent(in), dimension(:) :: db
    complex(DP), dimension(size(za),size(db)) :: c
    c = spread(za,dim=2,ncopies=size(db))+spread(db,dim=1,ncopies=size(za))
  end function outersum_dpzd

  pure function outersum_dpdz(da,zb) result(c)
    use constants, only : DP
    real(DP), intent(in), dimension(:) :: da
    complex(DP), intent(in), dimension(:) :: zb
    complex(DP), dimension(size(da),size(zb)) :: c
    c = spread(da,dim=2,ncopies=size(zb))+spread(zb,dim=1,ncopies=size(da))
  end function outersum_dpdz

  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  elemental function ctanh(z) result(f)
    use constants, only : DP
    complex(DP), intent(in) :: z
    complex(DP) :: f
    real(DP) :: x,y
    x = real(z)
    y = aimag(z)
    f = cmplx(tanh(2*x)/(1+cos(2*y)/cosh(2*x)), &
         & sin(2*y)/(cosh(2*x)+cos(2*y)),DP)
  end function ctanh

end module utility


