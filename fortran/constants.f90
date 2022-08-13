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

module constants

  ! real with range 300 orders of mag, 15 sig figs (8 on both g95 & ifort)
  integer, parameter :: DP = selected_real_kind(p=15,r=300)

  !! extended range internal variables (10 on g95, 10 on gfortran, 16 on ifort)
  !! on x86 and x86_64 this is implemented in hardware
  !!integer, parameter :: EP = selected_real_kind(r=3000)

  !! full quad precision (only on gfortran >= 4.6 and ifort)
  !! implemented in software and is ~100x slower
  integer, parameter :: EP = selected_real_kind(p=33,r=3000)

  !! 3.141592653589793238462643383279503_EP
  real(DP), parameter :: PI =    4.0_DP*atan(1.0_DP)
  real(DP), parameter :: PISQ = (4.0_DP*atan(1.0_DP))**2
  real(DP), parameter :: TWOPI = 8.0_DP*atan(1.0_DP)
  real(DP), parameter :: RT2PI = sqrt(8.0_DP*atan(1.0_DP)) 
  real(DP), parameter :: RT2 = sqrt(2.0_DP)
  real(DP), parameter :: PIOV2 =  2.0_DP*atan(1.0_DP)
  real(DP), parameter :: PIOV4 = atan(1.0_DP)
  real(DP), parameter :: PISQOV2 = (4.0_DP*atan(1.0_DP))**2/2.0_DP
  complex(DP), parameter :: EYE = cmplx(0.0_DP,1.0_DP,DP)
  real(DP), parameter :: ULP = 2.0_DP
  
end module constants

