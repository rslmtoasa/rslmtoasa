!> LMTO-ASA Structure Constants Library - Bessel Function Module
!> ================================================================
!>
!> This module provides spherical Bessel and Hankel functions used in
!> multiple scattering theory and LMTO calculations.
!>
!> Physical Significance:
!> - Bessel functions j_l(kr) are solutions of the free-particle Schrödinger equation
!> - Hankel functions h_l(kr) = j_l(kr) + i*n_l(kr) represent outgoing spherical waves
!> - Essential for describing electron wave functions in atomic sphere approximation
!> - Used in KKR Green's functions and LMTO basis construction
!>
!> Key Conventions (controlled by loka parameter):
!> - loka=0: Methfessel conventions (MSM) - real functions for bound states
!> - loka=1: Andersen conventions (OKA) - LMTO-specific normalization
!> - loka=2: Modified Andersen - alternative scaling for NMTO
!>
!> Author: Retained local Bessel/Hankel kernel implementation
!> Date: 2025-01-26

module lmto_strux_bessel
    use strux_errors, only: rx, rxi
    implicit none
    private
    
    ! Public interfaces
    public :: besslr      ! General Bessel/Hankel function calculator
    public :: bessl2      ! Andersen conventions wrapper
    public :: besslm      ! Methfessel conventions wrapper  
    public :: bessl       ! Methfessel conventions (alias)
    public :: bessjs      ! High-precision Bessel using dbesnu
    
contains

subroutine besslr(y, loka, lmin, lmax, fi, gi)
!> Purpose:
!>   Compute radial Bessel and Hankel functions for multiple-scattering theory.
!> Notes:
!>   The detailed legacy input/output remarks below are preserved because this
!>   routine defines the core normalization conventions used across the library.
! ----------------------------------------------------------------------
!i Inputs:
!i   y       :argument y = (κr)² where κ is wave number (complex in general)
!i           :y > 0: propagating waves (real κ), y < 0: evanescent waves (imaginary κ)
!i   loka    :convention selector for function definitions
!i           :1s digit: 0=Methfessel, 1=Andersen, 2=modified Andersen
!i           :10s digit: numerical method (0=besslr, 1=besnu for large y, 2=dbesnu)
!i   lmin    :minimum angular momentum quantum number
!i   lmax    :maximum angular momentum quantum number
!o Outputs
!o   fi(l)   :radial function proportional to spherical Bessel j_l(κr) / r^l
!o           :For y ≤ 0: fi(0) = sinh(x)/x where x = √(-y) = κr
!o           :Power series expansion used for numerical evaluation
!o   gi(l)   :radial function for outgoing/incoming waves
!o           :For y > 0: proportional to spherical Neumann n_l(κr) × r^(l+1)
!o           :For y < 0: proportional to spherical Hankel h_l(κr) × r^(l+1)
!o           :gi(0) = cos(√y) for y > 0, exp(-x) for y < 0
!r Remarks
!r   This subroutine computes the radial parts of spherical wave functions
!r   that are fundamental to quantum scattering theory and electronic
!r   structure calculations. The functions satisfy the radial Schrödinger
!r   equation for free particles or the Helmholtz equation.
!r
!r   **Mathematical Foundation:**
!r   Spherical Bessel functions j_l(z) and spherical Hankel functions h_l(z)
!r   are solutions to the radial equation:
!r   d²u/dr² + [k² - l(l+1)/r²]u = 0
!r
!r   **Physical Significance in Electronic Structure:**
!r   - **Free particle wave functions**: Solutions for electrons in constant potential
!r   - **Scattering states**: Asymptotic forms for electron scattering
!r   - **Multiple scattering theory**: KKR and LMTO basis functions
!r   - **Green's functions**: Resolvent operators for perturbation theory
!r   - **Partial wave expansion**: Angular momentum decomposition
!r
!r   **Multiple Scattering Applications:**
!r   - Wave functions: ψ(r) = ∑_l [A_l j_l(kr) + B_l h_l(kr)] Y_lm(θ,φ)
!r   - Structure constants: C_L(L') = ∑_R h_l(k|R|) e^{ik·R} Y_L^*(R) Y_L'(R)
!r   - T-matrices: Linear relations between regular/irregular solution coefficients
!r   - Scattering path operator: τ = m - m τ m where m is single-site t-matrix
!r
!r   **Key Relations:**
!r   - Wronskians: {h,j} = 1/(iκ), {n,j} = -1/κ (normalization conditions)
!r   - Asymptotic forms: h_l(z) ∼ e^{iz}/z, j_l(z) ∼ sin(z)/z, n_l(z) ∼ -cos(z)/z
!r   - Small z limits: j_l(z) ∼ z^l/(2l+1)!!, h_l(z) ∼ -i(2l-1)!!/z^{l+1}
!r   - Parity: j_l(-z) = (-1)^l j_l(z), h_l(-z) = (-1)^{l+1} h_l(z)
!u Updates
!u   03 Jan 18 Added loka=2 convention
!u   22 Aug 17 Optionally calls dbesnu for high accuracy.  Especially important for y>100.
!u   23 Jul 08 bug fix: besslr doesn't make fi/gi(lmax+1) when lmax=0
!u   19 May 04 Changed loka from logical to integer
! ----------------------------------------------------------------------
    implicit none
    ! Passed variables:
    integer, intent(in) :: loka, lmin, lmax
    double precision, intent(in) :: y
    double precision, intent(out) :: fi(lmin:lmax), gi(lmin:lmax)
    
    ! Local variables:
    integer :: i, isn, j1, j2, k, l, lmx, lmxp1, lmxp2, nf, tlp1, ll1, ll2, opt1, lmin1
    integer, parameter :: nlmax = 20
    double precision :: dt, dt2, exppr, my, srmy, g1, t
    double precision :: dum(nlmax*4 + 2), fac2l(-nlmax:nlmax*2 + 3)
    double precision :: fil(0:1), gil(0:1)
    real(8), parameter :: tol = 1.d-15
    logical, parameter :: lhank = .true.

    if (lmin > 0) call rx('BESSLR : lmin gt 0')
    if (lmax < 0) call rx('BESSLR : lmax lt 0')

    opt1 = mod(loka/10, 10)
    lmx = max(lmax, 2)
    if (lmx > nlmax + nlmax) call rxi(' BESSL : lmax gt nlmax*2, lmax=', lmx)

    ! A table of fac2l(l)=(2l-1)!!
    fac2l(0) = 1d0
    do l = 1, lmx + 1
        fac2l(l) = fac2l(l - 1)*(l + l - 1)
    end do
    do l = -1, lmin, -1
        fac2l(l) = fac2l(l + 1)/(l + l + 1)
    end do

    ! Case call bessjs
    if (y >= 90 .and. opt1 == 1 .or. opt1 == 2) then
        call bessjs(y, lmax, 13, fi(0), gi(0))
        if (lmin >= 0) goto 100
        lmin1 = -1
        if (lmax < 1) then
            call bessjs(y, 1, 13, fil, gil)
            fi(-1) = fil(0) - y*fil(1)
            gi(-1) = (gil(0) - gil(1))/y
            lmin1 = -2
        end if
        ! Warning ! if lmin is very negative and y large, recursion can be unstable
        do l = lmin1, lmin, -1
            dt = 2*l + 3
            fi(l) = dt*fi(l + 1) - y*fi(l + 2)
            gi(l) = (dt*gi(l + 1) - gi(l + 2))/y
        end do
        goto 100
    end if

    ! Case akap=0
    if (y == 0) then
        do l = lmin, lmax
            fi(l) = 1/fac2l(l + 1)
            gi(l) = fac2l(l)
        end do
        goto 100
    end if
    my = -y

    ! Get dum(1) = j_{lmx}(x)/x^{lmx} = fi(lmx)
    tlp1 = lmx + lmx + 1
    dt = 1d0
    t = 1d0
    i = 0
    do k = 1, 1000
        if (abs(dt) < tol) goto 21
        i = i + 2
        dt2 = i + tlp1
        dt = dt*my/(i*dt2)
        t = t + dt
    end do
    call rx('BESSLR: series not convergent')
21  continue
    dum(1) = t/fac2l(lmx + 1)

    ! Get dum(2) = j_{lmx-1}(x)/x^{lmx-1} = fi(lmx-1)
    tlp1 = tlp1 - 2
    dt = 1d0
    t = 1d0
    i = 0
    do k = 1, 1000
        if (abs(dt) < tol) goto 31
        i = i + 2
        dt2 = i + tlp1
        dt = dt*my/(i*dt2)
        t = t + dt
    end do
    call rx('BESSLR: series not convergent')
31  continue
    dum(2) = t/fac2l(lmx)

    ! Recursion for dum(k)=j_{lmx+1-k}(x)/x^{lmx+1-k}=fi(lmx+1-k)
    ll1 = lmx + lmx + 1
    ll2 = ll1 + 1
    nf = ll1
    do k = 3, ll2
        nf = nf - 2
        dum(k) = nf*dum(k - 1) - y*dum(k - 2)
    end do

    ! Get fi and gi from dum
    lmxp1 = lmx + 1
    lmxp2 = lmx + 2
    isn = (-1)**lmin
    do k = lmin, lmax
        j1 = lmxp1 - k
        j2 = lmxp2 + k
        fi(k) = dum(j1)
        ! n_l(x) = j_{-l-1}*(-1)^{l+1}
        gi(k) = dum(j2)*isn
        isn = -isn
    end do

    ! For E<0, use Hankel functions rather than Neumann functions
    if (lhank .and. y < 0d0) then
        srmy = sqrt(-y)
        gi(0) = 1d0
        g1 = 1d0 + srmy
        if (lmax >= 1) gi(1) = g1
        if (lmax >= 2) then
            tlp1 = 1
            do l = 2, lmax
                tlp1 = tlp1 + 2
                gi(l) = tlp1*gi(l - 1) - y*gi(l - 2)
            end do
        end if
        if (lmin <= -1) then
            gi(-1) = (gi(0) - g1)/y
            if (lmin <= -2) then
                do l = -2, lmin, -1
                    gi(l) = ((l + l + 3)*gi(l + 1) - gi(l + 2))/y
                end do
            end if
        end if
        exppr = 1d0/exp(srmy)
        forall (l=lmin:lmax) gi(l) = gi(l)*exppr
    end if

100 continue
    ! Scaling to Andersen's 2nd generation LMTO conventions
    if (mod(loka, 10) == 1) then
        do l = lmin, lmax
            fi(l) = fi(l)*fac2l(l)*0.5d0
            gi(l) = gi(l)/fac2l(l)
        end do
    end if

    ! Scaling to conventions mode 2
    if (mod(loka, 10) == 2) then
        do l = lmin, lmax
            fi(l) = fi(l)*fac2l(l)*(2*l + 1)
            gi(l) = gi(l)/fac2l(l)
        end do
    end if

end subroutine besslr


subroutine bessl2(y, lmin, lmax, fi, gi)
!> Purpose:
!>   Wrapper around `besslr` for Andersen-style normalization conventions.
!> Notes:
!>   See `besslr` for the full definitions. For `y=0`, `fi(l)=0.5/(2l+1)` and
!>   `gi(l)=1`.
    implicit none
    integer, intent(in) :: lmin, lmax
    double precision, intent(in) :: y
    double precision, intent(out) :: fi(lmin:lmax), gi(lmin:lmax)

    call besslr(y, 1, lmin, lmax, fi, gi)
end subroutine bessl2


subroutine besslm(y, lmax, fi, gi)
!> Purpose:
!>   Wrapper around `besslr` for Methfessel normalization conventions.
!> Notes:
!>   See `besslr` for the full definitions. For `y=0`, `fi(l)=(2l-1)!!` and
!>   `gi(l)=1/(2l+1)!!`.
    implicit none
    integer, intent(in) :: lmax
    double precision, intent(in) :: y
    double precision, intent(out) :: fi(0:lmax), gi(0:lmax)

    call besslr(y, 0, 0, lmax, fi, gi)
end subroutine besslm


subroutine bessl(y, lmax, fi, gi)
!> Purpose:
!>   Alias wrapper for the Methfessel-style `besslr` conventions.
!> Notes:
!>   See `besslr` for the full definitions. For `y=0`, `fi(l)=(2l-1)!!` and
!>   `gi(l)=1/(2l+1)!!`.
    implicit none
    integer, intent(in) :: lmax
    double precision, intent(in) :: y
    double precision, intent(out) :: fi(0:lmax), gi(0:lmax)

    call besslr(y, 0, 0, lmax, fi, gi)
end subroutine bessl


subroutine bessjs(y, lmax, opt, jr, hr)
!> Purpose:
!>   Return spherical Bessel and/or Hankel functions using the high-accuracy path.
!> Notes:
!>   The detailed remarks below explain the option flags and the crossover
!>   between the asymptotic `dbesnu` path and the lower-`y` `besslr` path.
! ----------------------------------------------------------------------
!i Inputs
!i   y     :y = x**2, where x is argument to Bessel
!i         :y > 0 => return spherical Bessel/Neumann function
!i         :y < 0 => return modified spherical Bessel/Hankel function
!i   lmax  :maximum l for a given site
!i   opt   :1s digit
!i         : 0 return nothing in jr
!i         : 1 return spherical Bessel in jr if y>0
!i         :   return spherical Bessel of Im x in jr if y<0
!i         : 2 return spherical Hankel in hr if y>0
!i         :   return spherical Hankel of Im x in hr if y<0
!i         : combination of 1+2 is allowed
!i         :10s digit
!i         : 0 return spherical Bessel function in jr and/or spherical Hankel in hr
!i         :   depending on 1s digit opt
!i         : 1 scale jr by 1/x**l and/or hr by x**(l+1)
!i         :100s digit
!i         : 0 always call dbesnu
!i         : 1 call besslr for y<90
!i         : 2 always call besslr
!o Outputs
!i   jr    : spherical Bessel and/or functions, possibly scaled by r^l
!r Remarks
!r   This routine and besslr generate the same results to a relative
!r   precision of better 10^-14 for y<90.  For y>90 besslr becomes unstable
!u Updates
!u   21 Aug 17 First created
! ----------------------------------------------------------------------
    implicit none
    ! Passed parameters
    integer, intent(in) :: lmax, opt
    double precision, intent(in) :: y
    double precision, intent(out) :: jr(0:lmax), hr(0:lmax)
    
    ! Local parameters
    logical :: limg
    integer :: ierr, i, l, opt0, opt1, opt2, k2
    integer, parameter :: ltop = 50
    real(8), parameter :: pi = 4*atan(1d0), tol = 1d-15
    double precision :: x, xx(1), dfac, fi(0:ltop), gi(0:ltop)
    double precision :: twolpk, tlp1, my

    opt0 = mod(opt, 10)
    opt1 = mod(opt/10, 10)
    opt2 = mod(opt/100, 10)

    if (y == 0) then
        if (mod(opt0, 2) == 1) then
            jr(0) = 1
            dfac = 1
            do l = 0, lmax
                dfac = dfac*(2*l + 1)
                jr(l) = 1/dfac
            end do
            if (opt1 == 0) jr(1:lmax) = 0
        end if
        if (opt0 >= 2) then
            hr(0) = 1
            dfac = 1
            do l = 0, lmax
                hr(l) = dfac
                dfac = dfac*(2*l + 1)
            end do
            if (opt1 == 0) call rx('bessls: dbesnu invalid argument for hankel, x=0')
        end if
        return

    elseif (y < 90 .and. opt2 == 1 .or. opt2 == 2) then
        if (opt0 > 2) then
            call besslr(y, 0, 0, lmax, jr, hr)
        else
            if (lmax > ltop) call rx('bessjs: increase ltop')
            call besslr(y, 0, 0, lmax, fi, gi)
            if (mod(opt0, 2) == 1) jr(0:lmax) = fi(0:lmax)
            if (opt0 >= 2) hr(0:lmax) = gi(0:lmax)
        end if
        if (opt1 == 0) then
            call rx('bessls: add scaling to branch calling besslr')
        end if
        return
    end if

    limg = y < 0
    x = sqrt(abs(y))
    l = 11; if (limg) l = 13
    ierr = 0
    
    if (mod(opt0, 2) == 1) then
        call dbesnu(x, lmax + 0.5d0, l, jr, xx, ierr)
        ! dbesnu has trouble for small x and high l.  Use polynomial expansion
        if (ierr > 0 .and. x < 1) then
            my = -y
            tlp1 = 2*lmax + 1
            ierr = 0
            jr(ierr:lmax) = 1
            fi(0:lmax + 1 - ierr) = 1
            do k2 = 2, 1000, 2
                do i = 0, lmax - ierr
                    twolpk = k2 + tlp1
                    fi(i) = fi(i)*my/(k2*(twolpk - 2*i))
                    jr(lmax - i) = jr(lmax - i) + fi(i)
                end do
                if (abs(fi(lmax)) < tol) exit
            end do
            dfac = 1
            do i = 0, lmax
                dfac = dfac*(2*i + 1)
                jr(i) = jr(i)/dfac*x**i
            end do
        else
            jr(0:lmax) = sqrt(pi/2/x)*jr(0:lmax)
        end if
        if (ierr /= 0) call rx('bessls: dbesnu returned with error')
    end if
    
    l = l + 1
    if (opt0 >= 2) then
        call dbesnu(x, lmax + 0.5d0, l, hr, xx, ierr)
        if (ierr /= 0) call rx('bessls: dbesnu returned with error')
        dfac = -sqrt(pi/2/x); if (limg) dfac = sqrt(2/pi/x)
        hr(0:lmax) = dfac*hr(0:lmax)
    end if

    if (opt1 /= 0) then
        if (mod(opt0, 2) == 1) then
            forall (l=1:lmax) jr(l) = jr(l)/x**l
        end if
        if (opt0 >= 2) then
            forall (l=0:lmax) hr(l) = hr(l)*x**(l + 1)
        end if
    end if

end subroutine bessjs


end module lmto_strux_bessel

! Plain external subroutines — no module mangling so all callers can resolve them.

subroutine dbesnu(x, nu, opt, fr, fi, ierr)
!> Purpose:
!>   Evaluate Bessel functions of fractional order for a real argument.
!> Notes:
!>   Converted from Questaal `provided/dbesnu.f` to F90 and kept close to the
!>   original because it underpins the large-argument branch in `bessjs`.
    use strux_errors, only: rx
!  Inputs:
!    x   : evaluate Bessel at x
!    nu  : Bessel order
!    opt : 1s digit: 1=J, 2=Y/Neumann, 3=I (modified 1st kind), 4=K (modified 2nd kind)
!          10s digit: 0=return only order nu, 1=return all orders 0..int(nu)
!          100s digit: 0=function, 1=1st deriv, 2=2nd deriv
!  Outputs:
!    fr  : function values (real part)
!    fi  : imaginary part (unused, reserved)
!    ierr: 0=success, <0=error, >0=only orders 0..ierr-1 accurate
    implicit none
    double precision, intent(in) :: x, nu
    integer, intent(in) :: opt
    double precision, intent(out) :: fr(*), fi(*)
    integer, intent(out) :: ierr

    integer :: nb, opt0, opt1, opt2, nextra, i
    double precision :: alpha, b(1000), fac, fac2

    nb = int(nu) + 1
    alpha = nu + 1 - nb
    if (alpha >= 1d0) then
        alpha = alpha - 1d0
        nb = nb + 1
    end if
    if (nb > 1000) call rx('dbesnu: increase size of b')
    ierr = 0
    opt0 = mod(opt, 10)
    opt1 = mod(opt/10, 10)
    opt2 = mod(opt/100, 10)

    nextra = 0
    if (opt2 /= 0) nextra = 1

    if (opt0 == 1) then
        call djbesl(x, alpha, nb+nextra, b, ierr)
        fac  = -1d0
        fac2 =  1d0
    else if (opt0 == 2) then
        call dybesl(x, alpha, nb+nextra, b, ierr)
        fac  = -1d0
        fac2 =  1d0
    else if (opt0 == 3) then
        call dibesl(x, alpha, nb+nextra, 1, b, ierr)
        fac  =  1d0
        fac2 = -1d0
    else if (opt0 == 4) then
        call dkbesl(x, alpha, nb+nextra, 1, b, ierr)
        fac  = -1d0
        fac2 = -1d0
    else
        call rx('dbesnu: bad opt')
        return
    end if

    if (ierr < 0) return
    if (ierr >= nb+nextra) ierr = 0

    if (opt1 > 0) then
        if (opt2 == 0) then
            fr(1:nb) = b(1:nb)
        else
            do i = 1, nb
                fr(i) = b(i)*(dble(i-1)+alpha)/x + fac*b(i+1)
            end do
            if (opt2 == 2) then
                do i = 1, nb
                    fr(i) = b(i)*((nu/x)**2 - fac2) - fr(i)/x
                end do
            end if
        end if
    else
        if (opt2 == 0) then
            fr(1) = b(nb)
        else
            fr(1) = b(nb)*nu/x + fac*b(nb+1)
            if (opt2 == 2) then
                fr(1) = b(nb)*((nu/x)**2 - fac2) - fr(1)/x
            end if
        end if
    end if

end subroutine dbesnu


subroutine djbesl(x, alpha, nb, b, ncalc)
!> Purpose:
!>   Evaluate `J_{n+alpha}(x)` for `n=0..nb-1`, with `alpha` in `[0,1)`.
!> Notes:
!>   Converted from Cody/Argonne TOMS F77 `djbesl.f` to F90.
    implicit none
    double precision, intent(in)  :: x, alpha
    integer,          intent(in)  :: nb
    double precision, intent(out) :: b(nb)
    integer,          intent(out) :: ncalc

    integer :: i, j, k, l, m, magx, n, nbmx, nend, nstart
    double precision :: alpem, alp2em, capp, capq, em, en, gnu
    double precision :: halfx, one, p, pi2, plast, pold, psave, psavel
    double precision :: s, sum, t, t1, tempa, tempb, tempc, test
    double precision :: tover, two, twofiv, twopi1, twopi2, vcos, vsin
    double precision :: xc, xin, xk, xm, z, zero, half, three, four
    double precision :: one30, three5, eighth, enmten, ensig, enten, rtnsig, xlarge
    double precision, parameter :: fact(25) = [ &
        1.0d0,1.0d0,2.0d0,6.0d0,24.0d0,1.2d2,7.2d2,5.04d3, &
        4.032d4,3.6288d5,3.6288d6,3.99168d7,4.790016d8,6.2270208d9, &
        8.71782912d10,1.307674368d12,2.0922789888d13,3.55687428096d14, &
        6.402373705728d15,1.21645100408832d17,2.43290200817664d18, &
        5.109094217170944d19,1.12400072777760768d21, &
        2.585201673888497664d22,6.2044840173323943936d23 ]

    pi2    = 0.636619772367581343075535d0
    twopi1 = 6.28125d0
    twopi2 = 1.935307179586476925286767d-3
    zero   = 0.0d0; eighth = 0.125d0; half = 0.5d0; one = 1.0d0
    two    = 2.0d0; three  = 3.0d0;   four = 4.0d0; twofiv = 25.0d0
    one30  = 130.0d0; three5 = 35.0d0
    enten  = 1.0d38; ensig  = 1.0d17; rtnsig = 1.0d-4
    enmten = 1.2d-37; xlarge = 1.0d4

    magx  = int(x)
    ncalc = nb
    b(1:nb) = zero

    if (.not. ((nb > 0) .and. (x >= zero) .and. (x <= xlarge) .and. &
               (alpha >= zero) .and. (alpha < one))) then
        b(1) = zero
        ncalc = min(nb,0) - 1
        return
    end if

    if (x < rtnsig) then
        ! Two-term ascending series for small X
        tempa = one
        alpem = one + alpha
        halfx = zero
        if (x > enmten) halfx = half*x
        if (alpha /= zero) tempa = halfx**alpha / (alpha * gamma(alpha))
        tempb = zero
        if ((x+one) > one) tempb = -halfx*halfx
        b(1) = tempa + tempa*tempb/alpem
        if ((x /= zero) .and. (b(1) == zero)) ncalc = 0
        if (nb /= 1) then
            if (x <= zero) then
                b(2:nb) = zero
            else
                tempc = halfx
                tover = (enmten+enmten)/x
                if (tempb /= zero) tover = enmten/abs(tempb)
                do n = 2, nb
                    tempa = tempa/alpem
                    alpem = alpem + one
                    tempa = tempa*tempc
                    if (tempa <= tover*alpem) tempa = zero
                    b(n) = tempa + tempa*tempb/alpem
                    if ((b(n) == zero) .and. (ncalc > n)) ncalc = n-1
                end do
            end if
        end if

    else if ((x > twofiv) .and. (nb <= magx+1)) then
        ! Asymptotic series for X > 25, NB small
        xc  = sqrt(pi2/x)
        xin = (eighth/x)**2
        m   = 11
        if (x >= three5) m = 8
        if (x >= one30) m = 4
        xm  = four*dble(m)
        t   = aint(x/(twopi1+twopi2)+half)
        z   = ((x-t*twopi1)-t*twopi2) - (alpha+half)/pi2
        vsin = sin(z)
        vcos = cos(z)
        gnu = alpha + alpha
        do i = 1, 2
            s    = ((xm-one)-gnu)*((xm-one)+gnu)*xin*half
            t    = (gnu-(xm-three))*(gnu+(xm-three))
            capp = s*t/fact(2*m+1)
            t1   = (gnu-(xm+one))*(gnu+(xm+one))
            capq = s*t1/fact(2*m+2)
            xk   = xm
            k    = m + m
            t1   = t
            do j = 2, m
                xk   = xk - four
                s    = ((xk-one)-gnu)*((xk-one)+gnu)
                t    = (gnu-(xk-three))*(gnu+(xk-three))
                capp = (capp+one/fact(k-1))*s*t*xin
                capq = (capq+one/fact(k))*s*t1*xin
                k    = k - 2
                t1   = t
            end do
            capp = capp + one
            capq = (capq+one)*(gnu*gnu-one)*(eighth/x)
            b(i) = xc*(capp*vcos - capq*vsin)
            if (nb == 1) return
            t    = vsin
            vsin = -vcos
            vcos = t
            gnu  = gnu + two
        end do
        if (nb > 2) then
            gnu = alpha + alpha + two
            do n = 3, nb
                b(n) = gnu*b(n-1)/x - b(n-2)
                gnu  = gnu + two
            end do
        end if

    else
        ! Miller's backward recurrence
        nbmx = nb - magx
        n    = magx + 1
        en   = dble(n+n) + (alpha+alpha)
        plast = one
        p     = en/x
        test  = ensig + ensig
        if (nbmx >= 3) then
            tover  = enten/ensig
            nstart = magx + 2
            nend   = nb - 1
            en     = dble(nstart+nstart) - two + (alpha+alpha)
            do k = nstart, nend
                n     = k
                en    = en + two
                pold  = plast
                plast = p
                p     = en*plast/x - pold
                if (p > tover) then
                    tover  = enten
                    p      = p/tover
                    plast  = plast/tover
                    psave  = p
                    psavel = plast
                    nstart = n + 1
                    do
                        n     = n + 1
                        en    = en + two
                        pold  = plast
                        plast = p
                        p     = en*plast/x - pold
                        if (p > one) exit
                    end do
                    tempb = en/x
                    test  = pold*plast*(half-half/(tempb*tempb))
                    test  = test/ensig
                    p     = plast*tover
                    n     = n - 1
                    en    = en - two
                    nend  = min(nb,n)
                    do l = nstart, nend
                        pold   = psavel
                        psavel = psave
                        psave  = en*psavel/x - pold
                        if (psave*psavel > test) then
                            ncalc = l - 1
                            goto 190
                        end if
                    end do
                    ncalc = nend
                    goto 190
                end if
            end do
            n    = nend
            en   = dble(n+n) + (alpha+alpha)
            test = max(test, sqrt(plast*ensig)*sqrt(p+p))
        end if
        do
            n     = n + 1
            en    = en + two
            pold  = plast
            plast = p
            p     = en*plast/x - pold
            if (p >= test) exit
        end do

190     continue
        n    = n + 1
        en   = en + two
        tempb = zero
        tempa = one/p
        m     = 2*n - 4*(n/2)
        sum   = zero
        em    = dble(n/2)
        alpem = (em-one) + alpha
        alp2em = (em+em) + alpha
        if (m /= 0) sum = tempa*alpem*alp2em/em
        nend = n - nb
        if (nend > 0) then
            do l = 1, nend
                n      = n - 1
                en     = en - two
                tempc  = tempb
                tempb  = tempa
                tempa  = (en*tempb)/x - tempc
                m      = 2 - m
                if (m /= 0) then
                    em     = em - one
                    alp2em = (em+em) + alpha
                    if (n == 1) goto 210
                    alpem  = (em-one) + alpha
                    if (alpem == zero) alpem = one
                    sum    = (sum+tempa*alp2em)*alpem/em
                end if
            end do
        end if
210     continue
        b(n) = tempa
        if (nend >= 0) then
            if (nb <= 1) then
                alp2em = alpha
                if ((alpha+one) == one) alp2em = one
                sum = sum + b(1)*alp2em
                goto 250
            else
                n  = n - 1
                en = en - two
                b(n) = (en*tempa)/x - tempb
                if (n == 1) goto 240
                m = 2 - m
                if (m /= 0) then
                    em     = em - one
                    alp2em = (em+em) + alpha
                    alpem  = (em-one) + alpha
                    if (alpem == zero) alpem = one
                    sum    = (sum+b(n)*alp2em)*alpem/em
                end if
            end if
        end if
        nend = n - 2
        if (nend /= 0) then
            do l = 1, nend
                n      = n - 1
                en     = en - two
                b(n)   = (en*b(n+1))/x - b(n+2)
                m      = 2 - m
                if (m /= 0) then
                    em     = em - one
                    alp2em = (em+em) + alpha
                    alpem  = (em-one) + alpha
                    if (alpem == zero) alpem = one
                    sum    = (sum+b(n)*alp2em)*alpem/em
                end if
            end do
        end if
        b(1) = two*(alpha+one)*b(2)/x - b(3)
240     continue
        em     = em - one
        alp2em = (em+em) + alpha
        if (alp2em == zero) alp2em = one
        sum = sum + b(1)*alp2em
250     continue
        if ((alpha+one) /= one) sum = sum*gamma(alpha)*(x*half)**(-alpha)
        tempa = enmten
        if (sum > one) tempa = tempa*sum
        do n = 1, nb
            if (abs(b(n)) < tempa) b(n) = zero
            b(n) = b(n)/sum
        end do
    end if

end subroutine djbesl


subroutine dybesl(x, alpha, nb, by, ncalc)
!> Purpose:
!>   Evaluate `Y_{n+alpha}(x)` for `n=0..nb-1`, with `alpha` in `[0,1)`.
!> Notes:
!>   Converted from Cody/Argonne TOMS F77 `dybesl.f` to F90.
    implicit none
    double precision, intent(in)  :: x, alpha
    integer,          intent(in)  :: nb
    double precision, intent(out) :: by(nb)
    integer,          intent(out) :: ncalc

    integer :: i, k, na
    double precision :: alfa, aye, c, cosmu, d, ddiv, div, dmu, d1, d2
    double precision :: e, en, enu, en1, even, ex, f, g, h, odd
    double precision :: p, pa, pa1, piby2, pim5, q, qa, qa1, q0, r, s, sinmu
    double precision :: sq2bpi, term, three, twobyx, ya, ya1, xna, x2
    double precision :: zero, half, one, two, eight, one5, ten9, five_pi
    double precision :: onbpi, pi_c, eps_c, del_c, xinf_c, xmin_c, thresh_c, xlarge_c
    double precision :: ch(21)

    zero   = 0.0d0; half = 0.5d0; one = 1.0d0; two = 2.0d0
    three  = 3.0d0; eight = 8.0d0; one5 = 15.0d0; ten9 = 19.0d0
    five_pi  = 1.5707963267948966192d1
    piby2    = 1.5707963267948966192d0
    pi_c     = 3.1415926535897932385d0
    sq2bpi   = 7.9788456080286535588d-1
    pim5     = 7.0796326794896619231d-1
    onbpi    = 3.1830988618379067154d-1
    del_c    = 1.0d-8
    xmin_c   = 4.46d-308
    xinf_c   = 1.79d308
    eps_c    = 1.11d-16
    thresh_c = 16.0d0
    xlarge_c = 1.0d8
    ch = [ -0.67735241822398840964d-23, -0.61455180116049879894d-22, &
             0.29017595056104745456d-20,  0.13639417919073099464d-18, &
             0.23826220476859635824d-17, -0.90642907957550702534d-17, &
            -0.14943667065169001769d-14, -0.33919078305362211264d-13, &
            -0.17023776642512729175d-12,  0.91609750938768647911d-11, &
             0.24230957900482704055d-09,  0.17451364971382984243d-08, &
            -0.33126119768180852711d-07, -0.86592079961391259661d-06, &
            -0.49717367041957398581d-05,  0.76309597585908126618d-04, &
             0.12719271366545622927d-02,  0.17063050710955562222d-02, &
            -0.76852840844786673690d-01, -0.28387654227602353814d+00, &
             0.92187029365045265648d+00 ]

    ex  = x
    enu = alpha
    ncalc = min(nb,0) - 1
    if (.not. ((nb > 0) .and. (ex >= xmin_c) .and. (ex < xlarge_c) .and. &
               (enu >= zero) .and. (enu < one))) then
        by(1) = zero
        return
    end if

    xna = aint(enu + half)
    na  = int(xna)
    if (na == 1) enu = enu - xna

    if (enu == -half) then
        p   = sq2bpi/sqrt(ex)
        ya  =  p*sin(ex)
        ya1 = -p*cos(ex)

    else if (ex < three) then
        ! Temme's scheme for small X
        d = half*ex
        d = -log(d)
        f = enu*d
        e = ex**(-enu)
        if (abs(enu) < del_c) then
            c = onbpi
        else
            c = enu/sin(enu*pi_c)
        end if
        ! sinh(f)/f
        if (abs(f) < one) then
            x2 = f*f
            en = ten9
            s  = one
            do i = 1, 9
                s  = s*x2/en/(en-one) + one
                en = en - two
            end do
        else
            s = (e - one/e)*half/f
        end if
        ! 1/gamma(1-alpha) via Chebyshev
        x2   = enu*enu*eight
        aye  = ch(1)
        even = zero
        alfa = ch(2)
        odd  = zero
        do i = 3, 19, 2
            even = -(aye+aye+even)
            aye  = -even*x2 - aye + ch(i)
            odd  = -(alfa+alfa+odd)
            alfa = -odd*x2 - alfa + ch(i+1)
        end do
        even  = (even*half+aye)*x2 - aye + ch(21)
        odd   = (odd+alfa)*two
        g     = odd*enu + even
        e     = (e + one/e)*half
        f     = two*c*(odd*e + even*s*d)
        e     = enu*enu
        p     = g*c
        q     = onbpi/g
        c     = enu*piby2
        if (abs(c) < del_c) then
            r = one
        else
            r = sin(c)/c
        end if
        r  = pi_c*c*r*r
        c  = one
        d  = -(half*ex)*(half*ex)
        h  = zero
        ya = f + r*q
        ya1 = p
        en = zero
        do
            en  = en + one
            f   = (f*en+p+q)/(en*en-e)
            c   = c*d/en
            p   = p/(en-enu)
            q   = q/(en+enu)
            g   = c*(f+r*q)
            h   = c*p - en*g
            ya  = ya + g
            ya1 = ya1 + h
            if (abs(g/(one+abs(ya))) + abs(h/(one+abs(ya1))) <= eps_c) exit
        end do
        ya  = -ya
        ya1 = -ya1/(half*ex)

    else if (ex < thresh_c) then
        ! Temme's scheme for moderate X
        c  = (half-enu)*(half+enu)
        d  = ex + ex   ! used as b in original
        e  = (ex*onbpi*cos(enu*pi_c)/eps_c)**2
        p  = one
        q  = -ex
        r  = one + ex*ex
        s  = r
        en = two
        do
            if (r*en*en >= e) exit
            en1 = en + one
            d2  = (en-one+c/en)/s
            p   = (en+en-p*d2)/en1
            q   = (-d+q*d2)/en1      ! note: d here is b=ex+ex
            s   = p*p + q*q
            r   = r*s
            en  = en1
        end do
        f  = p/s
        p  = f
        g  = -q/s
        q  = g
        en = en - one
        do while (en > zero)
            r   = en1*(two-p) - two     ! en1 still holds last en+1
            s   = d + en1*q             ! d = b = ex+ex
            d2  = (en-one+c/en)/(r*r+s*s)
            p   = d2*r
            q   = d2*s
            e   = f + one
            f   = p*e - g*q
            g   = q*e + p*g
            en1 = en
            en  = en - one
        end do
        f   = one + f
        d2  = f*f + g*g
        pa  = f/d2
        qa  = -g/d2
        d2  = enu + half - p
        q   = q + ex
        pa1 = (pa*q - qa*d2)/ex
        qa1 = (qa*q + pa*d2)/ex
        d   = ex - piby2*(enu+half)    ! reuse d for phase
        c   = cos(d)
        s   = sin(d)
        d2  = sq2bpi/sqrt(ex)
        ya  = d2*(pa*s + qa*c)
        ya1 = d2*(qa1*s - pa1*c)

    else
        ! Campbell's asymptotic scheme
        na  = 0
        d1  = aint(ex/five_pi)
        i   = int(d1)
        dmu = ((ex-15.0d0*d1) - d1*pim5) - (alpha+half)*piby2
        if (mod(i,2) == 0) then
            cosmu = cos(dmu)
            sinmu = sin(dmu)
        else
            cosmu = -cos(dmu)
            sinmu = -sin(dmu)
        end if
        ddiv = eight*ex
        dmu  = alpha
        d    = sqrt(ex)
        do k = 1, 2
            p  = cosmu
            cosmu = sinmu
            sinmu = -p
            d1   = (two*dmu-one)*(two*dmu+one)
            d2   = zero
            div  = ddiv
            p    = zero
            q    = zero
            q0   = d1/div
            term = q0
            do i = 2, 20
                d2   = d2 + eight
                d1   = d1 - d2
                div  = div + ddiv
                term = -term*d1/div
                p    = p + term
                d2   = d2 + eight
                d1   = d1 - d2
                div  = div + ddiv
                term = term*d1/div
                q    = q + term
                if (abs(term) <= eps_c) exit
            end do
            p    = p + one
            q    = q + q0
            if (k == 1) then
                ya  = sq2bpi*(p*cosmu - q*sinmu)/d
            else
                ya1 = sq2bpi*(p*cosmu - q*sinmu)/d
            end if
            dmu = dmu + one
        end do
    end if

    if (na == 1) then
        h = two*(enu+one)/ex
        if (h > one) then
            if (abs(ya1) > xinf_c/h) then
                h   = zero
                ya  = zero
            end if
        end if
        h   = h*ya1 - ya
        ya  = ya1
        ya1 = h
    end if

    by(1) = ya
    if (nb == 1) then
        ncalc = 1
        return
    end if
    by(2) = ya1
    if (ya1 == zero) then
        ncalc = 1
        return
    end if
    aye    = one + alpha
    twobyx = two/ex
    ncalc  = 2
    do i = 3, nb
        if (twobyx < one) then
            if (abs(by(i-1))*twobyx >= xinf_c/aye) goto 450
        else
            if (abs(by(i-1)) >= xinf_c/aye/twobyx) goto 450
        end if
        by(i) = twobyx*aye*by(i-1) - by(i-2)
        aye   = aye + one
        ncalc = ncalc + 1
    end do
450 continue
    by(ncalc+1:nb) = zero

end subroutine dybesl


subroutine dibesl(x, alpha, nb, ize, b, ncalc)
!> Purpose:
!>   Evaluate modified Bessel `I_{n+alpha}(x)` for `n=0..nb-1`.
!> Notes:
!>   Converted from Cody/Argonne TOMS F77 `dibesl.f` (`RIBESL`) to F90.
    implicit none
    double precision, intent(in)  :: x, alpha
    integer,          intent(in)  :: nb, ize
    double precision, intent(out) :: b(nb)
    integer,          intent(out) :: ncalc

    integer :: k, l, magx, n, nbmx, nend, nstart
    double precision :: const, em, empal, emp2al, en
    double precision :: halfx, one, p, plast, pold, psave, psavel
    double precision :: sum, tempa, tempb, tempc, test, tover, two, zero, half
    double precision :: enmten, ensig, enten, exparg, rtnsig, xlarge

    one    = 1.0d0; two = 2.0d0; zero = 0.0d0; half = 0.5d0; const = 1.585d0
    enten  = 1.0d308; ensig = 1.0d16; rtnsig = 1.0d-4
    enmten = 8.9d-308; xlarge = 1.0d4; exparg = 709.0d0

    ncalc = min(nb,0) - 1
    if (.not. ((nb > 0) .and. (x >= zero) .and. (alpha >= zero) .and. &
               (alpha < one) .and. &
               ((ize == 1 .and. x <= exparg) .or. (ize == 2 .and. x <= xlarge)))) return

    ncalc = nb
    magx  = int(x)

    if (x >= rtnsig) then
        nbmx   = nb - magx
        n      = magx + 1
        en     = dble(n+n) + (alpha+alpha)
        plast  = one
        p      = en/x
        test   = ensig + ensig
        if (2*magx > 5*16) then
            test = sqrt(test*p)
        else
            test = test/const**magx
        end if
        if (nbmx >= 3) then
            tover  = enten/ensig
            nstart = magx + 2
            nend   = nb - 1
            do k = nstart, nend
                n     = k
                en    = en + two
                pold  = plast
                plast = p
                p     = en*plast/x + pold
                if (p > tover) then
                    tover  = enten
                    p      = p/tover
                    plast  = plast/tover
                    psave  = p
                    psavel = plast
                    nstart = n + 1
                    do
                        n     = n + 1
                        en    = en + two
                        pold  = plast
                        plast = p
                        p     = en*plast/x + pold
                        if (p > one) exit
                    end do
                    tempb = en/x
                    test  = pold*plast/ensig
                    test  = test*(half-half/(tempb*tempb))
                    p     = plast*tover
                    n     = n - 1
                    en    = en - two
                    nend  = min(nb,n)
                    do l = nstart, nend
                        ncalc  = l
                        pold   = psavel
                        psavel = psave
                        psave  = en*psavel/x + pold
                        if (psave*psavel > test) then
                            ncalc = ncalc - 1
                            goto 120
                        end if
                    end do
                    ncalc = nend
                    goto 120
                end if
            end do
            n    = nend
            en   = dble(n+n) + (alpha+alpha)
            test = max(test, sqrt(plast*ensig)*sqrt(p+p))
        end if
        do
            n     = n + 1
            en    = en + two
            pold  = plast
            plast = p
            p     = en*plast/x + pold
            if (p >= test) exit
        end do

120     continue
        n      = n + 1
        en     = en + two
        tempb  = zero
        tempa  = one/p
        em     = dble(n) - one
        empal  = em + alpha
        emp2al = (em-one) + (alpha+alpha)
        sum    = tempa*empal*emp2al/em
        nend   = n - nb
        if (nend < 0) then
            b(n) = tempa
            b(n+1:nb) = zero
        else
            if (nend > 0) then
                do l = 1, nend
                    n      = n - 1
                    en     = en - two
                    tempc  = tempb
                    tempb  = tempa
                    tempa  = (en*tempb)/x + tempc
                    em     = em - one
                    emp2al = emp2al - one
                    if (n == 1) goto 150
                    if (n == 2) emp2al = one
                    empal  = empal - one
                    sum    = (sum+tempa*empal)*emp2al/em
                end do
            end if
150         continue
            b(n) = tempa
            if (nb <= 1) then
                sum = (sum+sum) + tempa
                goto 230
            end if
            n      = n - 1
            en     = en - two
            b(n)   = (en*tempa)/x + tempb
            if (n == 1) goto 220
            em     = em - one
            emp2al = emp2al - one
            if (n == 2) emp2al = one
            empal  = empal - one
            sum    = (sum+b(n)*empal)*emp2al/em
        end if
        nend = n - 2
        if (nend > 0) then
            do l = 1, nend
                n      = n - 1
                en     = en - two
                b(n)   = (en*b(n+1))/x + b(n+2)
                em     = em - one
                emp2al = emp2al - one
                if (n == 2) emp2al = one
                empal  = empal - one
                sum    = (sum+b(n)*empal)*emp2al/em
            end do
        end if
        b(1) = two*empal*b(2)/x + b(3)
220     continue
        sum = (sum+sum) + b(1)
230     continue
        if (alpha /= zero) sum = sum*gamma(one+alpha)*(x*half)**(-alpha)
        if (ize == 1) sum = sum*exp(-x)
        tempa = enmten
        if (sum > one) tempa = tempa*sum
        do n = 1, nb
            if (b(n) < tempa) b(n) = zero
            b(n) = b(n)/sum
        end do

    else
        ! Two-term ascending series for small X
        tempa  = one
        empal  = one + alpha
        halfx  = zero
        if (x > enmten) halfx = half*x
        if (alpha /= zero) tempa = halfx**alpha/gamma(empal)
        if (ize == 2) tempa = tempa*exp(-x)
        tempb = zero
        if ((x+one) > one) tempb = halfx*halfx
        b(1) = tempa + tempa*tempb/empal
        if ((x /= zero) .and. (b(1) == zero)) ncalc = 0
        if (nb > 1) then
            if (x == zero) then
                b(2:nb) = zero
            else
                tempc = halfx
                tover = (enmten+enmten)/x
                if (tempb /= zero) tover = enmten/tempb
                do n = 2, nb
                    tempa = tempa/empal
                    empal = empal + one
                    tempa = tempa*tempc
                    if (tempa <= tover*empal) tempa = zero
                    b(n) = tempa + tempa*tempb/empal
                    if ((b(n) == zero) .and. (ncalc > n)) ncalc = n-1
                end do
            end if
        end if
    end if

end subroutine dibesl


subroutine dkbesl(x, alpha, nb, ize, bk, ncalc)
!> Purpose:
!>   Evaluate modified Bessel `K_{n+alpha}(x)` for `n=0..nb-1`.
!> Notes:
!>   Converted from Cody/Argonne TOMS F77 `dkbesl.f` to F90.
    implicit none
    double precision, intent(in)  :: x, alpha
    integer,          intent(in)  :: nb, ize
    double precision, intent(out) :: bk(nb)
    integer,          intent(out) :: ncalc

    integer :: i, iend, itemp, j, k, m, mplus1
    double precision :: a, blpha, bk1, bk2, c, d, dm, d1, d2, d3
    double precision :: enu, ex, f0, f1, f2, four, half, one, p0, q0
    double precision :: r2, ratio, twonu, twox, t1, t2, wminf, zero, two, tinyx
    double precision :: eps_k, sqxmin, xinf_k, xmax_k, xmin_k, d_c
    double precision, parameter :: p(8) = [ &
        0.805629875690432845d0,  0.204045500205365151d2, &
        0.157705605106676174d3,  0.536671116469207504d3, &
        0.900382759291288778d3,  0.730923886650660393d3, &
        0.229299301509425145d3,  0.822467033424113231d0  ]
    double precision, parameter :: q_coef(7) = [ &
        0.294601986247850434d2,  0.277577868510221208d3, &
        0.120670325591027438d4,  0.276291444159791519d4, &
        0.344374050506564618d4,  0.221063190113378647d4, &
        0.572267338359892221d3  ]
    double precision, parameter :: r_coef(5) = [ &
        -0.48672575865218401848d0,  0.13079485869097804016d2, &
        -0.10196490580880537526d3,  0.34765409106507813131d3, &
         0.34958981245219347820d-3 ]
    double precision, parameter :: s_coef(4) = [ &
        -0.25579105509976461286d2,  0.21257260432226544008d3, &
        -0.61069018684944109624d3,  0.42269668805777760407d3  ]
    double precision, parameter :: t_coef(6) = [ &
        0.16125990452916363814d-9, 0.25051878502858255354d-7, &
        0.27557319615147964774d-5, 0.19841269840928373686d-3, &
        0.83333333333334751799d-2, 0.16666666666666666446d0  ]
    double precision, parameter :: estm(6) = [ &
        5.20583d1, 5.7607d0, 2.7782d0, 1.44303d1, 1.853004d2, 9.3715d0 ]
    double precision, parameter :: estf(7) = [ &
        4.18341d1, 7.1075d0, 6.4306d0, 4.25110d1, 1.35633d0, 8.45096d1, 2.0d1 ]
    ! LOG(2) - Euler's constant; D = sqrt(2/pi)
    double precision, parameter :: a_c = 0.11593151565841244881d0
    double precision, parameter :: d_sqrt2pi = 0.797884560802865364d0

    half   = 0.5d0; one = 1.0d0; two = 2.0d0; zero = 0.0d0
    four   = 4.0d0; tinyx = 1.0d-10
    eps_k  = 2.22d-16; sqxmin = 1.49d-154; xinf_k = 1.79d308
    xmin_k = 2.23d-308; xmax_k = 705.342d0

    ex  = x
    enu = alpha
    ncalc = min(nb,0) - 2
    if (.not. ((nb > 0) .and. (enu >= zero) .and. (enu < one) .and. &
               (ize >= 1) .and. (ize <= 2) .and. &
               ((ize /= 1) .or. (ex <= xmax_k)) .and. (ex > zero))) return

    k = 0
    if (enu < sqxmin) enu = zero
    if (enu > half) then
        k   = 1
        enu = enu - one
    end if
    twonu = enu + enu
    iend  = nb + k - 1
    c     = enu*enu
    d3    = -c

    if (ex <= one) then
        ! Calculate P0 = gamma(1+alpha)*(2/x)^alpha, Q0 = gamma(1-alpha)*(x/2)^alpha
        d1 = zero; d2 = p(1); t1 = one; t2 = q_coef(1)
        do i = 2, 7, 2
            d1 = c*d1 + p(i)
            d2 = c*d2 + p(i+1)
            t1 = c*t1 + q_coef(i)
            t2 = c*t2 + q_coef(i+1)
        end do
        d1 = enu*d1; t1 = enu*t1
        f1 = log(ex)
        f0 = a_c + enu*(p(8)-enu*(d1+d2)/(t1+t2)) - f1
        q0 = exp(-enu*(a_c - enu*(p(8)+enu*(d1-d2)/(t1-t2)) - f1))
        f1 = enu*f0
        p0 = exp(f1)
        ! Calculate F0
        d1 = r_coef(5); t1 = one
        do i = 1, 4
            d1 = c*d1 + r_coef(i)
            t1 = c*t1 + s_coef(i)
        end do
        if (abs(f1) <= half) then
            f1 = f1*f1
            d2 = zero
            do i = 1, 6
                d2 = f1*d2 + t_coef(i)
            end do
            d2 = f0 + f0*f1*d2
        else
            d2 = sinh(f1)/enu
        end if
        f0 = d2 - enu*d1/(t1*p0)

        if (ex <= tinyx) then
            bk(1) = f0 + ex*f0
            if (ize == 1) bk(1) = bk(1) - ex*bk(1)
            ratio = p0/f0
            c     = ex*xinf_k
            if (k /= 0) then
                ncalc = -1
                if (bk(1) >= c/ratio) goto 500
                bk(1) = ratio*bk(1)/ex
                twonu = twonu + two
                ratio = twonu
            end if
            ncalc = 1
            if (nb == 1) goto 500
            ncalc = -1
            do i = 2, nb
                if (ratio >= c) goto 500
                bk(i) = ratio/ex
                twonu = twonu + two
                ratio = twonu
            end do
            ncalc = 1
            goto 420
        else
            ! 1e-10 < X <= 1
            c      = one
            r2     = ex*ex/four       ! x2by4
            p0     = half*p0
            q0     = half*q0
            d1     = -one
            d2     = zero
            bk1    = zero
            bk2    = zero
            f1     = f0
            f2     = p0
            do
                d1   = d1 + two
                d2   = d2 + one
                d3   = d1 + d3
                c    = r2*c/d2
                f0   = (d2*f0 + p0 + q0)/d3
                p0   = p0/(d2-enu)
                q0   = q0/(d2+enu)
                t1   = c*f0
                t2   = c*(p0 - d2*f0)
                bk1  = bk1 + t1
                bk2  = bk2 + t2
                if (abs(t1/(f1+bk1)) <= eps_k .and. abs(t2/(f2+bk2)) <= eps_k) exit
            end do
            bk1 = f1 + bk1
            bk2 = two*(f2+bk2)/ex
            if (ize == 2) then
                d1  = exp(ex)
                bk1 = bk1*d1
                bk2 = bk2*d1
            end if
            wminf = estf(1)*ex + estf(2)
        end if

    else if (eps_k*ex > one) then
        ! X >> 1/eps
        ncalc = nb
        bk1 = one/(d_sqrt2pi*sqrt(ex))
        bk(1:nb) = bk1
        goto 500

    else
        ! X > 1
        twox  = ex + ex
        blpha = zero
        ratio = zero
        if (ex <= four) then
            ! 1 <= X <= 4
            d2 = aint(estm(1)/ex + estm(2))
            m  = int(d2)
            d1 = d2 + d2
            d2 = d2 - half; d2 = d2*d2
            do i = 2, m
                d1    = d1 - two
                d2    = d2 - d1
                ratio = (d3+d2)/(twox+d1-ratio)
            end do
            ! Compute I(|alpha|,x) and I(|alpha|+1,x) by backward recurrence
            d2 = aint(estm(3)*ex + estm(4))
            m  = int(d2)
            c  = abs(enu)
            d3 = c + c
            d1 = d3 - one
            f1 = xmin_k
            f0 = (two*(c+d2)/ex + half*ex/(c+d2+one))*xmin_k
            do i = 3, m
                d2    = d2 - one
                f2    = (d3+d2+d2)*f0
                blpha = (one+d1/d2)*(f2+blpha)
                f2    = f2/ex + f1
                f1    = f0
                f0    = f2
            end do
            f1 = (d3+two)*f0/ex + f1
            d1 = zero; t1 = one
            do i = 1, 7
                d1 = c*d1 + p(i)
                t1 = c*t1 + q_coef(i)
            end do
            p0  = exp(c*(a_c+c*(p(8)-c*d1/t1)-log(ex)))/ex
            f2  = (c+half-ratio)*f1/ex
            bk1 = p0 + (d3*f0-f2+f0+blpha)/(f2+f1+f0)*p0
            if (ize == 1) bk1 = bk1*exp(-ex)
            wminf = estf(3)*ex + estf(4)
        else
            ! X > 4
            dm  = aint(estm(5)/ex + estm(6))
            m   = int(dm)
            d2  = dm - half; d2 = d2*d2
            d1  = dm + dm
            do i = 2, m
                dm    = dm - one
                d1    = d1 - two
                d2    = d2 - d1
                ratio = (d3+d2)/(twox+d1-ratio)
                blpha = (ratio+ratio*blpha)/dm
            end do
            bk1 = one/((d_sqrt2pi+d_sqrt2pi*blpha)*sqrt(ex))
            if (ize == 1) bk1 = bk1*exp(-ex)
            wminf = estf(5)*(ex-abs(ex-estf(7))) + estf(6)
        end if
        bk2 = bk1 + bk1*(enu+half-ratio)/ex
    end if

    ! Recurrence to get all orders
    ncalc = nb
    bk(1) = bk1
    if (iend == 0) goto 500
    j = 2 - k
    if (j > 0) bk(j) = bk2
    if (iend == 1) goto 500
    m     = min(int(wminf-enu), iend)
    itemp = 1
    do i = 2, m
        t1    = bk1
        bk1   = bk2
        twonu = twonu + two
        if (ex < one) then
            if (bk1 >= (xinf_k/twonu)*ex) goto 195
        else
            if (bk1/ex >= xinf_k/twonu) goto 195
        end if
        bk2   = twonu/ex*bk1 + t1
        itemp = i
        j     = j + 1
        if (j > 0) bk(j) = bk2
    end do
195 continue
    m = itemp
    if (m == iend) goto 500
    ratio  = bk2/bk1
    mplus1 = m + 1
    ncalc  = -1
    do i = mplus1, iend
        twonu = twonu + two
        ratio = twonu/ex + one/ratio
        j     = j + 1
        if (j > 1) then
            bk(j) = ratio
        else
            if (bk2 >= xinf_k/ratio) goto 500
            bk2 = ratio*bk2
        end if
    end do
    ncalc = max(mplus1-k, 1)
    if (ncalc == 1) bk(1) = bk2
    if (nb == 1) goto 500
420 continue
    j = ncalc + 1
    do i = j, nb
        if (bk(ncalc) >= xinf_k/bk(i)) goto 500
        bk(i) = bk(ncalc)*bk(i)
        ncalc = i
    end do
500 continue

end subroutine dkbesl
