!> @anchor      nrlmsise00Class
!!
!> @brief       Object oriented NRLMSISE-00 based on C source code package - release 20041227
!! @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>CHK: 08.01.2018 (Initial (back-)port to FORTRAN from C)
!!              </ul>
!!
!> @details     The NRLMSISE-00 model was developed by Mike Picone, Alan Hedin, and
!!              Doug Drob. They also wrote a NRLMSISE-00 distribution package in
!!              FORTRAN which is available at
!!              http://uap-www.nrl.navy.mil/models_web/msis/msis_home.htm
!!
!!              Dominik Brodowski implemented and maintains this C version. You can
!!              reach him at mail@brodo.de. See the file "DOCUMENTATION" for details,
!!              and check http://www.brodo.de/english/pub/nrlmsise/index.html for
!!              updated releases of this package.
!!
!> @anchor      nrlmsise00Class
!!
!!------------------------------------------------------------------------------------------------
module nrlmsise00Class
    use slam_types,             only: sp
    use nrlmsise00data

    type, public :: ap_array
        real(sp),dimension(0:6)     :: a
    end type ap_array

    type, public :: nrlmsise_input
        integer         :: year     ! year, currently ignored
        integer         :: doy      ! day of year
        real(sp)        :: sec      ! seconds in day (UT)
        real(sp)        :: alt      ! altitude in kilometers
        real(sp)        :: g_lat    ! geodetic latitude
        real(sp)        :: g_long   ! geodetic longitude
        real(sp)        :: lst      ! local apparent solar time (hours), see note below
        real(sp)        :: f107A    ! 81 day average of F10.7 flux (centered on doy)
        real(sp)        :: f107     ! daily F10.7 flux for previous day
        real(sp)        :: ap       ! magnetic index(daily)
        type(ap_array)  :: ap_a     ! see above
    end type nrlmsise_input

    type, public :: nrlmsise_output
        real(sp),dimension(0:8)     :: d   ! densities
        real(sp),dimension(0:1)     :: t   ! temperatures    
    end type nrlmsise_output

    type, public :: nrlmsise_flags
        integer,dimension(0:23)     :: switches
        real(sp),dimension(0:23)    :: sw
        real(sp),dimension(0:23)    :: swc
    end type nrlmsise_flags

    type, public :: Nrlmsise00_class

        type(nrlmsise_flags)            :: flags
        
        ! PARMB
        real(sp)                        :: gsurf
        real(sp)                        :: re

        ! GTS3C
        real(sp)                        :: dd

        ! DMIX
        real(sp)                        :: dm04
        real(sp)                        :: dm16
        real(sp)                        :: dm28
        real(sp)                        :: dm32
        real(sp)                        :: dm40
        real(sp)                        :: dm01
        real(sp)                        :: dm14

        ! MESO7
        real(sp),dimension(0:4)         :: meso_tn1
        real(sp),dimension(0:3)         :: meso_tn2
        real(sp),dimension(0:4)         :: meso_tn3
        real(sp),dimension(0:1)         :: meso_tgn1
        real(sp),dimension(0:1)         :: meso_tgn2
        real(sp),dimension(0:1)         :: meso_tgn3

        ! LPOLY
        real(sp)                        :: dfa
        real(sp),dimension(0:3,0:8)     :: plg
        real(sp)                        :: ctloc
        real(sp)                        :: stloc
        real(sp)                        :: c2tloc
        real(sp)                        :: s2tloc
        real(sp)                        :: s3tloc
        real(sp)                        :: c3tloc
        real(sp)                        :: apdf
        real(sp),dimension(0:3)         :: apt

    contains

        procedure,public  :: tselec
        procedure,private :: glatf
        procedure,private :: ccor
        procedure,private :: ccor2
        procedure,private :: scalh
        procedure,private :: dnet
        procedure,private :: splini
        procedure,private :: splint
        procedure,private :: spline
        procedure,private :: zeta
        procedure,private :: densm
        procedure,private :: densu
        procedure,private :: g0
        procedure,private :: sumex
        procedure,private :: sg0
        procedure,private :: globe7
        procedure,private :: glob7s
        procedure,public  :: gtd7
        procedure,public  :: gtd7d
        procedure,public  :: ghp7
        procedure,public  :: gts7

    end type Nrlmsise00_class

    ! Constructor
    interface Nrlmsise00_class
        module procedure constructor
    end interface Nrlmsise00_class

contains

    ! ====================================================================
    !!
    !>  @brief      Constructor that should be used to initialize variables.
    !!
    !>  @author     Christopher Kebschull
    !>  @date       <ul>
    !!                  <li>ChK: 08.01.2018 (initial implementation)</li>
    !!              </ul>
    !>  @anchor     constructor
    !!
    ! --------------------------------------------------------------------
    type(Nrlmsise00_class) function constructor()

        integer :: i,k

        ! Copying arrays from C-style arrays to Fortran style - this can be optimized
        !  in the future by als back porting the array access in the routines.
        ! POWER7
        do i=1,150
            pt(i-1) = pt_data(i)
        end do
        do i=1,150
            do k=1,9
                pd(k-1,i-1) = pd_data(i,k)
            end do
        end do
        do i=1,150
            ps(i-1) = ps_data(i)
        end do
        do i=1,25
            do k=1,2
                pdl(k-1,i-1) = pdl_data(i,k)
            end do
        end do
        do i=1,100
            do k=1,4
                ptl(k-1,i-1) = ptl_data(i,k)
            end do
        end do
        do i=1,100
            do k=1,10
                pma(k-1,i-1) = pma_data(i,k)
            end do
        end do

        ! LOWER7
        do i=1,10
            ptm(i-1) = ptm_data(i)
        end do
        do i=1,10
            do k=1,8
                pdm(k-1,i-1) = pdm_data(i,k)
            end do
        end do
        do i=1,10
            pavgm(i-1) = pavgm_data(i)
        end do

        constructor%meso_tn1 = 0.0
        constructor%meso_tn2 = 0.0
        constructor%meso_tn3 = 0.0
        constructor%meso_tgn1 = 0.0
        constructor%meso_tgn2 = 0.0
        constructor%meso_tgn3 = 0.0
        constructor%plg = 0.0
        constructor%apt = 0.0

    end function constructor

    ! -------------------------------------------------------------------
    ! ------------------------------ TSELEC -----------------------------
    ! -------------------------------------------------------------------
    subroutine tselec(this, flags)

        implicit none

        class(Nrlmsise00_class),intent(inout)   :: this
        type(nrlmsise_flags),intent(inout)      :: flags
        integer                                 :: i

        do i=0,23
            if (i /= 9) then
                if (flags%switches(i)==1) then
                    flags%sw(i)=1.0
                else
                    flags%sw(i)=0.0
                end if
                if (flags%switches(i)>0) then
                    flags%swc(i)=1.0
                else
                    flags%swc(i)=0.0
                end if
            else
                flags%sw(i)=flags%switches(i)
                flags%swc(i)=flags%switches(i)
            end if
        end do
    end subroutine tselec

    ! -------------------------------------------------------------------
    ! ------------------------------ GLATF ------------------------------
    ! -------------------------------------------------------------------
    subroutine glatf(this, lat, gv, reff)

        implicit none

        class(Nrlmsise00_class),intent(inout)   :: this
        real(sp),intent(in)     :: lat
        real(sp),intent(inout)  :: gv
        real(sp),intent(inout)  :: reff
        real(sp),parameter      :: dgtr = 1.74533e-2
        real(sp)                :: c2

        c2 = cos(2.0 * dgtr * lat)
        gv = 980.616 * (1.0 - 0.0026373 * c2)
        reff = 2.0 * (gv) / (3.085462e-6 + 2.27e-9 * c2) * 1.e-5
    end subroutine glatf

    ! -------------------------------------------------------------------
    ! ------------------------------ CCOR -------------------------------
    ! -------------------------------------------------------------------
    !        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
    !        ALT - altitude
    !        R - target ratio
    !        H1 - transition scale length
    !        ZH - altitude of 1/2 R
    real(sp) function ccor(this, alt, r, h1, zh)

        implicit none

        class(Nrlmsise00_class),intent(inout)   :: this
        real(sp),intent(in)     :: alt
        real(sp),intent(in)     :: r
        real(sp),intent(in)     :: h1
        real(sp),intent(in)     :: zh

        real(sp)                :: e
        real(sp)                :: ex
        real(sp)                :: res

        e = (alt - zh) / h1
        if (e > 70.0) then
            ccor = exp(0.0)
            return
        end if
        if (e < -70.0) then
            ccor = exp(r)
            return
        end if
        ex = exp(e)
        e = r / (1.0 + ex)

        ccor = exp(e)

    end function ccor

    ! -------------------------------------------------------------------
    ! ------------------------------ CCOR -------------------------------
    ! -------------------------------------------------------------------
    !        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
    !          ALT - altitude
    !          R - target ratio
    !          H1 - transition scale length
    !          ZH - altitude of 1/2 R
    !          H2 - transition scale length #2 ?
    real(sp) function ccor2(this, alt, r, h1, zh, h2)

        implicit none

        class(Nrlmsise00_class),intent(inout)   :: this
        real(sp),intent(in)     :: alt
        real(sp),intent(in)     :: r
        real(sp),intent(in)     :: h1
        real(sp),intent(in)     :: zh
        real(sp),intent(in)     :: h2

        real(sp)                :: e1, e2
        real(sp)                :: ex1, ex2
        real(sp)                :: ccor2v

        e1 = (alt - zh) / h1
        e2 = (alt - zh) / h2
        if ((e1 > 70.0) &
            .or. (e2 > 70.0)) then
            ccor2 = exp(0.0)
            return
        end if
        if ((e1 < -70.0) &
            .and. (e2 < -70.0)) then
            ccor2 = exp(r)
            return
        end if
        ex1 = exp(e1)
        ex2 = exp(e2)
        ccor2v = r / (1.0 + 0.5 * (ex1 + ex2))

        ccor2 =  exp(ccor2v)
        return
    end function ccor2

    ! -------------------------------------------------------------------
    ! ------------------------------- SCALH -----------------------------
    ! -------------------------------------------------------------------
    real(sp) function scalh(this, alt, xm, temp)

        implicit none

        class(Nrlmsise00_class),intent(inout)   :: this
        real(sp),intent(in)   :: alt
        real(sp),intent(in)   :: xm
        real(sp),intent(in)   :: temp

        real(sp)              :: g
        real(sp),parameter    :: rgas=831.4

        g = this%gsurf / ((1.0 + alt/this%re)**2.0)
        g = rgas * temp / (g * xm)
        scalh = g
        return

    end function scalh

    ! -------------------------------------------------------------------
    ! -------------------------------- DNET -----------------------------
    ! -------------------------------------------------------------------
    !       TURBOPAUSE CORRECTION FOR MSIS MODELS
    !        Root mean density
    !        DD - diffusive density
    !        DM - full mixed density
    !        ZHM - transition scale length
    !        XMM - full mixed molecular weight
    !        XM  - species molecular weight
    !        DNET - combined density
    real(sp) function dnet(this, dd, dm, zhm, xmm, xm)

        implicit none

        class(Nrlmsise00_class),intent(inout)   :: this

        real(sp),intent(inout)  :: dd
        real(sp),intent(in)     :: dm
        real(sp),intent(in)     :: zhm
        real(sp),intent(in)     :: xmm
        real(sp),intent(in)     :: xm

        real(sp) a
        real(sp) ylog

        a  = zhm / (xmm-xm)

        if (.not. ((dm > 0.0) &
            .and. (dd > 0.0))) then
            write (*,*) "dnet log error", dm, dd, xm
            if ((dd ==0.0) .and. (dm == 0.0)) then
                dd = 1.0
            end if
            if (dm == 0.0) then
                dnet = dd
                return
            end if
            if (dd == 0.0) then
                dnet = dm
                return
            end if
        end if

        ylog = a * log(dm/dd)

        if (ylog < -10.0) then
            dnet = dd
            return
        end if
        if (ylog > 10.0) then
            dnet = dm
            return
        end if
        a = dd * (1.0 + exp(ylog))**(1.0/a)
        dnet = a
        return
    end function dnet

    ! -------------------------------------------------------------------
    ! ------------------------------- SPLINI ----------------------------
    ! -------------------------------------------------------------------
    !      INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
    !        XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
    !        Y2A: ARRAY OF SECOND DERIVATIVES
    !        N: SIZE OF ARRAYS XA,YA,Y2A
    !        X: ABSCISSA ENDPOINT FOR INTEGRATION
    !        Y: OUTPUT VALUE
    !
    subroutine splini(this, xa, ya, y2a, n, x, y)

        implicit none

        class(Nrlmsise00_class),intent(inout)   :: this
        real(sp),dimension(0:n-1),intent(in)    :: xa
        real(sp),dimension(0:n-1),intent(in)    :: ya
        real(sp),dimension(0:n-1),intent(in)    :: y2a
        integer,intent(in)                      :: n
        real(sp),intent(in)                     :: x
        real(sp),intent(inout)                  :: y

        real(sp)                    :: yi
        integer                     :: klo
        integer                     :: khi
        real(sp)                    :: xx, h, a, b, a2, b2

        yi=0
        klo=0
        khi=1

        do while ((x > xa(klo)) .and. (khi < n))
            xx=x
            if (khi < (n-1)) then
                if (x < xa(khi)) then
                    xx = x
                else 
                    xx = xa(khi)
                end if
            end if
            h = xa(khi) - xa(klo)
            a = (xa(khi) - xx)/h
            b = (xx - xa(klo))/h
            a2 = a*a
            b2 = b*b
            yi = yi + ((1.0 - a2) * ya(klo) / 2.0 + b2 * ya(khi) / 2.0 + ((-(1.0+a2*a2)/4.0 + a2/2.0) * y2a(klo) + (b2*b2/4.0 - b2/2.0) * y2a(khi)) * h * h / 6.0) * h
            klo = klo + 1
            khi = khi + 1
        end do
        y = yi

    end subroutine splini

    ! -------------------------------------------------------------------
    ! ------------------------------- SPLINT ----------------------------
    ! -------------------------------------------------------------------
    !      CALCULATE CUBIC SPLINE INTERP VALUE
    !       ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL.
    !       XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
    !       Y2A: ARRAY OF SECOND DERIVATIVES
    !       N: SIZE OF ARRAYS XA,YA,Y2A
    !       X: ABSCISSA FOR INTERPOLATION
    !       Y: OUTPUT VALUE
    !
    subroutine splint (this, xa, ya, y2a, n, x, y)

        implicit none

        class(Nrlmsise00_class),intent(inout)   :: this
        real(sp),dimension(0:n-1),intent(in)    :: xa
        real(sp),dimension(0:n-1),intent(in)    :: ya
        real(sp),dimension(0:n-1),intent(in)    :: y2a
        integer,intent(in)                      :: n
        real(sp),intent(in)                     :: x
        real(sp),intent(inout)                  :: y

        integer                                 :: klo
        integer                                 :: khi
        integer                                 :: k
        real(sp)                                :: h
        real(sp)                                :: a, b, yi

        klo=0
        khi=n-1

        do while ((khi - klo)>1)
            k=(khi+klo)/2
            if (xa(k)>x) then
                khi=k
            else
                klo=k
            end if
        end do

        h = xa(khi) - xa(klo)
        if (h == 0.0) then
            write (*,*) "bad XA input to splint"
        end if
        a = (xa(khi) - x)/h
        b = (x - xa(klo))/h
        yi = a * ya(klo) + b * ya(khi) + ((a*a*a - a) * y2a(klo) + (b*b*b - b) * y2a(khi)) * h * h/6.0
        y = yi
    end subroutine splint

    ! -------------------------------------------------------------------
    ! ------------------------------- SPLINE ----------------------------
    ! -------------------------------------------------------------------
    !       CALCULATE 2ND DERIVATIVES OF CUBIC SPLINE INTERP FUNCTION
    !       ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL
    !       X,Y: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
    !       N: SIZE OF ARRAYS X,Y
    !       YP1,YPN: SPECIFIED DERIVATIVES AT X(0) AND X(N-1) VALUES
    !                >= 1E30 SIGNAL SIGNAL SECOND DERIVATIVE ZERO
    !       Y2: OUTPUT ARRAY OF SECOND DERIVATIVES
    subroutine spline(this, x, y, n, yp1, ypn, y2)

        implicit none

        class(Nrlmsise00_class),intent(inout)       :: this
        real(sp),dimension(0:(n-1)),intent(in)      :: x
        real(sp),dimension(0:(n-1)),intent(in)      :: y
        integer,intent(in)                          :: n
        real(sp),intent(in)                         :: yp1
        real(sp),intent(in)                         :: ypn
        real(sp),dimension(0:(n-1)),intent(inout)   :: y2

        real(sp),dimension(:),allocatable           :: u
        real(sp)                                    :: sig, p, qn, un
        integer                                     :: i, k

        allocate(u(0:(n-1)))
        if (.not. allocated(u)) then
            write (*,*) "Out Of Memory in spline - ERROR"
            return
        end if

        if (yp1 > 0.99D30) then
            y2(0) = 0.0
            u(0) = 0.0
        else
            y2(0)=-0.5
            u(0)=(3.0/(x(1)-x(0)))*((y(1)-y(0))/(x(1)-x(0))-yp1)
        end if

        do i=1,(n-2)
            sig = (x(i)-x(i-1))/(x(i+1) - x(i-1))
            p = sig * y2(i-1) + 2.0
            y2(i) = (sig - 1.0) / p
            u(i) = (6.0 * ((y(i+1) - y(i))/(x(i+1) - x(i)) -(y(i) - y(i-1)) / (x(i) - x(i-1)))/(x(i+1) - x(i-1)) - sig * u(i-1))/p
        end do
        if (ypn > 0.99d30) then
            qn = 0.0
            un = 0.0
        else
            qn = 0.5
            un = (3.0 / (x(n-1) - x(n-2))) * (ypn - (y(n-1) - y(n-2))/(x(n-1) - x(n-2)))
        end if

        y2(n-1) = (un - qn * u(n-2)) / (qn * y2(n-2) + 1.0)
        do k=n-2, 0, -1
            y2(k) = y2(k) * y2(k+1) + u(k)
        end do

        deallocate(u)

    end subroutine spline

    ! -------------------------------------------------------------------
    ! ------------------------------- DENSM -----------------------------
    ! -------------------------------------------------------------------

    real(sp) pure function zeta(this, zz, zl)

        implicit none

        class(Nrlmsise00_class),intent(in)  :: this
        real(sp),intent(in)                 :: zz
        real(sp),intent(in)                 :: zl

        zeta = ((zz-zl)*(this%re+zl)/(this%re+zz))
        return

    end function zeta

    !      Calculate Temperature and Density Profiles for lower atmos.
    real(sp) function densm (this, alt, d0, xm, tz, mn3, zn3, tn3, tgn3, mn2, zn2, tn2, tgn2)

        implicit none

        class(Nrlmsise00_class),intent(inout)       :: this
        real(sp),intent(in)                         :: alt
        real(sp),intent(in)                         :: d0
        real(sp),intent(in)                         :: xm
        real(sp),intent(inout)                      :: tz
        integer,intent(in)                          :: mn3
        real(sp),dimension(0:(mn3-1)),intent(in)    :: zn3
        real(sp),dimension(0:(mn3-1)),intent(in)    :: tn3
        real(sp),dimension(0:1),intent(in)          :: tgn3
        integer,intent(in)                          :: mn2
        real(sp),dimension(0:(mn2-1)),intent(in)    :: zn2
        real(sp),dimension(0:(mn2-1)),intent(in)    :: tn2
        real(sp),dimension(0:1),intent(in)          :: tgn2

        real(sp),dimension(0:9)                     :: xs, ys, y2out
        real(sp),parameter                          :: rgas = 831.4
        real(sp)                                    :: z, z1, z2, t1, t2, zg, zgdif
        real(sp)                                    :: yd1, yd2
        real(sp)                                    :: x, y, yi
        real(sp)                                    :: expl, gamm, glb
        real(sp)                                    :: densm_tmp
        integer                                     :: mn
        integer                                     :: k

        densm_tmp = d0
        yd1 = 0.0
        yd2 = 0.0

        if (alt > zn2(0)) then
            if (xm == 0.0) then
                densm = tz
                return
            else
                densm = d0
                return
            end if
        end if

        ! STRATOSPHERE/MESOSPHERE TEMPERATURE
        if (alt > zn2(mn2-1)) then
            z=alt
        else
            z=zn2(mn2-1)
        end if

        mn = mn2
        z1 = zn2(0)
        z2 = zn2(mn-1)
        t1 = tn2(0)
        t2 = tn2(mn-1)
        zg = this%zeta(z, z1)
        zgdif = this%zeta(z2, z1)

        ! set up spline nodes
        do k=0,(mn-1)
            xs(k)=this%zeta(zn2(k),z1)/zgdif
            ys(k)=1.0 / tn2(k)
        end do
        yd1 = -tgn2(0) / (t1*t1) * zgdif
        yd2 = -tgn2(1) / (t2*t2) * zgdif * (((this%re+z2)/(this%re+z1))**2.0)

        ! calculate spline coefficients
        call this%spline (xs, ys, mn, yd1, yd2, y2out)
        x = zg/zgdif
        call this%splint (xs, ys, y2out, mn, x, y)

        ! temperature at altitude
        tz = 1.0 / y
        if (xm /= 0.0) then
            ! calaculate stratosphere / mesospehere density
            glb = this%gsurf / ((1.0 + z1/this%re)**2.0)
            gamm = xm * glb * zgdif / rgas

            ! Integrate temperature profile
            call this%splini(xs, ys, y2out, mn, x, yi)
            expl=gamm*yi
            if (expl>50.0) then
                expl=50.0
            end if
            ! Density at altitude
            densm_tmp = densm_tmp * (t1 / tz) * exp(-expl)
        end if

        if (alt>zn3(0)) then
            if (xm==0.0) then
                densm = tz
                return
            else
                densm = densm_tmp
                return
            end if
        end if

        ! troposhere / stratosphere temperature
        z = alt
        mn = mn3
        z1=zn3(0)
        z2=zn3(mn-1)
        t1=tn3(0)
        t2=tn3(mn-1)
        zg=this%zeta(z,z1)
        zgdif=this%zeta(z2,z1)

        ! set up spline nodes
        do k=0,mn-1
            xs(k) = this%zeta(zn3(k),z1) / zgdif
            ys(k) = 1.0 / tn3(k)
        end do
        yd1= yd1 - (tgn3(0) / (t1*t1) * zgdif)
        yd2= yd2 - (tgn3(1) / (t2*t2) * zgdif * (((this%re+z2)/(this%re+z1))**2.0))

        ! calculate spline coefficients
        call this%spline (xs, ys, mn, yd1, yd2, y2out)
        x = zg/zgdif
        call this%splint (xs, ys, y2out, mn, x, y)

        ! temperature at altitude
        tz = 1.0 / y
        if (xm /= 0.0) then
            ! calaculate tropospheric / stratosphere density
            glb = this%gsurf / ((1.0 + z1/this%re)**2.0)
            gamm = xm * glb * zgdif / rgas

            ! Integrate temperature profile
            call this%splini(xs, ys, y2out, mn, x, yi)
            expl=gamm*yi
            if (expl > 50.0) then
                expl = 50.0
            end if
            ! Density at altitude
            densm_tmp = densm_tmp * (t1 / tz) * exp(-expl)

        end if

        if (xm == 0.0) then
            densm = tz
            return
        else
            densm = densm_tmp
            return
        end if
    end function densm

    ! -------------------------------------------------------------------
    ! ------------------------------- DENSU -----------------------------
    ! -------------------------------------------------------------------
    !      Calculate Temperature and Density Profiles for MSIS models
    !       New lower thermo polynomial
    !
    real(sp) function densu (this, alt, dlb, tinf, tlb, xm, alpha, tz, zlb, s2, mn1, zn1, tn1, tgn1)

        implicit none

        class(Nrlmsise00_class),intent(inout)       :: this
        real(sp),intent(in)                         :: alt
        real(sp),intent(in)                         :: dlb
        real(sp),intent(in)                         :: tinf
        real(sp),intent(in)                         :: tlb
        real(sp),intent(in)                         :: xm
        real(sp),intent(in)                         :: alpha
        real(sp),intent(inout)                      :: tz
        real(sp),intent(in)                         :: zlb
        real(sp),intent(in)                         :: s2
        integer,intent(in)                          :: mn1
        real(sp),dimension(0:(mn1-1)),intent(inout) :: zn1
        real(sp),dimension(0:(mn1-1)),intent(inout) :: tn1
        real(sp),dimension(0:1),intent(inout)       :: tgn1

        real(sp)                :: yd2, yd1, x, y
        real(sp), parameter     :: rgas=831.4
        real(sp)                :: densu_temp
        real(sp)                :: za, z, zg2, tt, ta
        real(sp)                :: dta, z1, z2, t1, t2, zg, zgdif
        integer                 :: mn
        integer                 :: k
        real(sp)                :: glb
        real(sp)                :: expl
        real(sp)                :: yi
        real(sp)                :: densa
        real(sp)                :: gamma, gamm
        real(sp),dimension(0:4) :: xs, ys, y2out

        x = 0
        z1 = 0
        t1 = 0
        mn=0
        zgdif = 0
        densu_temp=1.0
        yd1 = 0.0
        yd2 = 0.0

        !* joining altitudes of Bates and spline
        za=zn1(0)
        if (alt > za) then
            z=alt
        else
            z=za
        end if

        ! geopotential altitude difference from ZLB
        zg2 = this%zeta(z, zlb)

        ! Bates temperature
        tt = tinf - (tinf - tlb) * exp(-s2*zg2)
        ta = tt
        tz = tt
        densu_temp = tz

        if (alt < za) then
            ! calculate temperature below ZA
            !temperature gradient at ZA from Bates profile
            dta = (tinf - ta) * s2 * ((this%re+zlb)/(this%re+za))**2.0
            tgn1(0)=dta
            tn1(0)=ta
            if (alt > zn1(mn1-1)) then
                z=alt
            else
                z=zn1(mn1-1)
            end if
            mn=mn1
            z1=zn1(0)
            z2=zn1(mn-1)
            t1=tn1(0)
            t2=tn1(mn-1)
            ! geopotental difference from z1
            zg = this%zeta (z, z1)
            zgdif = this%zeta(z2, z1)
            ! set up spline nodes
            do k=0,mn-1
                xs(k) = this%zeta(zn1(k), z1) / zgdif
                ys(k) = 1.0 / tn1(k)
            end do
            ! end node derivatives 
            yd1 = -tgn1(0) / (t1*t1) * zgdif
            yd2 = -tgn1(1) / (t2*t2) * zgdif * ((this%re+z2)/(this%re+z1))**2.0
            ! calculate spline coefficients 
            call this%spline (xs, ys, mn, yd1, yd2, y2out)
            x = zg / zgdif
            call this%splint (xs, ys, y2out, mn, x, y)
            ! temperature at altitude
            tz = 1.0 / y
            densu_temp = densu_temp * tz
        end if
        if (xm==0) then
            densu = densu_temp
            return
        end if
        ! calculate density above za
        glb = this%gsurf / (1.0 + zlb/this%re)**2.0
        gamma = xm * glb / (s2 * rgas * tinf)
        expl = exp(-s2 * gamma * zg2)
        if (expl > 50.0) then
            expl=50.0
        end if
        if (tt <= 0.0) then
            expl=50.0
        end if
        ! density at altitude
        densa = dlb * (tlb/tt)**((1.0+alpha+gamma)) * expl
        densu_temp=densa
        if (alt >= za) then
            densu = densu_temp
            return
        end if
        ! calculate density below za
        glb = this%gsurf / (1.0 + z1/this%re)**2.0
        gamm = xm * glb * zgdif / rgas

        ! integrate spline temperatures
        call this%splini (xs, ys, y2out, mn, x, yi)
        expl = gamm * yi
        if (expl > 50.0) then
            expl=50.0
        end if
        if (tz <= 0.0) then
            expl=50.0
        end if

        ! density at altitude
        densu_temp = densu_temp * (t1 / tz)**(1.0 + alpha) * exp(-expl)
        densu = densu_temp
        return
    end function densu

    ! -------------------------------------------------------------------
    ! ------------------------------- GLOBE7 ----------------------------
    ! -------------------------------------------------------------------
    !    3hr Magnetic activity functions
    !    Eq. A24d
    real(sp) pure function g0(this, a, p)

        implicit none

        class(Nrlmsise00_class),intent(in)      :: this
        real(sp),intent(in)                     :: a
        real(sp),dimension(0:149),intent(in)    :: p

        g0 = (a - 4.0 + (p(25) - 1.0) * (a - 4.0 + (exp(-sqrt(p(24)*p(24)) * (a - 4.0)) - 1.0) / sqrt(p(24)*p(24))))

        return

    end function g0

    !    Eq. A24c
    real(sp) pure function sumex(this, ex)

        implicit none

        class(Nrlmsise00_class),intent(in)  :: this
        real(sp),intent(in)                 :: ex

        sumex = (1.0 + (1.0 - ex**19.0) / (1.0 - ex) * ex**0.5)

        return
    end function sumex

    !    Eq. A24a
    real(sp) pure function sg0(this, ex, p, ap)

        implicit none

        class(Nrlmsise00_class),intent(in)      :: this
        real(sp),intent(in)                     :: ex
        real(sp),dimension(0:149),intent(in)    :: p
        real(sp),dimension(0:6),intent(in)      :: ap

        sg0 = (this%g0(ap(1),p) + (this%g0(ap(2),p)*ex + this%g0(ap(3),p)*ex*ex +   &
                    this%g0(ap(4),p)*ex**3.0 + (this%g0(ap(5),p)*ex**4.0 +      &
                    this%g0(ap(6),p)*ex**12.0)*(1.0-ex**8.0)/(1.0-ex)))         &
                / this%sumex(ex)
        return
    end function sg0

    !       CALCULATE G(L) FUNCTION 
    !       Upper Thermosphere Parameters
    real(sp) function globe7(this, p, input, flags)

        implicit none

        class(Nrlmsise00_class),intent(inout)   :: this
        real(sp),dimension(0:149),intent(inout) :: p
        type(nrlmsise_input),intent(inout)      :: input
        type(nrlmsise_flags),intent(inout)      :: flags

        real(sp)                :: t(0:14)
        integer                 :: i,j
        real(sp)                :: apd
        real(sp)                :: tloc
        real(sp)                :: c, s, c2, c4, s2
        real(sp),parameter      :: sr = 7.2722e-5
        real(sp),parameter      :: dgtr = 1.74533e-2
        real(sp),parameter      :: dr = 1.72142e-2
        real(sp),parameter      :: hr = 0.2618
        real(sp)                :: cd32, cd18, cd14, cd39
        real(sp)                :: df
        real(sp)                :: f1, f2
        real(sp)                :: tinf
        type(ap_array)          :: ap

        ! Aux variables 
        real(sp) p44, p45
        real(sp) t71, t72
        real(sp) t81, t82
        real(sp) exp1

        tloc=input%lst
        do j=0,13
            t(j)=0.0
        end do

        ! calculate legendre polynomials 
        c = sin(input%g_lat * dgtr)
        s = cos(input%g_lat * dgtr)
        c2 = c*c
        c4 = c2*c2
        s2 = s*s
        this%plg(0,1) = c
        this%plg(0,2) = 0.5*(3.0*c2 -1.0)
        this%plg(0,3) = 0.5*(5.0*c*c2-3.0*c)
        this%plg(0,4) = (35.0*c4 - 30.0*c2 + 3.0)/8.0
        this%plg(0,5) = (63.0*c2*c2*c - 70.0*c2*c + 15.0*c)/8.0
        this%plg(0,6) = (11.0*c*this%plg(0,5) - 5.0*this%plg(0,4))/6.0
    !      this%plg(0,7) = (13.0*c*this%plg(0,6) - 6.0*this%plg(0,5))/7.0 
        this%plg(1,1) = s
        this%plg(1,2) = 3.0*c*s
        this%plg(1,3) = 1.5*(5.0*c2-1.0)*s
        this%plg(1,4) = 2.5*(7.0*c2*c-3.0*c)*s
        this%plg(1,5) = 1.875*(21.0*c4 - 14.0*c2 +1.0)*s
        this%plg(1,6) = (11.0*c*this%plg(1,5)-6.0*this%plg(1,4))/5.0
    !      this%plg(1,7) = (13.0*c*this%plg(1,6)-7.0*this%plg(1,5))/6.0 
    !      this%plg(1,8) = (15.0*c*this%plg(1,7)-8.0*this%plg(1,6))/7.0 
        this%plg(2,2) = 3.0*s2
        this%plg(2,3) = 15.0*s2*c
        this%plg(2,4) = 7.5*(7.0*c2 -1.0)*s2
        this%plg(2,5) = 3.0*c*this%plg(2,4)-2.0*this%plg(2,3)
        this%plg(2,6) =(11.0*c*this%plg(2,5)-7.0*this%plg(2,4))/4.0
        this%plg(2,7) =(13.0*c*this%plg(2,6)-8.0*this%plg(2,5))/5.0
        this%plg(3,3) = 15.0*s2*s
        this%plg(3,4) = 105.0*s2*s*c 
        this%plg(3,5) =(9.0*c*this%plg(3,4)-7.*this%plg(3,3))/2.0
        this%plg(3,6) =(11.0*c*this%plg(3,5)-8.*this%plg(3,4))/3.0

        if (.not. (((flags%sw(7)==0) &
            .and. (flags%sw(8)==0)) &
            .and. (flags%sw(14)==0))) then
            this%stloc = sin(hr*tloc)
            this%ctloc = cos(hr*tloc)
            this%s2tloc = sin(2.0*hr*tloc)
            this%c2tloc = cos(2.0*hr*tloc)
            this%s3tloc = sin(3.0*hr*tloc)
            this%c3tloc = cos(3.0*hr*tloc)
        end if

        cd32 = cos(dr*(input%doy-p(31)))
        cd18 = cos(2.0*dr*(input%doy-p(17)))
        cd14 = cos(dr*(input%doy-p(13)))
        cd39 = cos(2.0*dr*(input%doy-p(38)))

        ! F10.7 EFFECT 
        df = input%f107 - input%f107A
        this%dfa = input%f107A - 150.0
        t(0) =  p(19)*df*(1.0+p(59)*this%dfa) + p(20)*df*df + p(21)*this%dfa + p(29)*this%dfa**2.0
        f1 = 1.0 + (p(47)*this%dfa +p(19)*df+p(20)*df*df)*flags%swc(1)
        f2 = 1.0 + (p(49)*this%dfa+p(19)*df+p(20)*df*df)*flags%swc(1)

        !  TIME INDEPENDENT 
        t(1) = (p(1)*this%plg(0,2)+ p(2)*this%plg(0,4)+p(22)*this%plg(0,6)) + &
              (p(14)*this%plg(0,2))*this%dfa*flags%swc(1) +p(26)*this%plg(0,1)

        !  SYMMETRICAL ANNUAL 
        t(2) = p(18)*cd32

        !  SYMMETRICAL SEMIANNUAL 
        t(3) = (p(15)+p(16)*this%plg(0,2))*cd18

        !  ASYMMETRICAL ANNUAL 
        t(4) =  f1*(p(9)*this%plg(0,1)+p(10)*this%plg(0,3))*cd14

        !  ASYMMETRICAL SEMIANNUAL 
        t(5) =    p(37)*this%plg(0,1)*cd39

            ! DIURNAL 
        if (flags%sw(7) /= 0) then

            t71 = (p(11)*this%plg(1,2))*cd14*flags%swc(5)
            t72 = (p(12)*this%plg(1,2))*cd14*flags%swc(5)
            t(6) = f2*((p(3)*this%plg(1,1) + p(4)*this%plg(1,3) + p(27)*this%plg(1,5) + t71) * &
                   this%ctloc + (p(6)*this%plg(1,1) + p(7)*this%plg(1,3) + p(28)*this%plg(1,5) &
                        + t72)*this%stloc)
        end if

        ! SEMIDIURNAL 
        if (flags%sw(8) /= 0) then

            t81 = (p(23)*this%plg(2,3)+p(35)*this%plg(2,5))*cd14*flags%swc(5)
            t82 = (p(33)*this%plg(2,3)+p(36)*this%plg(2,5))*cd14*flags%swc(5)
            t(7) = f2*((p(5)*this%plg(2,2)+ p(41)*this%plg(2,4) + t81)*this%c2tloc +(p(8)*this%plg(2,2) + p(42)*this%plg(2,4) + t82)*this%s2tloc)
        end if

        ! TERDIURNAL 
        if (flags%sw(14) /= 0) then
            t(13) = f2 * ((p(39)*this%plg(3,3)+(p(93)*this%plg(3,4)+p(46)*this%plg(3,6))*cd14*flags%swc(5))* this%s3tloc +(p(40)*this%plg(3,3)+(p(94)*this%plg(3,4)+p(48)*this%plg(3,6))*cd14*flags%swc(5))* this%c3tloc)
        end if

        ! magnetic activity based on daily ap 
        if (flags%sw(9) == -1) then
            ap = input%ap_a
            if (p(51) /= 0) then

                exp1 = exp(-10800.0*sqrt(p(51)*p(51))/(1.0+p(138)*(45.0-sqrt(input%g_lat*input%g_lat))))
                if (exp1 > 0.99999) then
                    exp1=0.99999
                end if
                if (p(24) < 1.e-4) then
                    p(24) = 1.e-4
                end if
                this%apt(0)=this%sg0(exp1,p,ap%a)
                !apt(1)=sg2(exp1,p,ap%a)
                !apt(2)=this%sg0(exp2,p,ap%a)
                !apt(3)=sg2(exp2,p,ap%a)
                
                if (flags%sw(9) /= 0) then
                    t(8) = this%apt(0)*(p(50)+p(96)*this%plg(0,2)+p(54)*this%plg(0,4)+ &
                            (p(125)*this%plg(0,1)+p(126)*this%plg(0,3)+p(127)*this%plg(0,5))*cd14*flags%swc(5)+ &
                            (p(128)*this%plg(1,1)+p(129)*this%plg(1,3)+p(130)*this%plg(1,5))*flags%swc(7)* &
                               cos(hr*(tloc-p(131))))
                end if
            end if
        else
            apd=input%ap-4.0
            p44=p(43)
            p45=p(44)
            if (p44 < 0.0) then
                p44 = 1.e-5
            end if
            this%apdf = apd + (p45-1.0)*(apd + (exp(-p44 * apd) - 1.0)/p44)
            if (flags%sw(9) /= 0) then
                t(8)=this%apdf*(p(32)+p(45)*this%plg(0,2)+p(34)*this%plg(0,4)+ &
                (p(100)*this%plg(0,1)+p(101)*this%plg(0,3)+p(102)*this%plg(0,5))*cd14*flags%swc(5)+ &
                (p(121)*this%plg(1,1)+p(122)*this%plg(1,3)+p(123)*this%plg(1,5))*flags%swc(7)* &
                        cos(hr*(tloc-p(124))))
            end if
        end if

        if ((flags%sw(10) /= 0) &
            .and. (input%g_long>-1000.0)) then

            ! longitudinal 
            if (flags%sw(11) /= 0) then
                t(10) = (1.0 + p(80)*this%dfa*flags%swc(1))* &
                        ((p(64)*this%plg(1,2)+p(65)*this%plg(1,4)+p(66)*this%plg(1,6)&
                        +p(103)*this%plg(1,1)+p(104)*this%plg(1,3)+p(105)*this%plg(1,5)&
                        +flags%swc(5)*(p(109)*this%plg(1,1)+p(110)*this%plg(1,3)+p(111)*this%plg(1,5))*cd14)* &
                        cos(dgtr*input%g_long) &
                        +(p(90)*this%plg(1,2)+p(91)*this%plg(1,4)+p(92)*this%plg(1,6)&
                        +p(106)*this%plg(1,1)+p(107)*this%plg(1,3)+p(108)*this%plg(1,5)&
                        +flags%swc(5)*(p(112)*this%plg(1,1)+p(113)*this%plg(1,3)+p(114)*this%plg(1,5))*cd14)* &
                        sin(dgtr*input%g_long))
            end if

            ! ut and mixed ut, longitude 
            if (flags%sw(12) /= 0) then
                t(11)=(1.0+p(95)*this%plg(0,1))*(1.0+p(81)*this%dfa*flags%swc(1))*&
                    (1.0+p(119)*this%plg(0,1)*flags%swc(5)*cd14)*&
                    ((p(68)*this%plg(0,1)+p(69)*this%plg(0,3)+p(70)*this%plg(0,5))*&
                    cos(sr*(input%sec-p(71))))
                t(11)= t(11) + flags%swc(11)*&
                    (p(76)*this%plg(2,3)+p(77)*this%plg(2,5)+p(78)*this%plg(2,7))*&
                    cos(sr*(input%sec-p(79))+2.0*dgtr*input%g_long)*(1.0+p(137)*this%dfa*flags%swc(1))
            end if

            ! ut, longitude magnetic activity 
            if (flags%sw(13) /= 0) then
                if (flags%sw(9) == -1) then
                    if (p(51) /= 0) then
                        t(12)=this%apt(0)*flags%swc(11)*(1.+p(132)*this%plg(0,1))*&
                            ((p(52)*this%plg(1,2)+p(98)*this%plg(1,4)+p(67)*this%plg(1,6))*&
                             cos(dgtr*(input%g_long-p(97))))&
                            +this%apt(0)*flags%swc(11)*flags%swc(5)*&
                            (p(133)*this%plg(1,1)+p(134)*this%plg(1,3)+p(135)*this%plg(1,5))*&
                            cd14*cos(dgtr*(input%g_long-p(136))) &
                            +this%apt(0)*flags%swc(12)* &
                            (p(55)*this%plg(0,1)+p(56)*this%plg(0,3)+p(57)*this%plg(0,5))*&
                            cos(sr*(input%sec-p(58)))
                    end if
                else
                    t(12) = this%apdf*flags%swc(11)*(1.0+p(120)*this%plg(0,1))*&
                        ((p(60)*this%plg(1,2)+p(61)*this%plg(1,4)+p(62)*this%plg(1,6))*&
                        cos(dgtr*(input%g_long-p(63))))&
                        +this%apdf*flags%swc(11)*flags%swc(5)* &
                        (p(115)*this%plg(1,1)+p(116)*this%plg(1,3)+p(117)*this%plg(1,5))* &
                        cd14*cos(dgtr*(input%g_long-p(118))) &
                        + this%apdf*flags%swc(12)* &
                        (p(83)*this%plg(0,1)+p(84)*this%plg(0,3)+p(85)*this%plg(0,5))* &
                        cos(sr*(input%sec-p(75)))
                end if
            end if
        end if

        ! parms not used: 82, 89, 99, 139-149 
        tinf = p(30)
        do i=0,13
            tinf = tinf + abs(flags%sw(i+1))*t(i)
        end do
        globe7 = tinf
        return

    end function globe7

    ! -------------------------------------------------------------------
    ! ------------------------------- GLOB7S ----------------------------
    ! -------------------------------------------------------------------
    !    VERSION OF GLOBE FOR LOWER ATMOSPHERE 10/26/99
    !
    real(sp) function glob7s(this, p, input, flags)

        implicit none

        class(Nrlmsise00_class),intent(inout)   :: this
        real(sp),dimension(0:99),intent(inout)  :: p
        type(nrlmsise_input),intent(inout)      :: input
        type(nrlmsise_flags),intent(inout)      :: flags

        real(sp),parameter          :: pset = 2.0
        real(sp),dimension(0:13)    :: t
        real(sp)                    :: tt
        real(sp)                    :: cd32, cd18, cd14, cd39
        integer                     :: i,j
        real(sp),parameter          :: dr=1.72142e-2
        real(sp),parameter          :: dgtr=1.74533e-2

        ! Auxiliary variable
        real(sp) t71, t72
        real(sp) t81, t82

        ! confirm parameter set
        if (p(99)==0) then
            p(99)=pset
        end if

        if (p(99) /= pset) then
            write (*,*) "Wrong parameter set for glob7s"
            glob7s = -1
            return
        end if

        do j=0,13
            t(j)=0.0
        end do

        cd32 = cos(dr*(input%doy-p(31)))
        cd18 = cos(2.0*dr*(input%doy-p(17)))
        cd14 = cos(dr*(input%doy-p(13)))
        cd39 = cos(2.0*dr*(input%doy-p(38)))

        ! F10.7 */
        t(0) = p(21)*this%dfa

        ! time independent
        t(1)=p(1)*this%plg(0,2) + p(2)*this%plg(0,4) + p(22)*this%plg(0,6) + p(26)*this%plg(0,1) + p(14)*this%plg(0,3) + p(59)*this%plg(0,5)

        ! SYMMETRICAL ANNUAL
        t(2)=(p(18)+p(47)*this%plg(0,2)+p(29)*this%plg(0,4))*cd32

        ! SYMMETRICAL SEMIANNUAL
        t(3)=(p(15)+p(16)*this%plg(0,2)+p(30)*this%plg(0,4))*cd18

        ! ASYMMETRICAL ANNUAL
        t(4)=(p(9)*this%plg(0,1)+p(10)*this%plg(0,3)+p(20)*this%plg(0,5))*cd14

        ! ASYMMETRICAL SEMIANNUAL
        t(5)=(p(37)*this%plg(0,1))*cd39

        ! DIURNAL
        if (flags%sw(7) /= 0) then
            t71 = p(11)*this%plg(1,2)*cd14*flags%swc(5)
            t72 = p(12)*this%plg(1,2)*cd14*flags%swc(5)
            t(6) = ((p(3)*this%plg(1,1) + p(4)*this%plg(1,3) + t71) * this%ctloc + (p(6)*this%plg(1,1) + p(7)*this%plg(1,3) + t72) * this%stloc)
        end if

        ! SEMIDIURNAL
        if (flags%sw(8) /= 0) then
            t81 = (p(23)*this%plg(2,3)+p(35)*this%plg(2,5))*cd14*flags%swc(5)
            t82 = (p(33)*this%plg(2,3)+p(36)*this%plg(2,5))*cd14*flags%swc(5)
            t(7) = ((p(5)*this%plg(2,2) + p(41)*this%plg(2,4) + t81) * this%c2tloc + (p(8)*this%plg(2,2) + p(42)*this%plg(2,4) + t82) * this%s2tloc)
        end if

        ! TERDIURNAL
        if (flags%sw(14) /= 0) then
            t(13) = p(39) * this%plg(3,3) * this%s3tloc + p(40) * this%plg(3,3) * this%c3tloc
        end if

        ! MAGNETIC ACTIVITY
        if (flags%sw(9) /= 0) then
            if (flags%sw(9) == 1) then
                t(8) = this%apdf * (p(32) + p(45) * this%plg(0,2) * flags%swc(2))
            end if
            if (flags%sw(9) == -1) then
                t(8)=(p(50)*this%apt(0) + p(96)*this%plg(0,2) * this%apt(0)*flags%swc(2))
            end if
        end if

        ! LONGITUDINAL
        if (.not. ((flags%sw(10) == 0) &
            .or. (flags%sw(11) == 0) &
            .or. (input%g_long <= -1000.0))) then
            t(10) = (1.0 + this%plg(0,1)*(p(80)*flags%swc(5)*cos(dr*(input%doy-p(81)))&
                    +p(85)*flags%swc(6)*cos(2.0*dr*(input%doy-p(86))))&
                +p(83)*flags%swc(3)*cos(dr*(input%doy-p(84)))&
                +p(87)*flags%swc(4)*cos(2.0*dr*(input%doy-p(88))))&
                *((p(64)*this%plg(1,2)+p(65)*this%plg(1,4)+p(66)*this%plg(1,6)&
                +p(74)*this%plg(1,1)+p(75)*this%plg(1,3)+p(76)*this%plg(1,5)&
                )*cos(dgtr*input%g_long)&
                +(p(90)*this%plg(1,2)+p(91)*this%plg(1,4)+p(92)*this%plg(1,6)&
                +p(77)*this%plg(1,1)+p(78)*this%plg(1,3)+p(79)*this%plg(1,5)&
                )*sin(dgtr*input%g_long))
        end if
        tt=0.0
        do i=0,13
            tt = tt + abs(flags%sw(i+1))*t(i)
        end do
        glob7s = tt
        return
    end function glob7s

    ! -------------------------------------------------------------------
    ! ------------------------------- GTD7 ------------------------------
    ! -------------------------------------------------------------------

    subroutine gtd7(this, input, flags, output)

        implicit none

        class(Nrlmsise00_class),intent(inout)   :: this
        type(nrlmsise_input),intent(inout)      :: input
        type(nrlmsise_flags),intent(inout)      :: flags
        type(nrlmsise_output),intent(inout)     :: output

        real(sp)                        :: xlat
        real(sp)                        :: xmm
        integer,parameter               :: mn3 = 5
        real(sp),dimension(0:mn3-1),parameter     :: zn3 = [32.50,20.0,15.0,10.0,0.0]
        integer,parameter               :: mn2 = 4
        real(sp),dimension(0:mn2-1),parameter     :: zn2 = [72.50,55.0,45.0,32.5]
        real(sp)                        :: altt
        real(sp),parameter              :: zmix = 62.5
        real(sp)                        :: tmp
        real(sp)                        :: dm28m
        real(sp)                        :: tz
        real(sp)                        :: dmc
        real(sp)                        :: dmr
        real(sp)                        :: dz28
        type(nrlmsise_output)           :: soutput
        integer                         :: i

        call this%tselec(flags)

        ! Latitude variation of gravity (none for sw(2)=0)
        xlat=input%g_lat
        if (flags%sw(2) == 0) then
            xlat=45.0
        end if
        call this%glatf(xlat, this%gsurf, this%re)

        xmm = pdm(2,4)

        ! THERMOSPHERE / MESOSPHERE (above zn2(0))
        if (input%alt>zn2(0)) then
            altt=input%alt
        else
            altt=zn2(0)
        end if

        tmp=input%alt
        input%alt=altt
        call this%gts7(input, flags, soutput)
        altt=input%alt
        input%alt=tmp
        if (flags%sw(0) /= 0) then  ! metric adjustment
            dm28m=this%dm28*1.0d6
        else
            dm28m=this%dm28
        end if
        output%t(0)=soutput%t(0)
        output%t(1)=soutput%t(1)
        if (input%alt>=zn2(0)) then
            do i=0,8
                output%d(i)=soutput%d(i)
            end do
            return
        end if

    !       LOWER MESOSPHERE/UPPER STRATOSPHERE (between zn3(0) and zn2(0))
    !          Temperature at nodes and gradients at end nodes
    !          Inverse temperature a linear function of spherical harmonics
        this%meso_tgn2(0)=this%meso_tgn1(1)
        this%meso_tn2(0)=this%meso_tn1(4)
        this%meso_tn2(1)=pma(0,0)*pavgm(0)/(1.0-flags%sw(20)*this%glob7s(pma(0,:), input, flags))
        this%meso_tn2(2)=pma(1,0)*pavgm(1)/(1.0-flags%sw(20)*this%glob7s(pma(1,:), input, flags))
        this%meso_tn2(3)=pma(2,0)*pavgm(2)/(1.0-flags%sw(20)*flags%sw(22)*this%glob7s(pma(2,:), input, flags))
        this%meso_tgn2(1)=pavgm(8)*pma(9,0)*(1.0+flags%sw(20)*flags%sw(22)*this%glob7s(pma(9,:), input, flags))*this%meso_tn2(3)*this%meso_tn2(3)/((pma(2,0)*pavgm(2))**2.0)
        this%meso_tn3(0)=this%meso_tn2(3)

        if (input%alt <= zn3(0)) then
    !       LOWER STRATOSPHERE AND TROPOSPHERE (below zn3(0))
    !          Temperature at nodes and gradients at end nodes
    !          Inverse temperature a linear function of spherical harmonics
            this%meso_tgn3(0)=this%meso_tgn2(1)
            this%meso_tn3(1)=pma(3,0)*pavgm(3)/(1.0-flags%sw(22)*this%glob7s(pma(3,:), input, flags))
            this%meso_tn3(2)=pma(4,0)*pavgm(4)/(1.0-flags%sw(22)*this%glob7s(pma(4,:), input, flags))
            this%meso_tn3(3)=pma(5,0)*pavgm(5)/(1.0-flags%sw(22)*this%glob7s(pma(5,:), input, flags))
            this%meso_tn3(4)=pma(6,0)*pavgm(6)/(1.0-flags%sw(22)*this%glob7s(pma(6,:), input, flags))
            this%meso_tgn3(1)=pma(7,0)*pavgm(7)*(1.0+flags%sw(22)*this%glob7s(pma(7,:), input, flags)) *this%meso_tn3(4)*this%meso_tn3(4)/((pma(6,0)*pavgm(6))**2.0)
        end if

        ! LINEAR TRANSITION TO FULL MIXING BELOW zn2(0)

        dmc=0
        if (input%alt>zmix) then
            dmc = 1.0 - (zn2(0)-input%alt)/(zn2(0) - zmix)
        end if
        dz28=soutput%d(2)
        
        !**** N2 density ****
        dmr=soutput%d(2) / dm28m - 1.0
        output%d(2)=this%densm(input%alt,dm28m,xmm, tz, mn3, zn3, this%meso_tn3, this%meso_tgn3, &
                            mn2, zn2, this%meso_tn2, this%meso_tgn2)
        output%d(2)=output%d(2) * (1.0 + dmr*dmc)

        !**** HE density ****
        dmr = soutput%d(0) / (dz28 * pdm(0,1)) - 1.0
        output%d(0) = output%d(2) * pdm(0,1) * (1.0 + dmr*dmc)

        !**** O density ****
        output%d(1) = 0
        output%d(8) = 0

        !**** O2 density ****
        dmr = soutput%d(3) / (dz28 * pdm(3,1)) - 1.0
        output%d(3) = output%d(2) * pdm(3,1) * (1.0 + dmr*dmc)

        !**** AR density ***
        dmr = soutput%d(4) / (dz28 * pdm(4,1)) - 1.0
        output%d(4) = output%d(2) * pdm(4,1) * (1.0 + dmr*dmc)

        !**** Hydrogen density ****
        output%d(6) = 0

        !**** Atomic nitrogen density ****
        output%d(7) = 0

        !**** Total mass density ****
        output%d(5) = 1.66e-24 * (4.0 * output%d(0) + 16.0 * output%d(1) + 28.0 &
                    * output%d(2) + 32.0 * output%d(3) + 40.0 * output%d(4) &
                    + output%d(6) + 14.0 * output%d(7))

        if (flags%sw(0) /= 0) then
            output%d(5)=output%d(5)/1000.0
        end if

        !**** temperature at altitude ****
        this%dd = this%densm(input%alt, 1.0, 0.0, tz, mn3, zn3,this%meso_tn3, this%meso_tgn3, mn2, zn2, this%meso_tn2, this%meso_tgn2)
        output%t(1)=tz

    end subroutine gtd7

    ! -------------------------------------------------------------------
    ! ------------------------------- GTD7D -----------------------------
    ! -------------------------------------------------------------------
    subroutine gtd7d(this, input, flags, output)

        implicit none

        class(Nrlmsise00_class),intent(inout)   :: this
        type(nrlmsise_input),intent(inout)      :: input
        type(nrlmsise_flags),intent(inout)      :: flags
        type(nrlmsise_output),intent(inout)     :: output

        call this%gtd7(input, flags, output)

        output%d(5) = 1.66e-24 * (4.0 * output%d(0) + 16.0 * output%d(1) &
                        + 28.0 * output%d(2) + 32.0 * output%d(3) + 40.0 &
                        * output%d(4) + output%d(6) + 14.0 * output%d(7) &
                        + 16.0 * output%d(8))

        if (flags%sw(0) /= 0) then
            output%d(5)=output%d(5)/1000.0
        end if

    end subroutine gtd7d

    ! -------------------------------------------------------------------
    ! -------------------------------- GHP7 -----------------------------
    ! -------------------------------------------------------------------
    subroutine ghp7(this, input, flags, output, press)

        implicit none

        class(Nrlmsise00_class),intent(inout)   :: this
        type(nrlmsise_input),intent(inout)      :: input
        type(nrlmsise_flags),intent(inout)      :: flags
        type(nrlmsise_output),intent(inout)     :: output
        real(sp),intent(in)                     :: press

        real(sp)                                :: bm = 1.3806e-19
        real(sp)                                :: rgas = 831.4
        real(sp)                                :: test = 0.00043
        real(sp)                                :: ltest = 12.0
        real(sp)                                :: pl, p
        real(sp)                                :: zi
        real(sp)                                :: z
        real(sp)                                :: cl, cl2
        real(sp)                                :: ca, cd
        real(sp)                                :: xn, xm, diff
        real(sp)                                :: g, sh
        integer                                 :: l

        pl = log10(press)

        if (pl >= -5.0) then
            if (pl>2.5) then
                zi = 18.06 * (3.0 - pl)
            else if ((pl>0.075) .and. (pl<=2.5)) then
                zi = 14.98 * (3.08 - pl)
            else if ((pl>-1.0) .and. (pl<=0.075)) then
                zi = 17.8 * (2.72 - pl)
            else if ((pl>-2.0) .and. (pl<=-1.0)) then
                zi = 14.28 * (3.64 - pl)
            else if ((pl>-4.0) .and. (pl<=-2.0)) then
                zi = 12.72 * (4.32 -pl)
            else
                zi = 25.3 * (0.11 - pl)
            end if
            cl = input%g_lat/90.0
            cl2 = cl*cl
            if (input%doy<182.0) then
                cd = (1.0 - dble (input%doy)) / 91.25
            else 
                cd = (dble(input%doy)) / 91.25 - 3.0
            end if
            ca = 0
            if ((pl > -1.11) .and. (pl<=-0.23)) then
                ca = 1.0
            end if
            if (pl > -0.23) then
                ca = (2.79 - pl) / (2.79 + 0.23)
            end if
            if ((pl <= -1.11) .and. (pl>-3.0)) then
                ca = (-2.93 - pl)/(-2.93 + 1.11)
            end if
            z = zi - 4.87 * cl * cd * ca - 1.64 * cl2 * ca + 0.31 * ca * cl
        else
            z = 22.0 * (pl + 4.0)**2.0 + 110.0
        end if
        ! iteration  loop */
        l = 0
        do
            l = l + 1
            input%alt = z
            call this%gtd7(input, flags, output)
            z = input%alt
            xn = output%d(0) + output%d(1) + output%d(2) + output%d(3) + output%d(4) + output%d(6) + output%d(7)
            p = bm * xn * output%t(1)
            if (flags%sw(0) /= 0) then
                p = p*1.e-6
            end if
            diff = pl - log10(p)
            if (sqrt(diff*diff)<test) then
                return
            end if
            if (l==ltest) then
                write(*,*) "ERROR: ghp7 not converging for press, diff:" ,press, diff
                return
            end if
            xm = output%d(5) / xn / 1.66e-24
            if (flags%sw(0) /= 0) then
                xm = xm * 1.3
            end if
            g = this%gsurf / ((1.0 + z/this%re)**2.0)
            sh = rgas * output%t(1) / (xm * g)

            ! new altitude estimate using scale height
            if (l <  6) then
                z = z - sh * diff * 2.302
            else
                z = z - sh * diff
            end if
        end do
    end subroutine ghp7

    ! -------------------------------------------------------------------
    ! ------------------------------- GTS7 ------------------------------
    ! -------------------------------------------------------------------
    !     Thermospheric portion of NRLMSISE-00
    !     See GTD7 for more extensive comments
    !      alt > 72.5 km!
    !
    subroutine gts7(this, input, flags, output)

        implicit none

        class(Nrlmsise00_class),intent(inout)   :: this
        type(nrlmsise_input),intent(inout)      :: input
        type(nrlmsise_flags),intent(inout)      :: flags
        type(nrlmsise_output),intent(inout)     :: output

        real(sp)                            :: za
        integer                             :: i, j
        real(sp)                            :: ddum, z
        real(sp),dimension(0:4)             :: zn1 = [120.0, 110.0, 100.0, 90.0, 72.50]
        real(sp)                            :: tinf
        integer,parameter                   :: mn1 = 5
        real(sp)                            :: g0
        real(sp)                            :: tlb
        real(sp)                            :: s
        real(sp)                            :: db01, db04, db14, db16, db28, db32, db40
        real(sp)                            :: zh28, zh04, zh16, zh32, zh40, zh01, zh14
        real(sp)                            :: zhm28, zhm04, zhm16, zhm32, zhm40, zhm01, zhm14
        real(sp)                            :: xmd
        real(sp)                            :: b28, b04, b16, b32, b40, b01, b14
        real(sp)                            :: tz
        real(sp)                            :: g28, g4, g16, g32, g40, g1, g14
        real(sp)                            :: zhf, xmm
        real(sp)                            :: zc04, zc16, zc32, zc40, zc01, zc14
        real(sp)                            :: hc04, hc16, hc32, hc40, hc01, hc14
        real(sp)                            :: hcc16, hcc32, hcc01, hcc14
        real(sp)                            :: zcc16, zcc32, zcc01, zcc14
        real(sp)                            :: rc16, rc32, rc01, rc14
        real(sp)                            :: rl
        real(sp)                            :: g16h, db16h, tho, zsht, zmho, zsho
        real(sp),parameter                  :: dgtr = 1.74533e-2
        real(sp),parameter                  :: dr = 1.72142e-2
        real(sp),dimension(0:8),parameter   :: alpha = [-0.380, 0.0, 0.0, 0.0, 0.170, 0.0, -0.380, 0.0, 0.0]
        real(sp),dimension(0:7),parameter   :: altl  = [200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0]
        real(sp)                            :: dd
        real(sp)                            :: hc216, hcc232

        za = pdl(1,15)
        zn1(0) = za
        do j=0,8
            output%d(j)=0.0
        end do

        ! TINF VARIATIONS NOT IMPORTANT BELOW ZA OR ZN1(1)
        if (input%alt>zn1(0)) then
            tinf = ptm(0)*pt(0) * &
                (1.0+flags%sw(16)*this%globe7(pt,input,flags))
        else
            tinf = ptm(0)*pt(0)
        end if
        output%t(0)=tinf

        !  GRADIENT VARIATIONS NOT IMPORTANT BELOW ZN1(5)
        if (input%alt>zn1(4)) then
            g0 = ptm(3)*ps(0) * &
                (1.0+flags%sw(19)*this%globe7(ps,input,flags))
        else
            g0 = ptm(3)*ps(0)
        end if
        tlb = ptm(1) * (1.0 + flags%sw(17)*this%globe7(pd(3,:),input,flags))*pd(3,0)
        s = g0 / (tinf - tlb)

    !      Lower thermosphere temp variations not significant for
    !      density above 300 km
        if (input%alt < 300.0) then
            this%meso_tn1(1)=ptm(6)*ptl(0,0)/(1.0-flags%sw(18)*this%glob7s(ptl(0,:), input, flags))
            this%meso_tn1(2)=ptm(2)*ptl(1,0)/(1.0-flags%sw(18)*this%glob7s(ptl(1,:), input, flags))
            this%meso_tn1(3)=ptm(7)*ptl(2,0)/(1.0-flags%sw(18)*this%glob7s(ptl(2,:), input, flags))
            this%meso_tn1(4)=ptm(4)*ptl(3,0)/(1.0-flags%sw(18)*flags%sw(20)*this%glob7s(ptl(3,:), input, flags))
            this%meso_tgn1(1)=ptm(8)*pma(8,0)*(1.0+flags%sw(18)*flags%sw(20)*this%glob7s(pma(8,:), input, flags))*this%meso_tn1(4)*this%meso_tn1(4)/((ptm(4)*ptl(3,0))**2.0)
        else
            this%meso_tn1(1)=ptm(6)*ptl(0,0)
            this%meso_tn1(2)=ptm(2)*ptl(1,0)
            this%meso_tn1(3)=ptm(7)*ptl(2,0)
            this%meso_tn1(4)=ptm(4)*ptl(3,0)
            this%meso_tgn1(1)=ptm(8)*pma(8,0)*this%meso_tn1(4)*this%meso_tn1(4)/((ptm(4)*ptl(3,0))**2.0)
        end if

        ! N2 variation factor at Zlb
        g28=flags%sw(21)*this%globe7(pd(2,:), input, flags)

        ! VARIATION OF TURBOPAUSE HEIGHT
        zhf=pdl(1,24)*(1.0+flags%sw(5)*pdl(0,24)*sin(dgtr*input%g_lat)*cos(dr*(input%doy-pt(13))))
        output%t(0)=tinf
        xmm = pdm(2,4)
        z = input%alt


        !**** N2 DENSITY ****
        ! Diffusive density at Zlb
        db28 = pdm(2,0)*exp(g28)*pd(2,0)
        ! Diffusive density at Alt */
        output%d(2)=this%densu(z,db28,tinf,tlb,28.0,alpha(2),output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
        dd=output%d(2)
        ! Turbopause
        zh28=pdm(2,2)*zhf
        zhm28=pdm(2,3)*pdl(1,5)
        xmd=28.0-xmm
        ! Mixed density at Zlb
        b28=this%densu(zh28,db28,tinf,tlb,xmd,(alpha(2)-1.0),tz,ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
        if ((flags%sw(15) /= 0) .and. (z<=altl(2))) then
            !  Mixed density at Alt
            this%dm28=this%densu(z,b28,tinf,tlb,xmm,alpha(2),tz,ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
            !  Net density at Alt
            output%d(2)=this%dnet(output%d(2),this%dm28,zhm28,xmm,28.0)
        end if


        !**** HE DENSITY ****
        !   Density variation factor at Zlb
        g4 = flags%sw(21)*this%globe7(pd(0,:), input, flags)
        !  Diffusive density at Zlb
        db04 = pdm(0,0)*exp(g4)*pd(0,0)
        ! Diffusive density at Alt
        output%d(0)=this%densu(z,db04,tinf,tlb,4.0,alpha(0),output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
        dd=output%d(0)
        if ((flags%sw(15) /= 0) .and. (z<altl(0))) then
            !  Turbopause
            zh04=pdm(0,2)
            !  Mixed density at Zlb
            b04=this%densu(zh04,db04,tinf,tlb,4.0-xmm,alpha(0)-1.0,output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
            !  Mixed density at Alt
            this%dm04=this%densu(z,b04,tinf,tlb,xmm,0.0,output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
            zhm04=zhm28
            !  Net density at Alt
            output%d(0)=this%dnet(output%d(0),this%dm04,zhm04,xmm,4.0)
            !  Correction to specified mixing ratio at ground
            rl=log(b28*pdm(0,1)/b04)
            zc04=pdm(0,4)*pdl(1,0)
            hc04=pdm(0,5)*pdl(1,1)
            !  Net density corrected at Alt
            output%d(0)=output%d(0)*this%ccor(z,rl,hc04,zc04)
        end if


        !**** O DENSITY ****

        !  Density variation factor at Zlb
        g16 = flags%sw(21)*this%globe7(pd(1,:),input,flags)
        !  Diffusive density at Zlb
        db16 = pdm(1,0)*exp(g16)*pd(1,0)
        ! Diffusive density at Alt
        output%d(1)=this%densu(z,db16,tinf,tlb,16.0,alpha(1),output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
        dd=output%d(1)
        if ((flags%sw(15) /= 0) .and. (z<=altl(1))) then
            !   Turbopause
            zh16=pdm(1,2)
            !  Mixed density at Zlb
            b16=this%densu(zh16,db16,tinf,tlb,16.0-xmm,(alpha(1)-1.0),output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
            !  Mixed density at Alt
            this%dm16=this%densu(z,b16,tinf,tlb,xmm,0.0,output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
            zhm16=zhm28
            !  Net density at Alt
            output%d(1)=this%dnet(output%d(1),this%dm16,zhm16,xmm,16.0)
            rl=pdm(1,1)*pdl(1,16)*(1.0+flags%sw(1)*pdl(0,23)*(input%f107A-150.0))
            hc16=pdm(1,5)*pdl(1,3)
            zc16=pdm(1,4)*pdl(1,2)
            hc216=pdm(1,5)*pdl(1,4)
            output%d(1)=output%d(1)*this%ccor2(z,rl,hc16,zc16,hc216)
            !   Chemistry correction
            hcc16=pdm(1,7)*pdl(1,13)
            zcc16=pdm(1,6)*pdl(1,12)
            rc16=pdm(1,3)*pdl(1,14)
            !  Net density corrected at Alt
            output%d(1)=output%d(1)*this%ccor(z,rc16,hcc16,zcc16)
        end if


        !**** O2 DENSITY ****
        ! Density variation factor at Zlb
        g32= flags%sw(21)*this%globe7(pd(4,:), input, flags)
        ! Diffusive density at Zlb
        db32 = pdm(3,0)*exp(g32)*pd(4,0)
        ! Diffusive density at Alt
        output%d(3)=this%densu(z,db32,tinf,tlb,32.0,alpha(3),output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
        dd=output%d(3)
        if (flags%sw(15) /= 0) then
            if (z<=altl(3)) then
                !   Turbopause
                zh32=pdm(3,2)
                !  Mixed density at Zlb
                b32=this%densu(zh32,db32,tinf,tlb,32.0-xmm,alpha(3)-1.0,output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
                !  Mixed density at Alt
                this%dm32=this%densu(z,b32,tinf,tlb,xmm,0.0,output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
                zhm32=zhm28
                !  Net density at Alt
                output%d(3)=this%dnet(output%d(3),this%dm32,zhm32,xmm,32.0)
                !   Correction to specified mixing ratio at ground
                rl=log(b28*pdm(3,1)/b32)
                hc32=pdm(3,5)*pdl(1,7)
                zc32=pdm(3,4)*pdl(1,6)
                output%d(3)=output%d(3)*this%ccor(z,rl,hc32,zc32)
            end if
            !  Correction for general departure from diffusive equilibrium above Zlb
            hcc32=pdm(3,7)*pdl(1,22)
            hcc232=pdm(3,7)*pdl(0,22)
            zcc32=pdm(3,6)*pdl(1,21)
            rc32=pdm(3,3)*pdl(1,23)*(1.0+flags%sw(1)*pdl(0,23)*(input%f107A-150.0))
            !  Net density corrected at Alt */
            output%d(3)=output%d(3)*this%ccor2(z,rc32,hcc32,zcc32,hcc232)
        end if


        !**** AR DENSITY ****
        ! Density variation factor at Zlb
        g40 = flags%sw(21)*this%globe7(pd(5,:),input,flags)
        ! Diffusive density at Zlb
        db40 = pdm(4,0)*exp(g40)*pd(5,0)
        ! Diffusive density at Alt
        output%d(4)=this%densu(z,db40,tinf,tlb,40.0,alpha(4),output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
        dd=output%d(4)
        if ((flags%sw(15) /= 0) .and. (z<=altl(4))) then
            !   Turbopause
            zh40=pdm(4,2)
            !  Mixed density at Zlb
            b40=this%densu(zh40,db40,tinf,tlb,40.0-xmm,alpha(4)-1.0,output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
            !  Mixed density at Alt
            this%dm40=this%densu(z,b40,tinf,tlb,xmm,0.0,output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
            zhm40=zhm28
            !  Net density at Alt
            output%d(4)=this%dnet(output%d(4),this%dm40,zhm40,xmm,40.0)
            !   Correction to specified mixing ratio at ground
            rl=log(b28*pdm(4,1)/b40)
            hc40=pdm(4,5)*pdl(1,9)
            zc40=pdm(4,4)*pdl(1,8)
            !  Net density corrected at Alt */
            output%d(4)=output%d(4)*this%ccor(z,rl,hc40,zc40)
        end if


        !**** HYDROGEN DENSITY ****

        ! Density variation factor at Zlb
        g1 = flags%sw(21)*this%globe7(pd(6,:), input, flags)
        ! Diffusive density at Zlb
        db01 = pdm(5,0)*exp(g1)*pd(6,0)
        ! Diffusive density at Alt
        output%d(6)=this%densu(z,db01,tinf,tlb,1.0,alpha(6),output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
        dd=output%d(6)
        if ((flags%sw(15) /= 0) .and. (z<=altl(6))) then
            !   Turbopause
            zh01=pdm(5,2)
            !  Mixed density at Zlb
            b01=this%densu(zh01,db01,tinf,tlb,1.0-xmm,alpha(6)-1.0,output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
            !  Mixed density at Alt
            this%dm01=this%densu(z,b01,tinf,tlb,xmm,0.0,output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
            zhm01=zhm28
            !  Net density at Alt
            output%d(6)=this%dnet(output%d(6),this%dm01,zhm01,xmm,1.0)
            !   Correction to specified mixing ratio at ground
            rl=log(b28*pdm(5,1)*sqrt(pdl(1,17)*pdl(1,17))/b01)
            hc01=pdm(5,5)*pdl(1,11)
            zc01=pdm(5,4)*pdl(1,10)
            output%d(6)=output%d(6)*this%ccor(z,rl,hc01,zc01)
            !   Chemistry correction */
            hcc01=pdm(5,7)*pdl(1,19)
            zcc01=pdm(5,6)*pdl(1,18)
            rc01=pdm(5,3)*pdl(1,20)
            !  Net density corrected at Alt
            output%d(6)=output%d(6)*this%ccor(z,rc01,hcc01,zcc01)
        end if


        !**** ATOMIC NITROGEN DENSITY ****/

        !   Density variation factor at Zlb
        g14 = flags%sw(21)*this%globe7(pd(7,:),input,flags)
        ! Diffusive density at Zlb
        db14 = pdm(6,0)*exp(g14)*pd(7,0)
        ! Diffusive density at Alt
        output%d(7)=this%densu(z,db14,tinf,tlb,14.0,alpha(7),output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
        dd=output%d(7)
        if ((flags%sw(15) /= 0) .and. (z<=altl(7))) then
            !   Turbopause
            zh14=pdm(6,2)
            !  Mixed density at Zlb
            b14=this%densu(zh14,db14,tinf,tlb,14.0-xmm,alpha(7)-1.0,output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
            !  Mixed density at Alt
            this%dm14=this%densu(z,b14,tinf,tlb,xmm,0.0,output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
            zhm14=zhm28
            !  Net density at Alt
            output%d(7)=this%dnet(output%d(7),this%dm14,zhm14,xmm,14.0)
            !   Correction to specified mixing ratio at ground */
            rl=log(b28*pdm(6,1)*sqrt(pdl(0,2)*pdl(0,2))/b14)
            hc14=pdm(6,5)*pdl(0,1)
            zc14=pdm(6,4)*pdl(0,0)
            output%d(7)=output%d(7)*this%ccor(z,rl,hc14,zc14)
            !   Chemistry correction */
            hcc14=pdm(6,7)*pdl(0,4)
            zcc14=pdm(6,6)*pdl(0,3)
            rc14=pdm(6,3)*pdl(0,5)
            !  Net density corrected at Alt */
            output%d(7)=output%d(7)*this%ccor(z,rc14,hcc14,zcc14)
        end if


        ! **** Anomalous OXYGEN DENSITY ****/
        g16h = flags%sw(21)*this%globe7(pd(8,:),input,flags)
        db16h = pdm(7,0)*exp(g16h)*pd(8,0)
        tho = pdm(7,9)*pdl(0,6)
        dd=this%densu(z,db16h,tho,tho,16.0,alpha(8),output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
        zsht=pdm(7,5)
        zmho=pdm(7,4)
        zsho=this%scalh(zmho,16.0,tho)
        output%d(8)=dd*exp(-zsht/zsho*(exp(-(z-zmho)/zsht)-1.0))


        ! total mass density
        output%d(5) = 1.66e-24*(4.0*output%d(0)+16.0*output%d(1)+28.0*output%d(2) &
                        +32.0*output%d(3)+40.0*output%d(4)+ output%d(6)+14.0*output%d(7))


        ! temperature
        z = sqrt(input%alt*input%alt)
        ddum = this%densu(z,1.0, tinf, tlb, 0.0, 0.0,output%t(1),ptm(5),s,mn1,zn1,this%meso_tn1,this%meso_tgn1)
        if (flags%sw(0) /= 0) then
            do i=0,8
                output%d(i)=output%d(i)*1.6
            end do
            output%d(5)=output%d(5)/1000.0
        end if

    end subroutine gts7


end module nrlmsise00Class
