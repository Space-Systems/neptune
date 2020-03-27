!===============================================================================================
!
!> @anchor      hwm07Class
!!
!> @brief       Horizontal Wind Model
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>ChK:  04.01.2018 (Added to class)</li>
!!              </ul>
!!
!> @details     This module contains parameters, subroutines and functions required for Earth's
!!              geopotential modeling, specifically in the context of numerical
!!              integration of a satellite trajectory.
!!
!> @copyright   A. D. Richmond (NCAR-HAO), Adapted for HWM07 by Douglas P. Drob, Used for NEPTUNE by Vitali Braun and adapted into the class structure by Christopher Kebschull
!!
!> @anchor      hwm07Class
!!
!!------------------------------------------------------------------------------------------------
module hwm07Class

    implicit none

    type vshengine_t
        integer                    :: mmax = 11
        integer                    :: nmax = 11
        real(8)                    :: p0(0:11,0:11),p1(0:11,0:11)
        real(8)                    :: sf(0:11,0:11)
        real(8)                    :: a(0:11,0:11),b(0:11,0:11)
        real(8)                    :: c(0:11,0:11),d(0:11,0:11)
    end type vshengine_t

    type :: fn_map                                                              ! Structure of VSH basis function map
      integer(4) :: realimag
      integer(4) :: irrotational
      integer(4) :: M
      integer(4) :: N
    end type fn_map

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !                           Module vsh_basis_module
    !
    ! Description: This is a common data module for model definition.  These
    !  parameters set by the first calling the subroutine loaddwm().
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    type vsh_basis_t

        type (fn_map), allocatable :: fn_map1(:)                                ! Working array of VSH map structure

        real(8), allocatable :: A(:,:)                                          !Derivative of the Associated Legendre Functions
        real(8), allocatable :: B(:,:)                                          !m/sin(theta) times the Associated Legendre Functions
        real(8), allocatable :: anm(:,:)
        real(8), allocatable :: bnm(:,:)
        real(8), allocatable :: cm(:)
        real(8), allocatable :: cn(:)
        real(8), allocatable :: e0n(:)
        real(8), allocatable :: fnm(:,:)
        real(8), allocatable :: m_arr(:)
        real(8), allocatable :: n_arr(:)
        real(8), allocatable :: cosmz(:)
        real(8), allocatable :: sinmz(:)

    end type vsh_basis_t

    !=========================================================================================
    !
    !
    !  Code used in HWM07 for computing quasi-dipole magnetic 
    !  coordinates. Modified from apex.zip code package obtained 
    !  from the NSF-CEDAR database. 
    !  See readme.txt file for detailed release notes.
    !
    !  AUTHOR
    !    A. D. Richmond (NCAR-HAO)
    !    Adapted for HWM07 by Douglas P. Drob
    !
    ! Point of Contact
    !    msishwmhelp@nrl.navy.mil
    !
    !  DATE
    !    19 August 2008
    !
    !  REFERENCE
    !    Richmond, A. D. (1995), Ionospheric electrodynamics using 
    !    Magnetic Apex Coordinates, J. Geomag. & Geoelec., 47, 
    !    191212.
    !
    ! ===============================================================
    type apexcord_t
        integer(4)               :: kgma
        real(4)                  :: glalmn,glalmx
        integer(4)               :: nla,nlo,nal
        integer(4)               :: lbx,lby,lbz,lbv
        integer(4)               :: lla,llo,lal

        integer(4)               :: lwk
        real(4),allocatable      :: wk(:)

        real(4)                  :: colat
        real(4)                  :: elon
        real(4)                  :: vp
        real(4)                  :: ctp
        real(4)                  :: stp
          
        real(4)                  :: rtod
        real(4)                  :: dtor
        real(4)                  :: re
        real(4)                  :: req
        real(4)                  :: msgu
        real(4)                  :: pola

        integer(4)               :: io
        integer(4)               :: jo
        integer(4)               :: ko

        logical                  :: loaddata

    end type apexcord_t

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !                           Module DWM
    !
    ! Description: This is a common data module for model definition.  These
    !  parameters are set by the first calling the subroutine loaddwm().
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type dwm_t

        integer(4)             :: nterm          ! Number of terms in the model
        integer(4)             :: lmax           ! Max latitudinal degree
        integer(4)             :: mmax           ! Max order of MLT var.

        integer(4),allocatable :: termarr(:,:)   ! 3 x nterm index of coupled terms
        real(4),allocatable    :: coeff(:)       ! Model coefficients
        real(4),allocatable    :: vsh_terms(:,:) ! VSH basis values
        integer(4)             :: nvshfn         ! # of VSH basis functions
        real(4)                :: twidth         ! Transition width of high-lat mask
        real(4),allocatable    :: termval(:,:)   ! Term values to which coefficients are applied

        ! Initialization flags and information

        logical                 :: modelinit
        character(128)          :: defaultdata

    end type dwm_t

    ! From former vshengine module
    type(vshengine_t)       :: vshdata
    ! From former apexcord module
    type(apexcord_t)        :: apexdata
    !
    type(dwm_t)             :: dwmdata

    type(vsh_basis_t)       :: vshbasisdata

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !                           Global Data Module NEWmodel
    !
    ! Description: This is a common data module for model definition.  These 
    !  parameters set by the first calling the subroutine loadmodel().
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    type Hwm07_class

        integer                 :: nbf              ! Count of basis terms per model level   
        integer                 :: maxs             ! s seasonal
        integer                 :: maxm             ! m stationary
        integer                 :: maxl             ! l migrating
        integer                 :: maxn             ! n latitude
        
        integer                 :: p                ! B-splines order, p=4 cubic, p=3 quadratic 
        integer                 :: nlev             ! e.g. Number of B-spline nodes
        integer                 :: nnode            ! nlev + p

        real(8)                 :: alttns           ! Transition 1
        real(8)                 :: altsym           ! Transition 2
        real(8)                 :: altiso           ! Constant Limit

        integer,allocatable     :: nb(:)            ! total number of basis functions @ level
        integer,allocatable     :: order(:,:)       ! spectral content @ level
        real(8),allocatable     :: vnode(:)         ! Vertical Altitude Nodes
        real(8),allocatable     :: mparm(:,:)       ! Model Parameters
        
        ! Global store for quasi-static model space parameters
        ! These will change internally depending on the input parameters
        
        real(8),allocatable     :: gfs(:,:),gfm(:,:),gfl(:,:)
        real(8),allocatable     :: gvbar(:,:),gwbar(:,:)
        real(8),allocatable     :: gbz(:,:),gbm(:,:)
        
        real(8),allocatable     :: gzwght(:)
        integer                 :: glev
         
        ! Miscellaneous flags and indicies
        
        integer                 :: maxo
        integer                 :: cseason = 0
        integer                 :: cwave = 0
        integer                 :: ctide = 0
        
        logical                 :: content(5)                                   ! Season/Waves/Tides
        logical                 :: component(0:1)                               ! Compute zonal/meridional

        real(8)                 :: last(5)
        
        ! Initialization flags and information
        
        logical                 :: modelinit
        logical                 :: reset
        character(128)          :: defaultdata

    contains

        procedure :: destroy
        procedure :: hwm07
        procedure :: hwmqt
        procedure :: HWMupdate
        procedure :: loadmodel
        procedure :: vertwght
        procedure :: vshengineinit
        procedure :: vshbasis
        procedure :: apxrda
        procedure :: apex
        procedure :: intrp
        procedure :: adpl
        procedure :: grapxyzv
        procedure :: gradxyzv
        procedure :: gradlpv
        procedure :: trilin
        procedure :: basvec
        procedure :: dwm07b_hwm_interface
        procedure :: dwm07b
        procedure :: loaddwm
        procedure :: vsh_basis_init
        procedure :: vsh_basis
        procedure :: dwm_kpspl3_calc
        procedure :: dwm_latwgt2
        procedure :: bspline_calc
        procedure :: pershift
        procedure :: ap_to_kp
        procedure :: dwm_altwgt
        procedure :: gd2qd
        
    end type Hwm07_class

        ! Constructor
    interface Hwm07_class
        module procedure constructor
    end interface Hwm07_class

contains

    ! ====================================================================
    !!
    !>  @brief      Constructor that should be used to initialize variables.
    !!
    !>  @author     Christopher Kebschull
    !>  @date       <ul>
    !!                  <li>ChK: 24.12.2017 (initial implementation)</li>
    !!              </ul>
    !>  @anchor     constructor
    !!
    ! --------------------------------------------------------------------
    type(Hwm07_class) function constructor()
        constructor%cseason = 0
        constructor%cwave = 0
        constructor%ctide = 0
        constructor%content = .true.          ! Season/Waves/Tides
        constructor%component = .true.      ! Compute zonal/meridional
        
        ! Initialization flags and information
        
        constructor%modelinit = .true.
        constructor%reset = .true.
        constructor%defaultdata = 'hwm071308e.dat'

        apexdata%kgma = 0
        apexdata%io = 1
        apexdata%jo = 1
        apexdata%ko = 1
        apexdata%loaddata = .true.

        dwmdata%modelinit = .true.
        dwmdata%defaultdata = 'dwm07b_104i.dat'
    end function constructor

    !===========================================================================
    !!
    !>  @anchor     destroy
    !!
    !>  @brief      Destroys all that needs destruction
    !>  @author     Christopher Kebschull
    !!
    !!
    !>  @date       <ul>
    !!                <li>03.01.2018 (initial design)</li>
    !!              </ul>
    !!
    !---------------------------------------------------------------------------
    subroutine destroy(this)
        class(Hwm07_class)    :: this

        if (allocated(this%nb)) deallocate(this%nb)
        if (allocated(this%order)) deallocate(this%order)
        if (allocated(this%vnode)) deallocate(this%vnode)
        if (allocated(this%mparm)) deallocate(this%mparm)
        if (allocated(this%gfs)) deallocate(this%gfs)
        if (allocated(this%gfm)) deallocate(this%gfm)
        if (allocated(this%gfl)) deallocate(this%gfl)
        if (allocated(this%gvbar)) deallocate(this%gvbar)
        if (allocated(this%gwbar)) deallocate(this%gwbar)
        if (allocated(this%gbz)) deallocate(this%gbz)
        if (allocated(this%gbm)) deallocate(this%gbm)
        if (allocated(this%gzwght)) deallocate(this%gzwght)

        if (allocated(vshbasisdata%fn_map1)) deallocate(vshbasisdata%fn_map1)
        if (allocated(vshbasisdata%A)) deallocate(vshbasisdata%A)
        if (allocated(vshbasisdata%B)) deallocate(vshbasisdata%B)
        if (allocated(vshbasisdata%anm)) deallocate(vshbasisdata%anm)
        if (allocated(vshbasisdata%bnm)) deallocate(vshbasisdata%bnm)
        if (allocated(vshbasisdata%cm)) deallocate(vshbasisdata%cm)
        if (allocated(vshbasisdata%cn)) deallocate(vshbasisdata%cn)
        if (allocated(vshbasisdata%e0n)) deallocate(vshbasisdata%e0n)
        if (allocated(vshbasisdata%fnm)) deallocate(vshbasisdata%fnm)
        if (allocated(vshbasisdata%m_arr)) deallocate(vshbasisdata%m_arr)
        if (allocated(vshbasisdata%n_arr)) deallocate(vshbasisdata%n_arr)
        if (allocated(vshbasisdata%cosmz)) deallocate(vshbasisdata%cosmz)
        if (allocated(vshbasisdata%sinmz)) deallocate(vshbasisdata%sinmz)

        if (allocated(dwmdata%termarr)) deallocate(dwmdata%termarr)
        if (allocated(dwmdata%coeff)) deallocate(dwmdata%coeff)
        if (allocated(dwmdata%vsh_terms)) deallocate(dwmdata%vsh_terms)
        if (allocated(dwmdata%termval)) deallocate(dwmdata%termval)

        if (allocated(apexdata%wk)) deallocate(apexdata%wk)

    end subroutine destroy

    !
    !  Horizontal Wind Model 07 (HWM07)
    !  Version HWM071308E, 13 July 2008
    !  See readme.txt file for detailed release notes.
    !
    !  AUTHOR
    !    Douglas P. Drob
    !    Space Science Division
    !    Naval Research Laboratory
    !    4555 Overlook Ave.
    !    Washington, DC 20375
    !
    !  Point of Contact
    !    msishwmhelp@nrl.navy.mil
    !
    !  DATE
    !    19 August 2008
    !
    !  REFERENCE
    !    Drob, D. P, J. T. Emmert, G. Crowley, J. M. Picone, G. G. Shepherd, 
    !      W. Skinner, Paul Hayes, R. J. Niciejewski, M. Larsen, C.Y. She, 
    !      J. W. Meriwether, G. Hernandez, M. J. Jarvis, D. P. Sipler, C. A. Tepley,
    !      M. S. OBrien, J. R. Bowman, Q. Wu, Y. Murayama, S. Kawamura, I.M. Reid,
    !      and R.A. Vincent (2008), An Empirical Model of the Earths Horizontal 
    !      Wind Fields: HWM07, J. Geophy. Res., doi:10.1029/2008JA013668.
    !
    !==================================================================================
    ! Input arguments:
    !        iyd - year and day as yyddd
    !        sec - ut(sec)
    !        alt - altitude(km)
    !        glat - geodetic latitude(deg)
    !        glon - geodetic longitude(deg)
    !        stl - not used
    !        f107a - not used
    !        f107 - not used
    !        ap - two element array with
    !             ap(1) = not used
    !             ap(2) = current 3hr ap index
    !
    ! Output argument:
    !        w(1) = meridional wind (m/sec + northward)
    !        w(2) = zonal wind (m/sec + eastward)
    !
    !================================================================================

    subroutine hwm07(this,iyd,sec,alt,glat,glon,stl,f107a,f107,ap,w)

        implicit none

        class(Hwm07_class)      :: this
        integer(4),intent(in)   :: iyd
        real(4),intent(in)      :: sec,alt,glat,glon,stl,f107a,f107
        real(4),intent(in)      :: ap(2)
        real(4),intent(out)     :: w(2)

        real(4)                 :: qw(2),dw(2)

        call this%hwmqt(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,qw)
        
        if (ap(2) .ge. 0.0) then
          call this%dwm07b_hwm_interface(iyd,sec,alt,glat,glon,ap,dw)
          w = qw + dw
        else
          w = qw
        endif
        
        return
        
    end subroutine HWM07

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !                       HWM-07 Legacy Wrapper
    !
    !  Description: Emulate HWM-93 subroutine calling convension
    !
    !  Programming Notes:
    !
    !
    !  Required Subroutines:
    !
    !    loadmodel()
    !    HWMupdate()
    !
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine HWMQT(this,IYD,SEC,ALT,GLAT,GLON,STL,F107A,F107,AP,W)
        
        implicit none

        class(Hwm07_class)      :: this
        INTEGER,intent(in)      :: IYD
        REAL(4),intent(in)      :: SEC,ALT,GLAT,GLON,STL,F107A,F107
        REAL(4),intent(in)      :: AP(2)
        REAL(4),intent(out)     :: W(2)

        real(8)                 :: input(5)
        real(8)                 :: u,v

        
        input(1) = dble(mod(IYD,1000))
        input(2) = dble(sec)
        input(3) = dble(glon)
        input(4) = dble(glat)
        input(5) = dble(alt)

        if (this%modelinit) then
            !call loadmodel(defaultdata)
            this%last = 1.0d-32                                                 ! in order to init that one... (VB)
            this%modelinit = .false.                                            !** initialisation performed by NEPTUNE
            call this%HWMupdate(input,this%last,this%gfs,this%gfl,this%gfm, &
                        this%gvbar,this%gwbar,this%gbz,this%gbm,this%gzwght,this%glev,u,v)
        endif
      
        if (this%reset) then
            this%last = 1.0d-32
            this%reset = .false.
        endif
      
        call this%HWMupdate(input,this%last,this%gfs,this%gfl,this%gfm, &
                    this%gvbar,this%gwbar,this%gbz,this%gbm,this%gzwght,this%glev,u,v)
        
        w(1) = real(-v)
        w(2) = real(u)

        return

    end subroutine HWMQT

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !   This subroutine is used to calculate the Vector Spherical
    !   Harmonic basis functions for a given observation.
    !
    !  Programming Notes:
    ! 
    !   This subroutine is only OPENMP/THREAD SAFE when no calls to
    !   loadmodel() are made.
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine HWMupdate(this,input,last,fs,fl,fm,vbar,wbar,ebz,ebm,zwght,lev,u,v)

        implicit none

        class(Hwm07_class)              :: this
        real(8),intent(in)              :: input(5) ! jday,utsec,glon,glat,alt
        real(8),intent(inout)           :: last(5)
        real(8),intent(inout)           :: fs(0:this%maxs,2)
        real(8),intent(inout)           :: fm(0:this%maxm,2)
        real(8),intent(inout)           :: fl(0:this%maxl,2)
        
        real(8),intent(inout)           :: vbar(0:this%maxn,0:this%maxo)
        real(8),intent(inout)           :: wbar(0:this%maxn,0:this%maxo)
        real(8),intent(inout),target    :: ebz(this%nbf,0:this%p)
        real(8),intent(inout),target    :: ebm(this%nbf,0:this%p)

        real(8),intent(inout)           :: zwght(0:this%p)
        integer,intent(inout)           :: lev
       
        real(8),intent(out)             :: u,v
     
        ! Local variables

        real(8),pointer           :: bz(:)
        real(8),pointer           :: bm(:)
        
        real(8)                   :: cs,ss,cm,sm,cl,sl
        real(8)                   :: cmcs,smcs,cmss,smss
        real(8)                   :: clcs,slcs,clss,slss
        real(8)                   :: AA,BB,CC,DD
        real(8)                   :: vb,wb
        
        integer                   :: b,c,d,m,n,s,l
        
        integer                   :: amaxs,amaxn
        integer                   :: pmaxm,pmaxs,pmaxn
        integer                   :: tmaxl,tmaxs,tmaxn

        logical                   :: refresh(5)
           
        real(8),parameter         :: twoPi = 2.0d0*3.1415926535897932384626433832795d0
        real(8),parameter         :: deg2rad = twoPi/360.0d0
        
    ! ====================================================================
    ! Update VSH model terms based on any change in the input parameters
    ! ====================================================================
        
        refresh = .false.

        ! Seasonal variations

        if (input(1) .ne. last(1)) then
            AA = input(1)*twoPi/365.25d0
            do s = 0,this%MAXS
                BB = dble(s)*AA
                fs(s,1) = dcos(BB)
                fs(s,2) = dsin(BB)
            enddo
            refresh(1:5) = .true.
            last(1) = input(1)
        endif

        ! Hourly time changes, tidal variations

        if (input(2) .ne. last(2) .or. input(3) .ne. last(3)) then

            AA = mod(input(2)/3600.d0 + input(3)/15.d0 + 48.d0,24.d0)
            BB = AA*twoPi/24.d0
            do l = 0,this%MAXL
                CC = dble(l)*BB
                fl(l,1) = dcos(CC)
                fl(l,2) = dsin(CC)
            enddo

            refresh(3) = .true.
            last(2) = input(2)
        endif

        ! Longitudinal variations, planetary waves
       
        if (input(3) .ne. last(3)) then
            AA = input(3)*deg2rad
            do m = 0,this%MAXM
                BB = dble(m)*AA
                fm(m,1) = dcos(BB)
                fm(m,2) = dsin(BB)
            enddo
            refresh(2) = .true.
            last(3) = input(3)
        endif

        ! Latitude

        if (input(4) .ne. last(4)) then
            AA = (90.0d0 - input(4))*deg2rad        ! theta = colatitude in radians
            call this%vshbasis(this%maxn,this%maxo,AA,vbar,wbar)
            refresh(1) = .true.
            refresh(2) = .true.
            refresh(3) = .true.
            refresh(4) = .true.
            last(4) = input(4)
        endif

        ! Altitude
        
        if (input(5) .ne. last(5)) then
            call this%vertwght(input(5),zwght,lev)
            last(5) = input(5)
        endif

        ! ====================================================================
        ! Linearize the VSH functions
        ! ====================================================================
           
        u = 0.0d0
        v = 0.0d0
           
        refresh = .true.
           
        ebz = 0.0
        ebm = 0.0
        
        do b = 0,this%p

            if (zwght(b) .eq. 0.d0) cycle
            
            d = b + lev
            
            bz => ebz(:,b)
            bm => ebm(:,b)

            amaxs = this%order(1,d)
            amaxn = this%order(2,d)
            pmaxm = this%order(3,d)
            pmaxs = this%order(4,d)
            pmaxn = this%order(5,d)
            tmaxl = this%order(6,d)
            tmaxs = this%order(7,d)
            tmaxn = this%order(8,d)

            c = 1
     
            ! ------------- Seasonal - Zonal average (m = 0) ----------------

            if (refresh(1) .and. this%content(1)) then

                c = 1                           ! global constants
                do n = 1,AMAXN
                    bz(c)   = -0.5d0*vbar(n,0)  ! Cr
                    bz(c+1) =  0.0d0            ! Br
                    bm(c)   =  0.0d0            ! Cr
                    bm(c+1) =  0.5d0*vbar(n,0)  ! Br
                    c = c + 2
                enddo
                    
                do s = 1,AMAXS                   ! seasonal variations
                    cs = fs(s,1)
                    ss = fs(s,2)
                    do n = s,AMAXN
                        vb = vbar(n,s)
                        wb = wbar(n,s)
                        AA =  vb*cs
                        BB =  vb*ss
                        CC = -wb*ss
                        DD = -wb*cs
                          bz(c) = -AA   ! Cr
                        bz(c+1) =  BB   ! Ci
                        bz(c+2) =  CC   ! Br
                        bz(c+3) =  DD   ! Bi
                          bm(c) =  CC   ! Cr
                        bm(c+1) =  DD   ! Ci
                        bm(c+2) =  AA   ! Br
                        bm(c+3) = -BB   ! Bi
                        c = c + 4
                    enddo
                enddo
                this%cseason = c
            else
                c = this%cseason
            endif
                
            ! ---------------- Stationary planetary waves --------------------

            if (refresh(2) .and. this%content(2)) then

                do m = 1,pmaxm

                   cm = fm(m,1)
                   sm = fm(m,2)

                   do n = m,pmaxn           ! s = 0
                   
                        vb = vbar(n,m)
                        wb = wbar(n,m)
                   
                        bz(c) =   -vb*cm    ! Cr * (cm) * -vb
                        bz(c+1) =  vb*sm    ! Ci * (sm) *  vb
                        bz(c+2) = -wb*sm    ! Br * (sm) * -wb
                        bz(c+3) = -wb*cm    ! Bi * (cm) * -wb
                        
                        bm(c) =   -wb*sm    ! Cr * (sm) * -wb
                        bm(c+1) = -wb*cm    ! Ci * (sm) * -wb
                        bm(c+2) =  vb*cm    ! Br * (cm) *  vb
                        bm(c+3) = -vb*sm    ! Bi * (sm) * -vb
                        
                        c = c + 4
                   
                   enddo
                   
                   do s = 1,pmaxs
                   
                      cs = fs(s,1)
                      ss = fs(s,2)
                   
                      do n = m,pmaxn
                         vb = vbar(n,m)
                         wb = wbar(n,m)
                         
                         bz(c) =   -vb*cm*cs    ! Crc * (cmcs) * -vb
                         bz(c+1) =  vb*sm*cs    ! Cic * (smcs) *  vb
                         bz(c+2) = -wb*sm*cs    ! Brc * (smcs) * -wb
                         bz(c+3) = -wb*cm*cs    ! Bic * (cmcs) * -wb
                         bz(c+4) = -vb*cm*ss    ! Crs * (cmss) * -vb
                         bz(c+5) =  vb*sm*ss    ! Cis * (smss) *  vb
                         bz(c+6) = -wb*sm*ss    ! Brs * (smss) * -wb
                         bz(c+7) = -wb*cm*ss    ! Bis * (cmss) * -wb
                         
                         bm(c) =   -wb*sm*cs    ! Crc * (smcs) * -wb
                         bm(c+1) = -wb*cm*cs    ! Cic * (smcs) * -wb
                         bm(c+2) =  vb*cm*cs    ! Brc * (cmcs) *  vb
                         bm(c+3) = -vb*sm*cs    ! Bic * (smcs) * -vb
                         bm(c+4) = -wb*sm*ss    ! Crs * (smss) * -wb
                         bm(c+5) = -wb*cm*ss    ! Cis * (smss) * -wb
                         bm(c+6) =  vb*cm*ss    ! Brs * (cmss) *  vb
                         bm(c+7) = -vb*sm*ss    ! Bis * (smss) * -vb
                         
                         c = c + 8
                      
                      enddo
                
                   enddo
                   this%cwave = c
                enddo        
            else
                c = this%cwave
            endif

            ! ---------------- Migrating Solar Tides ---------------------

            if (refresh(3) .and. this%content(3)) then
                 do l = 1,tmaxl
                
                   cl = fl(l,1)
                   sl = fl(l,2)

                   s = 0
                   do n = l,tmaxn
                   
                        vb = vbar(n,l)
                        wb = wbar(n,l)                                          
                   
                        bz(c) =   -vb*cl    ! Cr * (cl) * -vb
                        bz(c+1) =  vb*sl    ! Ci * (sl) *  vb
                        bz(c+2) = -wb*sl    ! Br * (sl) * -wb
                        bz(c+3) = -wb*cl    ! Bi * (cl) * -wb
                        
                        bm(c) =   -wb*sl    ! Cr * (sl) * -wb
                        bm(c+1) = -wb*cl    ! Ci * (sl) * -wb
                        bm(c+2) =  vb*cl    ! Br * (cl) *  vb
                        bm(c+3) = -vb*sl    ! Bi * (sl) * -vb
                        
                        c = c + 4
                   
                   enddo
                   
                   do s = 1,tmaxs
                      
                      cs = fs(s,1)
                      ss = fs(s,2)
                      
                      do n = l,tmaxn
                      
                         vb = vbar(n,l)
                         wb = wbar(n,l)
                                       
                         bz(c) =   -vb*cl*cs    ! Crc * (clcs) * -vb
                         bz(c+1) =  vb*sl*cs    ! Cic * (slcs) *  vb
                         bz(c+2) = -wb*sl*cs    ! Brc * (slcs) * -wb
                         bz(c+3) = -wb*cl*cs    ! Bic * (clcs) * -wb
                         
                         bz(c+4) = -vb*cl*ss    ! Crs * (clss) * -vb
                         bz(c+5) =  vb*sl*ss    ! Cis * (slss) *  vb
                         bz(c+6) = -wb*sl*ss    ! Brs * (slss) * -wb
                         bz(c+7) = -wb*cl*ss    ! Bis * (clss) * -wb
                                             
                         bm(c) =   -wb*sl*cs    ! Crc * (slcs) * -wb
                         bm(c+1) = -wb*cl*cs    ! Cic * (slcs) * -wb
                         bm(c+2) =  vb*cl*cs    ! Brc * (clcs) *  vb
                         bm(c+3) = -vb*sl*cs    ! Bic * (slcs) * -vb

                         bm(c+4) = -wb*sl*ss    ! Crs * (slss) * -wb
                         bm(c+5) = -wb*cl*ss    ! Cis * (slss) * -wb
                         bm(c+6) =  vb*cl*ss    ! Brs * (clss) *  vb
                         bm(c+7) = -vb*sl*ss    ! Bis * (slss) * -vb
                         
                         c = c + 8
                      enddo
                
                   enddo
                   this%ctide = c
                enddo
            else
                c = this%ctide
            endif
                
            ! ---------------- Non-Migrating Solar Tides ------------------
            
            ! TBD
                
            c = c - 1
            
            ! ====================================================================
            ! Calculate the wind components 
            ! ====================================================================

            if (this%component(0)) u = u + zwght(b)*dot_product(bz(1:c),this%mparm(1:c,d))
            if (this%component(1)) v = v + zwght(b)*dot_product(bm(1:c),this%mparm(1:c,d))

        enddo

        return 

    end subroutine HWMupdate

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! This subroutine initalizes the NEWmodel
    !
    !
    ! Required Subroutines:
    !
    !   vshengineinit()
    !
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine loadmodel(this, datafile)

        use slam_io, only: openFile, closeFile, DIRECT, IN_UNFORMATTED

        implicit none

        !character(128),intent(in)   :: datafile
        class(Hwm07_class)              :: this
        character(len=*),intent(in)     :: datafile

        integer                     :: i,j
        integer                     :: ich      ! file channel to read from (VB)
        integer                     :: ncomp

        if (allocated(this%vnode)) then
            deallocate(this%order,this%nb,this%vnode,this%mparm)
            deallocate(this%gfs,this%gfm,this%gfl,this%gvbar,this%gwbar, &
                        this%gzwght,this%gbz,this%gbm)
        endif

        !open(unit=23,file=trim(datafile),form='unformatted')
        ich = openFile(trim(adjustl(datafile)), DIRECT, IN_UNFORMATTED)

        read(ich) this%nbf,this%maxs,this%maxm,this%maxl,this%maxn,ncomp
        read(ich) this%nlev,this%p
        this%nnode = this%nlev + this%p
        allocate(this%nb(0:this%nnode))
        allocate(this%order(ncomp,0:this%nnode))
        allocate(this%vnode(0:this%nnode))
        read(ich) this%vnode
        this%vnode(3) = 0.0
        allocate(this%mparm(this%nbf,0:this%nlev))
        this%mparm = 0.0d0
        do i = 0,this%nlev-this%p+1-2
            read(ich) this%order(1:ncomp,i)
            read(ich) this%nb(i)
            read(ich) this%mparm(1:this%nbf,i)
        enddo
        ich = closeFile(ich)

        ! Set transition levels
        
        this%alttns = this%vnode(this%nlev-2)
        this%altsym = this%vnode(this%nlev-1)
        this%altiso = this%vnode(this%nlev)

        ! Initialize the vector spherical harmonics engine

        call this%vshengineinit()

        ! Allocate the global store of quasi-static parameters
        
        this%maxo = max(this%maxs,this%maxm,this%maxl)

        allocate(this%gfs(0:this%maxs,2),this%gfm(0:this%maxm,2),this%gfl(0:this%maxl,2))
        allocate(this%gvbar(0:this%maxn,0:this%maxo),this%gwbar(0:this%maxn,0:this%maxo))
        allocate(this%gbz(this%nbf,0:this%p),this%gbm(this%nbf,0:this%p))
        allocate(this%gzwght(0:this%p))
        
        this%gvbar = 0.0d0
        this%gwbar = 0.0d0
        this%gbz = 0.0d0
        this%gbm = 0.0d0
      
        ! Signal that the model has been initalized
        
        !modelinit = .false.
        
        ! Signal a reset of the input variable comparison flags
        
        this%reset = .true.

        return

    end subroutine loadmodel
     
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! This subroutine updates the vertical weights
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine vertwght(this,alt,wght,iz)

        implicit none

        class(Hwm07_class)      :: this
        real(8),intent(in)      :: alt
        real(8),intent(out)     :: wght(4)
        integer,intent(out)     :: iz
        
        real(8)             :: we(0:4)

        real(8)             :: e1(0:4) = &
            (/1.d0, 0.428251121076233d0,0.192825112107623d0,0.484304932735426d0,0.0d0/)
        real(8)             :: e2(0:4) = &
            (/0.d0, 0.571748878923767d0,0.807174887892377d0,-0.484304932735426d0,1.0d0/)

        real(8),parameter   :: H = 60.0d0
            
        iz = findspan(this%nnode-this%p-1,this%p,alt,this%vnode) - this%p

        iz = min(iz,26)
        
        wght(1) = bspline(this%p,this%nnode,this%vnode,iz,alt)
        wght(2) = bspline(this%p,this%nnode,this%vnode,iz+1,alt)
        if (iz .le. 25) then
            wght(3) = bspline(this%p,this%nnode,this%vnode,iz+2,alt)
            wght(4) = bspline(this%p,this%nnode,this%vnode,iz+3,alt)
            return
        endif
            if (alt .gt. 250.0d0) then
                we(0) = 0.0d0
                we(1) = 0.0d0
                we(2) = 0.0d0
                we(3) = exp(-(alt - 250.0d0)/H) 
                we(4) = 1.0d0
            else
                we(0) = bspline(this%p,this%nnode,this%vnode,iz+2,alt)
                we(1) = bspline(this%p,this%nnode,this%vnode,iz+3,alt)
                we(2) = bspline(this%p,this%nnode,this%vnode,iz+4,alt)
                we(3) = 0.0d0
                we(4) = 0.0d0
            endif
            wght(3) = dot_product(we,e1)
            wght(4) = dot_product(we,e2)
            
        return
        
    contains

        function bspline(p,m,V,i,u)

            implicit none
            
            real(8)     :: bspline
            integer     :: p,m
            real(8)     :: V(0:m)
            integer     :: i
            real(8)     :: u
            
            real(8)     :: N(0:p+1)
            real(8)     :: Vleft,Vright
            real(8)     :: saved,temp
            integer     :: j,k
                    
            if ((i .eq. 0) .and. (u .eq. V(0))) then
                bspline = 1.d0
                return
            endif
            
            if ((i .eq. (m-p-1)) .and. (u .eq. V(m))) then
                bspline = 1.d0
                return
            endif

            if (u .lt. V(i) .or. u .ge. V(i+p+1)) then
                bspline = 0.d0
                return
            endif
            
            N = 0.0d0
            do j = 0,p
                if (u .ge. V(i+j) .and. u .lt. V(i+j+1)) then
                    N(j) = 1.0d0
                else
                    N(j) = 0.0d0
                endif
            enddo
            
            do k = 1,p
                if (N(0) .eq. 0.d0) then
                    saved = 0.d0
                else
                    saved = ((u - V(i))*N(0))/(V(i+k) - V(i))
                endif
                do j = 0,p-k
                    Vleft = V(i+j+1)
                    Vright = V(i+j+k+1)
                    if (N(j+1) .eq. 0.d0) then
                        N(j) = saved
                        saved = 0.d0
                    else
                        temp = N(j+1)/(Vright - Vleft)
                        N(j) = saved + (Vright - u)*temp
                        saved = (u - Vleft)*temp
                    endif
                enddo
            enddo
            
            bspline = N(0)

            return

        end function bspline

        ! =====================================================
        ! Function to locate the knot span
        ! =====================================================

        integer function findspan(n,p,u,V)

            implicit none
            
            integer,intent(in)      :: n,p
            real(8),intent(in)      :: u
            real(8),intent(in)      :: V(0:n+1)
            integer                 :: low,mid,high
            
            if (u .ge. V(n+1)) then
                findspan = n
                return
            endif
            
            low = p
            high = n+1
            mid = (low + high)/2

            do while (u .lt. V(mid) .or. u .ge. V(mid + 1))
                if (u .lt. V(mid)) then
                    high = mid
                else
                    low = mid
                endif
                mid = (low + high)/2
            end do
        
            findspan = mid
            return

        end function findspan
        
    end subroutine vertwght

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !   Description:  Computational Engine for calculating Normalized Associated Vector
    !                 Spherical Harmonic basis functions Vbar and Wbar.
    !
    !   Author: Douglas P. Drob
    !           Space Science Division
    !           Naval Research Laboratory
    !           4555 Overlook Ave
    !           Washington, DC
    !
    !   Date: January 2007.
    !
    !   Notes:
    !
    !   The routines for the calculation Pbar have been adapted from ALFPACK developed 
    !   by Paul Swarztrauber at The National Center for Atmospheric Research, Boulder, 
    !   Colorado (80307) U.S.A. which is sponsored by The National Science Foundation.
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ==========================================================================
    ! Description: Subroutine to calculate the VSH coeffiecents @ theta
    !
    ! Notes: Before calling vshbasis() the static coeffiecents must calculated
    !        by calling vshengineinit().
    !
    ! ==========================================================================

    subroutine vshbasis(this,maxn,maxo,theta,vbar,wbar)

        implicit none

        class(Hwm07_class)  :: this
        integer,intent(in)  :: maxn,maxo
        real(8),intent(in)  :: theta
        real(8),intent(out) :: vbar(0:maxn,0:maxo)
        real(8),intent(out) :: wbar(0:maxn,0:maxo)

        real(8)             :: pbar(0:vshdata%nmax,0:maxo+1)
        real(8)             :: pb(0:vshdata%nmax),td(0:vshdata%nmax)
        real(8)             :: p0i,p1i,cost,r
        integer             :: i,m,mm1,mp1,n,nm1,nm2,nmm

        ! Normalized Associated Legendre Polynomials

        pbar = 0.0d0

        nm1 = vshdata%nmax - 1
        nm2 = vshdata%nmax - 2
        do m = 0,maxo+1
            nmm = vshdata%nmax - m
            p0i = lfpt(nm1,m,theta,vshdata%p0(:,m))
            p1i = lfpt(nm2,m,theta,vshdata%p1(:,m))
            pbar(nm1,m) = p0i
            if (nmm .le. 0) cycle
            pbar(nm2,m) = p1i
            if (nmm .eq. 1) cycle
            cost = dcos(theta)
            do n = 0,nmm-1
                pb(n) = -cost
                td(n) = vshdata%sf(n,m)
            enddo
            if (abs(p0i) .ge. abs(p1i)) then
                pb(0) = p0i
                r = -td(0)*pb(0)
                call tridiag(nmm-1,r,td(0),pb(1),td(1))
            else
                pb(0) = p0i
                pb(1) = p1i
                r = -td(1)*pb(1)
                call tridiag(nmm-2,r,td(1),pb(2),td(2))
            endif
            do n = m,vshdata%nmax-1
                i = vshdata%nmax-n-1
                pbar(n,m) = pb(i)
            enddo
        enddo

        ! Vector Spherical Harmonic Basis Functions from Pbar

        do n = 0,maxn
            vbar(n,0) = -pbar(n,1)
        enddo

        do m = 1,maxo
            mm1 = m - 1
            mp1 = m + 1
             do n = m,maxn
                nm1 = n-1
                vbar(n,m) = vshdata%a(n,m)*pbar(n,mp1) + vshdata%b(n,m)*pbar(n,mm1)
                wbar(n,m) = vshdata%c(n,m)*pbar(nm1,mm1) + vshdata%d(n,m)*pbar(nm1,mp1)
          enddo
        enddo

        return

    contains

        ! --------------------------------------------------------------
        ! Function for the normalized associated legendre polynomial 
        ! along the diagonal from a recursion relation.
        ! --------------------------------------------------------------

        real(8) function lfpt(n,m,theta,cp)
            implicit none
            integer,intent(in)      :: n,m
            real(8),intent(in)      :: theta
            real(8),intent(in)      :: cp(n)
            real(8)                 :: cdt,sdt
            real(8)                 :: ct,st,cth
            integer                 :: kdo,k
            lfpt = 0.0d0
            if (m .gt. n) return
            if (n .le. 0 .and. m .le. 0) then
                lfpt = dsqrt(0.5d0)
                return
            endif
            cdt = dcos(theta+theta)
            sdt = dsin(theta+theta)
            if (mod(n,2) .le. 0) then
                ct = 1.0d0
                st = 0.0d0
                if (mod(m,2) .le. 0) then
                    kdo = n/2+1
                    lfpt = 0.5d0*cp(1)
                    do k = 2,kdo
                        cth = cdt*ct - sdt*st
                        st = sdt*ct + cdt*st
                        ct = cth
                        lfpt = lfpt + cp(k)*ct
                    enddo
                else
                    kdo = n/2
                    do k = 1,kdo
                        cth = cdt*ct - sdt*st
                        st = sdt*ct + cdt*st
                        ct = cth
                        lfpt = lfpt + cp(k)*st
                    enddo
                endif
            else
                kdo = (n+1)/2
                ct = dcos(theta)
                st = -dsin(theta)
                if (mod(m,2) .le. 0) then
                    do k = 1,kdo
                        cth = cdt*ct - sdt*st
                        st = sdt*ct + cdt*st
                        ct = cth
                        lfpt = lfpt + cp(k)*ct
                    enddo
                else
                    do k = 1,kdo
                        cth = cdt*ct - sdt*st
                        st = sdt*ct + cdt*st
                        ct = cth
                        lfpt = lfpt + cp(k)*st
                    enddo
                endif
            endif
            return
        end function lfpt

        ! -----------------------------
        !  Tri-diagonal matrix solver 
        ! -----------------------------

        subroutine tridiag(n,r,a,b,c)
            implicit none
            integer,intent(in)      :: n
            real(8),intent(in)      :: r,a(n)
            real(8),intent(inout)   :: b(n),c(n)
            real(8)                 :: bih,bih1,b1
            real(8)                 :: qih,q2,ratio,rih
            integer                 :: i,j
            select case (n)
                case(:0)
                    return
                case(1)
                    b(1) = r/b(1)
                    return
                case(2)
                    qih = a(2)
                    bih = b(2)
                case(3:)
                    qih = a(n)
                    bih = b(n)
                    do j = 3,n
                        i = n - j + 2
                        if (abs(bih) .ge. abs(c(i))) then
                            ratio = c(i)/bih
                            c(i) = 0.0d0
                            b(i+1) = qih/bih
                            bih = b(i) - ratio*qih
                            qih = a(i)
                        else
                            b(i+1) = b(i)/c(i)
                            c(i) = a(i)/c(i)
                            bih1 = qih - bih*b(i+1)
                            qih = -bih*c(i)
                            bih = bih1
                        endif
                    enddo
            end select
            if (abs(bih) .ge. abs(c(1))) then
                q2 = qih/bih
                bih = b(1) - c(1)/bih*qih
                b(1) = r/bih
                b(2) = -q2*b(1)
            else
                ratio = bih/c(1)
                bih = qih - ratio*b(1)
                rih = -ratio*r
                b1 = rih/bih
                b(2) = (r - b(1)*b1)/c(1)
                b(1) = b1
            endif
            if (n-3 .ge. 0) then
                do i = 3,n
                    b(i) = -b(i)*b(i-1) - c(i-1)*b(i-2)
                enddo
            endif
            return
        end subroutine tridiag

    end subroutine vshbasis
        
    ! ========================================================================
    ! This subroutine intializes the scaling and recursion coeffiecents
    ! for normalized associated vector spherical harmonic basis funcions.
    ! ========================================================================

    subroutine vshengineinit(this)

        implicit none

        class(Hwm07_class)  :: this

        real(8)             :: fnmm,fnpn,fnpm
        real(8)             :: nfac1,nfac2
        integer             :: m,n,npm,nmm,nm1,nm2,i
        
        ! Factors for the normalized associated legendre polynomial

        vshdata%sf = 0.0d0
        nm1 = vshdata%nmax - 1
        nm2 = vshdata%nmax - 2
        do m = 0,vshdata%mmax
            call alfk(nm1,m,vshdata%p0(:,m))
            call alfk(nm2,m,vshdata%p1(:,m))
            fnmm = dble(vshdata%nmax - m)
            fnpn = dble(vshdata%nmax + nm1)
            fnpm = dble(vshdata%nmax + m)
            i = 0
            do n = m,nm1
                fnmm = fnmm - 1.d0
                fnpn = fnpn - 2.d0
                fnpm = fnpm - 1.d0
                vshdata%sf(i,m) = dsqrt(dble(fnmm*fnpm/(fnpn*(fnpn+2.d0))))
                i = i + 1
            end do
        enddo
            
        ! Theta indepedent scaling factors for calculation of Vbar and Wbar 
        ! from Pbar.
        
        do n = 1,vshdata%nmax
            nfac1 = 0.5d0/dsqrt(dble(n*(n+1)))
            nfac2 = nfac1*dsqrt(dble(2*n+1)/dble(2*n-1))
            do m = 1,n
                npm = n + m 
                nmm = n - m
                vshdata%a(n,m) = -nfac1*dsqrt(dble(nmm*(npm+1)))
                vshdata%b(n,m) =  nfac1*dsqrt(dble(npm*(nmm+1)))
                vshdata%c(n,m) =  nfac2*dsqrt(dble(npm*(npm-1)))
                vshdata%d(n,m) =  nfac2*dsqrt(dble(nmm*(nmm-1)))
            enddo
        enddo

        return

    contains

        ! -------------------------------------------------------------
        ! For a given M and N, the unnormalized Legendre polynomial 
        ! are calculated for all cases where (m .ge. 0).
        ! ------------------------------------------------------------

        subroutine alfk(n,m,cp)
            implicit none
            integer,intent(in)      :: n,m
            real(8),intent(out)     :: cp(n)
            real(8)                 :: fnum,fden,fnmh,fnnp1,fnmsq,fk
            real(8)                 :: pm1,t1,t2,cp2,a1,b1,c1
            integer                 :: nmms2,l,i
            if (m .gt. n) then
                cp(1) = 0.0d0
                return
            endif
            if (n .le. 0) then
                cp(1) = dsqrt(2.d0)
                return
            endif
            if (n .eq. 1) then
                if (m .eq. 0) then
                    cp(1) = dsqrt(1.5d0)
                else
                    cp(1) = dsqrt(0.75d0)
                endif 
                return
            endif
            if (mod(n+m,2) .eq. 0) then
                nmms2 = (n - m)/2
                fnum = dble(n + m + 1)
                fnmh = dble(n - m + 1)
                pm1 = 1.0d0
            else
                nmms2 = (n - m - 1)/2
                fnum = dble(n + m + 2)
                fnmh = dble(n - m + 2)
                pm1 = -1.0d0
            endif 
            t1 = 1.0d0
            t2 = 1.0d0
            if (nmms2 .ge. 1) then
                fden = 2.0d0
                do i = 1,nmms2
                    t1 = fnum*t1/fden
                    fnum = fnum + 2.0d0
                    fden = fden + 2.0d0
                enddo
            endif
            if (m .ne. 0) then
                do i = 1,m
                    t2 = fnmh*t2/(fnmh + pm1)
                    fnmh = fnmh + 2.0d0
                enddo
            endif
            if (mod(m/2,2) .ne. 0) t1 = -t1
            cp2 = t1 * dsqrt( (dble(n) + 0.5d0) *t2 ) / (2.0d0**(n-1))
            fnnp1 = dble(n*(n+1))
            fnmsq = fnnp1 - 2.0d0*dble(m*m)
            l = (n+1)/2
            if (mod(n,2) .eq. 0 .and. mod(m,2) .eq. 0) l = l + 1
            cp(l) = cp2
            if (l .le. 1) return
            fk = dble(n)
            a1 = (fk - 2.0d0)*(fk - 1.0d0) - fnnp1
            b1 = 2.0d0*(fk*fk - fnmsq)
            cp(l-1) = b1*cp(l)/a1
            l = l - 1
            do while (l .gt. 1)
                fk = fk - 2.0d0
                a1 = (fk - 2.0d0)*(fk - 1.0d0) - fnnp1
                b1 = -2.0d0*(fk*fk - fnmsq)
                c1 = (fk + 1.0d0)*(fk + 2.0d0) - fnnp1
                cp(l-1) = -(b1*cp(l) + c1*cp(l+1))/a1
                l = l - 1
            enddo
            return
         end subroutine alfk
        
    end subroutine vshengineinit

    ! ===============================================================
    ! Initalize the apex coordinate module variables
    ! ===============================================================
    subroutine apxrda(this, datafile)

        use slam_io, only: openFile, closeFile, DIRECT, IN_UNFORMATTED
        implicit none

        class(Hwm07_class)          :: this
        character(len=*),intent(in) :: datafile  ! datafile for QD interpolation grid VB
        integer(4) :: i
        integer    :: ich                     ! input file channel VB

        !open(unit=77,file='apexgrid.dat',form='unformatted')
        ich = openFile(datafile, DIRECT, IN_UNFORMATTED) 
        
        read(ich) apexdata%kgma,apexdata%glalmn,apexdata%glalmx,apexdata%nla,apexdata%nlo,apexdata%nal, &
            apexdata%lbx,apexdata%lby,apexdata%lbz,apexdata%lbv,apexdata%lla,apexdata%llo,apexdata%lal
            
        read(ich) apexdata%colat,apexdata%elon,apexdata%vp,apexdata%ctp,apexdata%stp
        read(ich) apexdata%rtod,apexdata%dtor,apexdata%re,apexdata%req,apexdata%msgu,apexdata%pola

        apexdata%lwk = apexdata%nla*apexdata%nlo*apexdata%nal*5 &
                                + apexdata%nla + apexdata%nlo + apexdata%nal

        if (allocated(apexdata%wk)) deallocate(apexdata%wk)
        allocate(apexdata%wk(apexdata%lwk))
        
        read(ich) (apexdata%wk(i),i=1,apexdata%lwk)
        ich = closeFile(ich)
        
        apexdata%loaddata = .false.
        
        return

    end subroutine apxrda

    ! =======================================================================
    ! Convert from (glat,glon) to apex coordinates
    ! =======================================================================
    subroutine apex(this,glat,glon,alt,hr,alon,xlatqd,f1,f2,ist)
         
        implicit none

        class(Hwm07_class)  :: this
        real(4),intent(in)  :: glat,glon,alt,hr

        real(4)        :: b(3),bhat(3)
        real(4)        :: bmag,si
        real(4)        :: alon
        real(4)        :: xlatm,vmp,w,d
        real(4)        :: be3,sim
        real(4)        :: xlatqd
        real(4)        :: d1(3),d2(3),d3(3)
        real(4)        :: e1(3),e2(3),e3(3)
        real(4)        :: f1(2),f2(2)
        
        integer(4)     :: ist
      
        real(4)        :: dfxdth,dfydth,dfzdth,dfvdth,dfxdln
        real(4)        :: dfydln,dfzdln,dfvdln,dfxdh,dfydh,dfzdh,dfvdh
        
        real(4)        :: fx,fy,fz,fv
        
        real(4)        :: cth,sth
        real(4)        :: gradx(3),grady(3),gradz(3),gradv(3)
        real(4)        :: glatx
        
        real(4)        :: fxdum,fydum,fzdum,fvdum
        real(4)        :: dmxdth,dmydth,dmzdth,dmvdth
        real(4)        :: dmxdh,dmydh,dmzdh,dmvdh
        
        real(4)        :: clm,r3_2
        real(4)        :: grclm(3),clmgrp(3),rgrlp(3)
            
        real(4)        :: f
        
    ! Initalize the module if needed (done by NEPTUNE) VB
        
        !if (loaddata) then
        !    call apxrda()
        !    loaddata = .false.
        !endif
                       
    ! Interpolate the 3d data cube
                       
        call this%intrp(glat,glon,alt,apexdata%wk(apexdata%lbx), &
                   apexdata%wk(apexdata%lby), &
                   apexdata%wk(apexdata%lbz), &
                   apexdata%wk(apexdata%lbv), &
                   apexdata%nla,apexdata%nlo, &
                   apexdata%nal,apexdata%wk(apexdata%lla), &
                   apexdata%wk(apexdata%llo), &
                   apexdata%wk(apexdata%lal),fx,fy,fz,fv, &
                   dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln, &
                   dfvdln,dfxdh,dfydh,dfzdh,dfvdh,ist)

        call this%adpl(glat,glon,cth,sth,fx,fy,fz,fv, &
             dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln)
         
        call this%gradxyzv(alt,cth,sth,  &
            dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln, &
            dfxdh,dfydh,dfzdh,dfvdh,gradx,grady,gradz,gradv)

        ! If the point is very close to either the North or South
        ! geographic pole, recompute the east-west gradients after
        ! stepping a small distance from the pole.
        
        if (glat .gt. apexdata%glalmx .or. glat .lt. apexdata%glalmn) then
                
            glatx = apexdata%glalmx
              
            if (glat .lt. 0.0) glatx = apexdata%glalmn
              
            call this%intrp(glatx,glon,alt,apexdata%wk(apexdata%lbx), &
                  apexdata%wk(apexdata%lby), &
                  apexdata%wk(apexdata%lbz), &
                  apexdata%wk(apexdata%lbv), &
                  apexdata%nla,apexdata%nlo,apexdata%nal, &
                  apexdata%wk(apexdata%lla), &
                  apexdata%wk(apexdata%llo), &
                  apexdata%wk(apexdata%lal),fxdum,fydum,fzdum,fvdum, &
                  dmxdth,dmydth,dmzdth,dmvdth,dfxdln,dfydln,dfzdln, &
                  dfvdln,dmxdh,dmydh,dmzdh,dmvdh,ist)
              
            call this%adpl(glatx,glon,cth,sth,fxdum,fydum,fzdum,fvdum,&
                dmxdth,dmydth,dmzdth,dmvdth,dfxdln,dfydln,dfzdln,dfvdln)
             
            call this%grapxyzv(alt,cth,sth,dfxdln, &
                dfydln,dfzdln,dfvdln,gradx,grady,gradz,gradv)
          
        endif

        call this%gradlpv(hr,alt,fx,fy,fz,fv,gradx,grady,gradz,gradv, &
            xlatm,alon,vmp,grclm,clmgrp,xlatqd,rgrlp,b,clm,r3_2)
         
        call this%basvec(hr,xlatm,grclm,clmgrp,rgrlp,b,clm,r3_2, &
            bmag,sim,si,f,d,w,bhat,d1,d2,d3,e1,e2,e3,f1,f2)

        be3 = bmag/d

        ist = 0

        return

    end subroutine apex

    ! =========================================================================
    ! Interpolation subroutine
    ! =========================================================================
    subroutine intrp(this,glat,glon,alt,x,y,z,v,nlat,nlon,nalt, &
                 gplat,gplon,gpalt,fx,fy,fz,fv, &
                 dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln, &
                 dfvdln,dfxdh,dfydh,dfzdh,dfvdh,ist)
         
        implicit none

        class(Hwm07_class)  :: this

        real(4)        :: glat,glon,alt
        integer(4)     :: nlat,nlon,nalt
        real(4)        :: x(nlat,nlon,nalt)
        real(4)        :: y(nlat,nlon,nalt)
        real(4)        :: z(nlat,nlon,nalt)
        real(4)        :: v(nlat,nlon,nalt)
        real(4)        :: gplat(nlat),gplon(nlon),gpalt(nalt)
        real(4)        :: fx,fy,fz,fv
        real(4)        :: dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln, &
                       dfvdln,dfxdh,dfydh,dfzdh,dfvdh
        integer(4)     :: ist

        
        real(4)        :: dlat,dlon
        real(4)        :: xi,yj
        real(4)        :: hti,diht
        real(4)        :: zk

        real(4)        :: dfxdn,dfxde,dfxdd
        real(4)        :: dfydn,dfyde,dfydd
        real(4)        :: dfzdn,dfzde,dfzdd
        real(4)        :: dfvdn,dfvde,dfvdd
        real(4)        :: dmf,dmdfdn,dmdfde,dmdfdd
        
        real(4)        :: fac
        real(4)        :: omfac

        integer(4)     :: i,j,k

        ist = 0

    ! input bounds checks

        if (glat .lt. gplat(1) .or. glat .gt. gplat(nlat)) stop 'glat out of bounds'
        if (alt  .lt. gpalt(1) .or. alt  .gt. gpalt(nalt)) stop 'alt out of bounds'
        if (glon .lt. gplon(1))    glon = glon + 360.0
        if (glon .gt. gplon(nlon)) glon = glon - 360.0
        if (glon .lt. gplon(1) .or. glon .gt. gplon(nlon)) stop 'glon out of bounds'

    ! Locate the closest indicies of the points in the array and
    ! and calculate the fractional offset for linear interpolation

        ! Find I
        
        i = apexdata%io
        if (gplat(i+1) .lt. glat) then
            do while (gplat(i+1) .lt. glat)
                i = i + 1
            enddo
        endif
        if (gplat(i) .gt. glat) then
            do while (gplat(i) .gt. glat)
                i = i - 1
            enddo
        endif
        apexdata%io = i
        dlat = gplat(i+1) - gplat(i)
        xi = (glat - gplat(i)) / dlat
            
        ! Find J

        j = apexdata%jo
        if (gplon(j+1) .lt. glon) then
            do while (gplon(j+1) .lt. glon)
                j = j + 1
            enddo
        endif
        if (gplon(j) .gt. glon) then
            do while (gplon(j) .gt. glon)
                j = j - 1
            enddo
        endif
        apexdata%jo = j
        dlon = gplon(j+1) - gplon(j)
        yj = (glon - gplon(j)) / dlon

        ! Find K

        k = apexdata%ko
        if (gpalt(k+1) .lt. alt) then
            do while (gpalt(k+1) .lt. alt)
                k = k + 1
            enddo
        endif
        if (gpalt(k) .gt. alt) then
            do while (gpalt(k) .gt. alt)
                k = k - 1
            enddo
        endif
        apexdata%ko = k
       
        hti  = apexdata%re/(apexdata%re+alt)
        diht = apexdata%re/(apexdata%re+gpalt(k+1)) &
                     - apexdata%re/(apexdata%re+gpalt(k))
        zk  = (hti - apexdata%re/(apexdata%re+gpalt(k))) / diht
      
    !  For intrp:

        call this%trilin(x(i,j,k),nlat,nlon,xi,yj,zk,fx,dfxdn,dfxde,dfxdd)
        dfxdth = -dfxdn*apexdata%rtod/dlat
        dfxdln =  dfxde*apexdata%rtod/dlon
        dfxdh  = -hti*hti*dfxdd/(apexdata%re*diht)

        call this%trilin(y(i,j,k),nlat,nlon,xi,yj,zk,fy,dfydn,dfyde,dfydd)
        dfydth = -dfydn*apexdata%rtod/dlat
        dfydln =  dfyde*apexdata%rtod/dlon
        dfydh  = -hti*hti*dfydd/(apexdata%re*diht)

        call this%trilin(z(i,j,k),nlat,nlon,xi,yj,zk,fz,dfzdn,dfzde,dfzdd)
        dfzdth = -dfzdn*apexdata%rtod/dlat
        dfzdln =  dfzde*apexdata%rtod/dlon
        dfzdh  = -hti*hti*dfzdd/(apexdata%re*diht)

        call this%trilin(v(i,j,k),nlat,nlon,xi,yj,zk,fv,dfvdn,dfvde,dfvdd)
        dfvdth = -dfvdn*apexdata%rtod/dlat
        dfvdln =  dfvde*apexdata%rtod/dlon
        dfvdh  = -hti*hti*dfvdd/(apexdata%re*diht)

    ! Improve calculation of longitudinal derivatives near poles

        if (glat .lt. dlat - 90.0) then
            fac = 0.5*xi
            omfac = 1.0 - fac
            xi = xi - 1.0
            i = i + 1
            call this%trilin(x(i,j,k),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
            dfxdln = dfxdln*omfac + fac*dmdfde*apexdata%rtod/dlon
            call this%trilin(y(i,j,k),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
            dfydln = dfydln*omfac + fac*dmdfde*apexdata%rtod/dlon
            call this%trilin(v(i,j,k),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
            dfvdln = dfvdln*omfac + fac*dmdfde*apexdata%rtod/dlon
            return
        endif

        if (glat .gt. 90.0-dlat) then
            fac = 0.5*(1.0- xi)
            omfac = 1.0 - fac
            xi = xi + 1.0
            i = i - 1
            call this%trilin(x(i,j,k),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
            dfxdln = dfxdln*omfac + fac*dmdfde*apexdata%rtod/dlon
            call this%trilin(y(i,j,k),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
            dfydln = dfydln*omfac + fac*dmdfde*apexdata%rtod/dlon
            call this%trilin(v(i,j,k),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
            dfvdln = dfvdln*omfac + fac*dmdfde*apexdata%rtod/dlon
        endif
          
        return

    end subroutine intrp

    ! ============================================================================
    !  Trilinear interpolation of u and its derivatives
    ! 940803 A. D. Richmond
    ! Inputs:
    !   u(1,1,1) = address of lower corner of interpolation box 
    !   nlat = first dimension of u from calling routine
    !   nlon = second dimension of u from calling routine
    !   xi = fractional distance across box in x direction 
    !   yj = fractional distance across box in y direction 
    !   zk = fractional distance across box in z direction 
    ! Outputs:
    !   fu = interpolated value of u
    !   dfudx = interpolated derivative of u with respect to i (x direction)
    !   dfudy = interpolated derivative of u with respect to j (y direction)
    !   dfudz = interpolated derivative of u with respect to k (z direction)
    ! ============================================================================

    subroutine trilin(this,u,nlat,nlon,xi,yj,zk,fu,dfudx,dfudy,dfudz)

        implicit none

        class(Hwm07_class)  :: this

        integer(4) :: nlat,nlon
        real(4)    :: u(nlat,nlon,2)

        real(4)    :: xi,yj,zk
        real(4)    :: fu
        real(4)    :: dfudx,dfudy,dfudz
        real(4)    :: omxi,omyj,omzk

        omxi = 1.0 - xi
        omyj = 1.0 - yj
        omzk = 1.0 - zk

        fu = u(1,1,1)*omxi*omyj*omzk &
            + u(2,1,1)*xi*omyj*omzk &
            + u(1,2,1)*omxi*yj*omzk &
            + u(1,1,2)*omxi*omyj*zk &
            + u(2,2,1)*xi*yj*omzk &
            + u(2,1,2)*xi*omyj*zk &
            + u(1,2,2)*omxi*yj*zk &
            + u(2,2,2)*xi*yj*zk

        dfudx = (u(2,1,1)-u(1,1,1))*omyj*omzk + (u(2,2,1)-u(1,2,1))*yj*omzk &
            + (u(2,1,2)-u(1,1,2))*omyj*zk &
            + (u(2,2,2)-u(1,2,2))*yj*zk
            
        dfudy = (u(1,2,1)-u(1,1,1))*omxi*omzk &
            + (u(2,2,1)-u(2,1,1))*xi*omzk &
            + (u(1,2,2)-u(1,1,2))*omxi*zk &
            + (u(2,2,2)-u(2,1,2))*xi*zk
            
        dfudz = (u(1,1,2)-u(1,1,1))*omxi*omyj &
            + (u(2,1,2)-u(2,1,1))*xi*omyj &
            + (u(1,2,2)-u(1,2,1))*omxi*yj &
            + (u(2,2,2)-u(2,2,1))*xi*yj
            
        return
          
    end subroutine trilin

    ! =============================================================================
    !C  v is used for vr2n
    !C  Add-back of pseudodipole component to x,y,z,v and their derivatives.
    !C  940715 A. D. Richmond
    !C Inputs:
    !C   glat = latitude, degrees
    !C   glon = longitude, degrees
    !C   fx = interpolated value of x
    !C   fy = interpolated value of y
    !C   fz = interpolated value of z
    !C   fv = interpolated value of v
    !C   dfxdth,dfydth,dfzdth,dfvdth = derivatives of x,y,z,v with respect to 
    !C	colatitude, in radians-1
    !C   dfxdln,dfydln,dfzdln,dfvdln = derivatives of x,y,z,v with respect to 
    !C	longitude, in radians-1
    !C Output:
    !C   cth,sth = cos(colatitude), sin(colatitude)
    !C   fx = interpolated value of x
    !C   fy = interpolated value of y
    !C   fz = interpolated value of z
    !C   fv = interpolated value of v
    !C   dfxdth,dfydth,dfzdth,dfvdth = derivatives of x,y,z,v with respect to 
    !C	colatitude, in radians-1
    !C   dfxdln,dfydln,dfzdln,dfvdln = derivatives of x,y,z,v with respect to 
    !C	longitude, in radians-1
    ! ==============================================================================

    subroutine adpl(this,glat,glon,cth,sth,fx,fy,fz,fv, &
        dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln)

        implicit none

        class(Hwm07_class)  :: this
        real(4)    :: glat,glon
        real(4)    :: cth,sth
        real(4)    :: fx,fy,fz,fv
        real(4)    :: dfxdth,dfydth,dfzdth,dfvdth
        real(4)    :: dfxdln,dfydln,dfzdln,dfvdln
        
        real(4)    :: cph,sph
        real(4)    :: ctm

        cph = cos((glon-apexdata%elon)*apexdata%dtor)
        sph = sin((glon-apexdata%elon)*apexdata%dtor)
        cth = sin(glat*apexdata%dtor)
        sth = cos(glat*apexdata%dtor)
        ctm = apexdata%ctp*cth + apexdata%stp*sth*cph

        fx = fx + sth*apexdata%ctp*cph - cth*apexdata%stp
        fy = fy + sth*sph
        fz = fz + ctm
        fv = fv - ctm

        dfxdth = dfxdth + apexdata%ctp*cth*cph + apexdata%stp*sth
        dfydth = dfydth + cth*sph
        dfzdth = dfzdth - apexdata%ctp*sth + apexdata%stp*cth*cph
        dfvdth = dfvdth + apexdata%ctp*sth - apexdata%stp*cth*cph
        dfxdln = dfxdln - apexdata%ctp*sth*sph
        dfydln = dfydln + sth*cph
        dfzdln = dfzdln - apexdata%stp*sth*sph
        dfvdln = dfvdln + apexdata%stp*sth*sph

        return
    end subroutine adpl

    !*******************************************************************************
    !          Calculates east, north, and up components of gradients of x,y,z,v in
    !          geodetic coordinates.  All gradients are in inverse km.  Assumes
    !          flatness of 1/298.25 and equatorial radius (REQ) of 6378.16 km.
    !         940803 A. D. Richmond
    !
    !
    !          40680925. = req**2 (rounded off)
    !          272340.   = req**2 * e2, where e2 = (2. - 1./298.25)/298.25
    !                      is square of eccentricity of ellipsoid.
    !
    ! ==============================================================================

    subroutine gradxyzv(this,alt,cth,sth, &
                dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln, &
                dfxdh,dfydh,dfzdh,dfvdh,gradx,grady,gradz,gradv)

        implicit none  

        class(Hwm07_class)  :: this
        real(4)    :: alt,cth,sth
        real(4)    :: dfxdth,dfydth,dfzdth,dfvdth
        real(4)    :: dfxdln,dfydln,dfzdln,dfvdln
        real(4)    :: dfxdh,dfydh,dfzdh,dfvdh
        real(4)    :: gradx(3),grady(3),gradz(3),gradv(3)

        real(4)    :: d2,d,rho
        real(4)    :: dddthod,drhodth,dzetdth,ddisdth

        d2 = 40680925.e0 - 272340.e0*cth*cth
        d = sqrt(d2)
        rho = sth*(alt + 40680925.e0/d)
        dddthod = 272340.e0*cth*sth/d2
        drhodth = alt*cth + (40680925.e0/d)*(cth-sth*dddthod)
        dzetdth =-alt*sth - (40408585.e0/d)*(sth+cth*dddthod)
        ddisdth = sqrt(drhodth*drhodth + dzetdth*dzetdth)
        
        gradx(1) = dfxdln/rho
        grady(1) = dfydln/rho
        gradz(1) = dfzdln/rho
        gradv(1) = dfvdln/rho

        gradx(2) = -dfxdth/ddisdth
        grady(2) = -dfydth/ddisdth
        gradz(2) = -dfzdth/ddisdth
        gradv(2) = -dfvdth/ddisdth
        
        gradx(3) = dfxdh
        grady(3) = dfydh
        gradz(3) = dfzdh
        gradv(3) = dfvdh

        return
        
    end subroutine gradxyzv

    !*******************************************************************************
    !          Calculates east, north, and up components of gradients of x,y,z,v in
    !          geodetic coordinates.  All gradients are in inverse km.  Assumes
    !          flatness of 1/298.25 and equatorial radius (REQ) of 6378.16 km.
    !         940803 A. D. Richmond
    !
    !          40680925. = req**2 (rounded off)
    !          272340.   = req**2 * e2, where e2 = (2. - 1./298.25)/298.25
    !                      is square of eccentricity of ellipsoid.
    !
    ! ==============================================================================

    subroutine grapxyzv(this,alt,cth,sth, &
            dfxdln,dfydln,dfzdln,dfvdln,gradx,grady,gradz,gradv)

        implicit none

        class(Hwm07_class)  :: this
        real(4)    :: alt,cth,sth
        real(4)    :: dfxdln,dfydln,dfzdln,dfvdln
        real(4)    :: gradx(3),grady(3),gradz(3),gradv(3)

        real(4)    :: d2,d,rho
        real(4)    :: dddthod,drhodth,dzetdth,ddisdth


        d2 = 40680925.e0 - 272340.e0*cth*cth  
        d = sqrt(d2)
        rho = sth*(alt + 40680925.e0/d)
        dddthod = 272340.e0*cth*sth/d2
        drhodth = alt*cth + (40680925.e0/d)*(cth-sth*dddthod)
        dzetdth =-alt*sth - (40408585.e0/d)*(sth+cth*dddthod)
        ddisdth = sqrt(drhodth*drhodth + dzetdth*dzetdth)
        
        gradx(1) = dfxdln/rho
        grady(1) = dfydln/rho
        gradz(1) = dfzdln/rho
        gradv(1) = dfvdln/rho

        return

    end subroutine grapxyzv
          
    !*******************************************************************************
    !          Uses gradients of x,y,z,v to compute geomagnetic field and
    !          gradients of apex latitude, longitude.
    !          940819 A. D. Richmond
    !          INPUT:
    !            HR     = reference altitude
    !            ALT    = altitude
    !            FX,FY,FZ,FV = interpolated values of x,y,z,v, plus pseudodipole
    !                     component
    !            GRADX,GRADY,GRADZ,GRADV = interpolated gradients of x,y,z,v,
    !                     including pseudodipole components (east,north,up)
    !          OUTPUT:
    !            XLATM  = modified apex latitude (lambda_m), degrees
    !            XLONM  = apex longitude (phi_a), degrees
    !            VMP    = magnetic potential, in T.m.
    !            GRCLM  = grad(cos(lambda_m)), in km-1
    !            CLMGRP = cos(lambda_m)*grad(phi_a), in km-1
    !            QDLAT  = quasi-dipole latitude, degrees
    !            RGRLP  = (re + alt)*grad(lambda')
    !            B      = magnetic field, in nT
    !            CLM    = cos(lambda_m)
    !            R3_2   = ((re + alt)/(re + hr))**(3/2)
    ! ===========================================================================
          
    subroutine gradlpv(this,hr,alt,fx,fy,fz,fv,gradx,grady,gradz,gradv, &
                xlatm,xlonm,vmp,grclm,clmgrp,qdlat,rgrlp,b,clm,r3_2)

        implicit none

        class(Hwm07_class)  :: this
        real(4)    :: hr,alt
        real(4)    :: fx,fy,fz,fv
        real(4)    :: gradx(3),grady(3),gradz(3),gradv(3)
        real(4)    :: xlatm,xlonm,vmp
        real(4)    :: grclm(3),clmgrp(3),rgrlp(3)
        real(4)    :: b(3)
        real(4)    :: r3_2

        real(4)    :: rr,r,rn
        real(4)    :: sqrror
        
        real(4)    :: cpm,spm
        real(4)    :: bo
        real(4)    :: rn2
        real(4)    :: x2py2
        real(4)    :: xnorm
        real(4)    :: xlp,slp,clp
        real(4)    :: qdlat,clm
        real(4)    :: grclp

        integer(4) :: i

        !integer, save :: ncalled = 0
        !ncalled = ncalled + 1
        !write(*,*) ncalled, "hr = ", hr
        rr = apexdata%re + hr
        r  = apexdata%re + alt
        !write(*,*) "hr = ", hr
        !write(*,*) "alt = ", alt
        rn = r/apexdata%re
        sqrror = sqrt(rr/r)
        r3_2 = 1.0/sqrror/sqrror/sqrror
        xlonm = atan2(fy,fx)
        cpm = cos(xlonm)
        spm = sin(xlonm)
        xlonm = apexdata%rtod*xlonm
        bo = apexdata%vp*1.e6

    ! 1.E6 converts T to nT and km-1 to m-1.

        rn2 = rn*rn
        vmp = apexdata%vp*fv/rn2
        b(1) = -bo*gradv(1)/rn2
        b(2) = -bo*gradv(2)/rn2
        b(3) = -bo*(gradv(3) - 2.0*fv/r)/rn2

        x2py2 = fx*fx + fy*fy
        xnorm = sqrt(x2py2 + fz*fz)
        xlp = atan2(fz,sqrt(x2py2))
        slp = sin(xlp)
        clp = cos(xlp)
        qdlat = xlp*apexdata%rtod
        clm = sqrror*clp
        !write(*,*) "sq = ", sqrror, clm, clp

        if (clm .gt. 1.0) then
            !write(*,*) "sqrror = ", sqrror
            !write(*,*) "clp    = ", clp
            !write(*,*) "clm    = ", clm
            stop &
            'gradlpv, point lies below field line that peaks at reference height.'
                   
        end if

        xlatm = apexdata%rtod*acos(clm)
        
    !  If southern magnetic hemisphere, reverse sign of xlatm

        if (slp .lt. 0.0) xlatm = - xlatm

        do i = 1,3
            grclp = cpm*gradx(i) + spm*grady(i)
            rgrlp(i) = r*(clp*gradz(i) - slp*grclp)
            grclm(i) = sqrror*grclp
            clmgrp(i) = sqrror*(cpm*grady(i)-spm*gradx(i))
        enddo
          
        grclm(3) = grclm(3) - sqrror*clp/(2.0*r)

        return

    end subroutine gradlpv

    !*******************************************************************************
    !          Computes base vectors and other parameters for apex coordinates.
    !          Vector components:  east, north, up
    !          940801 A. D. Richmond
    !          Reference:
    !            Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex
    !            Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
    !          INPUTS:
    !            HR     = reference altitude
    !            XLATM  = modified apex latitude, degrees
    !            GRCLM  = grad(cos(lambda_m)), in km-1
    !            CLMGRP = cos(lambda_m)*grad(phi_a), in km-1
    !            RGRLP  = (re + altitude)*grad(lambda')
    !            B      = magnetic field, in nT
    !            CLM    = cos(lambda_m)
    !            R3_2   = ((re + altitude)/(re + hr))**(3/2)
    !          RETURNS:
    !            BMAG    = magnitude of magnetic field, in nT
    !            SIM     = sin(I_m) of article
    !            SI      = sin(I)
    !            F       = F of article
    !            D       = D of article
    !            W       = W of article
    !            BHAT    = unit vector along geomagnetic field direction
    !            D1...F2 = base vectors of article
    !
    !=============================================================================

    subroutine basvec(this,hr,xlatm,grclm,clmgrp,rgrlp,b,clm,r3_2, &
                bmag,sim,si,f,d,w,bhat,d1,d2,d3,e1,e2,e3,f1,f2)

        implicit none

        class(Hwm07_class)  :: this
        real(4)    :: hr,xlatm
        real(4)    :: grclm(3),clmgrp(3),rgrlp(3),b(3)
        real(4)    :: clm,r3_2
        real(4)    :: bmag
        real(4)    :: sim
        real(4)    :: si
        real(4)    :: f
        real(4)    :: d
        real(4)    :: w
        real(4)    :: bhat(3)
        real(4)    :: d1(3),d2(3),d3(3)
        real(4)    :: e1(3),e2(3),e3(3)
        real(4)    :: f1(2),f2(2)
        
        real(4)    :: rr
        real(4)    :: simoslm
        
        integer(4) :: i
        real(4)    :: d1db,d2db
        
        rr = apexdata%re + hr
        simoslm = 2.0/sqrt(4.0 - 3.0*clm*clm)
        sim = simoslm*sin(xlatm*apexdata%dtor)
        bmag = sqrt(b(1)*b(1) + b(2)*b(2) + b(3)*b(3))
        d1db = 0.0
        d2db = 0.0
        do i = 1,3
            bhat(i) = b(i)/bmag
            d1(i) = rr*clmgrp(i)
            d1db = d1db + d1(i)*bhat(i)
            d2(i) = rr*simoslm*grclm(i)
            d2db = d2db + d2(i)*bhat(i)
        enddo

    !   Ensure that d1,d2 are exactly perpendicular to B:

        do i = 1,3
            d1(i) = d1(i) - d1db*bhat(i)
            d2(i) = d2(i) - d2db*bhat(i)
        enddo  

        e3(1) = d1(2)*d2(3) - d1(3)*d2(2)
        e3(2) = d1(3)*d2(1) - d1(1)*d2(3)
        e3(3) = d1(1)*d2(2) - d1(2)*d2(1)
        d = bhat(1)*e3(1) + bhat(2)*e3(2) + bhat(3)*e3(3)

    ! Following step may be unnecessary, but it ensures that e3 lies along bhat.
        
        do i = 1,3
            d3(i) = bhat(i)/d
            e3(i) = bhat(i)*d
        enddo
       
        e1(1) = d2(2)*d3(3) - d2(3)*d3(2)
        e1(2) = d2(3)*d3(1) - d2(1)*d3(3)
        e1(3) = d2(1)*d3(2) - d2(2)*d3(1)
        e2(1) = d3(2)*d1(3) - d3(3)*d1(2)
        e2(2) = d3(3)*d1(1) - d3(1)*d1(3)
        e2(3) = d3(1)*d1(2) - d3(2)*d1(1)

        w = rr*rr*clm*abs(sim)/(bmag*d)

        si = -bhat(3)

        f1(1) =  rgrlp(2) 
        f1(2) = -rgrlp(1)
        f2(1) = -d1(2)*r3_2
        f2(2) =  d1(1)*r3_2
        f = f1(1)*f2(2) - f1(2)*f2(1)

        return

    end subroutine basvec


    !=========================================================================================
    !
    !
    !
    !
    !
    !=========================================================================================

    !
    !  Disturbance wind part of Horizontal Wind Model HWM07
    !  Version DWM07B104i
    !  See readme.txt file for detailed release notes.
    !
    !  AUTHOR
    !    John Emmert
    !    Space Science Division
    !    Naval Research Laboratory
    !    4555 Overlook Ave.
    !    Washington, DC 20375
    !
    !  Point of Contact
    !    msishwmhelp@nrl.navy.mil
    !
    !  DATE
    !    19 August 2008
    !
    !  REFERENCE
    !    Emmert, J. T., D. P. Drob, G. G. Shepherd, G. Hernandez, M. J. Jarvis, J. W. 
    !      Meriwether, R. J. Niciejewski, D. P. Sipler, and C. A. Tepley (2008),
    !      DWM07 global empirical model of upper thermospheric storm-induced 
    !      disturbance winds, J. Geophys Res., 113, doi:10.1029/2008JA013541.
    !


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !                       DWM-HWM Interface
    !
    !  Author: John Emmert
    !  Date: 7/10/2007
    !  Description: Using HWM inputs, computes Quasi-dipole latitude and local time,
    !               and Kp. Retrieves DWM results for these conditions, converts
    !               to geographic directions, and applies artificial height profile.
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine dwm07b_hwm_interface(this,IYD,SEC,ALT,GLAT,GLON,AP,DW)

        implicit none

        class(Hwm07_class)      :: this
        INTEGER,intent(in)      :: IYD
        REAL(4),intent(in)      :: SEC,ALT,GLAT,GLON
        REAL(4),intent(in)      :: AP(2)
        REAL(4),intent(out)     :: DW(2)

        real(4)                 :: ut, mlat, mlt, kp, mmpwind, mzpwind
        real(4)                 :: f1e, f1n, f2e, f2n
        real(4)                 :: dummy
        real(4)                 :: day, mlon, asun_glat, asun_glon, asun_mlat, asun_mlon

        real(8), parameter      :: pi=3.141592653590, dtor=pi/180D0, sin_eps=0.39781868

        !CONVERT AP TO KP
        kp = this%ap_to_kp(ap(2))

        !CONVERT GEO LAT/LON TO QD LAT/LON

        call this%gd2qd(glat,glon,mlat,mlon,f1e,f1n,f2e,f2n)

        !COMPUTE QD MAGNETIC LOCAL TIME (LOW-PRECISION)
        day = real(mod(iyd,1000))
        ut = sec / 3600.0
        asun_glat = -asin(sin((day-80.0)*real(dtor)) * sin_eps) / dtor
        asun_glon = -ut * 15.0
        call this%gd2qd(asun_glat, asun_glon, asun_mlat, asun_mlon, &
                   dummy,dummy,dummy,dummy)
        mlt = (mlon - asun_mlon) / 15.0

        !RETRIEVE DWM WINDS
        call this%dwm07b(mlt, mlat, kp, mmpwind, mzpwind)

        !CONVERT TO GEOGRAPHIC COORDINATES
        dw(1) = f2n*mmpwind + f1n*mzpwind
        dw(2) = f2e*mmpwind + f1e*mzpwind

        !APPLY HEIGHT PROFILE
        dw = dw * this%dwm_altwgt(alt)

        return

    end subroutine dwm07b_hwm_interface



    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !   This subroutine is used to evaluate DWM
    !
    !  Programming Notes:
    !
    !   This subroutine is only OPENMP/THREAD SAFE when no calls to
    !   loaddwm() are made.
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine dwm07b(this, mlt, mlat, kp, mmpwind, mzpwind)

        implicit none

        class(Hwm07_class)        :: this
        real(4),intent(in)        :: mlt            !Magnetic local time (hours)
        real(4),intent(in)        :: mlat           !Magnetic latitude (degrees)
        real(4),intent(in)        :: kp             !3-hour Kp

        real(4),intent(out)       :: mmpwind        !Mer. disturbance wind (+north, QD coordinates)
        real(4),intent(out)       :: mzpwind        !Zon. disturbance wind (+east, QD coordinates)

        ! Local variables
        real(4)                   :: mltdeg
        real(4)                   :: kp_terms(0:2)
        real(4)                   :: latwgt_terms
        real(4)                   :: termval_temp(0:1)
        integer(4)                :: iterm

        real(4), external         :: dwm_latwgt2

        
        !LOAD MODEL PARAMETERS IF NECESSARY (done by NEPTUNE!) VB
        !if (modelinit) call loaddwm(defaultdata)

        !COMPUTE VSH TERMS
        mltdeg = 15.0*mlt
        call this%vsh_basis(mlat, mltdeg)

        !COMPUTE KP TERMS
        call this%dwm_kpspl3_calc(kp, kp_terms)

        !COMPUTE LATITUDINAL WEIGHTING TERMS
        latwgt_terms = this%dwm_latwgt2(mlat, mlt, kp, dwmdata%twidth)

        !GENERATE COUPLED TERMS
        do iterm = 0, dwmdata%nterm-1
          termval_temp = 1
          if (dwmdata%termarr(0,iterm) .ne. 999) termval_temp = termval_temp * dwmdata%vsh_terms(0:1,dwmdata%termarr(0,iterm))
          if (dwmdata%termarr(1,iterm) .ne. 999) termval_temp = termval_temp * kp_terms(dwmdata%termarr(1,iterm))
          if (dwmdata%termarr(2,iterm) .ne. 999) termval_temp = termval_temp * latwgt_terms
          dwmdata%termval(0:1,iterm) = termval_temp
        end do

        !APPLY COEFFICIENTS
        mmpwind = dot_product(dwmdata%coeff, dwmdata%termval(0,0:dwmdata%nterm-1))
        mzpwind = dot_product(dwmdata%coeff, dwmdata%termval(1,0:dwmdata%nterm-1))

        return

    end subroutine dwm07b



    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! This subroutine loads the disturbance wind model parameters
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine loaddwm(this,datafile)
        use slam_io, only: openFile, closeFile, DIRECT, IN_UNFORMATTED

        implicit none

        class(Hwm07_class)           :: this
        character(len=*), intent(in) :: datafile

        integer :: ich    ! input channel to read from

        ich = openFile(datafile, DIRECT, IN_UNFORMATTED)

        if (allocated(dwmdata%termarr)) deallocate(dwmdata%termarr)
        if (allocated(dwmdata%termval)) deallocate(dwmdata%termval)
        if (allocated(dwmdata%coeff)) deallocate(dwmdata%coeff)
        read(ich) dwmdata%nterm, dwmdata%mmax, dwmdata%lmax
        allocate(dwmdata%termarr(0:2, 0:dwmdata%nterm-1))
        read(ich) dwmdata%termarr
        allocate(dwmdata%coeff(0:dwmdata%nterm-1))
        allocate(dwmdata%termval(0:1, 0:dwmdata%nterm-1))
        read(ich) dwmdata%coeff
        read(ich) dwmdata%twidth

        ich = closeFile(ich)

        call this%vsh_basis_init

        dwmdata%modelinit = .false.

        return

    end subroutine loaddwm

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! This subroutine initializes the VSH basis recursion coefficients
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine vsh_basis_init(this)

        !use dwm_module, only: vsh_terms
        !use vsh_basis_module

        implicit none

        class(Hwm07_class)  :: this
        integer(4)       :: i, j, n, m, ifn
        type (fn_map)    :: fn_map0(0:2*2*(dwmdata%mmax+1)*(dwmdata%lmax+1))


        !CREATE VSH FUNCTION MAP
        ifn = -1
        do n = 0, dwmdata%lmax
          do m = 0, dwmdata%mmax
            if ((m .eq. 0) .and. (n .eq. 0)) cycle
            if (m .gt. n) cycle
            do j = 0, 1
              do i = 0, 1
                if ((m .eq. 0) .and. (i .eq. 1)) cycle
                ifn = ifn + 1
                fn_map0(ifn)%realimag     = i
                fn_map0(ifn)%irrotational = j
                fn_map0(ifn)%M            = m
                fn_map0(ifn)%N            = n
              enddo
            enddo
          enddo
        enddo
        dwmdata%nvshfn = ifn + 1
        if (allocated(vshbasisdata%fn_map1)) deallocate(vshbasisdata%fn_map1)
        allocate(vshbasisdata%fn_map1(0:dwmdata%nvshfn-1))
        do ifn = 0, dwmdata%nvshfn-1
          vshbasisdata%fn_map1(ifn)%realimag     = fn_map0(ifn)%realimag
          vshbasisdata%fn_map1(ifn)%irrotational = fn_map0(ifn)%irrotational
          vshbasisdata%fn_map1(ifn)%M            = fn_map0(ifn)%M
          vshbasisdata%fn_map1(ifn)%N            = fn_map0(ifn)%N
        end do

        !CREATE ARRAY THAT WILL CONTAIN VSH BASIS VALUES
        if (allocated(dwmdata%vsh_terms)) deallocate(dwmdata%vsh_terms)
        allocate(dwmdata%vsh_terms(0:1, 0:dwmdata%nvshfn-1))

        !CREATE RECURSION COEFFICIENT ARRAYS
        if (allocated(vshbasisdata%A)) deallocate(vshbasisdata%A)
        if (allocated(vshbasisdata%B)) deallocate(vshbasisdata%B)
        if (allocated(vshbasisdata%anm)) deallocate(vshbasisdata%anm)
        if (allocated(vshbasisdata%bnm)) deallocate(vshbasisdata%bnm)
        if (allocated(vshbasisdata%fnm)) deallocate(vshbasisdata%fnm)
        if (allocated(vshbasisdata%cm)) deallocate(vshbasisdata%cm)
        if (allocated(vshbasisdata%cn)) deallocate(vshbasisdata%cn)
        if (allocated(vshbasisdata%e0n)) deallocate(vshbasisdata%e0n)
        if (allocated(vshbasisdata%m_arr)) deallocate(vshbasisdata%m_arr)
        if (allocated(vshbasisdata%n_arr)) deallocate(vshbasisdata%n_arr)
        if (allocated(vshbasisdata%cosmz)) deallocate(vshbasisdata%cosmz)
        if (allocated(vshbasisdata%sinmz)) deallocate(vshbasisdata%sinmz)

        allocate( vshbasisdata%A(0:dwmdata%mmax, 0:dwmdata%lmax) )
        allocate( vshbasisdata%B(0:dwmdata%mmax, 0:dwmdata%lmax) )
        allocate( vshbasisdata%anm(0:dwmdata%mmax, 0:dwmdata%lmax) )
        allocate( vshbasisdata%bnm(0:dwmdata%mmax, 0:dwmdata%lmax) )
        allocate( vshbasisdata%fnm(0:dwmdata%mmax, 0:dwmdata%lmax) )
        allocate( vshbasisdata%cm(0:dwmdata%mmax) )
        allocate( vshbasisdata%cn(0:dwmdata%lmax) )
        allocate( vshbasisdata%e0n(0:dwmdata%lmax) )
        allocate( vshbasisdata%m_arr(0:dwmdata%mmax) )
        allocate( vshbasisdata%n_arr(0:dwmdata%lmax) )
        allocate( vshbasisdata%cosmz(0:dwmdata%mmax) )
        allocate( vshbasisdata%sinmz(0:dwmdata%mmax) )
        vshbasisdata%A = 0
        vshbasisdata%B = 0

        !COMPUTE RECURSION COEFFICIENTS
        do m = 0, dwmdata%mmax
          vshbasisdata%cm(m) =  dsqrt(1 + 0.5/dble(max(m,1)))
          vshbasisdata%m_arr(m) = dble(m)
        end do
        do n = 0, dwmdata%lmax
          vshbasisdata%n_arr(n) = dble(n)
          vshbasisdata%cn(n) = 1 / dsqrt(dble(max(n,1)) * dble(n+1))
          vshbasisdata%e0n(n) = dsqrt(dble(n*(n+1)) / 2.0)
        end do
        vshbasisdata%anm = 0
        vshbasisdata%bnm = 0
        vshbasisdata%fnm = 0
        do m = 1, dwmdata%mmax
          if (m .eq. dwmdata%lmax) cycle
          do n = m+1, dwmdata%lmax
            vshbasisdata%anm(m,n) = dsqrt( dble((2*n-1)*(2*n+1)) / dble((n-m)*(n+m)) )
            vshbasisdata%bnm(m,n) = dsqrt( dble((2*n+1)*(n+m-1)*(n-m-1)) / dble((n-m)*(n+m)*(2*n-3)) )
            vshbasisdata%fnm(m,n) = dsqrt( dble((n-m)*(n+m)*(2*n+1)) / dble(2*n-1) )
          end do
        enddo

        return

    end subroutine vsh_basis_init

    !***************************************************************************************************
    !
    !JOHN EMMERT  3/26/07
    !EVALUATES A SERIES OF 2-COMPONENT (HORIZONTAL) VECTOR SPHERICAL HARMONIC FUNCTIONS, UP TO SPECIFIED
    !DEGREE AND ORDER. THE CONVENTION OF KILLEEN ET AL. [ADV. SPACE RES., 1987, P. 207]
    !IS USED, EXCEPT THAT NORMALIZED ASSOCIATED LEGENDRE FUNCTIONS ARE USED.
    !
    !7/7/07 TRANSLATED TO F90 FROM IDL VERSION. TRANSLATION IS LIMITED, AND ASSUMES THAT THE INPUT
    !ARE IN DEGREES, THE THETA ARGUMENT IS A LATITUDE, AND ONLY VALID (NON-ZERO) FUNCTIONS ARE TO BE 
    !RETURNED. MUST BE INITIALIZED WITH vsh_basis_init
    !
    !CALLING SEQUENCE
    !        Result = vsh_basis(theta,phi)
    !
    !ARGUMENTS
    !        theta        The latitudinal coordinate, specified in degrees.
    !        phi            The azimuthal coordinate, specified in degrees.
    !
    !RESULT
    !       A 2 x nfn array, where nfn is the number of valid functions in the basis.
    !        The two elements of the first dimension correspond to the meridional and zonal components
    !       (in that order) of the horizontal spherical harmonic vectors.
    !
    !ROUTINES USED
    !        None.
    subroutine vsh_basis(this,theta, phi)

        !use vsh_basis_module

        implicit none

        class(Hwm07_class) :: this
        real(4)            :: theta, phi
        real(8)            :: x, y, z, mz, norm_m
        integer(4)         :: n, m, ifn
        integer(4)         :: fn_id

        real(8), parameter :: pi=3.141592653590, dtor=pi/180D0

        !PROCESS INPUT ARGUMENTS

        !write(*,*) "THETA = "
        !write(*,*) theta
        !write(*,*) "DTOR = "
        !write(*,*) dtor
        x = cos((90. - dble(theta)) * dtor)
        y = dsqrt(1-x*x)                        !y = sin(theta)
        z = dble(phi) * dtor

        !CALCULATE Pnm / sin(theta)
        if (dwmdata%mmax .ge. 1) vshbasisdata%B(1,1) = dsqrt(3D0)
        do m = 2, dwmdata%mmax
          vshbasisdata%B(m,m) = y * vshbasisdata%cm(m) * vshbasisdata%B(m-1,m-1)
        enddo
        do m = 1, dwmdata%mmax
          do n = m+1, dwmdata%lmax
            vshbasisdata%B(m,n) = vshbasisdata%anm(m,n) * x * vshbasisdata%B(m,n-1) - vshbasisdata%bnm(m,n) * vshbasisdata%B(m,n-2)
          end do
        end do

        !CALCULATE d(Pnm) / d(theta)
        do m = 1, dwmdata%mmax
          do n = m, dwmdata%lmax
            vshbasisdata%A(m,n) = vshbasisdata%n_arr(n) * x * vshbasisdata%B(m,n) - vshbasisdata%fnm(m,n) * vshbasisdata%B(m,n-1)
          end do
        end do
        do n = 1, dwmdata%lmax
          vshbasisdata%A(0,n) = -vshbasisdata%e0n(n) * y * vshbasisdata%B(1,n)
        end do

        !CALCULATE m(Pnm) / sin(theta) AND APPLY SECTORAL NORMALIZATION
        do m = 0, dwmdata%mmax
          if (m .eq. 0) then
            norm_m = 1D0/sqrt(2D0)  !Holmes and Featherstone norm factor adjusted to match Swartztrauber
          else
            norm_m = 0.5D0          !Holmes and Featherstone norm factor adjusted to match Swartztrauber
          end if
          do n = m, dwmdata%lmax
            vshbasisdata%B(m,n) = vshbasisdata%B(m,n) * vshbasisdata%m_arr(m) * vshbasisdata%cn(n) * norm_m
            vshbasisdata%A(m,n) = vshbasisdata%A(m,n) *            vshbasisdata%cn(n) * norm_m
          end do
        end do

        !CALCULATE VECTOR SPHERICAL HARMONIC FUNCTIONS

        do m = 0, dwmdata%mmax
          mz = dble(m)*z
          vshbasisdata%cosmz(m) = cos(mz)
          vshbasisdata%sinmz(m) = sin(mz)
        end do
        do ifn = 0, dwmdata%nvshfn-1
          m = vshbasisdata%fn_map1(ifn)%M
          n = vshbasisdata%fn_map1(ifn)%N
          fn_id = vshbasisdata%fn_map1(ifn)%realimag + 2*vshbasisdata%fn_map1(ifn)%irrotational
          select case (fn_id)
          case (0)
            dwmdata%vsh_terms(0,ifn) = real( -vshbasisdata%A(m,n) * vshbasisdata%cosmz(m) )  !Multiplies real  pt of irrot. coeff. (b)
            dwmdata%vsh_terms(1,ifn) = real( -vshbasisdata%B(m,n) * vshbasisdata%sinmz(m) )  !Multiplies real  pt of irrot. coeff. (b)
          case (1)
            dwmdata%vsh_terms(0,ifn) = real(  vshbasisdata%A(m,n) * vshbasisdata%sinmz(m) )  !Multiplies imag. pt of irrot. coeff. (b)
            dwmdata%vsh_terms(1,ifn) = real( -vshbasisdata%B(m,n) * vshbasisdata%cosmz(m) )  !Multiplies imag. pt of irrot. coeff. (b)
          case (2)
            dwmdata%vsh_terms(0,ifn) = real(  vshbasisdata%B(m,n) * vshbasisdata%sinmz(m) )  !Multiplies real  pt of solen. coeff. (c)
            dwmdata%vsh_terms(1,ifn) = real( -vshbasisdata%A(m,n) * vshbasisdata%cosmz(m) )  !Multiplies real  pt of solen. coeff. (c)
          case (3)
            dwmdata%vsh_terms(0,ifn) = real(  vshbasisdata%B(m,n) * vshbasisdata%cosmz(m) )  !Multiplies imag. pt of solen. coeff. (c)
            dwmdata%vsh_terms(1,ifn) = real(  vshbasisdata%A(m,n) * vshbasisdata%sinmz(m) )  !Multiplies imag. pt of solen. coeff. (c)
          end select
        end do

        return

    end subroutine vsh_basis



    !***************************************************************************************************
    !
    !JOHN EMMERT 5/3/07
    !EVALUATES A BASIS OF QUADRATIC KP SPLINES WITH NODES AT [-10,-8,0,2,5,8,18,20]
    !AND SUMS THE FIRST TWO AND LAST TWO BASIS FUNCTIONS. THIS IS EQUIVALENT TO IMPOSING
    !A CONSTRAINT OF ZERO SLOPE AT Kp=0 and Kp=8. THE FUNCTION RETURNS THREE TERMS. INPUT
    !Kp VALUES GREATER THAN 8 ARE TRUNCATED TO 8.
    !
    !TRANSLATED TO F90 FROM IDL VERSION, 7/7/07
    !
    !CALLING SEQUENCE
    !   Result = CALL DWM_KPSPL3_CALC( kp, dwm_kpspl3 )
    !
    !ARGUMENTS
    !   kp          Kp index (0-8)
    !   dwm_kpspl3    A 3-element real array containing the resulting basis functions.
    !
    !ROUTINES USED
    !    bspline_calc

    subroutine dwm_kpspl3_calc(this, kp0, dwm_kpspl3)

        implicit none

        class(Hwm07_class)        :: this
        real(4), intent(in)       :: kp0
        real(4), intent(out)      :: dwm_kpspl3(0:2)

        real(4)                   :: kp
        real(4)                   :: kpspl(0:4)

        kp = max(kp0, 0.)
        kp = min(kp,  8.)
        call this%bspline_calc(8,  kp, (/-10., -8., 0., 2., 5., 8., 18., 20./), kpspl, 2, 0)
        dwm_kpspl3(0) = kpspl(0) + kpspl(1)
        dwm_kpspl3(1) = kpspl(2)
        dwm_kpspl3(2) = kpspl(3) + kpspl(4)

        return

    end subroutine dwm_kpspl3_calc

    !***************************************************************************************************

    !JOHN EMMERT 5/4/07
    !COMPUTES A LATITUDE DEPENDENT WEIGHTING FUNCTION THAT GOES TO
    !ZERO AT LOW LATITUDES AND ONE AT HIGH LATITUDES. THE TRANSITION
    !IS AN EXPONENTIAL S-CURVE, WITH THE TRANSITION LATITUDE DETERMINED
    !BY THE MLT/Kp MODEL GENERATED BY dwm_maxgrad1_02.pro, AND A
    !TRANSITION WITH SPECIFIED BY THE CALLING ROUTINE.
    !
    !TRANSLATED TO F90 FROM IDL VERSION, 7/7/07
    !
    !
    !CALLING SEQUENCE
    !   Result = DWM_LATWGT2(mlat, mlt, kp, twidth)
    !
    !ARGUMENTS
    !     mlt:       An n-element array of magnetic local times (hours)
    !     mlat:      Magnetic latitude (scalar value)
    !     kp:        3-hour Kp index
    !     twidth:    Latitude transition width
    !
    !ROUTINES USED
    !    None


    function dwm_latwgt2(this, mlat, mlt, kp0, twidth)

        implicit none

        class(Hwm07_class)        :: this
        real(4)                   :: dwm_latwgt2
        real(4)                   :: mlat, mlt, kp0, kp, twidth
        real(4)                   :: mltrad, sinmlt, cosmlt, tlat

        real(4), parameter :: coeff(0:5) = (/ 65.7633,  -4.60256,  -3.53915,  &
                                             -1.99971,  -0.752193,  0.972388 /)
        real(4), parameter :: pi=3.141592653590, dtor=pi/180.0


        mltrad = mlt * 15.0 * dtor
        sinmlt = sin(mltrad)
        cosmlt = cos(mltrad)
        kp = max(kp0, 0.)
        kp = min(kp,  8.)
        tlat = coeff(0) + coeff(1)*cosmlt + coeff(2)*sinmlt +   &
               kp*(coeff(3) + coeff(4)*cosmlt + coeff(5)*sinmlt)
        dwm_latwgt2 = 1.0 / ( 1 + exp(-(abs(mlat)-tlat)/twidth) )

        return

    end function dwm_latwgt2



    !******************************************************************************
    !******************************************************************************
    !
    !BSPLINE
    !JOHN EMMERT 3/31/05
    !TRANSLATED TO FORTRAN-90 10/4/06. FORTRAN VERSION ONLY ALLOWS SCALAR X ARGUMENT
    !CALCULATES A BASIS OF B-SPLINES OF SPECIFIED ORDER
    !USES DE BOOR'S RECURSION RELATION
    !DeBoor, C. A. (1978), A practical guide to splines, Appl. Math. Sci., 27,
    !  Springer, New York.
    !
    !CALLING SEQUENCE:   CALL BSPLINE_CALC(m, x, nodes, bspline, order, periodic)
    !
    !ARGUMENTS
    !     m:         The number of elements in the 'nodes' argument array (integer).
    !     x:         The dependent variable location at which the spline basis is
    !                to be evaluated (real).
    !     nodes:     Vector of ordered node positions (real array).
    !
    !                For a non-periodic basis:
    !                  Add k-1 (=the value of the order keyword, below) nodes to
    !                  either end of the range (domain) of the data. The extra
    !                  nodes on the left are the starting points of splines that
    !                  end inside the domain.  The extra nodes on the right are the
    !                  end points of splines that start inside the domain. The
    !                  splines actually used in the basis are those that begin at
    !                  the first m-k nodes. The sum of splines is then one (i.e.,
    !                  normalized) over the domain.
    !
    !                For a periodic basis:
    !                  The last node is identified with the first, so that there are
    !                  m-1 splines in the basis. For example, for a periodic basis
    !                  of 6 evenly spaced splines over the interval [0,24], with the
    !                  first node at 2, use the following array of 7 nodes:
    !                  [2,6,10,14,18,22,26]
    !
    !     order:     Set this integer argument to the desired spline order (=k-1
    !                in De Boor's notation).
    !                  Simple bins:           order=0
    !                  Linear interpolation:  order=1
    !                  Quadratic splines:     order=2
    !                  Cubic splines:         order=3 (default)
    !                  etc.
    !                Note that the number of nodes must be greater than k.
    !
    !     periodic:  Set this integer argument to 1 if a basis of periodic
    !                splines is desired; 0 otherwise.

    !RESULT
    !     The function returns bspline, a m1-element array, where m1 = m-k in the
    !     case of a non-periodic basis and m-1 in the case of a periodic basis. The
    !     elements of the array represent each of the basis splines.

    !EXAMPLE
    !     x = 18.0
    !     nodes = (/1.0, 3.5, 6.0, 18.0, 20.5, 23.0, 25.0/)
    !     m = size(nodes)
    !     order = 3
    !     periodic = 1
    !     call bspline_calc(m, x, nodes, bspline, order, periodic)
    !     print *, bspline
    !       0.025355    0.390467    0.584179    0.000000    0.000000    0.000000
    !
    !ROUTINES USED THAT ARE NOT IN THE STANDARD FORTRAN-90 LIBRARY
    !     pershift
    !
    subroutine bspline_calc(this, nnode0, x0, node0, bspline, order, periodic)

      implicit none

      class(Hwm07_class)     :: this
      integer(4), intent(in) :: nnode0, order, periodic
      real(4), intent(in)    :: x0, node0(0:nnode0-1)
      integer(4)             :: i, j, k, nnode, nspl
      real(4)                :: x, perint(0:1), perspan, dx1, dx2
      real(4), allocatable   :: node(:), node1(:)
      real(4), allocatable   :: bspline0(:)
      real(4), intent(out)   :: bspline(0:nnode0-1-periodic-(1-periodic)*(order+1))

      !PREPARE WORKING ARRAYS
      k = order + 1
      nnode=nnode0 + periodic*order
      nspl = nnode - k
      x = x0
      allocate(node(0:nnode-1), node1(0:nnode-1))
      allocate(bspline0(0:nnode-2))
      if (periodic .eq. 1) then
        perint = (/ node0(0), node0(nnode0-1) /)
        perspan = perint(1) - perint(0)
        x = this%pershift(x, perint)
        node = (/ node0, node0(1:order) + perspan /)
        do i=0, nnode0+order-1
          node1(i) = this%pershift(node(i), perint)
        end do
      else
        node = node0
        node1 = node
      end if

      !COMPUTE SPLINES
      do i = 0, nnode-2
        bspline0(i) = 0.0
        if (node1(i+1) .gt. node1(i)) then
          if ((x .ge. node1(i)) .and. (x .lt. node1(i+1))) bspline0(i) = 1.0
        else
          if ((x .ge. node1(i)) .or.  (x .lt. node1(i+1))) bspline0(i) = 1.0
        end if
      end do
      do j = 2, k
        do i = 0, nnode-j-1
          dx1 = x - node1(i)
          dx2 = node1(i+j) - x
          if (periodic .eq. 1) then
            if (dx1 .lt. 0) dx1 = dx1 + perspan
            if (dx2 .lt. 0) dx2 = dx2 + perspan
          end if
          bspline0(i) =    bspline0(i)   * dx1 / (node(i+j-1) - node(i)) &
                         + bspline0(i+1) * dx2 / (node(i+j)   - node(i+1))
        end do
      end do
      bspline = bspline0(0:nspl-1)
      deallocate(node, node1, bspline0)

      return

    end subroutine bspline_calc



    !******************************************************************************
    !******************************************************************************
    !
    !PERSHIFT
    !JOHN EMMERT   9/12/03
    !TRANSLATED TO FORTRAN-90 10/4/06. FORTRAN VERSION ONLY ALLOWS SCALAR INPUTS
    !SHIFTS INPUT VALUES INTO A SPECIFIED PERIODIC INTERVAL
    !
    !CALLING SEQUENCE:   Result = PERSHIFT(x, range)
    !
    !ARGUMENTS
    !      x:        The value to be shifted
    !      perint:   2-element vector containing the start and end values
    !                of the desired periodic interval.  The periodicity is
    !                determined by the span of the range.
    !
    !ROUTINES USED THAT ARE NOT IN THE STANDARD FORTRAN-90 LIBRARY
    !      None
    function pershift(this,x, perint)

      class(Hwm07_class) :: this
      real(4), parameter :: tol=1e-4
      real(4)            :: x, perint(0:1)
      real(4)            :: a, span, offset, offset1, pershift

      pershift = x
      a = perint(0)
      span = perint(1) - perint(0)
      if (span .ne. 0) then
        offset = x-a
        offset1 = mod(offset,span)
        if (abs(offset1) .lt. tol) offset1 = 0
      endif
      pershift = a + offset1
      if ((offset .lt. 0) .and. (offset1 .ne. 0)) pershift = pershift + span

      return

    end function pershift


    !******************************************************************************
    !******************************************************************************
    !
    !AP_TO_KP
    !JOHN EMMERT   7/10/07
    !CONVERTS AP VALUES TO KP VALUES, VIA LINEAR INTERPOLATION ON THE LOOKUP
    !TABLE
    !
    !CALLING SEQUENCE:   Result = AP_TO_KP(ap)
    !
    !ARGUMENTS
    !      ap:       The ap index. Values <0 or >400 are truncated
    !                to 0 and 400, respectively
    !
    !ROUTINES USED THAT ARE NOT IN THE STANDARD FORTRAN-90 LIBRARY
    !      None
    function ap_to_kp(this, ap0)

      class(Hwm07_class) :: this

      real(4), parameter :: apgrid(0:27) = (/0.,2.,3.,4.,5.,6.,7.,9.,12.,15.,18., &
                                             22.,27.,32.,39.,48.,56.,67.,80.,94., &
                                           111.,132.,154.,179.,207.,236.,300.,400./)
      real(4), parameter :: kpgrid(0:27) = (/0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11., &
                                             12.,13.,14.,15.,16.,17.,18.,19.,20.,21., &
                                             22.,23.,24.,25.,26.,27./) / 3.0
      real(4)            :: ap0, ap, ap_to_kp
      integer(4)         :: i


      ap = ap0
      if (ap .lt. 0) ap=0
      if (ap .gt. 400) ap=400

      i = 1
      do while (ap .gt. apgrid(i))
        i = i + 1
      end do
      if (ap .eq. apgrid(i)) then
        ap_to_kp = kpgrid(i)
      else
        ap_to_kp = kpgrid(i-1) + (ap - apgrid(i-1)) / (3.0 * (apgrid(i) - apgrid(i-1)))
      end if

      return

    end function ap_to_kp

    !******************************************************************************
    !******************************************************************************
    !
    !DWM_ALTWGT
    !JOHN EMMERT   7/10/07
    !COMPUTES AN EXPONENTIAL STEP FUNCTION IN HEIGHT TO BE APPLIED TO DWM WINDS (WHICH
    !HAVE NO HEIGHT DEPENDENCE). THE FUNCTION CUTS OFF DISTURBANCE WINDS BELOW 125 KM.
    !
    !CALLING SEQUENCE:   Result = DWM_ALTWGT(alt)
    !
    !ARGUMENTS
    !      alt:      Height in km
    !
    !ROUTINES USED THAT ARE NOT IN THE STANDARD FORTRAN-90 LIBRARY
    !      None

    function dwm_altwgt(this, alt)

      class(Hwm07_class) :: this
      real(4), parameter :: talt=125.0, twidth=5.0
      real(4)            :: alt, dwm_altwgt

      dwm_altwgt = 1.0 / (1 + exp(-(alt - talt)/twidth))

      return

    end function dwm_altwgt

    !*********************************************************************

    subroutine gd2qd(this,glat,glon,qdlat,qdlon,f1e,f1n,f2e,f2n)

        implicit none

        class(Hwm07_class)       :: this

        real(4), intent(in)      :: glat,glon
        real(4), intent(out)     :: qdlat, qdlon
        real(4), intent(out)     :: f1e,f1n,f2e,f2n
        real(4)                  :: f1(2),f2(2)
        integer(4)               :: ist
        real(4), parameter       :: alt = 250.0
        real(4)                  :: hr  = 0.0

        call this%apex(glat,glon,alt,hr,qdlon,qdlat,f1,f2,ist)
                        
        if (ist .gt. 0) stop

        f1e = f1(1)
        f1n = f1(2)
        f2e = f2(1)
        f2n = f2(2)

        return

    end subroutine gd2qd

end module hwm07Class
