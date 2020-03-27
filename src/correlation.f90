!>-----------------------------------------------------------------------------------------------
!> @brief       Correlation matrix computation module
!!
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>VB: 22.07.2013 (initial design) </li>
!!                <li>CHK 16.11.2015 (updated to use libslam) </li>
!!                <li>VB: 18.11.2015 (fixed partial derivatives for jacobian computation) </li>
!!                <li>VB: 15.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>CHK: 11.07.2018 (updated to class)</li>
!!              </ul>
!!
!> @details     This module contains all routines and parameters for the computation
!!              of the correlation matrix.
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      correlation
!!
!!------------------------------------------------------------------------------------------------
module correlation

    use slam_astro,               only: EGM96, getEarthGravity, getEarthRadius, getEarthRotation, getOrbitalPeriod
    use slam_astro_conversions,   only: true2mean
    use averaging,                only: getMeanElements
    use neptune_error_handling,   only: E_GEO_INIT, E_CORR_INIT, E_CORR_SINGULAR, E_CORR_INTERPOLATION, setNeptuneError
    use slam_error_handling,      only: isControlled, hasToReturn, hasFailed, FATAL, WARNING, checkIn, checkOut
    use gravity,                  only: Gravity_class, COEFF_DC, COEFF_DS
    use slam_math,                only: deg2rad, rad2deg, redang, mag, pi, halfPi, twopi, eps15, eps6
    use slam_interpolation,       only: lagrange_interpolation
    use slam_orbit_types,         only: kepler_t, nonSingOrbEl_t, kepler2NonSing
    use slam_types,               only: dp

    real(dp), parameter :: ckaula = 1.2d-5                                      ! Kaula's constant (Kaula, 1966)
    real(dp), parameter :: step_degrees = 1.d0                                  ! step in uj and uk for the integral evaluation
    real(dp), parameter :: step_rad = step_degrees*deg2rad                      ! step size of bins in radiant

    integer, parameter  :: ncor = 180                                           ! number of auto-correlation function bins
    integer, parameter  :: maxDegree = 1000                                     ! max. degree of geopotential to take into account for the estimation of the correlation function
    integer,  parameter :: nbins = int(360.d0/step_degrees)                     ! number of bins for double integral evaluation

    ! The correlation model class (type defintion)
    type, public :: Correlation_class

        integer             :: nrevs                                            ! number of revolutions

        real(dp)            :: cosi                                             ! cos of inclination
        real(dp)            :: ctgi                                             ! co-tangens of inclination
        real(dp)            :: ome2                                             ! 1-e**2
        real(dp)            :: orbitPeriod                                      ! orbital period in days
        real(dp)            :: sigmaR                                           ! uncertainty in radial direction due to geopotential (sigma_r)
        real(dp)            :: sini                                             ! sin of inclination
        real(dp)            :: trackShift                                       ! groundtrack shift of subsequent orbits
        real(dp)            :: u0                                               ! initial argument of true latitude

        real(dp), dimension(0:ncor)             :: Kr                           ! correlation function array
        real(dp), dimension(:), allocatable     :: cSigmaSq                     ! array containing squared sum: [Sum_k(delta_C**2 + delta_S**2)]/(2n+1)
        real(dp), dimension(:,:,:), allocatable :: corrMat                      ! final correlation matrix containing the value of matrix Q for each revolution

        type(kepler_t) :: meanEls                                               ! mean orbital elements

        logical :: corrInitialized                                              ! initialization flag
        logical :: propagateNoise                                               ! correlation matrix integration flag

    contains

        !============================================
        !
        ! public methods
        !
        !------------------------
        procedure :: getCorrelationMatrix
        procedure :: getNoisePropagationFlag
        procedure :: initCorrelation
        procedure :: updateCorrelationMatrix
        procedure :: setNoisePropagationFlag
        procedure :: getUconv
        procedure :: getPhixB
        procedure :: getIntersections
        procedure :: getCombMatrix
        procedure :: getCorrFunction
        procedure :: getCorrelationTerm
        procedure :: getRevTransition
        procedure :: jacobi_ke2rv
        procedure :: jacobi_oe2rv

    end type

    ! Constructor
    interface Correlation_class
        module procedure constructor
    end interface Correlation_class

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
    type(Correlation_class) function constructor()

        constructor%corrInitialized = .false.                                   ! initialization flag
        constructor%propagateNoise  = .false.                                   ! correlation matrix integration flag
    end function constructor

!--------------------------------------------------------------------------------------------------
!
!> @anchor      initCorrelation
!!
!> @brief       Initialization of correlation module
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 29.07.2013 (initial design)    </li>
!!              </ul>
!!
!> @param[in]   gravity_model - the gravity model
!!
!> @details     This routine initializes the correlation module. In its current version, the SigmaN
!!              array is computed for later use in the computation of the auto-correlation function.
!!              A prior initialization of the geopotential (gravity) module is required, which
!!              reads the normalized coefficients of the geopotential.
!!
!!------------------------------------------------------------------------------------------------
    subroutine initCorrelation(this,gravity_model)

        class(Correlation_class)            :: this
        type(Gravity_class),intent(inout)   :: gravity_model
        character(len=*), parameter :: csubid = 'initCorrelation'
        integer  :: i,j     ! loop counter
        real(dp) :: dc2     ! sigma of harmonic c squared
        real(dp) :: ds2     ! sigma of harmonic s squared
        real(dp) :: summ    ! auxiliary

        if(isControlled()) then
          if(hasToReturn()) return
          call checkIn(csubid)
        end if

        !** check if geopotential has already been initialized
        !------------------------------------------------------
        if(.not. gravity_model%isInitializedGeopotential()) then
          call setNeptuneError(E_GEO_INIT, FATAL)
          return
        end if
        !-------------------------------------------------------

        !** allocate array which holds the sigmaN values for different n
        if (allocated(this%cSigmaSq)) deallocate(this%cSigmaSq)
        allocate(this%cSigmaSq(0:gravity_model%getGeopotentialDegree()))

        !** compute sigma_n squared
        !---------------------------------------------
        do i = 2, gravity_model%getGeopotentialDegree()
          summ = 0.d0
          do j = 0,i
            dc2  = gravity_model%getHarmonicCoefficient(COEFF_DC, i, j)
            ds2  = gravity_model%getHarmonicCoefficient(COEFF_DS, i, j)
            !dc2  = gravity_model%getCoefficientSigma('C',i,j)
            !ds2  = gravity_model%getCoefficientSigma('S',i,j)
            summ = summ + dc2*dc2 + ds2*ds2
          end do
          this%cSigmaSq(i) = summ/(2.d0*i + 1.d0)
          !write(*,*) sqrt(cSigmaSq(i))
        end do
        !-----------------------------------------------

        this%corrInitialized = .true.
        
        !** done!
        if(isControlled()) then
          call checkOut(csubid)
        end if
        
        return

    end subroutine initCorrelation

    !========================================================================
    !
    !> @anchor      getCorrFunction
    !!
    !> @brief       Returns the auto-correlation function in Earth fixed spherical coordinates
    !> @author      Vitali Braun
    !!
    !> @param[in]   gravity_model - the gravity model
    !> @param[in]   r     3-dim radius vector in GCRF
    !!
    !> @return      an array
    !!
    !> @date        <ul>
    !!                <li> 29.07.2013 (initial design)</li>
    !!              </ul>
    !!
    !> @details     The auto-correlation function is an array in 1 degree steps, its size being
    !!              dependent on the module parameter 'ncorDegrees'. 
    !!
    !> @todo        get correlation function also for eccentric orbits
    !!
    !------------------------------------------------------------------------
    function getCorrFunction(this,gravity_model,r)

        class(Correlation_class)            :: this
        real(dp), dimension(0:ncor)         :: getCorrFunction
        type(Gravity_class),intent(inout)   :: gravity_model
        real(dp), dimension(3), intent(in)  :: r

        character(len=*), parameter :: csubid = 'getCorrFunction'

        integer  :: i             ! loop counter
        integer  :: n             ! loop counter

        real(dp) :: cospsi        ! cos of angular distance between points r1 and r2
        real(dp) :: eta           ! auxiliary
        real(dp) :: naux          ! auxiliary
        real(dp) :: pl1,pl2,pl3   ! legendre polynomials required for the iteration procedure 
        real(dp) :: Reor2         ! (RE/r)**2
        real(dp) :: sumCorr       ! sum of auto-correlation terms for each n (degree)
        real(dp) :: sumSigmaR     ! auxiliary for the computation of sigma_r
        real(dp) :: xn            ! auxiliary variable to compute the (RE/r)**2n terms

        if(isControlled()) then
          if(hasToReturn()) return
          call checkIn(csubid)
        end if

        !** check if already initialized
        if(.not. this%corrInitialized) then
          call this%initCorrelation(gravity_model)
        end if

        Reor2 = (getEarthRadius()/mag(r))**2.d0 

        !================================================
        !
        ! Start main loop over 0...ncorDegrees degrees
        !
        !-----------------------------------------------

        sumSigmaR = 0.d0

        do i = 0, ncor

          cospsi  = cos(i*step_rad)   ! angular distance of points r1 and r2: psi
          sumCorr = 0.d0
          xn      = Reor2
          
          !** Compute first two Legendre polynomials Pn0
          pl1 = 1.5d0*cospsi*cospsi - 0.5d0
          pl2 = (2.5d0*cospsi*cospsi - 1.5d0)*cospsi

          !============================================
          !
          ! Start loop over degree n
          !
          !-------------------------------------
          do n = 2, maxDegree

            if(n == 2) then
              pl3 = pl1
            else if(n == 3) then
              pl3 = pl2
            else
              pl3 = gravity_model%getLegendrePolynomial(pl1, pl2, cospsi, n)
            end if

            xn = Reor2*xn   ! build powers of Reor2

            !** estimate errors in harmonic coefficients
            !--------------------------------------------
            if(n <= gravity_model%getGeopotentialDegree()) then
              eta = 2.d0*this%cSigmaSq(n)
            else
              eta = (ckaula/n/n)**2.d0
            end if

            naux = xn*(n + 1.d0)*(n + 1.d0)*(n + 0.5d0)

            if(i == 0) then ! sum up terms for sigma_r
              sumSigmaR = sumSigmaR + naux*eta
            end if

            sumCorr = sumCorr + naux*eta*pl3

            !** set new legendre polynomials for recursive computation if n > 3
            if(n > 3) then
              pl1 = pl2
              pl2 = pl3
            end if

          end do

          !if(i == 0) this%sigmaR = sqrt(sumSigmaR)

          getCorrFunction(i) = sumCorr/sumSigmaR

        end do

        this%sigmaR = sqrt(sumSigmaR)*getEarthGravity()/mag(r)**2.d0
        !write(*,*) "sigmaR = ", this%sigmaR

        !====================================================
        !
        ! FEATURE/REMOVE: Auto-correlation function output
        !
        !-------------------------------------------------
        !open(unit=111, file='corr.out')
        !do i=0,ncor
        !  write(111,*) i*step_degrees, getCorrFunction(i)
        !end do
        !close(111)

        !stop
        !--------------------------------------------------

        !** done!
        if(isControlled()) then
          call checkOut(csubid)
        end if

        return

    end function getCorrFunction

!========================================================================
!
!> @anchor      getCorrelationTerm
!!
!> @brief       Returns a term 'j-k' of the correlation double integral
!> @author      Vitali Braun
!!
!> @param[in]   revDiff   Difference j-k (j==k if '0') in revolution numbers
!!
!> @return      correlation terms as 6x6 matrix
!!
!> @date        <ul>
!!                <li> 30.07.2013 (initial design)</li>
!!              </ul>
!!
!> @details     The correlation matrix is computed via the evaluation of
!!              a double integral. This routine returns single terms of
!!              the double integral which has been transformed to a sum
!!              of discrete terms, containing the correlation for uj = uk
!!              as well as the correlation of intersection points in sub-
!!              sequent orbits (j .ne. k).
!!
!------------------------------------------------------------------------
  function getCorrelationTerm(this,revDiff)

    class(Correlation_class)            :: this
    integer, intent(in)                 :: revDiff

    real(dp), dimension(1:6,1:6) :: getCorrelationTerm
    character(len=*), parameter :: csubid = 'getCorrelationTerm'

    integer :: bidx      ! bin index for correlation function
    integer :: ij, ik    ! loop counter
    integer :: jidx,kidx ! indices to find arg of latitude in conversion array
    
    !real(dp), dimension(2,2) :: ArgLatIntersect   ! containing the argument of latitude for the 
                                                  ! two intersection points of subsequent orbits
    real(dp), dimension(6,6) :: combMat           ! combined matrix: Phi x B x P x B^T x Phi^T
    real(dp) :: coslat_j                          ! cos of latitude for index j
    real(dp) :: coslat_k                          ! cos of latitude for index k
    real(dp) :: cospsi                            ! cos of separation angle psi
    real(dp) :: dc2                               ! combination of sigmaR and integration variables duj and duk
    real(dp) :: deltaLamb0                        ! ground track shift for given revolution number difference
    real(dp) :: dtemp                             ! auxiliary
    real(dp), dimension(6,3) :: PhixB_j           ! multiplication of state error transition matrix and force matrix: Phi X B for index j
    real(dp), dimension(6,3) :: PhixB_k           ! multiplication of state error transition matrix and force matrix: Phi X B for index k
    real(dp) :: psi                               ! separation angle psi
    real(dp) :: sinlat_j                          ! sin of latitude for index j
    real(dp) :: sinlat_k                          ! sin of latitude for index k
    real(dp) :: sinpsi                            ! sin of separation angle psi
    real(dp), dimension(1:2,0:nbins-1,1:2) :: uconv ! conversion array to obtain geocentric longitude and
                                                  ! latitude for given argument of true latitude (initialized by routine getUconv() 
                                                  ! the size of the second dimension depends on the binning and has to account for
                                                  ! a complete revolution. First dimension is for uj (=1) and uk(=2), while
                                                  ! the third contains latitude (=1) and longitude (=2)
    real(dp) :: xx                                ! auxiliary (tells whether there is a correlation (<>0) or not (=0))
    
    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** initialize
    getCorrelationTerm = 0.d0
    
    deltaLamb0 = this%trackShift*revDiff
    deltaLamb0 = redang(deltaLamb0, 2, 1, .false.)  ! reduce to 0..2pi

    xx  = deltaLamb0*rad2deg 
    dc2 = (step_degrees*step_degrees*this%sigmaR*this%orbitPeriod/dble(nbins))**2.d0  ! sigmaR*duj/duj_dot*duk/duk_dot, in degrees**2
                                                                                ! orbitPeriod/'360' (-> nbins) is the time step related to step_degrees

    
    !==========================================================================
    !
    ! Compute intersection points for subsequent orbits
    ! based on idiff
    !
    !--------------------------------------------------------------------------
    !write(*,*) "getting intersection..."
    !ArgLatIntersect = getIntersections(deltaLamb0)
    ! 
    !write(*,*) ArgLatIntersect*rad2deg
    !--------------------------------------------------------------------------

    !==========================================================================
    !
    ! Get conversion array to transform longitude and latitude to argument
    ! of true latitude
    !
    !--------------------------------------------------------------------------
    !write(*,*) "getting uconv..."
    uconv = this%getUconv(revDiff)

!    do i=0,nbins
!      write(*,'(4(1x,f12.6))') uconv(1,i,1)*rad2deg, uconv(1,i,2)*rad2deg, uconv(2,i,1)*rad2deg, uconv(2,i,2)*rad2deg
!    end do

    !--------------------------------------------------------------------------

    !==========================================================================
    !
    ! Now start integration by summing up double integral
    !
    !--------------------------------------------------------------------------
    do ij = 0, nbins-1

      jidx = nint(this%u0) + ij
      if(jidx > (nbins - 1)) jidx = jidx - nbins

      coslat_j = cos(uconv(1,jidx,1))
      sinlat_j = sin(uconv(1,jidx,1))

      PhixB_j  = this%getPhixB(jidx*step_rad)

      do ik = 0, nbins-1

        kidx = nint(this%u0) + ik
        if(kidx > (nbins-1)) kidx = kidx - nbins
       
        coslat_k = cos(uconv(2,kidx,1))
        sinlat_k = sin(uconv(2,kidx,1))

        !** separation angle
        !----------------------------------------------------------------------------------
        cospsi   = coslat_j*coslat_k*cos(uconv(2,kidx,2)-uconv(1,jidx,2)) + sinlat_j*sinlat_k

        !if(revDiff == 1) then
        !  write(*,*) "uj,uk,cospsi = ", ij, ik, cospsi
        !  read(*,*)
        !end if

        if(cospsi > 0.d0) then

          dtemp = 1.d0 - cospsi*cospsi
          
          if(dtemp > 0.d0) then
            sinpsi = sqrt(dtemp)
          else 
            sinpsi = 0.d0
          end if

          if(abs(cospsi) > eps15) then
            psi = atan(sinpsi/cospsi)
          else
            psi = halfPi
          end if

          if(psi < 0.d0) psi = -psi

          !** get bin index
          bidx = nint(psi/step_rad)

          if(bidx <= ncor) then
            xx = this%Kr(bidx)
          else
            xx = 0.d0
          end if

        else
          xx = 0.d0
        end if
        !----------------------------------------------------------------------------------

        !==================================================================================
        !
        ! FEATURE/REMOVE: Output of correlation function for different revolutions
        !
        !-----------------------------------------------------------------------------
        !write(112,'(2(f7.1,x),e12.5e2)') dble(ij + twopi*rad2deg*revdiff), dble(ik), xx

        !** perform integration, if xx is not zero
        !----------------------------------------------------------
        if(xx /= 0.d0) then
          PhixB_k                 = this%getPhixB(kidx*step_rad)
          combMat                 = this%getCombMatrix(PhixB_j, PhixB_k)
          getCorrelationTerm(:,:) = getCorrelationTerm(:,:) + combMat(:,:)*xx*dc2
        end if

        !----------------------------------------------------------
      end do
      !write(112,*) ! REMOVE
    end do
    !--------------------------------------------------------------------------
    
    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if
 
    return

  end function getCorrelationTerm

    !====================================================================================================
    !
    !> @anchor      getCombMatrix
    !!
    !> @brief       Returns the matrix product of (SET x B)_j x Kr x (B^T x SET^T)_k.
    !> @author      Vitali Braun
    !!
    !> @param[in]   pj    SET x B for j component
    !> @param[in]   pk    SET x B for k component
    !!
    !> @return      a 6x6 matrix
    !!
    !! @date        <ul> 
    !!                <li> 01.08.2013 (initial design)</li>
    !!              </ul>
    !!
    !! @details     This function computes the combined matrix product from the already preprocessed (
    !!              'getPhixB') product of state error transition matrix Phi and force matrix B. 
    !!              Two of these products have to be provided for the components j and k within the
    !!              double integral.
    !!
    !----------------------------------------------------------------------------------------------------

    function getCombMatrix(this, pj, pk)

        class(Correlation_class)             :: this
        real(dp), dimension(6,3), intent(in) :: pj
        real(dp), dimension(6,3), intent(in) :: pk
        real(dp), dimension(6,6)             :: getCombMatrix

        getCombMatrix = 0.d0

        getCombMatrix(1,1) = pj(1,1)*pk(1,1) + 0.5d0*pj(1,2)*pk(1,2) + 0.5d0*pj(1,3)*pk(1,3)
        getCombMatrix(1,2) = pj(1,1)*pk(2,1) + 0.5d0*pj(1,2)*pk(2,2)
        getCombMatrix(1,3) = 0.5d0*pj(1,3)*pk(3,3)
        getCombMatrix(1,4) = 0.5d0*pj(1,3)*pk(4,3)
        getCombMatrix(1,5) = pj(1,1)*pk(5,1) + 0.5d0*pj(1,2)*pk(5,2) + 0.5d0*pj(1,3)*pk(5,3)
        getCombMatrix(1,6) = pj(1,1)*pk(6,1) + 0.5d0*pj(1,2)*pk(6,2) + 0.5d0*pj(1,3)*pk(6,3)
        getCombMatrix(2,1) = pj(2,1)*pk(1,1) + 0.5d0*pj(2,2)*pk(1,2)
        getCombMatrix(2,2) = pj(2,1)*pk(2,1) + 0.5d0*pj(2,2)*pk(2,2)
        getCombMatrix(2,5) = pj(2,1)*pk(5,1) + 0.5d0*pj(2,2)*pk(5,2)
        getCombMatrix(2,6) = pj(2,1)*pk(6,1) + 0.5d0*pj(2,2)*pk(6,2)
        getCombMatrix(3,1) = 0.5d0*pj(3,3)*pk(1,3)
        getCombMatrix(3,3) = 0.5d0*pj(3,3)*pk(3,3)
        getCombMatrix(3,4) = 0.5d0*pj(3,3)*pk(4,3)
        getCombMatrix(3,5) = 0.5d0*pj(3,3)*pk(5,3)
        getCombMatrix(3,6) = 0.5d0*pj(3,3)*pk(6,3)
        getCombMatrix(4,1) = 0.5d0*pj(4,3)*pk(1,3)
        getCombMatrix(4,3) = 0.5d0*pj(4,3)*pk(3,3)
        getCombMatrix(4,4) = 0.5d0*pj(4,3)*pk(4,3)
        getCombMatrix(4,5) = 0.5d0*pj(4,3)*pk(5,3)
        getCombMatrix(4,6) = 0.5d0*pj(4,3)*pk(6,3)
        getCombMatrix(5,1) = pj(5,1)*pk(1,1) + 0.5d0*pj(5,2)*pk(1,2) + 0.5d0*pj(5,3)*pk(1,3)
        getCombMatrix(5,2) = pj(5,1)*pk(2,1) + 0.5d0*pj(5,2)*pk(2,2) 
        getCombMatrix(5,3) = 0.5d0*pj(5,3)*pk(3,3)
        getCombMatrix(5,4) = 0.5d0*pj(5,3)*pk(4,3)
        getCombMatrix(5,5) = pj(5,1)*pk(5,1) + 0.5d0*pj(5,2)*pk(5,2) + 0.5d0*pj(5,3)*pk(5,3)
        getCombMatrix(5,6) = pj(5,1)*pk(6,1) + 0.5d0*pj(5,2)*pk(6,2) + 0.5d0*pj(5,3)*pk(6,3)
        getCombMatrix(6,1) = pj(6,1)*pk(1,1) + 0.5d0*pj(6,2)*pk(1,2) + 0.5d0*pj(6,3)*pk(1,3)
        getCombMatrix(6,2) = pj(6,1)*pk(2,1) + 0.5d0*pj(6,2)*pk(2,2) 
        getCombMatrix(6,3) = 0.5d0*pj(6,3)*pk(3,3)
        getCombMatrix(6,4) = 0.5d0*pj(6,3)*pk(4,3)
        getCombMatrix(6,5) = pj(6,1)*pk(5,1) + 0.5d0*pj(6,2)*pk(5,2) + 0.5d0*pj(6,3)*pk(5,3)
        getCombMatrix(6,6) = pj(6,1)*pk(6,1) + 0.5d0*pj(6,2)*pk(6,2) + 0.5d0*pj(6,3)*pk(6,3)

    end function getCombMatrix


    !====================================================================================================
    !
    !> @anchor      getIntersections
    !!
    !> @brief       Returns the intersection points (argument of true latitude) for subsequent orbits
    !> @author      Vitali Braun
    !!
    !> @param[in]   deltaLamb0    ground track shift in rad (w/o taking into account Earth's rotation)
    !!
    !> @return      2x2 matrix
    !!
    !> @date        <ul>
    !!                <li> 30.07.2013 (initial design)</li>
    !!              </ul>
    !!
    !> @details     This function computes the intersection points of subsequent orbits
    !!              with respect to the ground track. The result are two intersection points
    !!              which contain the argument of true latitude for each orbit. The
    !!              parameter revDifference gives the number of revolutions in between the two
    !!              orbits under consideration. If revDifference=0 it means that the same orbit is
    !!              considered, for revDifference=1 the intersection of the current orbit with the
    !!              following orbit is analysed, etc.
    !!
    !> @todo        compute intersections also for eccentric orbits
    !!
    !----------------------------------------------------------------------------------------------------
    function getIntersections(this,deltaLamb0)

        real(dp), dimension(1:2,1:2) :: getIntersections

        class(Correlation_class)    :: this
        real(dp), intent(in)        :: deltaLamb0

        real(dp) :: cosdL2       ! cos of deltaLambda/2
        real(dp) :: deltaLambda  ! delta in longitude
        real(dp) :: del          ! accounting for Earth's rotation
        real(dp) :: sindL2       ! sin of deltaLambda/2
        real(dp) :: u1a          ! arg of true latitude for first orbit in point 'a'
        real(dp) :: u1a_old = 10.d0 ! auxiliary for the iteration of u1a (arbitrary initial value of 10.0)

        del = 0.d0

        !===================================================================
        !
        ! Search for 'u1a' which is the argument of true latitude for the
        ! intersection along the first orbit. An iteration is required to
        ! account for the time the object travels from equator to the 
        ! the intersection point, which is a function of the intersection.
        !
        !-------------------------------------------------------------------
        do

          deltaLambda = deltaLamb0 + del

          sindL2 = sin(deltaLambda*0.5d0)
          cosdL2 = cos(deltaLambda*0.5d0)

          if(sindL2 < eps15) sindL2 = eps15 ! prevent division by zero

          u1a = atan(cosdL2/sindL2/this%cosi)

          if(cosdL2/this%cosi < 0.d0) u1a = u1a + pi
          
          del = -this%trackShift*(u1a/pi - 0.5d0)    ! <-- ADAPT here for eccentric orbits

          if(abs(u1a - u1a_old) < eps6) exit
          
          u1a_old = u1a

        end do
        !-------------------------------------------------------------------

        getIntersections(1,1) = u1a
        getIntersections(1,2) = pi - u1a ! u2a
        getIntersections(2,1) = u1a + pi ! u1b
        getIntersections(2,2) = getIntersections(1,2) + pi    ! u2b

        return

    end function getIntersections

    !====================================================================================================
    !
    !> @anchor      getPhixB
    !!
    !> @brief       Returns the matrix product of SET (Phi) and force matrix (B).
    !> @author      Vitali Braun
    !!
    !> @param[in]   argLat    argument of true latitude in radians
    !!
    !> @return      6x3 matrix
    !!
    !> @date        <ul>
    !!                <li> 01.08.2013 (initial design)</li>
    !!              </ul>
    !!
    !> @details     This function computes the matrix product of the state error transition matrix Phi
    !!              and the force matrix B. The coordinates used are non-singular for circular orbits,
    !!              however, they do not work for equatorial orbit in the current implementation. The
    !!              coordinates are (u,a,i,Om,l,h) with l = e*cos(om) and h = e*sin(om), according to
    !!              Nazarenko (2010).
    !!
    !> @todo        Check derivatives for state transition matrix
    !!
    !----------------------------------------------------------------------------------------------------
    function getPhixB(this,argLat)

        real(dp), dimension(6,3) :: getPhixB
        class(Correlation_class) :: this
        real(dp), intent(in)     :: argLat

        real(dp) :: b13     ! b(1,3)
        real(dp) :: b21     ! b(2,1)
        real(dp) :: b22     ! b(2,2)
        real(dp) :: b33     ! b(3,3)
        real(dp) :: b43     ! b(4,3)
        real(dp) :: b51     ! b(5,1)
        real(dp) :: b52     ! b(5,2)
        real(dp) :: b53     ! b(5,3)
        real(dp) :: b61     ! b(6,1)
        real(dp) :: b62     ! b(6,3)
        real(dp) :: b63     ! b(6,3)
        real(dp) :: cosu    ! cos of arg of true latitude
        real(dp) :: cosu0   ! cos of initial arg of true latitude
        real(dp) :: dtemp   ! auxiliary
        real(dp) :: h       ! variable e*sin(om)
        real(dp) :: l       ! variable e*cos(om)
        real(dp) :: M       ! mean anomaly
        real(dp) :: p       ! semi-latus rectum
        real(dp) :: r0      ! initial radius vector
        real(dp) :: rabs    ! orbit radius (km)
        real(dp) :: rop     ! r/p
        real(dp) :: sinu    ! sin of arg of true latitude
        real(dp) :: sinu0   ! sin of initial arg of true latitude
        real(dp) :: some2   ! sqrt(ome2)
        real(dp) :: u12     ! phi(1,2)
        real(dp) :: u15     ! phi(1,5)
        real(dp) :: u16     ! phi(1,6)
        real(dp) :: xx      ! auxiliary

        sinu  = sin(argLat)
        sinu0 = sin(this%meanEls%aop + this%meanEls%tran)
        cosu  = cos(argLat)
        cosu0 = cos(this%meanEls%aop + this%meanEls%tran)

        some2 = sqrt(this%ome2)
        h     = this%meanEls%ecc*sin(this%meanEls%aop)
        l     = this%meanEls%ecc*cos(this%meanEls%aop)
        p     = this%meanEls%sma*this%ome2
        r0    = p/(1.d0 + this%meanEls%ecc*cos(this%meanEls%tran))

        !** compute mean anomaly for parameter u12
        call true2mean(this%meanEls%ecc, argLat-this%meanEls%aop, dtemp, M)

        !** compute radius
        rabs = p/(1.d0 + this%meanEls%ecc*cos(argLat-this%meanEls%aop))
        rop  = rabs/p

        !======================================================
        !
        ! State transition matrix
        !
        !---------------------------------------
        u12  = -1.5d0*(this%meanEls%sma/rabs)**2.d0*some2/this%meanEls%sma*(M - this%meanEls%man)
        u15  = (sinu + h - (sinu0 + h)*(r0/rabs)**2.d0 + sinu/rop - sinu0*r0/(rop*rabs))/this%ome2
        u16  = (cosu + l - (cosu0 + l)*(r0/rabs)**2.d0 + cosu/rop - cosu0*r0/(rop*rabs))/this%ome2

        !=======================================================
        !
        ! Force matrix
        !
        !-----------------------------------------
        xx = this%orbitPeriod*some2/twopi/this%meanEls%sma

        b13  = -xx*this%ctgi*sinu*rop
        b21  = -this%orbitPeriod/pi/some2*(l*sinu - h*cosu)
        b22  = this%orbitPeriod/pi/some2*(1.d0 + l*cosu + h*sinu)
        b33  = xx*rop*cosu
        b43  = xx*rop*sinu/this%sini
        b51  = xx*sinu
        b52  = xx*((l + cosu)*rop + cosu)
        b53  = xx*h*this%ctgi*rop*sinu
        b61  = -xx*cosu
        b62  = xx*((h + sinu)*rop + sinu)
        b63  = -xx*l*this%ctgi*rop*sinu

        !========================================================
        !
        ! Multiplication
        !
        !------------------------------------
        getPhixB = 0.d0

        getPhixB(1,1) = u15*b51 + u16*b61
        getPhixB(1,2) = u12*b22 + u15*b52 + u16*b62
        getPhixB(1,3) = b13     + u15*b53 + u16*b63 
        getPhixB(2,1) = b21
        getPhixB(2,2) = b22
        getPhixB(3,3) = b33
        getPhixB(4,3) = b43
        getPhixB(5,1) = b51
        getPhixB(5,2) = b52
        getPhixB(5,3) = b53
        getPhixB(6,1) = b61
        getPhixB(6,2) = b62
        getPhixb(6,3) = b63

    end function getPhixB

    !====================================================================================================
    !
    !> @anchor      getUconv
    !!
    !> @brief       Returns the conversion between argument of true latitude and longitude/latitude for a given revolution number difference
    !> @author      Vitali Braun
    !!
    !> @param[in]   revDiff    Difference j-k (j==k if '0') in revolution numbers
    !!
    !> @return      multi-dimensional matrix
    !!
    !> @date        <ul>
    !!                <li> 31.07.2013 (initial design)</li>
    !!              </ul>
    !!
    !> @details     This function converts from argument of true latitude to longitude and latitude,
    !!              as the integration variable is the argument of true latitude but the correlation function
    !!              requires longitude and latitude in order to compute the angular separation of two points.
    !!
    !> @todo        compute conversion for non-spherical Earth
    !> @todo        do something for equatorial orbits
    !!
    !----------------------------------------------------------------------------------------------------
    function getUconv(this,revDiff)

        class(Correlation_class)            :: this
        integer, intent(in)                 :: revDiff

        real(dp), dimension(1:2,0:nbins-1,1:2) :: getUconv

        integer  :: i                                ! loop counter
        integer  :: p                                ! loop counter

        real(dp) :: cosdL                            ! cos of dL
        real(dp) :: cosL                             ! cos of longitude (w/o Earth's rotation)
        real(dp) :: coslam                           ! cos of longitude
        real(dp) :: coslat                           ! cos of latitude
        real(dp) :: dL                               ! corrected longitude, considering Earth's rotation
        real(dp) :: lat                              ! latitude in rad
        real(dp) :: lon                              ! longitude in rad
        real(dp) :: revs                             ! parameter to account for 'revDiff' in computation of groundtrack shift
        real(dp) :: sindL                            ! sin of dL
        real(dp) :: sinL                             ! sin of longitude (w/o Earth's rotation)
        real(dp) :: sinlam                           ! sin of longitude
        real(dp) :: sinlat                           ! sin of latitude
        real(dp) :: uj                               ! argument of true latitude 'j'

        !========================================
        !
        ! compute values for uj(p=2) and uk(p=1)
        !
        !----------------------------------------
        do p = 1,2

          if(p == 2) then
            revs = 0
          else
            revs = twopi*revDiff
          end if

          do i = 0,nbins-1

            uj     = step_rad*i
            sinlat = this%sini*sin(uj)                 ! <-- ADAPT here for non-spherical Earth
            coslat = sqrt(1.d0 - sinlat*sinlat)
            lat    = asin(sinlat)
            dL     = this%trackShift*(revs +uj)/twopi  ! <-- ADAPT here for eccentric orbits

            sindL = sin(dL)
            cosdL = cos(dL)
            cosL  = cos(uj)/coslat   ! right spherical triangle
            sinL  = sqrt(1.d0 - cosL*cosL)

            if(sinlat < 0.d0) sinL = -sinL

            sinlam = sinL*cosdL - cosL*sindL   ! sin(lam) = sin(lam_0 - deltaLam), accounting for Earth's rotation
            coslam = cosL*cosdL + sinL*sindL

            if(abs(coslam) > eps15) then

              lon  = atan(sinlam/coslam)
              if(coslam < 0.d0) lon = lon + pi

            else

              if(sinlam > 0.d0) then
                lon = halfPi
              else if(sinlam < 0.d0) then
                lon = pi + halfPi ! 3/2 pi
              end if

            end if

            if(lon < 0.d0) lon = lon + twopi

            getUconv(p,i,1) = lat
            getUconv(p,i,2) = lon

          end do

        end do

    end function getUconv
 
    !========================================================================
    !
    !> @anchor      getNoisePropagationFlag
    !!
    !> @brief       Get the correlation matrix integration flag
    !> @author      Vitali Braun
    !!
    !> @return      .true. when noise propagation is enabled
    !!
    !> @date        <ul>
    !!                <li> 22.07.2013 (initial design)</li>
    !!              </ul>
    !!
    !------------------------------------------------------------------------
    logical function getNoisePropagationFlag(this)

        class(Correlation_class)    :: this
        getNoisePropagationFlag = this%propagateNoise
        return

    end function getNoisePropagationFlag

    !========================================================================
    !
    !> @anchor      getRevTransition
    !!
    !> @brief       Returns the correlation matrix term for one combination
    !!              of coefficients j and k
    !!
    !> @author      Vitali Braun
    !!
    !> @param[in]   irev    number of integer revolutions
    !> @param[in]   j       first correlation variable
    !> @param[in]   k       second correlation variable
    !> @param[in]   imat    matrix I
    !!
    !> @return      6x6 matrix
    !!
    !> @date        <ul>
    !!                <li> 03.08.2013 (initial design)</li>
    !!              </ul>
    !!
    !------------------------------------------------------------------------
    function getRevTransition(this,irev,j,l,imat) result(k)

        class(Correlation_class)            :: this

        integer, intent(in) :: irev
        integer, intent(in) :: j
        integer, intent(in) :: l
        real(dp), dimension(6,6), intent(in) :: imat
        real(dp), dimension(6,6)             :: k

        real(dp) :: revj    ! revolution number for variable j
        real(dp) :: revk    ! revolution number for variable k
        real(dp) :: temp    ! auxiliary

        !** initialise
        k = imat

        !** compute single secular term in matrix
        temp = -3.d0*pi*(1.d0 + this%meanEls%ecc*cos(this%meanEls%aop))**2.d0/(this%ome2**1.5d0*this%meanEls%sma)
        revj = temp*(irev - j)
        revk = temp*(irev - l)

        k(1,1) = imat(1,1) + revj*imat(2,1) + revk*imat(1,2) + revj*revk*imat(2,2)
        k(1,2) = imat(1,2) + revj*imat(2,2)
        k(1,3) = imat(1,3) + revj*imat(2,3)
        k(1,4) = imat(1,4) + revj*imat(2,4)
        k(1,5) = imat(1,5) + revj*imat(2,5)
        k(1,6) = imat(1,6) + revj*imat(2,6)
        k(2,1) = imat(2,1) + revk*imat(2,2)
        k(3,1) = imat(3,1) + revk*imat(3,2)
        k(4,1) = imat(4,1) + revk*imat(4,2)
        k(5,1) = imat(5,1) + revk*imat(5,2)
        k(6,1) = imat(6,1) + revk*imat(6,2)

        !    do p=1,6
        !      write(*,'(6(1x,e12.4e2))') (k(p,q) - imat(p,q), q=1,6)
        !    end do

        return

    end function getRevTransition

    !====================================================================================
    !
    !> @anchor      getCorrelationMatrix
    !!
    !> @brief       Returns the correlation matrix for a given time
    !> @author      Vitali Braun
    !!
    !> @param[in]   rqtime    requested time in seconds since t0
    !!
    !> @return      6x6 matrix
    !!
    !> @date        <ul>
    !!                <li> 07.08.2013 (initial design)</li>
    !!                <li> 18.11.2015 (fixed issue where interpolation did not work for less than five revs.)</li>
    !!              </ul>
    !!
    !-------------------------------------------------------------------------------------
    function getCorrelationMatrix(this,rqtime) result(cm)

        class(Correlation_class)    :: this
        real(dp), intent(in)        :: rqtime
        real(dp), dimension(6,6)    :: cm

        character(len=*), parameter :: csubid = 'getCorrelationMatrix'
        integer :: end_idx                   ! upper index for interpolation in correlation matrix array
        integer :: ilow                      ! index in correlation matrix index which corresponds to current revolution number
        integer :: iup                       ! ilow + 1
        integer :: i,k,l                     ! loop counter
        integer :: npoints                   ! interpolation polynomial degree 
        integer :: start_idx                 ! lower index for interpolation in correlation matrix array

        integer, parameter :: MAX_POINTS = 5    ! maximum number of points used in interpolation
        real(dp), dimension(MAX_POINTS) :: xi   ! x-values for interpolation
        real(dp), dimension(MAX_POINTS) :: yi   ! y-values for interpolation

        if(isControlled()) then
          if(hasToReturn()) return
          call checkIn(csubid)
        end if

        if(allocated(this%corrMat)) then ! start interpolation procedure...

          !** find out interpolation degree, which is between 1 (min. - 2 points) and 4 (max. - 5 points)
          npoints = max(2,min(size(this%corrMat,1),MAX_POINTS))
          npoints = 2

          !** get bounding array elements
          ilow = int(rqtime/this%orbitPeriod) 
          iup  = ilow + 1

          if((ilow < 0) .or. (iup > (size(this%corrMat,1) - 1))) then
            call setNeptuneError(E_CORR_INTERPOLATION, FATAL)
            return
          end if

          !===============================================
          !
          ! Start interpolation
          !
          !------------------------------------------
          start_idx = ilow - (npoints+1)/2 + 1
          end_idx   = ilow +  npoints/2

          if(start_idx < 0) then
            start_idx = 0
            end_idx   = npoints-1
          else if(end_idx > (size(this%corrMat,1)-1)) then
            end_idx   = size(this%corrMat,1) - 1   ! because array starts at index '0'
            start_idx = end_idx - npoints + 1
          end if

          do i=0,npoints-1
            xi(i+1) = start_idx + i
          end do

          do i=1,6
            do k=1,6
              do l=0,npoints-1
                yi(l+1) = this%corrMat(start_idx + l, i, k)
                !if(i==1 .and. k==1) then
                !  write(*,*) this%corrMat(start_idx+l,i,k)
                !  write(*,*) "xi = ", xi
                !  write(*,*) "yi = ", yi(l+1)
                !  write(*,*) "x = ", ilow + mod(rqtime,this%orbitPeriod)/this%orbitPeriod 
                !  write(*,*)
                !  !read(*,*)
                !end if
              end do

              call lagrange_interpolation(                                                &
                                           xi(1:npoints),                                 &  ! <-- DBL() data points array - x-values
                                           yi(1:npoints),                                 &  ! <-- DBL() data points array - y-values
                                           npoints,                                       &  ! <-- INT   number of data points
                                           ilow + mod(rqtime,this%orbitPeriod)/this%orbitPeriod,    &  ! <-- DBL   searched-for x value, which is the fraction of the current revolution
                                           cm(i,k)                                        &  ! --> DBL   resulting y value
                                         )
            end do
          end do
        else  ! error
          call setNeptuneError(E_CORR_INIT, FATAL)
          return
        end if

        if(isControlled()) then
          call checkOut(csubid)
        end if
        return

    end function getCorrelationMatrix

    !====================================================================================
    !
    !> @anchor      updateCorrelationMatrix
    !!
    !> @brief       Update the correlation matrix for given time step
    !> @author      Vitali Braun
    !!
    !> @param[in]   r         radius vector in GCRF
    !> @param[in]   dt        update time span (t-t0) in seconds
    !!
    !> @date        <ul>
    !!                <li> 22.07.2013 (initial design)</li>
    !!              </ul>
    !!
    !-------------------------------------------------------------------------------------
    subroutine updateCorrelationMatrix(                 &
                                      this,             &
                                      gravity_model,    &
                                      r,                & ! <-- DBL()   radius vector in GCRF
                                      span              & ! <-- DBL     time span in seconds
                                    )

        class(Correlation_class)                :: this
        type(Gravity_class),intent(inout)       :: gravity_model
        real(dp), dimension(3),   intent(in)    :: r
        real(dp),                 intent(in)    :: span


        character(len=*), parameter :: csubid = 'updateCorrelationMatrix'
        integer :: j,k,l                             ! loop counter
        real(dp), dimension(6,6)    :: corrMat0    ! correlation matrix for i=j
        real(dp), dimension(6,6)    :: corrMatS    ! correlation matrix for i/=j
        real(dp), dimension(6,6)    :: corrMat_jj  ! matrix I_0    which reflects the correlation for uj = uk
        real(dp), dimension(6,6)    :: corrMat_jk  ! matrix I_jk   which reflects correlation for uj /= uk
        real(dp), dimension(6,6)    :: corrMatAux  ! auxiliary matrix to sum up single I terms for the correlation matrix
        real(dp), dimension(6,6)    :: jac         ! jacobian to convert correlation matrix from non-singular orbit elements to GCRF
        real(dp)                    :: dlambda     ! track shift after a given number of revolutions

        type(nonSingOrbEl_t) :: nonSingEls

        !sig0 = 0.d0
        !sigs = 0.d0

        if(isControlled()) then
          if(hasToReturn()) return
          call checkIn(csubid)
        end if

        !==============================================================================
        !
        ! Compute required orbit parameters
        !
        !-------------------------------------------------

        this%meanEls = getMeanElements()
        if(hasFailed()) return

        !** check inclination first, as there exists a singularity for equatorial orbits
        if(this%meanEls%inc < eps6) then
          call setNeptuneError(E_CORR_SINGULAR, WARNING)
          call this%setNoisePropagationFlag(.false.)
          return
        end if

        !** convert to non-singular element set
        nonSingEls = kepler2nonSing(this%meanEls)

        !** store initial argument of true latitude as integer index
        this%u0 = int(nonSingEls%arglat)

        !** orbit period in seconds
        this%orbitPeriod = getOrbitalPeriod(this%meanEls%sma)
        if(hasFailed()) return

        !** groundtrack shift (account for siderial rotation of the Earth)
        this%trackShift = getEarthRotation(0.d0, EGM96)*this%orbitPeriod

        !** number of revolutions
        this%nrevs = int(span/this%orbitPeriod) + 1   ! integer number of evolutions where t(nrev) > span

        !** inclination
        this%cosi       = cos(this%meanEls%inc)
        this%sini       = sin(this%meanEls%inc)
        this%ctgi       = this%cosi/this%sini

        !** eccentricity
        this%ome2  = 1.d0 - this%meanEls%ecc*this%meanEls%ecc

        !==============================================================================
        !
        ! Compute correlation function Kr
        !
        !------------------------------------------------------------------------
        this%Kr = this%getCorrFunction(gravity_model, r)

        !------------------------------------------------------------------------
            
            !==================================================================================
            !
            ! FEATURE/REMOVE: Output of correlation function for different revolutions
            !
            !-----------------------------------------------------------------------------
            !open(unit=112, file="corrRev.out")

        !==============================================================================
        !
        ! Compute Q_0 which corresponds to u_j = u_k (k = j)
        !
        !--------------------------------------------------------------
        corrMat_jj = this%getCorrelationTerm(0)
        !--------------------------------------------------------------

        !===============================================================================
        !
        ! Start loop over area defined by parameters j and k. Only a triangle area is
        ! analysed due to symmetry.
        !
        !------------------------------------------------------------

        !** allocate correlation matrix array for nrevs revolutions
        if (allocated(this%corrMat)) deallocate(this%corrMat)
        allocate(this%corrMat(0:this%nrevs,6,6))

        this%corrMat  = 0.d0
        corrMat0 = 0.d0
        corrMatS = 0.d0

        !** do square areas where j=k
        !--------------------------------------
        do j = 1,this%nrevs
          do k=1,j
            corrMatAux    = this%getRevTransition(j, k, k, corrMat_jj)
            corrMat0(:,:) = corrMat0(:,:) + corrMatAux(:,:)
          end do

          this%corrMat(j,:,:) = this%corrMat(j,:,:) + corrMat0

          ! do k=1,6
          !   sig0(j,k) = sqrt(corrMat0(k,k))
          ! end do
          ! write(*,'(i2.2,6(1x,e12.5e2))') j, (sig0(j,l), l=1,6)
        end do
        !--------------------------------------

         !write(*,*) "CorrMat0 = "
         !write(*,*) "-----------------------"
         !do l=1,6
         !  write(*,'(6(1x,e12.5e2))') (corrMat0(l,j), j=1,6)
         !end do
         !write(*,*) "-----------------------"

         !write(*,*) "CorrMat_first = ", sqrt(this%corrMat(1,1,1) )
         !read(*,*)
         !write(*,*) "Off-diagonals:"
         !write(*,*) "-----------------------------"

        !** now do off-diagonal squares
        !--------------------------------------
        do l = 2, this%nrevs
          do j = 1,l-1
            dlambda = this%trackShift*j
            dlambda = redang(dlambda, 2, 1, .false.)  ! reduce to 0..2pi

            !if(dlambda < (pi - this%trackShift*0.5d0) .or. dlambda > (pi + this%trackShift*0.5d0)) then ! I_(j-k) exists
            corrMat_jk = this%getCorrelationTerm(j)

            do k = 1, l - j
              corrMatAux   = this%getRevTransition(l, j, k+j, corrMat_jk)
              corrMatS(:,:) = corrMatS(:,:) + corrMatAux(:,:)
            end do
            !end if
          end do
          
          this%corrMat(l,:,:) = this%corrMat(l,:,:) + corrMatS + transpose(corrMatS)
        !       do j = 1,6
        !         sigs(l,j) = sig0(l,j)*sig0(l,j) + 2.d0*corrMatS(j,j)
        !         if(sigs(l,j) > 0.d0) then
        !           sigs(l,j) = sqrt(sigs(l,j))
        !         else
        !           sigs(l,j) = 0.d0
        !         end if
        !       end do

        !       write(*,'(i2.2,6(1x,e12.5e2))') l, (sigs(l,j), j=1,6)
        end do
        !--------------------------------------
         !write(*,*) "CorrMatS = "
         !write(*,*) "-----------------------"
        !     do l=0,nrevs
        !       write(112,'(6(1x,e12.5e2))') (corrMat(l,k,k), k=1,6)
        !     end do
         !write(*,*) "-----------------------"
         !do l=1,6
         !  write(*,'(6(1x,e12.5e2))') (corrMat(nrevs,l,j), j=1,6)
         !end do


        !** convert to GCRF
        !---------------------------------------------------------
        jac     = this%jacobi_oe2rv(nonSingEls)

        !do j=1,6
        !  write(*,'(6(x,e12.4e2))') (jac(j,k), k=1,6)
        !end do
        !read(*,*)

        do j=0,this%nrevs
          this%corrMat(j,:,:) = matmul(matmul(jac,this%corrMat(j,:,:)),transpose(jac))
         ! write(*,*) "j", j, (corrMat(j,k,k), k=1,6)
        end do
        !read(*,*)
        !--------------------------------------------------------
            
            !==================================================================================
            !
            ! FEATURE/REMOVE: Output of correlation function for different revolutions
            !
            !-----------------------------------------------------------------------------
            !close(112)
            !stop


        !** done!
        if(isControlled()) then
          call checkOut(csubid)
        end if

        return

    end subroutine updateCorrelationMatrix

    !========================================================================
    !
    !> @anchor      setNoisePropagationFlag
    !!
    !> @brief       Set the correlation matrix integration flag
    !> @author      Vitali Braun
    !!
    !> param[in]    indicates whether noise propagation is enabled
    !!
    !> @date        <ul>
    !!                <li> 22.07.2013 (initial design)</li>
    !!              </ul>
    !!
    !------------------------------------------------------------------------
    subroutine setNoisePropagationFlag(this,flag)

        class(Correlation_class)    :: this
        logical, intent(in)         :: flag

        this%propagateNoise = flag
        return

    end subroutine setNoisePropagationFlag

    !=====================================================================================================
    !
    !> @anchor      jacobi_oe2rv
    !!
    !> @brief       Computes jacobian matrix for the conversion of orbital elements to cartesian state
    !> @author      Vitali Braun
    !!
    !> @param[in]   oel   non-singular (for circular orbits) orbit elements as given by Nazarenko(2010)
    !!                    or Alfriend et al.(2010)
    !!
    !> @return      6x6 matrix
    !!
    !> @date        <ul>
    !!                <li> 07.08.2013 (initial design)</li>
    !!                <li> 15.11.2015 (corrected several partial derivatives)</li>
    !!              </ul>
    !!
    !-----------------------------------------------------------------------------------------------------
    function jacobi_oe2rv(this,oel) result(j)

        implicit none

        class(Correlation_class)            :: this
        type(nonSingOrbEl_t), intent(in)    :: oel

        real(dp), dimension(6,6) :: j

        real(dp) :: cosi          ! cos of inclination
        real(dp) :: cosra         ! cos of raan
        real(dp) :: cosu          ! cos of arg of latitude
        real(dp) :: c1            ! auxiliary
        real(dp) :: c2            ! auxiliary
        real(dp) :: c3            ! auxiliary
        real(dp) :: c4            ! auxiliary
        real(dp) :: c5            ! auxiliary
        real(dp) :: c6            ! auxiliary
        real(dp) :: c7            ! auxiliary
        real(dp) :: c8            ! auxiliary
        real(dp) :: c9            ! auxiliary
        real(dp) :: p             ! semilatus rectum (km)
        real(dp) :: rabs          ! radius magnitude
        real(dp) :: rabs2         ! rabs**2
        real(dp) :: r2op          ! rabs2/p
        real(dp) :: rop           ! rabs/p
        real(dp) :: rcosu         ! rabs*cosu
        real(dp) :: rsinu         ! rabs*sinu
        real(dp) :: sini          ! sin of inclination
        real(dp) :: sinra         ! sin of raan
        real(dp) :: sinu          ! sin of arg of latitude
        real(dp) :: LsuHcu        ! auxiliary

        real(dp) :: p11,p12   ! first line of rotation matrix ECI/PQW
        real(dp) :: p21,p22   ! second line of rotation matrix ECI/PQW
        real(dp) :: p32       ! third line of rotation matrix ECI/PQW

        real(dp) :: dp12di ! inclination deriv. of first line of rotation matrix ECI/PQW
        real(dp) :: dp22di ! inclination deriv. of second line of rotation matrix ECI/PQW
        real(dp) :: dp32di ! inclination deriv. of third line of rotation matrix ECI/PQW

        real(dp) :: dp11dr,dp12dr ! raan deriv. of first line of rotation matrix ECI/PQW
        real(dp) :: dp21dr,dp22dr ! raan deriv. of second line of rotation matrix ECI/PQW

        !** compute additional parameters and auxiliaries
        sinu  = sin(oel%arglat)
        cosu  = cos(oel%arglat)
        p     = oel%sma*(1.d0 - oel%q1**2.d0 - oel%q2**2.d0)
        rabs  = p/(1.d0 + oel%q1*cosu + oel%q2*sinu)
        LsuHcu = oel%q1*sinu - oel%q2*cosu
        rabs2 = rabs*rabs
        r2op  = rabs2/p
        rop   = rabs/p
        rsinu = rabs*sinu
        rcosu = rabs*cosu

        !** calculate trigonometric terms
        cosra = cos(oel%raan)
        cosi  = cos(oel%inc)

        sinra = sin(oel%raan)
        sini  = sin(oel%inc)

        !** process constants
        c1 = -rabs*sinu - r2op*cosu*(oel%q2*cosu - oel%q1*sinu)
        c2 =  rabs*cosu - r2op*sinu*(oel%q2*cosu - oel%q1*sinu)
        c3 = rabs/oel%sma
        c4 = -r2op*cosu*cosu
        c5 = -r2op*sinu*cosu
        c6 = -r2op*sinu*sinu

        c7 = -sqrt(getEarthGravity()/p)
        c8 =  c7*(oel%q2 + sinu)
        c9 = -c7*(oel%q1 + cosu)

        p11 =  cosra
        p12 = -sinra*cosi
        p21 =  sinra
        p22 =  cosra*cosi
        !p31 =  0.d0
        p32 =  sini

        !dp11di = 0.d0
        dp12di = sinra*sini
        !dp21di =  0.d0
        dp22di = -cosra*sini
        !dp31di = 0.d0
        dp32di = cosi

        dp11dr = -p21
        dp12dr = -p22
        dp21dr =  p11
        dp22dr =  p12

        !** process jacobian matrix
        ! j(1,1) = p11*c1 + p12*c2                       ! d(r_i)/d_u
        j(1,1) = r2op*LsuHcu*(p11*cosu + p12*sinu) - rabs*(p11*sinu - p12*cosu)    ! d(r_i)/d_u
        j(1,2) = c3*(cosu*p11 + sinu*p12)              ! d(r_i)/d_a
        j(1,3) = dp12di*rsinu                          ! d(r_i)/d_inc
        j(1,4) = dp11dr*rcosu + dp12dr*rsinu           ! d(r_i)/d_raan
        !j(1,5) = p11*c4 + p12*c5                       ! d(r_i)/d_q1
        j(1,5) = -rop*(2.d0*oel%sma*oel%q1 + rabs*cosu)*(p11*cosu + p12*sinu)    ! d(r_i)/d_q1
        !j(1,6) = p11*c5 + p12*c6                       ! d(r_i)/d_q2
        j(1,6) = -rop*(2.d0*oel%sma*oel%q2 + rabs*sinu)*(p11*cosu + p12*sinu)    ! d(r_i)/d_q2

        !j(2,1) = p21*c1 + p22*c2                       ! d(r_j)/d_u
        j(2,1) = r2op*LsuHcu*(p21*cosu + p22*sinu) - rabs*(p21*sinu - p22*cosu)  ! d(r_j)/d_u
        j(2,2) = c3*(cosu*p21 + sinu*p22)              ! d(r_j)/d_a
        j(2,3) = dp22di*rsinu                          ! d(r_j)/d_inc
        j(2,4) = dp21dr*rcosu + dp22dr*rsinu           ! d(r_j)/d_raan
        !j(2,5) = p21*c4 + p22*c5                       ! d(r_j)/d_q1
        j(2,5) = -rop*(2.d0*oel%sma*oel%q1 + rabs*cosu)*(p21*cosu + p22*sinu)    ! d(r_j)/d_q1
        j(2,6) = -rop*(2.d0*oel%sma*oel%q2 + rabs*sinu)*(p21*cosu + p22*sinu) ! d(r_j)/d_q2
        !j(2,6) = p21*c5 + p22*c6                       ! d(r_j)/d_q2

        !j(3,1) = p32*c2                                ! d(r_k)/d_u
        j(3,1) = r2op*LsuHcu*sinu*p32 - rabs*cosu*p32        ! d(r_k)/d_u
        j(3,2) = c3*sinu*p32                           ! d(r_k)/d_a
        j(3,3) = dp32di*rsinu                          ! d(r_k)/d_inc
        j(3,4) = 0.d0                                  ! d(r_k)/d_raan
        j(3,5) = -rop*(2.d0*oel%sma*oel%q1 + rabs*cosu)*p32*sinu          ! d(r_k)/d_q1
        j(3,6) = -rop*(2.d0*oel%sma*oel%q2 + rabs*sinu)*p32*sinu          ! d(r_k)/d_q2
        !j(3,6) = p32*c6                                ! d(r_k)/d_q2

        j(4,1) = c7*(p11*cosu + p12*sinu)              ! d(v_i)/d_u
        j(4,2) = 0.5d0*c7/oel%sma*(p11*(sinu+oel%q2) - p12*(cosu+oel%q1))    ! d(v_i)/d_a
        j(4,3) = c9*dp12di                             ! d(v_i)/d_inc
        !j(4,4) = c8*dp11dr                             ! d(v_i)/d_raa
        j(4,4) = c7*((sinu+oel%q2)*dp11dr - (cosu+oel%q1)*dp12dr)       ! d(v_i)/d_raa
        !j(4,5) = -c7*p12                              ! d(v_i)/d_q1
        j(4,5) = c7*oel%q1*oel%sma/p*(p11*(sinu+oel%q2) + p12*(cosu+oel%q1)) - c7*p12       ! d(v_i)/d_q1
        !j(4,6) =  c7*p11                               ! d(v_i)/d_q2
        j(4,6) = c7*oel%q2*oel%sma/p*(p11*(sinu+oel%q2) + p12*(cosu+oel%q1)) + c7*p11       ! d(v_i)/d_q2
                                                                          
        j(5,1) = c7*(p21*cosu + p22*sinu)              ! d(v_j)/d_u
        !j(5,2) = 0.5d0*c8/oel%sma*(c8*p21 + c9*p22)    ! d(v_j)/d_a
        j(5,2) = 0.5d0*c7/oel%sma*(p21*(sinu+oel%q2) - p22*(cosu+oel%q1))    ! d(v_j)/d_a
        j(5,3) = c9*dp22di                             ! d(v_j)/d_inc
        j(5,4) = c7*((sinu+oel%q2)*dp21dr - (cosu+oel%q1)*dp22dr)    ! d(v_j)/d_raa
        j(5,5) = -c7*oel%q1*oel%sma/p*(p21*(sinu+oel%q2) + p22*(cosu+oel%q1)) - c7*p22  ! d(v_j)/d_q1
        !j(5,6) =  c7*p21                               ! d(v_j)/d_q2
        j(5,6) =  -c7*oel%q2*oel%sma/p*(p21*(sinu+oel%q2) + p22*(cosu+oel%q1)) + c7*p21                            ! d(v_j)/d_q2
                                                                          
        j(6,1) = c7*p32*sinu                           ! d(v_k)/d_u
        !j(6,2) = 0.5d0*c8/oel%sma*p32*c9               ! d(v_k)/d_a
        j(6,2) = 0.5d0*c7/oel%sma*p32*(oel%q1+cosu)   ! d(v_k)/d_a
        j(6,3) = c9*dp32di                             ! d(v_k)/d_inc
        j(6,4) = 0.d0                                  ! d(v_k)/d_raa
        j(6,5) = (c9*oel%q1*oel%sma/p - c7)*p32        ! d(v_k)/d_q1
        j(6,6) = c9*oel%q2*oel%sma/p*p32               ! d(v_k)/d_q2

        return

    end function jacobi_oe2rv

    !=====================================================================================================
    !
    !> @anchor      jacobi_ke2rv
    !!
    !> @brief       Computes jacobian matrix for the conversion of keplerian elements to cartesian state
    !> @author      Vitali Braun
    !!
    !> @param[in]   el  keplerian elements
    !!
    !> @return      6x6 matrix
    !!
    !> @date        <ul>
    !!                <li> 07.08.2013 (initial design)</li>
    !!                <li> 15.11.2015 (corrected several partial derivatives)</li>
    !!              </ul>
    !!
    !-----------------------------------------------------------------------------------------------------
    subroutine jacobi_ke2rv(          &  ! computes jacobian matrix for covariance transformation
                          this,   &
                          el,     &  ! <-- DBL() keplerian elements
                                     !           el(1) = semi-major axis  [km]
                                     !           el(2) = eccentricity     [-]
                                     !           el(3) = inclination      [rad]
                                     !           el(4) = raan             [rad]
                                     !           el(5) = arg. of perigee  [rad]
                                     !           el(6) = true anomaly     [rad]
                          j_ke2rv &  ! --> DBL() jacobian matrix 
                       )
        !-------------------------------------------------------------------------------
        !  * Programming Language * Fortran 95/2003
        !-------------------------------------------------------------------------------
        !  * Purpose *
        !  processing the jacobian matrix for covariance matrix transformation from 
        !  keplerian to cartesian coordinates.
        !-------------------------------------------------------------------------------

        !  use astro
        !  use types

        implicit none

        class(Correlation_class)            :: this
        real(dp), dimension(6), intent(in)  :: el

        real(dp), dimension(6,6), intent(out) :: j_ke2rv

        real(dp) :: a             ! semi-major axis
        real(dp) :: aop           ! arg. of perigee
        real(dp) :: cosap         ! cos of aop
        real(dp) :: cosi          ! cos of inclination
        real(dp) :: cosra         ! cos of raan
        real(dp) :: costa         ! cos of true anomaly
        real(dp) :: c1            ! auxiliary
        real(dp) :: c2            ! auxiliary
        real(dp) :: c3            ! auxiliary
        real(dp) :: c4            ! auxiliary
        real(dp) :: c5            ! auxiliary
        real(dp) :: c6            ! auxiliary
        real(dp) :: c7            ! auxiliary
        real(dp) :: eps           ! eccentricity
        real(dp) :: eps2          ! eccentricity^2
        real(dp) :: inc           ! inclination
        real(dp) :: raan          ! right asc. of asc. node
        real(dp) :: sinap         ! sin of aop
        real(dp) :: sini          ! sin of inclination
        real(dp) :: sinra         ! sin of raan
        real(dp) :: sinta         ! sin of true anomaly
        real(dp) :: tran          ! true anomaly

        real(dp) :: p11,p12   ! first line of rotation matrix ECI/PQW
        real(dp) :: p21,p22   ! second line of rotation matrix ECI/PQW
        real(dp) :: p31,p32   ! third line of rotation matrix ECI/PQW

        real(dp) :: dp11di,dp12di ! inclination deriv. of first line of rotation matrix ECI/PQW
        real(dp) :: dp21di,dp22di ! inclination deriv. of second line of rotation matrix ECI/PQW
        real(dp) :: dp31di,dp32di ! inclination deriv. of third line of rotation matrix ECI/PQW

        real(dp) :: dp11dr,dp12dr ! raan deriv. of first line of rotation matrix ECI/PQW
        real(dp) :: dp21dr,dp22dr ! raan deriv. of second line of rotation matrix ECI/PQW

        real(dp) :: dp11do,dp12do ! arg. of perigee deriv. of first line of rotation matrix ECI/PQW
        real(dp) :: dp21do,dp22do ! arg. of perigee deriv. of second line of rotation matrix ECI/PQW
        real(dp) :: dp31do,dp32do ! arg. of perigee deriv. of third line of rotation matrix ECI/PQW

        !** assign keplerian elements
        a    = el(1)
        eps  = el(2)
        inc  = el(3)
        raan = el(4)
        aop  = el(5)
        tran = el(6)

        eps2 = eps*eps

        !** calculate trigonometric terms
        costa = cos(tran)
        cosap = cos(aop)
        cosra = cos(raan)
        cosi  = cos(inc)

        sinta = sin(tran)
        sinap = sin(aop)
        sinra = sin(raan)
        sini  = sin(inc)

        !** process constants
        c1 = (1.d0 - eps2)/(1.d0 + eps*costa)
        c2 = -(2.d0*a*eps + a*costa + a*eps2*costa)/(1.d0+eps*costa)**2.d0
        c3 = a*c1/(1.d0 + eps*costa)
        c4 = sqrt(getEarthGravity()/(a*(1.d0 - eps2)))  
        c5 = (1.d0 + eps*cos(tran))**2.d0/(1.d0 - eps2)**1.5d0
        c6 = c4*0.5d0/a
        c7 = c4/(1.d0 - eps2)

        p11 =  cosra*cosap - sinra*sinap*cosi
        p12 = -cosra*sinap - sinra*cosap*cosi
        p21 =  sinra*cosap + cosra*sinap*cosi
        p22 = -sinra*sinap + cosra*cosap*cosi
        p31 =  sinap*sini
        p32 =  cosap*sini

        dp11di = sinra*sinap*sini
        dp12di = sinra*cosap*sini
        dp21di = -cosra*sinap*sini
        dp22di = -cosra*cosap*sini
        dp31di = sinap*cosi
        dp32di = cosap*cosi

        dp11dr = -p21
        dp12dr = -p22
        dp21dr =  p11
        dp22dr =  p12

        dp11do =  p12 
        dp12do = -p11
        dp21do =  p22
        dp22do = -p21
        dp31do =  p32
        dp32do = -p31

        !** process jacobian matrix
        j_ke2rv(1,1) = c1*(costa*p11 + sinta*p12)
        j_ke2rv(1,2) = c2*(costa*p11 + sinta*p12)
        j_ke2rv(1,3) = a*c1*(costa*dp11di + sinta*dp12di)
        j_ke2rv(1,4) = a*c1*(costa*dp11dr + sinta*dp12dr)
        j_ke2rv(1,5) = a*c1*(costa*dp11do + sinta*dp12do)
        j_ke2rv(1,6) = c3*c5*((eps + costa)*p12 - sinta*p11)

        j_ke2rv(2,1) = c1*(costa*p21 + sinta*p22)
        j_ke2rv(2,2) = c2*(costa*p21 + sinta*p22)
        j_ke2rv(2,3) = a*c1*(costa*dp21di + sinta*dp22di)
        j_ke2rv(2,4) = a*c1*(costa*dp21dr + sinta*dp22dr)
        j_ke2rv(2,5) = a*c1*(costa*dp21do + sinta*dp22do)
        j_ke2rv(2,6) = c3*c5*((eps + costa)*p22 - sinta*p21)

        j_ke2rv(3,1) = c1*(costa*p31 + sinta*p32)
        j_ke2rv(3,2) = c2*(costa*p31 + sinta*p32)
        j_ke2rv(3,3) = a*c1*(costa*dp31di + sinta*dp32di)
        j_ke2rv(3,4) = 0.d0
        j_ke2rv(3,5) = a*c1*(costa*dp31do + sinta*dp32do)
        j_ke2rv(3,6) = c3*c5*((eps + costa)*p32 - sinta*p31)

        j_ke2rv(4,1) = c6*(sinta*p11 - (eps + costa)*p12)
        j_ke2rv(4,2) = c7*((1.d0 + eps*costa)*p12 - eps*sinta*p11)
        j_ke2rv(4,3) = c4*((eps + costa)*dp12di - sinta*dp11di)
        j_ke2rv(4,4) = c4*((eps + costa)*dp12dr - sinta*dp11dr)
        j_ke2rv(4,5) = c4*((eps + costa)*dp12do - sinta*dp11do)
        j_ke2rv(4,6) = c4*c5*(-costa*p11 - sinta*p12)

        j_ke2rv(5,1) = c6*(sinta*p21 - (eps + costa)*p22)
        j_ke2rv(5,2) = c7*((1.d0 + eps*costa)*p22 - eps*sinta*p21)
        j_ke2rv(5,3) = c4*((eps + costa)*dp22di - sinta*dp21di) 
        j_ke2rv(5,4) = c4*((eps + costa)*dp22dr - sinta*dp21dr)
        j_ke2rv(5,5) = c4*((eps + costa)*dp22do - sinta*dp21do)
        j_ke2rv(5,6) = c4*c5*(-costa*p21 - sinta*p22)

        j_ke2rv(6,1) = c6*(sinta*p31 - (eps + costa)*p32)
        j_ke2rv(6,2) = c7*((1.d0 + eps*costa)*p32 - eps*sinta*p31)
        j_ke2rv(6,3) = c4*((eps + costa)*dp32di - sinta*dp31di)
        j_ke2rv(6,4) = 0.d0
        j_ke2rv(6,5) = c4*((eps + costa)*dp32do - sinta*dp31do)
        j_ke2rv(6,6) = c4*c5*(-costa*p31 - sinta*p32)
        
    end subroutine jacobi_ke2rv
end module correlation

