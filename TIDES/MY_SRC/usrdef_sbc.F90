MODULE usrdef_sbc
   !!======================================================================
   !!                     ***  MODULE  usrdef_sbc  ***
   !!
   !!                     ===  EQ-WAVES configuration  ===
   !!
   !! User defined :   surface forcing of a user configuration
   !!======================================================================
   !! History :  1.0   ! 2023-03  (P. Oddo)  user defined interface
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usrdef_sbc    : user defined surface bounday conditions in GYRE case
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran    !

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usrdef_sbc_oce       ! routine called in sbcmod module
   PUBLIC   usrdef_sbc_ice_tau   ! routine called by icestp.F90 for ice dynamics
   PUBLIC   usrdef_sbc_ice_flx   ! routine called by icestp.F90 for ice thermo

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_sbc.F90 15145 2021-07-26 16:16:45Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usrdef_sbc_oce( kt, Kbb )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE usrdef_sbc  ***
      !!              
      !! ** Purpose :   provide at each time-step the GYRE surface boundary
      !!              condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   analytical seasonal cycle for GYRE configuration.
      !!                CAUTION : never mask the surface stress field !
      !!
      !! ** Action  : - set the ocean surface boundary condition, i.e.   
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, sfx
      !!
      !! Reference : Hazeleger, W., and S. Drijfhout, JPO, 30, 677-695, 2000.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      INTEGER, INTENT(in) ::   Kbb  ! ocean time index
      !!
      INTEGER  ::   ji, jj                 ! dummy loop indices
      REAL(wp) ::   zsumemp, zsurf
      REAL(wp) ::   ztx, zty, zmod, zcoef ! temporary variables
      REAL(wp) ::   zrhoa  = 1.22         ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3       ! drag coefficient
      REAL(wp) ::   Ampl_p  = 15.0        ! Ampl of the gaussian pertub
      REAL(wp) ::   sigma_p = 5.0         ! Sigma of the gaussian pertub
      REAL(wp) ::   lonc_p  = 200.0       ! Center lon of the gaussian pertub
      REAL(wp) ::   latc_p  = 0.0         ! Center lat of the gaussian pertub
      REAL(wp) ::   wind_1  = 7.5         ! Initial Wind speed 
      REAL(wp) ::   wind_2  = 7.5         ! Wind speed after the pertubation
      REAL(wp) ::   wind_a                ! Wind speed at kt time-step
      REAL(wp) ::   gauss                 ! pertubation
      !!---------------------------------------------------------------------

      ! --------------------------------------------------------------
      ! Null heat and freshwater fluxes  
      ! Do not vary in time
      ! --------------------------------------------------------------
      qsr(:,:) = 0.0_wp
      qns(:,:) = 0.0_wp
      emp(:,:) = 0.0_wp
      utau(:,:) = 0.0_wp
      vtau(:,:) = 0.0_wp
      taum(:,:) = 0.0_wp
      wndm(:,:) = 0.0_wp

!      zsumemp = GLOB_SUM( 'usrdef_sbc', emp  (:,:)   ) 
!      zsurf   = GLOB_SUM( 'usrdef_sbc', tmask(:,:,1) ) 
!      zsumemp = zsumemp / zsurf         ! Default GYRE configuration

      ! freshwater (mass flux) and update of qns with heat content of emp
!      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )   ! emp used in sshwzv over the whole domain
!         emp (ji,jj) = emp(ji,jj) - zsumemp * tmask(ji,jj,1)          ! freshwater flux (=0 in domain average)
!         sfx (ji,jj) = 0.0_wp                                         ! no salt flux
!!         qns (ji,jj) = qns(ji,jj) - emp(ji,jj) * sst_m(ji,jj) * rcp   ! evap and precip are at SST
!      END_2D


      ! --------------------------------------------------------------
      ! Momentum fluxes
      ! At the beginning apply constant/homogeneous wind in negative
      ! zonal direction
      ! Once the upwelling is well defined, apply a gaussian pertubation
      ! for approx a day
      ! --------------------------------------------------------------
      !wind_a = wind_1
      !IF (kt.GT.2572) wind_a = wind_2
      
      !zcoef = 1. / ( zrhoa * zcdrag ) 
      !utau(:,:) = -1.0_wp * ( zrhoa * zcdrag ) * 2.0_wp * (wind_a*wind_a)
    !   utau(:,:) = 0.0_wp
    !   vtau(:,:) = 0.0_wp
!      DO_2D( 0, 0, 0, 0 )
         !ztx = utau(ji-1,jj  ) + utau(ji,jj) 
         !zty = vtau(ji  ,jj-1) + vtau(ji,jj) 
         !zmod = 0.5 * SQRT( ztx * ztx + zty * zty )
         !taum(ji,jj) = zmod
         !wndm(ji,jj) = SQRT( zmod * zcoef )
!          taum(ji,jj) = 0.0_wp
!          wndm(ji,jj) = 0.0_wp
!      END_2D

      ! --------------------------------------------------------------
      ! Apply a Gaussian pertubation to the wind stress
      ! Here I am assuming time-step is 1800s and I apply the
      ! pertubation for approx 1 day. Clearly changing the time-step
      ! these numbers should be reviewed
      ! --------------------------------------------------------------
      !IF (kt.GE.2500.AND.kt.LE.2572) THEN
      !DO_2D( 1, 1, 1, 1 )
      !  gauss=Ampl_p*exp(-1/(sigma_p*sigma_p)*((glamt(ji,jj)-lonc_p)**2+(gphit(ji,jj)-latc_p)**2))
      !  utau(ji,jj) = utau(ji,jj) + gauss
      !END_2D
      !ENDIF

      ! ---------------------------------- !
      !  control print at first time-step  !
      ! ---------------------------------- !
      IF( kt == nit000 .AND. lwp ) THEN 
         WRITE(numout,*)
         WRITE(numout,*)'usrdef_sbc_oce : analytical surface fluxes for EQ-WAVE configuration'               
         WRITE(numout,*)'~~~~~~~~~~~ ' 
         WRITE(numout,*)'           nyear      = ', nyear
         WRITE(numout,*)'           nmonth     = ', nmonth
         WRITE(numout,*)'           nday       = ', nday
         WRITE(numout,*)'           nday_year  = ', nday_year
         WRITE(numout,*)'           zsumemp    = ', zsumemp
         WRITE(numout,*)'           zsurf      = ', zsurf
         WRITE(numout,*)'           ndastp     = ', ndastp
      ENDIF
      !
   END SUBROUTINE usrdef_sbc_oce


   SUBROUTINE usrdef_sbc_ice_tau( kt )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
   END SUBROUTINE usrdef_sbc_ice_tau


   SUBROUTINE usrdef_sbc_ice_flx( kt, phs, phi )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phs    ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phi    ! ice thickness
   END SUBROUTINE usrdef_sbc_ice_flx

   !!======================================================================
END MODULE usrdef_sbc
