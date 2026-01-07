MODULE usrdef_istate
   !!======================================================================
   !!                     ***  MODULE usrdef_istate   ***
   !!
   !!                      ===  TSUNAMI configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  NEMO !
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE dom_oce        
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lbclnk          ! lateral boundary conditions - mpp exchanges
   !   
   USE usrdef_nam
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate       ! called by istate.F90
   PUBLIC   usr_def_istate_ssh   ! called by domqco.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2016)
   !! $Id$ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here TSUNAMI configuration 
      !!
      !! ** Method  :   Set a gaussian anomaly of pressure and associated
      !!                geostrophic velocities
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s] 
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d   ! 2D workspace
      REAL(wp) ::   zfact
      INTEGER  ::   ji, jj, jk
      INTEGER  ::   igloi, igloj   ! to be removed in the future, see comment bellow
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : TSUNAMI configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   '
      !
      pts(:,:,:,jp_tem) = 20._wp
      pts(:,:,:,jp_sal) = 30._wp
      pu( :,:,:       ) =  0._wp
      pv( :,:,:       ) =  0._wp
      
   END SUBROUTINE usr_def_istate


   SUBROUTINE usr_def_istate_ssh( ptmask, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate_ssh  ***
      !! 
      !! ** Purpose :   Initialization of ssh
      !!                Here TSUNAMI configuration 
      !!
      !! ** Method  :   Set ssh
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask   [m]
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height   [m]
      !
      INTEGER  ::   ji, jj
      REAL(wp) :: zmax, zdist, zlon0, zlat0, zsigma, zdx, zdy
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate_ssh : TSUNAMI configuration, analytical definition of initial ssh'
      !
      ! Parametri della perturbazione (Angolo Sud-Est)
      zlon0 = 19.0_wp    ! Longitudine centro (Albania/Grecia)
      zlat0 = 40.5_wp    ! Latitudine centro (Canale d'Otranto)
      zsigma = 20000._wp ! Raggio 20 km
      zmax   = 0.1_wp    ! Altezza 5 cm
      pssh(:,:) = 0._wp
      
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         ! Calcola distanza approssimata in metri (1 grado ~ 111km)
         zdx = ( glamt(ji,jj) - zlon0 ) * 111000._wp * COS( rad * zlat0 )
         zdy = ( gphit(ji,jj) - zlat0 ) * 111000._wp
         zdist = SQRT( zdx * zdx + zdy * zdy )
         
         ! Applica gaussiana se dentro il raggio d'azione
         IF( zdist <= (4._wp * zsigma) ) THEN
             pssh(ji,jj) = zmax * EXP( - (zdist * zdist) / (zsigma * zsigma) ) * ptmask(ji,jj,1)
         ENDIF
      END_2D

      !
      CALL lbc_lnk('usrdef_istate', pssh, 'T',  1. )            ! apply boundary conditions
      !
   END SUBROUTINE usr_def_istate_ssh
   
   !!======================================================================
END MODULE usrdef_istate
