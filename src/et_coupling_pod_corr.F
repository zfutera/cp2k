!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief post-processing corrections of POD states
!> \author Z. Futera (07.2026)
! **************************************************************************************************
MODULE et_coupling_pod_corr

   USE cp_control_types,                ONLY: dft_control_type
   USE et_coupling_pod_types,           ONLY: et_cpl
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling_pod_corr'

   PUBLIC :: pod_corr, pod_corr_sigma

CONTAINS

! **************************************************************************************************
!> \brief apply various corrections to POD states 
!> \param output_unit unit number of the open output stream
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \author Z. Futera (07.2026)
! **************************************************************************************************
   SUBROUTINE pod_corr(output_unit, qs_env, ec)

      ! Routine arguments
      INTEGER                                            :: output_unit
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec

      ! DFT+Sigma
      CALL pod_corr_sigma(output_unit, qs_env, ec)

   END SUBROUTINE pod_corr

! **************************************************************************************************
!> \brief apply DFT+Sigma corrections to POD states 
!> \param output_unit unit number of the open output stream
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \author Z. Futera (07.2026)
! **************************************************************************************************
   SUBROUTINE pod_corr_sigma(output_unit, qs_env, ec)

      ! Routine arguments
      INTEGER                                            :: output_unit
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec

      ! Local variables
      INTEGER                                            :: i, j, k
      LOGICAL, DIMENSION(2)                              :: has_occ, has_unocc
      REAL(KIND=dp), DIMENSION(2)                        :: occ, unocc
      TYPE(dft_control_type), POINTER                    :: dft_cntrl
      TYPE(section_vals_type), POINTER                   :: inp_et_pod, inp_pod_block, inp_pod_sigma

      ! DFT setup
      CALL get_qs_env(qs_env, dft_control=dft_cntrl)

      ! POD setup
      inp_et_pod => section_vals_get_subs_vals(qs_env%input, 'PROPERTIES%ET_COUPLING%PROJECTION')
      inp_pod_block => section_vals_get_subs_vals(inp_et_pod, 'BLOCK')

      ! Blocks
      DO i = 1, ec%n_blocks

        inp_pod_sigma => section_vals_get_subs_vals(inp_pod_block, 'SIGMA', i_rep_section=i)

        CALL section_vals_val_get(inp_pod_sigma, 'OCCUPIED_ALPHA', &
                                  r_val=occ(1), explicit=has_occ(1))
        CALL section_vals_val_get(inp_pod_sigma, 'OCCUPIED_BETA', &
                                  r_val=occ(2), explicit=has_occ(2))
        CALL section_vals_val_get(inp_pod_sigma, 'UNOCCUPIED_ALPHA', &
                                  r_val=unocc(1), explicit=has_unocc(1))
        CALL section_vals_val_get(inp_pod_sigma, 'UNOCCUPIED_BETA', &
                                  r_val=unocc(2), explicit=has_unocc(2))

        ! Spin components
        DO j = 1, dft_cntrl%nspins

          IF (has_occ(j) .OR. has_unocc(j)) THEN

            ! States
            DO k = 1, ec%block(i)%mo(j)%nmo
              IF (has_occ(j) .AND. ec%block(i)%mo(j)%eigenvalues(k) <= ec%fermi(j)) &
                ec%block(i)%mo(j)%eigenvalues(k) = ec%block(i)%mo(j)%eigenvalues(k) + occ(j)
              IF (has_unocc(j) .AND. ec%block(i)%mo(j)%eigenvalues(k) > ec%fermi(j)) &
                ec%block(i)%mo(j)%eigenvalues(k) = ec%block(i)%mo(j)%eigenvalues(k) + unocc(j)
            END DO

          END IF

        END DO

      END DO

   END SUBROUTINE pod_corr_sigma

END MODULE et_coupling_pod_corr
