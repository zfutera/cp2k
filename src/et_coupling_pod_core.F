!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief core routines performing projection-operator based diabatization (POD)
!> \author Z. Futera (10.2021)
! **************************************************************************************************
MODULE et_coupling_pod_core

   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_equivalent,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_read_unformatted,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type,&
                                              cp_fm_write_unformatted
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_generate_filename,&
                                              cp_print_key_unit_nr
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE et_coupling_pod_mo,              ONLY: pod_mo_set
   USE et_coupling_pod_save,            ONLY: pod_save_coupl_close,&
                                              pod_save_coupl_elements,&
                                              pod_save_coupl_header,&
                                              pod_save_coupl_states,&
                                              pod_save_state_close,&
                                              pod_save_state_coeffs,&
                                              pod_save_state_header
   USE et_coupling_pod_types,           ONLY: et_cpl
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling_pod_core'

   PUBLIC :: pod_block_diag, &
             pod_block_hamiltonian

CONTAINS

! **************************************************************************************************
!> \brief reorder Hamiltonian matrix according to defined atomic blocks
!> \param ec electronic coupling data structure
!> \param mat_h the Hamiltonian matrix
!> \param mat_w working matrix of the same dimension
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_block_projection(ec, mat_h, mat_w)

      IMPLICIT NONE

      ! Routine arguments
      TYPE(et_cpl), POINTER                              :: ec
      TYPE(cp_fm_type)                                   :: mat_h
      TYPE(cp_fm_type)                                   :: mat_w

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_block_projection', &
                                     routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: ia1, ia2, ib1, ib2, id1, id2, handle
      INTEGER                                            :: is_r0, is_c0, is_r1, is_c1
      INTEGER                                            :: it_r0, it_c0, it_r1, it_c1

      ! Time mark
      CALL timeset(routineN, handle)

      ! Matrix-element reordering
      id1 = 1
      DO ib1 = 1, ec%n_blocks
        ! Atoms
        DO ia1 = 1, ec%block(ib1)%n_atoms
          ! AOs - diagonal block
          is_r0 = ec%block(ib1)%atom(ia1)%ao_pos
          is_c0 = ec%block(ib1)%atom(ia1)%ao_pos
          is_r1 = is_r0 + ec%block(ib1)%atom(ia1)%n_ao - 1
          is_c1 = is_c0 + ec%block(ib1)%atom(ia1)%n_ao - 1
          it_r0 = id1
          it_c0 = id1
          it_r1 = it_r0 + ec%block(ib1)%atom(ia1)%n_ao - 1
          it_c1 = it_c0 + ec%block(ib1)%atom(ia1)%n_ao - 1
          CALL cp_fm_to_fm_submat(mat_h,mat_w, &
                 ec%block(ib1)%atom(ia1)%n_ao,ec%block(ib1)%atom(ia1)%n_ao, &
                 is_r0,is_c0,it_r0,it_c0)
          ! Off-diagonal atom-atom AO blocks
          id2 = 1
          DO ib2 = 1, ib1
            DO ia2 = 1, ec%block(ib2)%n_atoms
              IF (ib2.LT.ib1 .OR. ia2.LT.ia1) THEN
                is_r0 = ec%block(ib2)%atom(ia2)%ao_pos
                is_c0 = ec%block(ib1)%atom(ia1)%ao_pos
                is_r1 = is_r0 + ec%block(ib2)%atom(ia2)%n_ao - 1
                is_c1 = is_c0 + ec%block(ib1)%atom(ia1)%n_ao - 1
                it_r0 = id2
                it_c0 = id1
                it_r1 = it_r0 + ec%block(ib2)%atom(ia2)%n_ao - 1
                it_c1 = it_c0 + ec%block(ib1)%atom(ia1)%n_ao - 1
                CALL cp_fm_to_fm_submat(mat_h,mat_w, &
                       ec%block(ib2)%atom(ia2)%n_ao,ec%block(ib1)%atom(ia1)%n_ao, &
                       is_r0,is_c0,it_r0,it_c0)
              END IF
              id2 = id2 + ec%block(ib2)%atom(ia2)%n_ao
            END DO
          END DO
          ! Updated target-matrix position
          id1 = id1 + ec%block(ib1)%atom(ia1)%n_ao
        END DO
      END DO

      ! Copy the reordered matrix to original data array
      CALL cp_fm_to_fm(mat_w, mat_h)

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_block_projection

! **************************************************************************************************
!> \brief transform KS hamiltonian to orthogonalized block-separated basis set
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param fm_s full-matrix structure used for allocation of KS matrices
!> \param mat_t storage for pointers to the transformed KS matrices
!> \param mat_w working matrix of the same dimension
!> \param n_ao total number of AO basis functions
!> \param n_spins number of spin components
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_block_hamiltonian(qs_env, ec, fm_s, mat_t, mat_w, n_ao, n_spins)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      TYPE(cp_fm_struct_type), POINTER                   :: fm_s
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mat_t
      TYPE(cp_fm_type)                                   :: mat_w
      INTEGER                                            :: n_ao, n_spins

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_block_hamiltonian', &
         routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: i, handle
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_h

      ! Time mark
      CALL timeset(routineN, handle)

      ! Memory allocation
      ALLOCATE (mat_t(n_spins))
      CPASSERT(ASSOCIATED(mat_t))

      ! KS Hamiltonian
      CALL get_qs_env(qs_env, matrix_ks=mat_h)

      DO i = 1, n_spins

         ! Full-matrix format
         CALL cp_fm_create(matrix=mat_t(i), matrix_struct=fm_s, &
                           name='KS HAMILTONIAN IN SEPARATED ORTHOGONALIZED BASIS SET')
         CALL copy_dbcsr_to_fm(mat_h(i)%matrix, mat_t(i))

         ! Transform KS Hamiltonian to the orthogonalized AO basis set
         CALL parallel_gemm("N", "N", n_ao, n_ao, n_ao, 1.0_dp, &
                            ec%tmatf, mat_t(i), 0.0_dp, mat_w)
         CALL parallel_gemm("N", "N", n_ao, n_ao, n_ao, 1.0_dp, &
                            mat_w, ec%tmatf, 0.0_dp, mat_t(i))

         ! Reorder KS Hamiltonian elements to defined block structure
         CALL pod_block_projection(ec, mat_t(i), mat_w)

      END DO

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_block_hamiltonian

! **************************************************************************************************
!> \brief Diagonalize diagonal blocks of the KS hamiltonian in separated orthogonalized basis set
!> \param output_unit unit number of the open output stream
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param mat_h Hamiltonian in separated orthogonalized basis set
!> \param n_spins Number of spin components
!> \param restart Read data from restart files
!> \param blocks Save Hamiltonian blocks to files, not to memory
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_block_diag(output_unit, qs_env, ec, mat_h, n_spins, restart, blocks)

      ! Routine arguments
      INTEGER                                            :: output_unit
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mat_h
      INTEGER                                            :: n_spins
      LOGICAL                                            :: restart
      LOGICAL                                            :: blocks

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_block_diag', routineP = moduleN//':'//routineN

      ! Local variables
      CHARACTER(LEN=default_path_length)                 :: filename, file_hdg
      LOGICAL                                            :: skip_i, skip_j, skip_off, skip_calc, read_hdg, ok
      INTEGER                                            :: i, id, j, k, l, spin, unit_hab, &
                                                            unit_wfn, unit_hbk, handle
      REAL(KIND=dp), DIMENSION(:), POINTER               :: vec_e
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_type), DIMENSION(:), ALLOCATABLE        :: dat
      TYPE(cp_fm_struct_type), POINTER                   :: fm_s
      TYPE(cp_fm_type)                                   :: mat_u
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: block_sec, print_sec, print_key

      ! Pointer initialization
      NULLIFY (vec_e)
      NULLIFY (blacs_env)
      NULLIFY (block_sec)
      NULLIFY (fm_s)
      NULLIFY (logger)
      NULLIFY (para_env)
      NULLIFY (print_sec)

      ! Time mark
      CALL timeset(routineN, handle)

      ! Settings
      read_hdg = .FALSE.
      IF (blocks) THEN

         ! Log file
         logger => cp_get_default_logger()
         ! Input file
         print_sec => section_vals_get_subs_vals(qs_env%input, &
                         'PROPERTIES%ET_COUPLING%PROJECTION%PRINT')
         block_sec => section_vals_get_subs_vals(qs_env%input, &
                         'PROPERTIES%ET_COUPLING%PROJECTION%BLOCK')

         ! Check whether the POD blocks are available
         read_hdg = .TRUE.
         print_key => section_vals_get_subs_vals(print_sec, 'HAMILTONIAN')
         DO i = 1, ec%n_blocks
            CALL section_vals_val_get(block_sec, keyword_name='SKIP', l_val=skip_i, i_rep_section=i)
            IF (skip_i) THEN
               read_hdg = .FALSE.
               EXIT
            ELSE
               DO j = 1, n_spins
                  WRITE (filename, '(A8,I1.1,A1,I1.1)') 'RESTART_', i, '_', j
                  file_hdg = cp_print_key_generate_filename(logger, print_key, middle_name=TRIM(filename), &
                                extension='.hdg', my_local=.FALSE.)
                  INQUIRE (file=file_hdg, exist=ok)
                  IF (.NOT. ok) THEN
                     read_hdg = .FALSE.
                     EXIT
                  END IF
               END DO
            END IF
         END DO
         ! Check whether the off-diagonal blocks are available
         IF (read_hdg) THEN
            DO i = 1, ec%n_blocks
               DO j = i+1, ec%n_blocks
                  DO k = 1, n_spins
                     WRITE (filename, '(A8,I1.1,A1,I1.1,A1,I1.1)') 'RESTART_', i, '-', j, '_', k
                     file_hdg = cp_print_key_generate_filename(logger, print_key, middle_name=TRIM(filename), &
                                   extension='.hbk', my_local=.FALSE.)
                     INQUIRE (file=file_hdg, exist=ok)
                     IF (.NOT. ok) THEN
                        read_hdg = .FALSE.
                        EXIT
                     END IF
                  END DO
               END DO
            END DO
         END IF

      END IF

      ! Parallel environment
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)

      ! Storage for block sub-matrices
      ALLOCATE (dat(ec%n_blocks))
      CPASSERT(ALLOCATED(dat))

      ! Storage for couplings and wavefunctions
      DO i = 1, ec%n_blocks
         ! Memory for states
         CPASSERT(.NOT. ASSOCIATED(ec%block(i)%mo))
         ALLOCATE (ec%block(i)%mo(n_spins))
         CPASSERT(ASSOCIATED(ec%block(i)%mo))
         ! Memory for couplings
         ec%block(i)%n_hab_blocks = ec%n_blocks - i
         CPASSERT(.NOT. ASSOCIATED(ec%block(i)%hab))
         IF (ec%block(i)%n_hab_blocks > 0) THEN
            ALLOCATE (ec%block(i)%hab(n_spins, ec%block(i)%n_hab_blocks))
            CPASSERT(ASSOCIATED(ec%block(i)%hab))
         END IF
      END DO

      ! Whole-system restart files - header
      IF (.NOT. blocks .OR. read_hdg) THEN
         CALL pod_save_coupl_header(qs_env, ec, n_spins, unit_hab); 
         CALL pod_save_state_header(qs_env, ec, n_spins, unit_wfn); 
      END IF

      ! Spin components
      DO spin = 1, n_spins

         skip_off = .FALSE.

         ! Diagonal blocks
         j = 1
         DO i = 1, ec%n_blocks

            skip_calc = .FALSE.

            ! Memory allocation
            CALL cp_fm_struct_create(fmstruct=fm_s, para_env=para_env, context=blacs_env, &
                                     nrow_global=ec%block(i)%n_ao, ncol_global=ec%block(i)%n_ao)
            CALL cp_fm_create(matrix=dat(i), matrix_struct=fm_s, &
                              name='H_KS DIAGONAL BLOCK')

            ! Save block to file only
            IF (blocks) THEN

               CALL section_vals_val_get(block_sec, keyword_name='SKIP', l_val=skip_i, i_rep_section=i)

               IF (.NOT. skip_i) THEN
               
                  WRITE (filename, '(A8,I1.1,A1,I1.1,A1,I1.1)') 'RESTART_', i, '-', i, '_', spin

                  ! Read block data from the file
                  IF (restart) THEN

                     ! Diagonalized block
                     IF (read_hdg) THEN

                        skip_calc = .TRUE.

                        IF (output_unit > 0) &
                           WRITE (output_unit, '(/,T3,A)') 'Reading block '//&
                           TRIM(ADJUSTL(cp_to_string(i)))//' energies and MOs from restart file'

                        WRITE (filename, '(A8,I1.1,A1,I1.1)') 'RESTART_', i, '_', spin
      
                        unit_hbk = cp_print_key_unit_nr(logger, print_sec, 'HAMILTONIAN', &
                                      extension='.hdg', middle_name=TRIM(filename), file_form='UNFORMATTED', &
                                      file_position='REWIND', file_action='READ', log_filename=.FALSE.)

                        CPASSERT(.NOT. ASSOCIATED(vec_e))
                        ALLOCATE (vec_e(ec%block(i)%n_ao))
                        CPASSERT(ASSOCIATED(vec_e))
                        IF (unit_hbk > 0) THEN
                           READ (UNIT=unit_hbk) (vec_e(k), k=1, ec%block(i)%n_ao)
                        END IF

                        CALL cp_fm_read_unformatted(dat(i), unit_hbk)

                        CALL pod_mo_set(qs_env, ec, i, spin, dat(i), vec_e)

                        CALL pod_save_coupl_states(ec, i, spin, unit_hab)
                        CALL pod_save_state_coeffs(ec, i, spin, unit_wfn)

                        IF (ASSOCIATED(vec_e)) &
                           DEALLOCATE (vec_e)
                        NULLIFY (vec_e)

                     ! Non-diagonalized block
                     ELSE

                        skip_off = .TRUE.
                  
                        IF (output_unit > 0) &
                           WRITE (output_unit, '(/,T3,A)') 'Reading Hamiltonian block '//&
                           TRIM(ADJUSTL(cp_to_string(i)))//' data from restart file'
                  
                        unit_hbk = cp_print_key_unit_nr(logger, print_sec, 'HAMILTONIAN', &
                                      extension='.hbk', middle_name=TRIM(filename), file_form='UNFORMATTED', &
                                      file_position='REWIND', file_action='READ', log_filename=.FALSE.)
                        CALL cp_fm_read_unformatted(dat(i), unit_hbk)

                     END IF
               
                  ! Save Hamiltonian block to file
                  ELSE

                     skip_calc = .TRUE.
               
                     unit_hbk = cp_print_key_unit_nr(logger, print_sec, 'HAMILTONIAN', &
                                   extension='.hbk', middle_name=TRIM(filename), file_form='UNFORMATTED', &
                                   file_position='REWIND', file_action='WRITE', file_status='REPLACE', &
                                   log_filename=.FALSE.)
                     CALL cp_fm_to_fm_submat(mat_h(spin), dat(i), ec%block(i)%n_ao, &
                             ec%block(i)%n_ao, j, j, 1, 1)
                     CALL cp_fm_write_unformatted(dat(i), unit_hbk)
               
                  END IF
               
                  CALL cp_print_key_finished_output(unit_hbk, logger, print_sec, 'HAMILTONIAN')

               ELSE
           
                 skip_calc = .TRUE.
               
               END IF

            END IF

            ! Energy calculation
            IF (.NOT. skip_calc) THEN

               ! Block elements
               IF (.NOT. blocks) THEN
                  CALL cp_fm_to_fm_submat(mat_h(spin), dat(i), ec%block(i)%n_ao, &
                          ec%block(i)%n_ao, j, j, 1, 1)
               END IF

               ! Diagonalization
               CPASSERT(.NOT. ASSOCIATED(vec_e))
               ALLOCATE (vec_e(ec%block(i)%n_ao))
               CPASSERT(ASSOCIATED(vec_e))

               CALL cp_fm_create(matrix=mat_u, matrix_struct=fm_s, name='UNITARY MATRIX')
               CALL choose_eigv_solver(dat(i), mat_u, vec_e)
               CALL cp_fm_to_fm(mat_u, dat(i))

               ! Save state energies / vectors
               CALL pod_mo_set(qs_env, ec, i, spin, mat_u, vec_e)

               ! Save energies and MO coefficients to restart file
               IF (.NOT. blocks) THEN
                  CALL pod_save_coupl_states(ec, i, spin, unit_hab)
                  CALL pod_save_state_coeffs(ec, i, spin, unit_wfn)
               ELSE

                  IF (output_unit > 0) &
                     WRITE (output_unit, '(T3,A)') 'Saving block '//&
                     TRIM(ADJUSTL(cp_to_string(i)))//' energies and MOs to file'

                  WRITE (filename, '(A8,I1.1,A1,I1.1)') 'RESTART_', i, '_', spin

                  unit_hbk = cp_print_key_unit_nr(logger, print_sec, 'HAMILTONIAN', &
                                extension='.hdg', middle_name=TRIM(filename), file_form='UNFORMATTED', &
                                file_position='REWIND', file_action='WRITE', file_status='REPLACE', &
                                log_filename=.FALSE.)
                  IF (unit_hbk > 0) then
                     WRITE (UNIT=unit_hbk) (vec_e(k), k=1, ec%block(i)%n_ao)
                  END IF
                  CALL cp_fm_write_unformatted(mat_u, unit_hbk)
                  CALL cp_print_key_finished_output(unit_hbk, logger, print_sec, 'HAMILTONIAN')

               END IF

               ! Clean memory
               CALL cp_fm_release(matrix=mat_u)
               IF (ASSOCIATED(vec_e)) &
                  DEALLOCATE (vec_e)
               NULLIFY (vec_e)

            END IF

            ! Clean memory
            CALL cp_fm_struct_release(fmstruct=fm_s)
            NULLIFY (fm_s)

            ! Off-set for next block
            j = j + ec%block(i)%n_ao

         END DO

         ! Off-diagonal blocks
         IF (.NOT. skip_off) THEN
            k = 1
            DO i = 1, ec%n_blocks
               l = 1
               id = 1
               DO j = 1, ec%n_blocks
                  IF (j > i) THEN
  
                     ! Memory allocation
                     CALL cp_fm_struct_create(fmstruct=fm_s, para_env=para_env, context=blacs_env, &
                                              nrow_global=ec%block(i)%n_ao, ncol_global=ec%block(j)%n_ao)
                     CALL cp_fm_create(matrix=ec%block(i)%hab(spin, id), matrix_struct=fm_s, &
                                       name='H_KS OFF-DIAGONAL BLOCK')
  
                     ! Read block elements from restart file
                     IF (read_hdg) THEN

                        WRITE (filename, '(A8,I1.1,A1,I1.1,A1,I1.1)') 'RESTART_', i, '-', j, '_', spin
                        unit_hbk = cp_print_key_unit_nr(logger, print_sec, 'HAMILTONIAN', &
                                      extension='.hbk', middle_name=TRIM(filename), file_form='UNFORMATTED', &
                                      file_position='REWIND', file_action='READ', log_filename=.FALSE.)
                        CALL cp_fm_read_unformatted(ec%block(i)%hab(spin, id), unit_hbk)
                        CALL cp_print_key_finished_output(unit_hbk, logger, print_sec, 'HAMILTONIAN')

                     ! Copy block data from Hamiltonian
                     ELSE

                        CALL cp_fm_to_fm_submat(mat_h(spin), &
                                                ec%block(i)%hab(spin, id), ec%block(i)%n_ao, &
                                                ec%block(j)%n_ao, k, l, 1, 1)

                        ! Save block to file
                        IF (blocks) THEN
    
                           CALL section_vals_val_get(block_sec, keyword_name='SKIP', l_val=skip_i, i_rep_section=i)
                           CALL section_vals_val_get(block_sec, keyword_name='SKIP', l_val=skip_j, i_rep_section=j)
    
                           IF (.NOT. skip_i .AND. .NOT. skip_j) THEN
                              ! open the file
                              WRITE (filename, '(A8,I1.1,A1,I1.1,A1,I1.1)') 'RESTART_', i, '-', j, '_', spin
                              unit_hbk = cp_print_key_unit_nr(logger, print_sec, 'HAMILTONIAN', &
                                            extension='.hbk', middle_name=TRIM(filename), file_form='UNFORMATTED', &
                                            file_position='REWIND', file_action='WRITE', file_status='REPLACE', &
                                            log_filename=.FALSE.)
                              CALL cp_fm_write_unformatted(ec%block(i)%hab(spin, id), unit_hbk)
                              CALL cp_print_key_finished_output(unit_hbk, logger, print_sec, 'HAMILTONIAN')
                              CALL cp_fm_release(matrix=ec%block(i)%hab(spin, id))
                           END IF
    
                        END IF

                     END IF
 
                     ! Transformation
                     IF (.NOT. blocks .OR. read_hdg) THEN
 
                        CALL cp_fm_create(matrix=mat_u, matrix_struct=fm_s, name='FULL WORK MATRIX')
                        CALL parallel_gemm("T", "N", ec%block(i)%n_ao, ec%block(j)%n_ao, ec%block(i)%n_ao, &
                                           1.0_dp, dat(i), ec%block(i)%hab(spin, id), 0.0_dp, mat_u)
                        CALL parallel_gemm("N", "N", ec%block(i)%n_ao, ec%block(j)%n_ao, ec%block(j)%n_ao, &
                                           1.0_dp, mat_u, dat(j), 0.0_dp, ec%block(i)%hab(spin, id))
 
                        ! Save coupling elements to binary data file
                        CALL pod_save_coupl_elements(ec, i, id, spin, unit_hab); 
 
                        ! Clean memory
                        CALL cp_fm_release(matrix=mat_u)
 
                     END IF
  
                     ! Clean memory
                     CALL cp_fm_struct_release(fmstruct=fm_s)
                     NULLIFY (fm_s)
  
                     id = id + 1
  
                  END IF
                  ! Off-set for next block
                  l = l + ec%block(j)%n_ao
               END DO
               ! Off-set for next block
               k = k + ec%block(i)%n_ao
            END DO
         END IF

         ! Clean memory
         IF (ALLOCATED(dat)) THEN
            DO i = 1, SIZE(dat)
               CALL cp_fm_release(matrix=dat(i))
            END DO
         END IF
      END DO

      ! Clean memory
      IF (ALLOCATED(dat)) &
         DEALLOCATE (dat)

      ! Close output streams
      IF (.NOT. blocks .OR. read_hdg) THEN
         CALL pod_save_coupl_close(qs_env, unit_hab)
         CALL pod_save_state_close(qs_env, unit_wfn)
      END IF

      ! Clean memory
      IF (ASSOCIATED(mat_h)) THEN
         DO i = 1, SIZE(mat_h)
            CALL cp_fm_release(matrix=mat_h(i))
         END DO
         IF (ASSOCIATED(mat_h)) &
            DEALLOCATE (mat_h)
      END IF

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_block_diag

END MODULE et_coupling_pod_core
