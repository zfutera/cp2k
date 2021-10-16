!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief calculates the electron transfer coupling elements by projection-operator approach
!>        Kondov et al. J.Phys.Chem.C 2007, 111, 11970-11981
!> \author Z. Futera (02.2017)
! **************************************************************************************************
MODULE et_coupling_pod

   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: Futera2017,&
                                              cite_reference
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE et_coupling_pod_core,            ONLY: pod_block_diag,&
                                              pod_block_hamiltonian
   USE et_coupling_pod_print,           ONLY: pod_print_couplings,&
                                              pod_print_state_fracs,&
                                              pod_print_states
   USE et_coupling_pod_read,            ONLY: pod_read_restart_data
   USE et_coupling_pod_save,            ONLY: pod_save_state_cube
   USE et_coupling_pod_types,           ONLY: et_cpl,&
                                              pod_data_release
   USE et_coupling_pod_utils,           ONLY: pod_check_restart,&
                                              pod_get_transf_mat
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE orbital_pointers,                ONLY: nso
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: mo_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling_pod'

   PUBLIC :: pod_coupling_calc

CONTAINS

! **************************************************************************************************
!> \brief check the electronic-coupling input section and set the atomic block data
!> \param ec electronic coupling data structure
!> \param fa system Fermi level (alpha spin)
!> \param fb system Fermi level (beta spin)
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_set_fermi(ec, fa, fb)

      ! Routine arguments
      TYPE(et_cpl), POINTER                              :: ec
      REAL(KIND=dp)                                      :: fa
      REAL(KIND=dp), OPTIONAL                            :: fb

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_set_fermi', routineP = moduleN//':'//routineN

      NULLIFY (ec%fermi)

      IF (PRESENT(fb)) THEN

         ALLOCATE (ec%fermi(2))
         CPASSERT(ASSOCIATED(ec%fermi))
         ec%fermi(1) = fa
         ec%fermi(2) = fb

      ELSE

         ALLOCATE (ec%fermi(1))
         CPASSERT(ASSOCIATED(ec%fermi))
         ec%fermi(1) = fa

      END IF

   END SUBROUTINE pod_set_fermi

! **************************************************************************************************
!> \brief check the electronic-coupling input section and set the atomic block data
!> \param qs_env QuickStep environment containing all system data
!> \param et_proj_sec the electronic-coupling input section
!> \param ec electronic coupling data structure
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_block_init(qs_env, et_proj_sec, ec)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: et_proj_sec
      TYPE(et_cpl), POINTER                              :: ec

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_block_init', routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: i, j, k, l, n, n_ao, n_atoms, n_set, handle
      INTEGER, DIMENSION(:), POINTER                     :: atom_id, atom_nf, atom_ps, n_shell, t
      INTEGER, DIMENSION(:, :), POINTER                  :: ang_mom_id
      TYPE(gto_basis_set_type), POINTER                  :: ao_basis_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: block_sec

      ! Time mark
      CALL timeset(routineN, handle)

      ! Initialization
      NULLIFY (ang_mom_id)
      NULLIFY (ao_basis_set)
      NULLIFY (atom_id)
      NULLIFY (atom_nf)
      NULLIFY (atom_ps)
      NULLIFY (block_sec)
      NULLIFY (n_shell)
      NULLIFY (particle_set)
      NULLIFY (qs_kind_set)
      NULLIFY (t)

      ! Initialization
      ec%n_atoms = 0
      ec%n_blocks = 0
      NULLIFY (ec%fermi)
      NULLIFY (ec%tmatf)
      NULLIFY (ec%tmatr)
      NULLIFY (ec%block)

      ! Number of atoms / atom types
      CALL get_qs_env(qs_env, particle_set=particle_set, qs_kind_set=qs_kind_set, natom=n_atoms)
      ! Number of AO basis functions
      CALL get_qs_kind_set(qs_kind_set, nsgf=n_ao)

      ! Number of AO functions per atom
      ALLOCATE (atom_nf(n_atoms))
      IF (.NOT. ASSOCIATED(atom_nf)) &
         CPABORT('Cannot allocate working "atom_nf" array. ')

      atom_nf = 0
      DO i = 1, n_atoms
         CALL get_atomic_kind(particle_set(i)%atomic_kind, kind_number=j)
         CALL get_qs_kind(qs_kind_set(j), basis_set=ao_basis_set)
         IF (.NOT. ASSOCIATED(ao_basis_set)) &
            CPABORT('Unsupported basis set type. ')
         CALL get_gto_basis_set(gto_basis_set=ao_basis_set, &
                                nset=n_set, nshell=n_shell, l=ang_mom_id)
         DO j = 1, n_set
            DO k = 1, n_shell(j)
               atom_nf(i) = atom_nf(i) + nso(ang_mom_id(k, j))
            END DO
         END DO
      END DO

      ! Sanity check
      n = 0
      DO i = 1, n_atoms
         n = n + atom_nf(i)
      END DO
      IF (n /= n_ao) &
         CPABORT('Inconsistent number of atomic orbitals. ')

      ! Atom position in AO array
      ALLOCATE (atom_ps(n_atoms))
      IF (.NOT. ASSOCIATED(atom_ps)) &
         CPABORT('Cannot allocate working "atom_ps" array. ')
      atom_ps = 1
      DO i = 1, n_atoms - 1
         atom_ps(i + 1) = atom_ps(i) + atom_nf(i)
      END DO

      ! Number of blocks
      block_sec => section_vals_get_subs_vals(et_proj_sec, 'BLOCK')
      CALL section_vals_get(block_sec, n_repetition=ec%n_blocks)
      ALLOCATE (ec%block(ec%n_blocks))
      IF (.NOT. ASSOCIATED(ec%block)) &
         CPABORT('Cannot allocate memory for POD block data. ')

      ! Block data
      ALLOCATE (t(n_atoms))
      IF (.NOT. ASSOCIATED(t)) &
         CPABORT('Cannot allocate memory for working "t" array. ')

      ec%n_atoms = 0
      DO i = 1, ec%n_blocks

         ! Initialization
         ec%block(i)%n_atoms = 0
         ec%block(i)%n_electrons = 0
         ec%block(i)%n_ao = 0
         NULLIFY (ec%block(i)%atom)
         NULLIFY (ec%block(i)%mo)
         NULLIFY (ec%block(i)%hab)

         ! Number of electrons
         CALL section_vals_val_get(block_sec, i_rep_section=i, &
                                   keyword_name='NELECTRON', i_val=ec%block(i)%n_electrons)

         ! User-defined atom array
         CALL section_vals_val_get(block_sec, i_rep_section=i, &
                                   keyword_name='ATOMS', i_vals=atom_id)

         ! Count unique atoms
         DO j = 1, SIZE(atom_id)
            ! Check atom ID validity
            IF (atom_id(j) < 1 .OR. atom_id(j) > n_atoms) &
               CPABORT('invalid fragment atom ID ('//TRIM(ADJUSTL(cp_to_string(atom_id(j))))//')')
            ! Check if the atom is not in previously-defined blocks
            DO k = 1, i - 1
               DO l = 1, ec%block(k)%n_atoms
                  IF (ec%block(k)%atom(l)%id == atom_id(j)) &
                     CPABORT('multiple assignment of atom'//TRIM(ADJUSTL(cp_to_string(atom_id(j)))))
               END DO
            END DO
            ! Check if the atom is not duplicated in the present block
            DO k = 1, ec%block(i)%n_atoms
               IF (t(k) == atom_id(j)) &
                  CPABORT('multiple assignment of atom'//TRIM(ADJUSTL(cp_to_string(atom_id(j)))))
            END DO
            ! Save the atom
            ec%block(i)%n_atoms = ec%block(i)%n_atoms + 1
            t(ec%block(i)%n_atoms) = atom_id(j)
         END DO

         ! Memory allocation
         ALLOCATE (ec%block(i)%atom(ec%block(i)%n_atoms))
         IF (.NOT. ASSOCIATED(ec%block(i)%atom)) &
            CPABORT('Cannot allocate memory for atom IDs')

         ! Save atom IDs and number of AOs
         DO j = 1, ec%block(i)%n_atoms
            ec%block(i)%atom(j)%id = t(j)
            ec%block(i)%atom(j)%n_ao = atom_nf(ec%block(i)%atom(j)%id)
            ec%block(i)%atom(j)%ao_pos = atom_ps(ec%block(i)%atom(j)%id)
            ec%block(i)%n_ao = ec%block(i)%n_ao + ec%block(i)%atom(j)%n_ao
         END DO

         ec%n_atoms = ec%n_atoms + ec%block(i)%n_atoms
      END DO

      ! Clean memory
      IF (ASSOCIATED(atom_nf)) &
         DEALLOCATE (atom_nf)
      IF (ASSOCIATED(atom_ps)) &
         DEALLOCATE (atom_ps)
      IF (ASSOCIATED(t)) &
         DEALLOCATE (t)

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_block_init

! **************************************************************************************************
!> \brief calculates the electron transfer coupling elements by projection-operator approach
!>        Kondov et al. J.Phys.Chem.C 2007, 111, 11970-11981
!> \param qs_env QuickStep environment containing all system data
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_coupling_calc(qs_env)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_coupling_calc', &
         routineP = moduleN//':'//routineN

      ! Local variables
      CHARACTER(LEN=default_path_length)                 :: file_hab, file_tmf, file_tmr, file_wfn
      INTEGER                                            :: i, j, k, n_ao, n_atoms, output_unit, handle
      LOGICAL                                            :: do_kp, restart, blocks_only
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mat_h
      TYPE(cp_fm_struct_type), POINTER                   :: fm_s
      TYPE(cp_fm_type), POINTER                          :: mat_w
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_cntrl
      TYPE(et_cpl), POINTER                              :: ec
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: et_proj_sec

      ! Pointer initialization
      NULLIFY (logger)

      NULLIFY (blacs_env)
      NULLIFY (para_env)
      NULLIFY (dft_cntrl)
      NULLIFY (kpoints)
      NULLIFY (mo)
      NULLIFY (qs_kind_set)
      NULLIFY (et_proj_sec)

      NULLIFY (fm_s)
      NULLIFY (mat_h)
      NULLIFY (mat_w)

      NULLIFY (ec)

      ! Time mark
      CALL timeset(routineN, handle)

      ! Reference
      CALL cite_reference(Futera2017)

      ! Stream for output to LOG file
      logger => cp_get_default_logger()

      et_proj_sec => section_vals_get_subs_vals(qs_env%input, 'PROPERTIES%ET_COUPLING%PROJECTION')

      output_unit = cp_print_key_unit_nr(logger, et_proj_sec, &
                                         'PROGRAM_RUN_INFO', extension='.log')

      ! Restart from data files
      restart = pod_check_restart(logger, et_proj_sec, &
                   file_tmf, file_tmr, file_hab, file_wfn, blocks_only)

      ! Header
      IF (output_unit > 0) THEN
         WRITE (output_unit, '(/,T2,A)') &
            '!-----------------------------------------------------------------------------!'
         WRITE (output_unit, '(T17,A)') &
            'Electronic coupling - Projection-operator method'
      END IF

      ! Main data structure
      ALLOCATE (ec)
      CPASSERT(ASSOCIATED(ec))
      CALL pod_block_init(qs_env, et_proj_sec, ec)

      ! Number of atoms and AO functions
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, natom=n_atoms)
      CALL get_qs_kind_set(qs_kind_set, nsgf=n_ao)

      ! Print out info about system partitioning
      IF (output_unit > 0) THEN

         WRITE (output_unit, '(/,T3,A,I10)') &
            'Number of atoms                    = ', n_atoms
         WRITE (output_unit, '(T3,A,I10)') &
            'Number of fragments                = ', ec%n_blocks
         WRITE (output_unit, '(T3,A,I10)') &
            'Number of fragment atoms           = ', ec%n_atoms
         WRITE (output_unit, '(T3,A,I10)') &
            'Number of unassigned atoms         = ', n_atoms - ec%n_atoms
         WRITE (output_unit, '(T3,A,I10)') &
            'Number of AO basis functions       = ', n_ao

         DO i = 1, ec%n_blocks

            WRITE (output_unit, '(/,T3,A,I0,A)') &
               'Block ', i, ':'
            WRITE (output_unit, '(T3,A,I10)') &
               'Number of block atoms              = ', ec%block(i)%n_atoms
            WRITE (output_unit, '(T3,A,I10)') &
               'Number of block electrons          = ', ec%block(i)%n_electrons
            WRITE (output_unit, '(T3,A,I10)') &
               'Number of block AO functions       = ', ec%block(i)%n_ao

            IF (ec%block(i)%n_atoms < 10) THEN

               WRITE (output_unit, '(T3,A,10I6)') &
                  'Block atom IDs                     =     ', &
                  (ec%block(i)%atom(j)%id, j=1, ec%block(i)%n_atoms)

            ELSE

               WRITE (output_unit, '(T3,A)') 'Block atom IDs                     ='
               DO j = 1, ec%block(i)%n_atoms/10
                  WRITE (output_unit, '(T3,A,10I6)') '      ', &
                     (ec%block(i)%atom((j - 1)*10 + k)%id, k=1, 10)
               END DO
               IF (MOD(ec%block(i)%n_atoms, 10) /= 0) THEN
                  WRITE (output_unit, '(T3,A,10I6)') '      ', &
                     (ec%block(i)%atom(k + 10*(ec%block(i)%n_atoms/10))%id, &
                      k=1, MOD(ec%block(i)%n_atoms, 10))
               END IF

            END IF

         END DO

         IF (restart) THEN
            WRITE (output_unit, '(/,T3,A)') &
               'POD data are read from restart files'
            IF (file_tmf /= '') &
               WRITE (output_unit, '(T3,A,A)') &
                  'Forward transformation matrix      = ', TRIM(ADJUSTL(file_tmf))
            IF (file_tmr /= '') &
               WRITE (output_unit, '(T3,A,A)') &
                  'Reverse transformation matrix      = ', TRIM(ADJUSTL(file_tmr))
            IF (file_hab /= '') &
               WRITE (output_unit, '(T3,A,A)') &
                  'Hamiltonian data file              = ', TRIM(ADJUSTL(file_hab))
            IF (file_wfn /= '') &
               WRITE (output_unit, '(T3,A,A)') &
                  'Wavefunction data file             = ', TRIM(ADJUSTL(file_wfn))
         END IF

      END IF

      ! Full matrix data structure
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      CALL cp_fm_struct_create(fmstruct=fm_s, para_env=para_env, context=blacs_env, &
                               nrow_global=n_ao, ncol_global=n_ao)
      CALL cp_fm_create(matrix=mat_w, matrix_struct=fm_s, name='FULL WORK MATRIX')

      ! Spin polarization / K-point sampling
      CALL get_qs_env(qs_env, dft_control=dft_cntrl, do_kpoints=do_kp)

      ! Spin polarized wavefunction
      IF (dft_cntrl%nspins == 2) THEN

         IF (output_unit > 0) &
            WRITE (output_unit, '(/,T3,A)') 'Spin-polarized calculation'

         IF (do_kp) THEN

            !<--- Open shell / K points --------------------------------------------------->!

            CALL get_qs_env(qs_env, kpoints=kpoints)
            IF (output_unit > 0) &
               WRITE (output_unit, '(T3,A,I10)') 'Number of K-points                 = ', &
               kpoints%nkp

         ELSE

            !<--- Open shell / No K-points ------------------------------------------------>!

            ! Open shell, no K-points
            IF (output_unit > 0) &
               WRITE (output_unit, '(T3,A)') 'No K-point sampling (Gamma point only)'

            ! State eneries of the whole system
            CALL get_qs_env(qs_env, mos=mo)
            IF (mo(1)%mo_set%nao /= mo(2)%mo_set%nao) &
               CPABORT('different number of alpha/beta AO basis functions')
            IF (output_unit > 0) THEN
               WRITE (output_unit, '(/,T3,A,I10)') &
                  'Number of AO basis funtions        = ', mo(1)%mo_set%nao
               WRITE (output_unit, '(T3,A,I10)') &
                  'Number of alpha states             = ', mo(1)%mo_set%nmo
               WRITE (output_unit, '(T3,A,I10)') &
                  'Number of beta states              = ', mo(2)%mo_set%nmo
            END IF
            CALL pod_print_states(output_unit, mo, dft_cntrl%nspins, &
                                  'the whole system', fermi=.TRUE.)
            CALL pod_set_fermi(ec, mo(1)%mo_set%mu, mo(2)%mo_set%mu)

            ! Read data from files
            IF (restart) THEN

               ! Hamiltonian & wavefunction
               CALL pod_read_restart_data(output_unit, qs_env, ec, dft_cntrl%nspins, fm_s, &
                                          mat_w, n_ao, file_tmf, file_tmr, file_hab, file_wfn, &
                                          blocks_only)

               ! POD calculation
            ELSE

               ! Lowdin transformation matrices
               CALL pod_get_transf_mat(output_unit, qs_env, ec, fm_s, mat_w, n_ao)
               ! KS Hamiltonian
               CALL pod_block_hamiltonian(qs_env, ec, fm_s, mat_h, mat_w, &
                                          n_ao, dft_cntrl%nspins)
               ! Block diagonization
               CALL pod_block_diag(output_unit, qs_env, ec, mat_h, dft_cntrl%nspins, &
                                   restart, blocks_only)

            END IF

            ! Print out energies and couplings
            DO i = 1, ec%n_blocks
               IF (output_unit > 0) THEN
                  CALL pod_print_states(output_unit, ec%block(i)%mo, dft_cntrl%nspins, &
                                        'block '//TRIM(ADJUSTL(cp_to_string(i)))//' states', &
                                        mx_mo_a=mo(1)%mo_set%nmo, mx_mo_b=mo(2)%mo_set%nmo, &
                                        fermi=.TRUE.)
               END IF
               CALL pod_print_state_fracs(output_unit, qs_env, ec, i, dft_cntrl%nspins)
            END DO

            CALL pod_print_couplings(output_unit, qs_env, ec, &
                                     mo(1)%mo_set%nmo, mo(2)%mo_set%nmo)

            ! Save electronic states
            CALL pod_save_state_cube(qs_env, ec, dft_cntrl%nspins)

         END IF

         ! Spin non-polarized wavefunction
      ELSE

         IF (output_unit > 0) &
            WRITE (output_unit, '(/,T3,A)') 'Spin-restricted calculation'

         IF (do_kp) THEN

            !<--- Close shell / K-points -------------------------------------------------->!

            CALL get_qs_env(qs_env, kpoints=kpoints)
            IF (output_unit > 0) &
               WRITE (output_unit, '(T3,A,I10)') 'Number of K-points                 = ', &
               kpoints%nkp

         ELSE

            !<--- Closed shell / No K-points ---------------------------------------------->!

            ! Closed shell, no K-points
            IF (output_unit > 0) &
               WRITE (output_unit, '(T3,A)') 'No K-point sampling (Gamma point only)'

            ! State eneries of the whole system
            CALL get_qs_env(qs_env, mos=mo)
            IF (output_unit > 0) THEN
               WRITE (output_unit, '(/,T3,A,I10)') &
                  'Number of AO basis funtions        = ', mo(1)%mo_set%nao
               WRITE (output_unit, '(T3,A,I10)') &
                  'Number of states                   = ', mo(1)%mo_set%nmo
            END IF
            CALL pod_print_states(output_unit, mo, dft_cntrl%nspins, &
                                  'the whole system', fermi=.TRUE.)
            CALL pod_set_fermi(ec, mo(1)%mo_set%mu)

            ! Read data from files
            IF (restart) THEN

               ! Hamiltonian & wavefunction
               CALL pod_read_restart_data(output_unit, qs_env, ec, dft_cntrl%nspins, fm_s, &
                                          mat_w, n_ao, file_tmf, file_tmr, file_hab, file_wfn, &
                                          blocks_only)

               ! POD calculation
            ELSE

               ! Lowdin transformation matrices
               CALL pod_get_transf_mat(output_unit, qs_env, ec, fm_s, mat_w, n_ao)
               ! KS Hamiltonian
               CALL pod_block_hamiltonian(qs_env, ec, fm_s, mat_h, mat_w, &
                                          n_ao, dft_cntrl%nspins)
               ! Block diagonization
               CALL pod_block_diag(output_unit, qs_env, ec, mat_h, dft_cntrl%nspins, &
                                   restart, blocks_only)

            END IF

            ! Print out energies 
            DO i = 1, ec%n_blocks
               IF (ASSOCIATED(ec%block(i)%mo(1)%mo_set) .AND. output_unit > 0) THEN
                  CALL pod_print_states(output_unit, ec%block(i)%mo, dft_cntrl%nspins, &
                                        'block '//TRIM(ADJUSTL(cp_to_string(i)))//' states', &
                                        mx_mo_a=mo(1)%mo_set%nmo, fermi=.TRUE.)
               END IF
               CALL pod_print_state_fracs(output_unit, qs_env, ec, i, dft_cntrl%nspins)
            END DO

            ! Print out coupling elements
            CALL pod_print_couplings(output_unit, qs_env, ec, mo(1)%mo_set%nmo)

            ! Save electronic states
            CALL pod_save_state_cube(qs_env, ec, dft_cntrl%nspins)

            ! Save electronic states
            CALL pod_save_state_cube(qs_env, ec, dft_cntrl%nspins)

            ! Save electronic states
            CALL pod_save_state_cube(qs_env, ec, dft_cntrl%nspins)

            !<----------------------------------------------------------------------------->!

         END IF

      END IF

      ! Footer
      IF (output_unit > 0) &
         WRITE (output_unit, '(/,T2,A)') &
         '!-----------------------------------------------------------------------------!'

      ! Clean memory
      CALL cp_fm_struct_release(fmstruct=fm_s)
      CALL cp_fm_release(matrix=mat_w)
      CALL pod_data_release(ec)

      ! Close output stream
      CALL cp_print_key_finished_output(output_unit, logger, et_proj_sec, 'PROGRAM_RUN_INFO')

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_coupling_calc

END MODULE et_coupling_pod
