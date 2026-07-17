!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines saving to files various data related to POD method
!> \author Z. Futera (03.2019)
! **************************************************************************************************
MODULE et_coupling_pod_save

   USE cp_control_types,                ONLY: dft_control_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_types,                     ONLY: cp_fm_write_unformatted
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_generate_filename,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE et_coupling_pod_mo,              ONLY: pod_mo_save
   USE et_coupling_pod_types,           ONLY: et_cpl
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: mo_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling_pod_save'

   PUBLIC :: pod_save_coupl_close, &
             pod_save_coupl_elements, &
             pod_save_coupl_header, &
             pod_save_coupl_states, &
             pod_save_state_close, &
             pod_save_state_coeffs, &
             pod_save_state_cube, &
             pod_save_state_header, &
             pod_save_states, &
             pod_save_transf_mat

CONTAINS

! **************************************************************************************************
!> \brief save header of binary POD coupling file
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param n_spins number of spin components
!> \param unit_nr binary-file output unit
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_save_coupl_header(qs_env, ec, n_spins, unit_nr)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: n_spins
      INTEGER, INTENT(OUT)                               :: unit_nr

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_coupl_header', &
         routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: i, handle
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key


      ! Initialization
      NULLIFY (logger)
      NULLIFY (print_key)

      ! Time mark
      CALL timeset(routineN, handle)

      ! Output settings
      logger => cp_get_default_logger()
      print_key => section_vals_get_subs_vals(qs_env%input, &
                                              'PROPERTIES%ET_COUPLING%PROJECTION%PRINT')

      ! Open the file
      unit_nr = cp_print_key_unit_nr(logger, print_key, 'HAMILTONIAN', &
                                     middle_name='RESTART', extension='.hab', file_form='UNFORMATTED', &
                                     file_position='REWIND', file_action='WRITE', file_status='REPLACE')

      IF (unit_nr > 0) THEN

         ! Number of spin components
         WRITE (UNIT=unit_nr) n_spins
         ! Number of blocks
         WRITE (UNIT=unit_nr) ec%n_blocks
         ! Block dimensions
         WRITE (UNIT=unit_nr) (ec%block(i)%n_ao, i=1, ec%n_blocks)

         FLUSH (unit_nr)

      END IF

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_save_coupl_header

! **************************************************************************************************
!> \brief save POD block energies to binary data file
!> \param ec electronic coupling data structure
!> \param ib POD block ID
!> \param is spin component
!> \param unit_nr binary-file output unit
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_save_coupl_states(ec, ib, is, unit_nr)

      ! Routine arguments
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: ib, is, unit_nr

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_coupl_states', &
         routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: i, handle

      ! Time mark
      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN

         ! Save block energies
         WRITE (UNIT=unit_nr) (ec%block(ib)%mo(is)%eigenvalues(i), &
                               i=1, ec%block(ib)%mo(is)%nmo)

         FLUSH (unit_nr)

      END IF

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_save_coupl_states

! **************************************************************************************************
!> \brief save POD coupling elements from one off-diagonal block to binary data file
!> \param ec electronic coupling data structure
!> \param i1 ID of the first POD block
!> \param i2 ID of the second POD block
!> \param is spin component
!> \param unit_nr binary-file output unit
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_save_coupl_elements(ec, i1, i2, is, unit_nr)

      ! Routine arguments
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: i1, i2, is, unit_nr

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_coupl_elements', &
         routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: handle

      ! Time mark
      CALL timeset(routineN, handle)

      ! Save data to file
      CALL cp_fm_write_unformatted(ec%block(i1)%hab(is, i2), unit_nr)

      IF (unit_nr > 0) &
         FLUSH (unit_nr)

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_save_coupl_elements

! **************************************************************************************************
!> \brief close output stream for saving the couplings
!> \param qs_env QuickStep environment containing all system data
!> \param unit_nr binary-file output unit
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_save_coupl_close(qs_env, unit_nr)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(OUT)                               :: unit_nr

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_coupl_close', &
         routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: handle
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      ! Initialization
      NULLIFY (logger)
      NULLIFY (print_key)

      ! Time mark
      CALL timeset(routineN, handle)

      ! Output settings
      logger => cp_get_default_logger()
      print_key => section_vals_get_subs_vals(qs_env%input, &
                                              'PROPERTIES%ET_COUPLING%PROJECTION%PRINT')

      ! Close the file
      CALL cp_print_key_finished_output(unit_nr, logger, print_key, 'HAMILTONIAN')

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_save_coupl_close

! **************************************************************************************************
!> \brief save specified electronic states to cube files
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param n_spins number of spin states
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_save_state_cube(qs_env, ec, n_spins)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: n_spins

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_state_cube', &
         routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: i, j, k, l, n, handle
      INTEGER, DIMENSION(:), POINTER                     :: list
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mo_set_type), POINTER                         :: mo
      TYPE(section_vals_type), POINTER                   :: block_sec, mo_sec, print_sec

      ! Initialization
      NULLIFY (logger)
      NULLIFY (block_sec)
      NULLIFY (print_sec)
      NULLIFY (mo_sec)

      ! Time mark
      CALL timeset(routineN, handle)

      ! Output logger
      logger => cp_get_default_logger()
      block_sec => section_vals_get_subs_vals(qs_env%input, &
                                              'PROPERTIES%ET_COUPLING%PROJECTION%BLOCK')

      ! Print states of all blocks
      DO i = 1, ec%n_blocks

         IF (ASSOCIATED(ec%block(i)%mo)) THEN

            print_sec => section_vals_get_subs_vals(block_sec, 'PRINT', i_rep_section=i)
 
            ! Check if the print input section is active
            IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                                 print_sec, 'MO_CUBES'), cp_p_file)) THEN
 
               mo_sec => section_vals_get_subs_vals(print_sec, 'MO_CUBES')
 
               ! Spin states
               DO j = 1, n_spins
 
                  mo => ec%block(i)%mo(j)
 
                  CALL section_vals_val_get(mo_sec, keyword_name='MO_LIST', n_rep_val=n)
 
                  ! List of specific MOs
                  IF (n > 0) THEN
 
                     DO k = 1, n
                        NULLIFY (list)
                        CALL section_vals_val_get(mo_sec, keyword_name='MO_LIST', &
                                                  i_rep_val=k, i_vals=list)
                        IF (ASSOCIATED(list)) THEN
                           DO l = 1, SIZE(list)
                              CALL pod_mo_save(qs_env, logger, print_sec, mo, i, list(l), j)
                           END DO
                        END IF
                     END DO
 
                     ! Frontier MOs
                  ELSE
 
                     ! Occupied states
                     CALL section_vals_val_get(mo_sec, keyword_name='NHOMO', i_val=n)
 
                     IF (n > 0) THEN
                        DO k = MAX(1, mo%homo - n + 1), mo%homo
                           CALL pod_mo_save(qs_env, logger, print_sec, mo, i, k, j)
                        END DO
                     END IF
 
                     ! Unoccupied states
                     CALL section_vals_val_get(mo_sec, keyword_name='NLUMO', i_val=n)
 
                     IF (n > 0) THEN
                        DO k = mo%lfomo, MIN(mo%lfomo + n - 1, mo%nmo)
                           CALL pod_mo_save(qs_env, logger, print_sec, mo, i, k, j)
                        END DO
                     END IF
 
                  END IF
 
               END DO
 
            END IF

         END IF

      END DO

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_save_state_cube

! **************************************************************************************************
!> \brief save header of binary POD wavefunction file
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param n_spins number of spin components
!> \param unit_nr binary-file output unit
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_save_state_header(qs_env, ec, n_spins, unit_nr)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: n_spins
      INTEGER, INTENT(OUT)                               :: unit_nr

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_state_header', &
         routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: i, handle
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      ! Initialization
      NULLIFY (logger)
      NULLIFY (print_key)

      ! Time mark
      CALL timeset(routineN, handle)

      ! Output settings
      logger => cp_get_default_logger()
      print_key => section_vals_get_subs_vals(qs_env%input, &
                                              'PROPERTIES%ET_COUPLING%PROJECTION%PRINT')

      ! Open the file
      unit_nr = cp_print_key_unit_nr(logger, print_key, 'WAVEFUNCTION', &
                                     middle_name='RESTART', extension='.wfn', file_form='UNFORMATTED', &
                                     file_position='REWIND', file_action='WRITE', file_status='REPLACE')

      IF (unit_nr > 0) THEN

         ! Number of spin components
         WRITE (UNIT=unit_nr) n_spins
         ! Number of blocks
         WRITE (UNIT=unit_nr) ec%n_blocks
         ! Block dimensions
         WRITE (UNIT=unit_nr) (ec%block(i)%n_ao, i=1, ec%n_blocks)

         FLUSH (unit_nr)

      END IF

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_save_state_header

! **************************************************************************************************
!> \brief save POD orbital coefficients to binary data file
!> \param ec electronic coupling data structure
!> \param ib POD block ID
!> \param is spin component
!> \param unit_nr binary-file output unit
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_save_state_coeffs(ec, ib, is, unit_nr)

      ! Routine arguments
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: ib, is, unit_nr

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_state_coeffs', &
         routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: handle

      ! Time mark
      CALL timeset(routineN, handle)

      ! Save data to file
      CALL cp_fm_write_unformatted(ec%block(ib)%mo(is)%mo_coeff, unit_nr)

      IF (unit_nr > 0) &
         FLUSH (unit_nr)

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_save_state_coeffs

! **************************************************************************************************
!> \brief close output stream for saving the localized wavefunctions
!> \param qs_env QuickStep environment containing all system data
!> \param unit_nr binary-file output unit
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_save_state_close(qs_env, unit_nr)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(OUT)                               :: unit_nr

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_state_close', &
         routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: handle
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      ! Initialization
      NULLIFY (logger)
      NULLIFY (print_key)

      ! Time mark
      CALL timeset(routineN, handle)

      ! Output settings
      logger => cp_get_default_logger()
      print_key => section_vals_get_subs_vals(qs_env%input, &
                                              'PROPERTIES%ET_COUPLING%PROJECTION%PRINT')

      ! Close the file
      CALL cp_print_key_finished_output(unit_nr, logger, print_key, 'WAVEFUNCTION')

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_save_state_close

! **************************************************************************************************
!> \brief save visualization script for plotting the POD states
!> \param name_gpt name of the Gnuplot script file
!> \param name_dat name of the data file
!> \param ns number of spin components
!> \param nb number of POD blocks
!> \param e0 lower boundary of the energy range
!> \param e1 upper boundary of the energy range
!> \param rel indicates whether the states are relative or absolute
!> \author Z. Futera (07.2026)
! **************************************************************************************************
   SUBROUTINE pod_save_states_gnuplot(name_gpt,name_dat,ns,nb,e0,e1,rel)

      ! Routine arguments
      CHARACTER(LEN=default_path_length)                 :: name_gpt
      CHARACTER(LEN=default_path_length)                 :: name_dat
      INTEGER                                            :: ns
      INTEGER                                            :: nb
      REAL(KIND=dp)                                      :: e0
      REAL(KIND=dp)                                      :: e1
      LOGICAL                                            :: rel

      ! Local variables
      INTEGER                                            :: unit_nr, i, n
      REAL(KIND=dp)                                      :: mx_occ
      CHARACTER(LEN=default_path_length)                 :: base

      CALL open_file(file_name=name_gpt, file_action="WRITE", &
                     unit_number=unit_nr, file_status='REPLACE')

      ! File base name
      base = name_gpt
      n = SCAN(TRIM(name_gpt),'.',BACK=.true.)
      IF (n > 0) &
        base = name_gpt(1:n)

      ! Full occupancy
      mx_occ = 2.0
      IF (ns > 1) &
        mx_occ = 1.0

      WRITE(unit_nr, '(A)') '#!gnuplot'
      WRITE(unit_nr, '(A)') &
        '################################################################################'
      WRITE(unit_nr, '(A)') 'fl = "'//TRIM(name_dat)//'"'
      WRITE(unit_nr, '(A)') &
        '################################################################################'
      WRITE(unit_nr, '(A,F8.3)') 'e_min = ', cp_unit_from_cp2k(e0,"ev")
      WRITE(unit_nr, '(A,F8.3)') 'e_max = ', cp_unit_from_cp2k(e1,"ev")
      WRITE(unit_nr, '(A)') &
        '################################################################################'
      WRITE(unit_nr, '(A,F4.1)') 'fig_size_x = ',4.0
      WRITE(unit_nr, '(A,F4.1)') 'fig_size_y = ',4.0
      WRITE(unit_nr, '(A,I4)')   'plt_num_x  = ',nb
      WRITE(unit_nr, '(A,I4)')   'plt_num_y  = ',1
      WRITE(unit_nr, '(A,F4.1)') 'plt_off_x0 = ',10.0
      WRITE(unit_nr, '(A,F4.1)') 'plt_off_x1 = ',1.0
      WRITE(unit_nr, '(A,F4.1)') 'plt_off_y0 = ',1.0
      WRITE(unit_nr, '(A,F4.1)') 'plt_off_y1 = ',7.0
      WRITE(unit_nr, '(A,F4.1)') 'plt_sep_x  = ',1.0
      WRITE(unit_nr, '(A,F4.1)') 'plt_sep_y  = ',1.0
      WRITE(unit_nr, '(A)') &
        'plt_size_x = (100.0-plt_off_x0-plt_off_x1-((plt_num_x-1)*plt_sep_x))/plt_num_x'
      WRITE(unit_nr, '(A)') &
        'plt_size_y = (100.0-plt_off_y0-plt_off_y1-((plt_num_y-1)*plt_sep_y))/plt_num_y'
      WRITE(unit_nr, '(A)') &
        '################################################################################'
      WRITE(unit_nr, '(A)') 'cl_oc = "#0055ff"'
      WRITE(unit_nr, '(A)') 'cl_un = "#ff5500"'
      WRITE(unit_nr, '(A)') &
        '################################################################################'
      WRITE(unit_nr, '(A)') 'set terminal pdfcairo transparent truecolor enhanced \'
      WRITE(unit_nr, '(A)') 'background rgb "white" dashed font "Helvetica, 12" \'
      WRITE(unit_nr, '(A)') 'size fig_size_x,fig_size_y'
      WRITE(unit_nr, '(A)') 'set encoding iso_8859_1'
      WRITE(unit_nr, '(A)') 'set output "'//TRIM(base)//'pdf"'
      WRITE(unit_nr, '(A)') 'set multiplot'
      WRITE(unit_nr, '(A)') 'set size plt_size_x/100,plt_size_y/100'
      WRITE(unit_nr, '(A)') &
        '################################################################################'
      WRITE(unit_nr, '(A,F3.1)') 'x_min = ',0.0
      WRITE(unit_nr, '(A,F3.1)') 'x_max = ',1.0
      WRITE(unit_nr, '(A)') 'y_off = 0.05*(e_max - e_min)'
      WRITE(unit_nr, '(A)') 'set xrange [x_min:x_max]'
      WRITE(unit_nr, '(A)') 'set yrange [e_min-y_off:e_max+y_off]'
      WRITE(unit_nr, '(A)') &
        '################################################################################'
      WRITE(unit_nr, '(A)') 'set tmargin 0'
      WRITE(unit_nr, '(A)') 'set bmargin 0'
      WRITE(unit_nr, '(A)') 'set lmargin 0'
      WRITE(unit_nr, '(A)') 'set rmargin 0'
      WRITE(unit_nr, '(A)') 'set nokey'
      WRITE(unit_nr, '(A)') 'set noxtics'
      WRITE(unit_nr, '(A)') 'set format x ""'
      WRITE(unit_nr, '(A,F3.1,A)') 'set cbrange [0.0:',mx_occ,']'
      WRITE(unit_nr, '(A,F3.1,A)') 'set palette defined ( 0.0 cl_un, ',mx_occ,' cl_oc )'
      WRITE(unit_nr, '(A)') 'unset colorbox'
      WRITE(unit_nr, '(A)') &
        '################################################################################'

      DO i = 1, nb
        WRITE(unit_nr, '(A,I3,A)') &
          'set origin (plt_off_x0+',i-1,'*(plt_size_x+plt_sep_x))/100,\'
        WRITE(unit_nr, '(A)') &
          '  (plt_off_y0+0*(plt_size_y+plt_sep_y))/100'
        WRITE(unit_nr, '(A,I3,A)') &
          'set title "Block ',i,'" offset 0.0,-0.5'
        IF (i == 1) THEN
          IF (rel) THEN
            WRITE(unit_nr, '(A)') 'set ylabel "E - E_F [eV]"'
          ELSE
            WRITE(unit_nr, '(A)') 'set ylabel "E [eV]"'
          END IF
          WRITE(unit_nr, '(A)') 'set format y "%3.0f"'
        ELSE
          WRITE(unit_nr, '(A)') 'set ylabel ""'
          WRITE(unit_nr, '(A)') 'set format y ""'
        END IF
        WRITE(unit_nr, '(A)') 'plot \'
        WRITE(unit_nr, '(A,I3,A)') &
          '  fl index ',i-1,' u (0.1):3:(0.8):(0):2 w vectors nohead \'
        WRITE(unit_nr, '(A)') &
          '     lt 1 lw 1.25 lc palette'
        WRITE(unit_nr, '(A)') &
          '################################################################################'
      END DO

      WRITE(unit_nr, '(A)') 'set nomultiplot'
      WRITE(unit_nr, '(A)') 'quit'

      CALL close_file(unit_nr)

   END SUBROUTINE pod_save_states_gnuplot

! **************************************************************************************************
!> \brief save the POD electronic states to external file and create a visualization script
!> \param output_unit unit number of the open output stream
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \author Z. Futera (07.2026)
! **************************************************************************************************
   SUBROUTINE pod_save_states(output_unit, qs_env, ec)

      ! Routine arguments
      INTEGER                                            :: output_unit
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec

      ! Local variables
      LOGICAL                                            :: do_save, relative
      LOGICAL                                            :: has_range_min, has_range_max
      INTEGER                                            :: i, j, n, na, nb, unit_nr
      CHARACTER(LEN=default_path_length)                 :: file_dat, file_gpt
      REAL(KIND=dp)                                      :: range_min, range_max
      REAL(KIND=dp)                                      :: e, ea, eb, e0, e1, e_min, e_max
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_cntrl
      TYPE(section_vals_type), POINTER                   :: inp_et_print, inp_et_states
      TYPE(section_vals_type), POINTER                   :: inp_et_gnuplot

      ! Check whether to save states to file or not
      logger => cp_get_default_logger()
      inp_et_print => section_vals_get_subs_vals(qs_env%input, &
                                                'PROPERTIES%ET_COUPLING%PROJECTION%PRINT')
      do_save = BTEST(cp_print_key_should_output(logger%iter_info, &
                      inp_et_print, 'STATES'), cp_p_file)

      ! Save the states
      IF (do_save) THEN

        ! Options

        inp_et_states => section_vals_get_subs_vals(inp_et_print, 'STATES')
        file_dat = cp_print_key_generate_filename(logger, inp_et_states, middle_name='STATES', &
                                                  extension='.dat', my_local=.FALSE.)

        inp_et_gnuplot => section_vals_get_subs_vals(inp_et_print, 'STATES%GNUPLOT')
        file_gpt = cp_print_key_generate_filename(logger, inp_et_gnuplot, middle_name="STATES", &
                                                  extension='.gpt', my_local=.FALSE.)

        CALL section_vals_val_get(inp_et_states,'RELATIVE',l_val=relative)
        CALL section_vals_val_get(inp_et_states,'RANGE_MIN',r_val=range_min,explicit=has_range_min)
        CALL section_vals_val_get(inp_et_states,'RANGE_MAX',r_val=range_max,explicit=has_range_max)

        ! Data file

        IF (file_dat /= '') THEN
           WRITE (output_unit, '(/,T3,A,A)') &
                 'Storing electronic states to ', TRIM(ADJUSTL(file_dat))

           IF (relative) THEN
              WRITE (output_unit, '(T3,A,F8.3,A)') &
                    ' -- Relative energies wrt to Fermi level (', &
                    cp_unit_from_cp2k(ec%fermi(1),"ev"),' eV)'
           ELSE
              WRITE (output_unit, '(T3,A)') &
                    ' -- Absolute energies'
           END IF

           e0 = 0.0d0
           e1 = 0.0d0

           IF (has_range_min .OR. has_range_max) THEN
             IF (has_range_min) THEN
                WRITE (output_unit, '(T3,A,F8.3,A)') &
                      ' -- Energy range lower boundary: ', &
                      cp_unit_from_cp2k(range_min,"ev"),' eV'
                e0 = range_min
             END IF
             IF (has_range_max) THEN
                WRITE (output_unit, '(T3,A,F8.3,A)') &
                      ' -- Energy range upper boundary: ', &
                      cp_unit_from_cp2k(range_max,"ev"),' eV'
                e1 = range_max
             END IF
           ELSE
              WRITE (output_unit, '(T3,A)') &
                    ' -- All states are saved'
           END IF

           ! Save states

           CALL get_qs_env(qs_env, dft_control=dft_cntrl)

           CALL open_file(file_name=file_dat, file_action="WRITE", &
                          unit_number=unit_nr, file_status='REPLACE')

           IF (.NOT.has_range_min) THEN
             IF (ec%n_blocks > 0) THEN
               e0 = ec%block(1)%mo(1)%eigenvalues(1)
               e1 = ec%block(1)%mo(1)%eigenvalues(1)
             END IF
           END IF

           DO i = 1, ec%n_blocks
             WRITE(unit_nr, '(A,I5)') &
               '# POD block ', i

             ! Spin polarized
             IF (dft_cntrl%nspins > 1) THEN

                na = ec%block(i)%mo(1)%nmo
                nb = ec%block(i)%mo(1)%nmo
                n = MAX(na, nb)

                ea = 0.0d0
                eb = 0.0d0
                IF (relative) THEN
                  ea = ec%fermi(1);
                  eb = ec%fermi(2);
                END IF
                IF (.NOT.has_range_min .AND. ec%block(i)%mo(1)%nmo > 0) THEN
                   e0 = MIN(e0,ec%block(i)%mo(1)%eigenvalues(1)-ea)
                   e0 = MIN(e0,ec%block(i)%mo(2)%eigenvalues(1)-eb)
                END IF
                IF (.NOT.has_range_max .AND. ec%block(i)%mo(1)%nmo > 0) THEN
                   e1 = MAX(e1,ec%block(i)%mo(1)%eigenvalues(ec%block(i)%mo(1)%nmo)-ea)
                   e1 = MAX(e1,ec%block(i)%mo(2)%eigenvalues(ec%block(i)%mo(1)%nmo)-eb)
                END IF

                DO j = 1, n
                   ea = ec%block(i)%mo(1)%eigenvalues(j)
                   eb = ec%block(i)%mo(2)%eigenvalues(j)
                   IF (relative) THEN
                     ea = ea - ec%fermi(1)
                     eb = eb - ec%fermi(2)
                   END IF
                   IF (has_range_min .AND. ea < range_min .AND. eb < range_min) &
                     CYCLE
                   IF (has_range_max .AND. ea > range_max .AND. eb > range_max) &
                     EXIT
                   IF (j > na) THEN
                      WRITE (unit_nr, '(I10,A20,F8.3,F12.6)') j, ' ', &
                        ec%block(i)%mo(2)%occupation_numbers(j), cp_unit_from_cp2k(eb,"ev")
                   ELSEIF (j > nb) THEN
                      WRITE (unit_nr, '(I10,F8.3,F12.6,A20)') j, &
                        ec%block(i)%mo(1)%occupation_numbers(j), cp_unit_from_cp2k(ea,"ev"), ' '
                   ELSE 
                      WRITE (unit_nr, '(I10,F8.3,F12.6,F8.3,F12.6)') j, &
                        ec%block(i)%mo(1)%occupation_numbers(j), cp_unit_from_cp2k(ea,"ev"), &
                        ec%block(i)%mo(2)%occupation_numbers(j), cp_unit_from_cp2k(eb,"ev")
                   END IF
                END DO

             ! Spin non-polarized
             ELSE

                e = 0.0d0
                IF (relative) &
                  e = ec%fermi(1);
                IF (.NOT.has_range_min .AND. ec%block(i)%mo(1)%nmo > 0) &
                   e0 = MIN(e0,ec%block(i)%mo(1)%eigenvalues(1)-e)
                IF (.NOT.has_range_max .AND. ec%block(i)%mo(1)%nmo > 0) &
                   e1 = MAX(e1,ec%block(i)%mo(1)%eigenvalues(ec%block(i)%mo(1)%nmo)-e)

                DO j = 1, ec%block(i)%mo(1)%nmo
                   e = ec%block(i)%mo(1)%eigenvalues(j)
                   IF (relative) &
                     e = e - ec%fermi(1)
                   IF (has_range_min .AND. e < range_min) &
                     CYCLE
                   IF (has_range_max .AND. e > range_max) &
                     EXIT
                   WRITE (unit_nr, '(I10,F8.3,F12.6)') &
                     j, ec%block(i)%mo(1)%occupation_numbers(j), &
                     cp_unit_from_cp2k(e,"ev")
                END DO

             END IF

             WRITE(unit_nr, '(/)') 
           END DO

           CALL close_file(unit_nr)

        END IF

        ! Gnuplot script

        IF (file_gpt /= '') THEN
           WRITE (output_unit, '(T3,A,A)') &
                 'Creating Gnuplot visualization script ', TRIM(ADJUSTL(file_gpt))

           CALL  pod_save_states_gnuplot(file_gpt,file_dat, &
                   dft_cntrl%nspins,ec%n_blocks,e0,e1,relative)

        END IF

      END IF

   END SUBROUTINE pod_save_states

! **************************************************************************************************
!> \brief save transformation matrix S^{+/-1/2} to binary data file
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param key 'F' for forward, 'B' for backward
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   SUBROUTINE pod_save_transf_mat(qs_env, ec, key)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      CHARACTER                                          :: key

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_save_transf_mat', &
         routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: unit_nr, handle
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      ! Initialization
      NULLIFY (logger)
      NULLIFY (print_key)

      ! Time mark
      CALL timeset(routineN, handle)

      ! Output settings
      logger => cp_get_default_logger()
      print_key => section_vals_get_subs_vals(qs_env%input, &
                                              'PROPERTIES%ET_COUPLING%PROJECTION%PRINT')

      ! Forward transformation matrix
      IF (key == 'F') THEN

         unit_nr = cp_print_key_unit_nr(logger, print_key, 'TRANSF_MATRIX_F', &
                                        middle_name='RESTART', extension='.tmf', file_form='UNFORMATTED', &
                                        file_position='REWIND', file_action='WRITE', file_status='REPLACE')

         CALL cp_fm_write_unformatted(ec%tmatf, unit_nr)

         CALL cp_print_key_finished_output(unit_nr, logger, print_key, 'TRANSF_MATRIX_F')

         ! Backward transformation matrix
      ELSE IF (key == 'B') THEN

         unit_nr = cp_print_key_unit_nr(logger, print_key, 'TRANSF_MATRIX_R', &
                                        middle_name='RESTART', extension='.tmr', file_form='UNFORMATTED', &
                                        file_position='REWIND', file_action='WRITE', file_status='REPLACE')

         CALL cp_fm_write_unformatted(ec%tmatr, unit_nr)

         CALL cp_print_key_finished_output(unit_nr, logger, print_key, 'TRANSF_MATRIX_R')

      END IF

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_save_transf_mat

END MODULE et_coupling_pod_save
