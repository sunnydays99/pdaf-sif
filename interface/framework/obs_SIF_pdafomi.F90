#ifdef CLMFIVE
MODULE obs_SIF_pdafomi

    USE mod_parallel_pdaf, ONLY: mype_filter
    USE PDAFomi, ONLY: obs_f, obs_l

    IMPLICIT NONE
    SAVE
    PUBLIC

    ! Variables which are inputs to the module (usually set in init_pdaf)
    LOGICAL :: assim_SIF      !< Whether to assimilate this data type
    REAL    :: rms_obs_SIF      !< Observation error standard deviation (for constant errors)

    ! longitude and latitude of grid cells and observation cells (same pattern as SM)
    INTEGER, ALLOCATABLE :: longxy(:), latixy(:), longxy_obs(:), latixy_obs(:)

    ! One can declare further variables, e.g. for file names which can
    ! be use-included in init_pdaf() and initialized there.

  ! Declare instances of observation data types used here
  ! We use generic names here, but one could rename the variables
    TYPE(obs_f), TARGET, PUBLIC :: thisobs      ! full observation
    TYPE(obs_l), TARGET, PUBLIC :: thisobs_l    ! local observation

  !$OMP THREADPRIVATE(thisobs_l)

CONTAINS

    SUBROUTINE init_dim_obs_SIF(step, dim_obs)

      USE mpi
      USE mod_parallel_pdaf, ONLY: mype_filter, comm_filter, npes_filter, &
                                  abort_parallel, mype_world
      USE mod_assimilation,  ONLY: obs_index_p, obs_filename, &
                                  obs_pdaf2nc, obs_nc2pdaf, &
                                  local_dims_obs, local_disp_obs, &
                                  longxy_obs_floor, latixy_obs_floor, &
                                  screen, cradius_SIF
      USE PDAFomi,           ONLY: PDAFomi_gather_obs, pi
      USE mod_read_obs,      ONLY: multierr, read_obs_nc_type
      USE enkf_clm_mod,      ONLY: clmstatevec_allcol, clmstatevec_only_active, &
                                  clmstatevec_max_layer, state_clm2pdaf_p, &
                                  domain_def_clm
      USE mod_parallel_pdaf, ONLY: mype_world
      USE shr_kind_mod,      ONLY: r8 => shr_kind_r8
      USE GridcellType,      ONLY: grc
      USE clm_varcon,        ONLY: spval, ispval
      USE decompMod,         ONLY: get_proc_bounds, get_proc_global
      USE ColumnType,        ONLY: col
      USE PatchType,         ONLY: patch          ! *** KEY: patch level for FSIF ***
      USE mod_tsmp,          ONLY: obs_interp_switch, point_obs, da_print_obs_index
      USE enkf_clm_mod,      ONLY: get_interp_idx

      IMPLICIT NONE

      INTEGER, INTENT(in)    :: step
      INTEGER, INTENT(inout) :: dim_obs

      ! Local variables — same pattern as obs_SM
      INTEGER :: i, c, g, p            ! counters
      INTEGER :: cnt, dim_obs_p
      REAL, ALLOCATABLE :: obs_p(:), obs_g(:), ivar_obs_p(:)
      REAL, ALLOCATABLE :: ocoord_p(:,:)
      REAL, ALLOCATABLE :: lon_obs(:), lat_obs(:), dr_obs(:), obserr(:)
      INTEGER, ALLOCATABLE :: layer_obs(:)
      REAL, ALLOCATABLE :: obscov(:,:)
      CHARACTER(len=110) :: current_observation_filename
      CHARACTER(len=20)  :: obs_type_name

      ! Patch-level variables for SIF
      ! NOTE: FSIF is ptr_patch, so we need to average over patches per gridcell
      REAL, ALLOCATABLE :: fsif_gridcell(:)   ! gridcell-average FSIF from history
      REAL :: fsif_patch_sum, coszen_g, fsno_g, elai_g
      INTEGER :: n_patch_active

      INTEGER :: begp, endp, begc, endc, begl, endl, begg, endg
      INTEGER :: numg, numl, numc, nump
      REAL(r8), POINTER :: lon(:), lat(:)
      INTEGER, POINTER :: mycgridcell(:)
      REAL :: deltax, deltay
      LOGICAL :: obs_snapped, newgridcell
      INTEGER :: ierror, cnt_p, sum_dim_obs_p
      INTEGER :: pe

      IF (mype_filter==0) &
          WRITE(*,*) 'Assimilate observations - obs type SIF (FSIF at 740nm)'

      IF (assim_SIF) thisobs%doassim = 1
      thisobs%disttype = 3   ! geographic haversine
      thisobs%ncoord   = 2

      obs_type_name = 'SIF'

      WRITE(current_observation_filename, '(a, i5.5)') TRIM(obs_filename)//'.', step

      IF (mype_filter == 0) THEN
          CALL read_obs_nc_type(current_observation_filename, obs_type_name, &
                              dim_obs, obs_g, lon_obs, lat_obs, layer_obs, &
                              dr_obs, obserr, obscov)
      END IF

      CALL mpi_bcast(dim_obs, 1, MPI_INTEGER, 0, comm_filter, ierror)

      ! Handle zero-obs case (same boilerplate as SM)
      IF (dim_obs == 0) THEN
          dim_obs_p = 0
          ALLOCATE(obs_p(1), ivar_obs_p(1), ocoord_p(2,1), thisobs%id_obs_p(1,1))
          thisobs%infile = 0
          CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
                                thisobs%ncoord, cradius_SIF, dim_obs)
          IF (mype_filter==0) DEALLOCATE(obs_g)
          DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)
          RETURN
      END IF

      ! Broadcast obs arrays (same pattern as SM)
      CALL mpi_bcast(multierr,   1,       MPI_INTEGER,          0, comm_filter, ierror)
      IF (mype_filter /= 0) THEN
          IF (ALLOCATED(obs_g))    DEALLOCATE(obs_g);    ALLOCATE(obs_g(dim_obs))
          IF (ALLOCATED(lon_obs))  DEALLOCATE(lon_obs);  ALLOCATE(lon_obs(dim_obs))
          IF (ALLOCATED(lat_obs))  DEALLOCATE(lat_obs);  ALLOCATE(lat_obs(dim_obs))
          IF (ALLOCATED(dr_obs))   DEALLOCATE(dr_obs);   ALLOCATE(dr_obs(2))
          IF (ALLOCATED(layer_obs))DEALLOCATE(layer_obs);ALLOCATE(layer_obs(dim_obs))
          IF (multierr==1) THEN
              IF (ALLOCATED(obserr)) DEALLOCATE(obserr); ALLOCATE(obserr(dim_obs))
          END IF
      END IF
      CALL mpi_bcast(obs_g,    dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
      IF (multierr==1) &
          CALL mpi_bcast(obserr, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
      CALL mpi_bcast(lon_obs,  dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
      CALL mpi_bcast(lat_obs,  dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
      CALL mpi_bcast(dr_obs,   2,       MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
      CALL mpi_bcast(layer_obs,dim_obs, MPI_INTEGER,          0, comm_filter, ierror)

      thisobs%infile = 1
      CALL domain_def_clm(lon_obs, lat_obs, dim_obs, longxy, latixy, longxy_obs, latixy_obs)

      lon => grc%londeg
      lat => grc%latdeg
      mycgridcell => col%gridcell

      CALL get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
      CALL get_proc_global(numg, numl, numc, nump)

      ! -----------------------------------------------------------------------
      ! SIF-SPECIFIC: Read gridcell-average FSIF from eCLM history
      !
      ! FSIF is a patch-level variable. We compute a weighted average over
      ! all vegetated patches in each gridcell. This matches what a satellite
      ! observes: the area-weighted canopy SIF from the gridcell footprint.
      !
      ! We also read COSZEN and FSNO for QC masking.
      ! -----------------------------------------------------------------------
      CALL read_fsif_from_history(step, begg, endg, begp, endp, &
                                 fsif_gridcell, coszen_gridcell, fsno_gridcell)
      ! (Implementation of read_fsif_from_history: see Section 1.2)

      ! -----------------------------------------------------------------------
      ! Snapping loop: obs → gridcells (same structure as SM but at gridcell level)
      ! SIF obs map to gridcells directly (no column/layer depth dimension)
      ! -----------------------------------------------------------------------
      IF (ALLOCATED(thisobs%id_obs_p)) DEALLOCATE(thisobs%id_obs_p)
      ALLOCATE(thisobs%id_obs_p(1, endg-begg+1))
      thisobs%id_obs_p(1,:) = 0

      dim_obs_p = 0
      cnt = 1

      DO i = 1, dim_obs
          obs_snapped = .FALSE.
          DO g = begg, endg

              deltax = ABS(lon(g) - lon_obs(i))
              IF (deltax > 180.0) deltax = 360.0 - deltax
              deltay = ABS(lat(g) - lat_obs(i))

              IF ((deltax <= dr_obs(1)) .AND. (deltay <= dr_obs(2))) THEN

                ! -------------------------------------------------------
                ! SIF QC MASKS — applied at gridcell level before counting
                ! -------------------------------------------------------
                ! 1. Nighttime: COSZEN < 0.01 (no SIF without sunlight)
                IF (coszen_gridcell(g-begg+1) < 0.01) CYCLE
                ! 2. Snow-covered: FSNO > 0.5 (snow confounds SIF signal)
                IF (fsno_gridcell(g-begg+1) > 0.5) CYCLE
                ! 3. No valid model SIF (spval check, covers lake/urban)
                IF (fsif_gridcell(g-begg+1) <= 0.0 .OR. &
                    fsif_gridcell(g-begg+1) == spval) CYCLE
                ! 4. Negative satellite obs (physically impossible)
                IF (obs_g(i) <= 0.0) CYCLE

                dim_obs_p = dim_obs_p + 1
                ! id_obs_p points to gridcell offset in state (set=1, not used
                ! directly in obs_op_SIF which reads from history array)
                thisobs%id_obs_p(1, g-begg+1) = 1

                IF (obs_snapped) THEN
                    WRITE(*,*) "ERROR: SIF obs snapped at multiple gridcells, i=", i
                    CALL abort_parallel()
                END IF
                obs_snapped = .TRUE.
                cnt = cnt + 1
              END IF
          END DO
      END DO

      ! --- Allocate PE-local obs arrays ---
      IF (ALLOCATED(ivar_obs_p)) DEALLOCATE(ivar_obs_p)
      ALLOCATE(ivar_obs_p(dim_obs_p))
      IF (ALLOCATED(ocoord_p)) DEALLOCATE(ocoord_p)
      ALLOCATE(ocoord_p(2, dim_obs_p))

      ! (Same MPI allreduce pattern as SM)
      ! Dimension of full observation vector
      ! ------------------------------------

      ! add and broadcast size of PE-local observation dimensions using mpi_allreduce
      call mpi_allreduce(dim_obs_p, sum_dim_obs_p, 1, MPI_INTEGER, MPI_SUM, &
        comm_filter, ierror)

      ! Check sum of dimensions of PE-local observation vectors against
      ! dimension of full observation vector
      if (.not. sum_dim_obs_p == dim_obs) then
        print *, "TSMP-PDAF mype(w)=", mype_world, ": ERROR Sum of PE-local observation dimensions"
        print *, "sum_dim_obs_p=", sum_dim_obs_p
        print *, "dim_obs=", dim_obs
        call abort_parallel()
      end if

      !  Gather PE-local observation dimensions and displacements in arrays
      ! ----------------------------------------------------------------

       ! Allocate array of PE-local observation dimensions
      IF (ALLOCATED(local_dims_obs)) DEALLOCATE(local_dims_obs)
      ALLOCATE(local_dims_obs(npes_filter))

      ! Gather array of PE-local observation dimensions
      call mpi_allgather(dim_obs_p, 1, MPI_INTEGER, local_dims_obs, 1, MPI_INTEGER, &
        comm_filter, ierror)

      ! Allocate observation displacement array local_disp_obs
      IF (ALLOCATED(local_disp_obs)) DEALLOCATE(local_disp_obs)
      ALLOCATE(local_disp_obs(npes_filter))

      ! Set observation displacement array local_disp_obs
      local_disp_obs(1) = 0
      do i = 2, npes_filter
        local_disp_obs(i) = local_disp_obs(i-1) + local_dims_obs(i-1)
      end do

      if (mype_filter==0 .and. screen > 2) then
          print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: local_disp_obs=", local_disp_obs
      end if


      ! Write index mapping array NetCDF->PDAF
      ! --------------------------------------
      ! Set index mapping `obs_pdaf2nc` between observation order in
      ! NetCDF input and observation order in pdaf as determined by domain
      ! decomposition.

      ! Use-case: Correct index order in loops over NetCDF-observation
      ! file input arrays.

      ! Trivial example: The order in the NetCDF file corresponds exactly
      ! to the order in the domain decomposition in PDAF, e.g. for a
      ! single PE per component model run.

      ! Non-trivial example: The first observation in the NetCDF file is
      ! not located in the domain/subgrid of the first PE. Rather, the
      ! second observation in the NetCDF file (`i=2`) is the only
      ! observation (`cnt = 1`) in the subgrid of the first PE
      ! (`mype_filter = 0`). This leads to a non-trivial index mapping,
      ! e.g. `obs_pdaf2nc(1)==2`:
      !
      ! i = 2
      ! cnt = 1
      ! mype_filter = 0
      !
      ! obs_pdaf2nc(local_disp_obs(mype_filter+1)+cnt) = i
      !-> obs_pdaf2nc(local_disp_obs(1)+1) = 2
      !-> obs_pdaf2nc(1) = 2


      if (allocated(obs_pdaf2nc)) deallocate(obs_pdaf2nc)
      allocate(obs_pdaf2nc(dim_obs))
      obs_pdaf2nc = 0
      if (allocated(obs_nc2pdaf)) deallocate(obs_nc2pdaf)
      allocate(obs_nc2pdaf(dim_obs))
      obs_nc2pdaf = 0


      if (point_obs==1) then

        cnt = 1
        do i = 1, dim_obs
          ! Many processes may not contain the observation / do not need
          ! to snap it, so default true
          obs_snapped = .true.

          do g = begg,endg
            newgridcell = .true.
            do c = begc,endc
              cg =   mycgridcell(c)
              if(cg == g) then
                if(newgridcell) then

                  if(is_use_dr) then
                    deltax = abs(lon(g)-lon_obs(i))
                    if (deltax > 180.0) then
                      deltax = 360.0 - deltax
                    end if
                    deltay = abs(lat(g)-lat_obs(i))
                  end if

                  if(((is_use_dr).and.(deltax<=dr_obs(1)).and.(deltay<=dr_obs(2))).or. &
                     ((.not. is_use_dr).and.(longxy_obs(i) == longxy(g-begg+1)) .and. &
                      (latixy_obs(i) == latixy(g-begg+1)))) then
#ifdef CLMFIVE
                    if(state_clm2pdaf_p(c,1)==ispval) then
                      ! `ispval`: column not in state vector, most likely
                      ! because it is hydrologically inactive

                      ! Observation not snapped, even though location is
                      ! right!
                      obs_snapped = .false.

                      ! Do not use this column for snapping an
                      ! observation, instead cycle to next column
                      cycle
                    end if
#endif
                    obs_pdaf2nc(local_disp_obs(mype_filter+1)+cnt) = i
                    obs_nc2pdaf(i) = local_disp_obs(mype_filter+1)+cnt
                    cnt = cnt + 1

                    ! Observation snapped at location (possibly
                    ! overwriting a false from inactive column before)
                    obs_snapped = .true.
                  end if

                  newgridcell = .false.

                end if
              end if
            end do
          end do

          ! Warning, when an observation has not been snapped to any
          ! active gridcell.
          if(.not. obs_snapped) then
            print *, "TSMP-PDAF mype(w)=", mype_world, ": ERROR observations exist at non-active gridcells."
            print *, "Consider removing the following observation from the observation files."
            print *, "Observation-index in NetCDF-file: i=", i
            call abort_parallel()
          end if
        end do

      end if

      ! collect values from all PEs, by adding all PE-local arrays (works
      ! since only the subsection belonging to a specific PE is non-zero)
      call mpi_allreduce(MPI_IN_PLACE,obs_pdaf2nc,dim_obs,MPI_INTEGER,MPI_SUM,comm_filter,ierror)
      call mpi_allreduce(MPI_IN_PLACE,obs_nc2pdaf,dim_obs,MPI_INTEGER,MPI_SUM,comm_filter,ierror)

      if (mype_filter==0 .and. screen > 2) then
          print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: obs_pdaf2nc=", obs_pdaf2nc
      end if

      ! Write process-local observation arrays

      ! --- Fill obs_p, ivar_obs_p, ocoord_p ---
      IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)
      ALLOCATE(obs_p(dim_obs_p))
      IF (ALLOCATED(obs_index_p)) DEALLOCATE(obs_index_p)
      ALLOCATE(obs_index_p(dim_obs_p))

      cnt = 1
      DO i = 1, dim_obs
          DO g = begg, endg
              deltax = ABS(lon(g) - lon_obs(i))
              IF (deltax > 180.0) deltax = 360.0 - deltax
              deltay = ABS(lat(g) - lat_obs(i))
              IF ((deltax <= dr_obs(1)) .AND. (deltay <= dr_obs(2))) THEN
                  ! Apply same QC masks as above
                  IF (coszen_gridcell(g-begg+1) < 0.01)    CYCLE
                  IF (fsno_gridcell(g-begg+1) > 0.5)       CYCLE
                  IF (fsif_gridcell(g-begg+1) <= 0.0 .OR. &
                      fsif_gridcell(g-begg+1) == spval)     CYCLE
                  IF (obs_g(i) <= 0.0)                     CYCLE

                  ! Haversine: coords in radians
                  ocoord_p(1, cnt) = lon_obs(i) * pi / 180.0
                  ocoord_p(2, cnt) = lat_obs(i) * pi / 180.0

                  ! obs_index_p for SIF stores the gridcell offset (g-begg+1)
                  ! This is used in obs_op_SIF to look up fsif_gridcell
                  obs_index_p(cnt) = g - begg + 1

                  obs_p(cnt) = obs_g(i)
                  IF (multierr==1) THEN
                      ivar_obs_p(cnt) = 1.0 / (obserr(i) * obserr(i))
                  ELSE
                      ivar_obs_p(cnt) = 1.0 / (rms_obs_SIF * rms_obs_SIF)
                  END IF
                  cnt = cnt + 1
              END IF
          END DO
      END DO

      CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
                            thisobs%ncoord, cradius_SIF, dim_obs)

      DEALLOCATE(obs_g, obs_p, ocoord_p, ivar_obs_p)
 
    END SUBROUTINE init_dim_obs_SIF    


    SUBROUTINE obs_op_SIF(dim_p, dim_obs, state_p, ostate)

      USE mod_assimilation,  ONLY: obs_index_p
      USE PDAFomi_obs_f,     ONLY: PDAFomi_gather_obsstate
      USE enkf_clm_mod,      ONLY: fsif_current   ! gridcell array — see Section 1.2

      IMPLICIT NONE
      INTEGER, INTENT(in)    :: dim_p, dim_obs
      REAL,    INTENT(in)    :: state_p(dim_p)
      REAL,    INTENT(inout) :: ostate(dim_obs)

      REAL, ALLOCATABLE :: ostate_p(:)
      INTEGER :: i

      IF (thisobs%dim_obs_p > 0) THEN
          ALLOCATE(ostate_p(thisobs%dim_obs_p))
      ELSE
          ALLOCATE(ostate_p(1))
      END IF

    ! The obs operator for SIF is a DIRECT READ from the eCLM history array.
    ! obs_index_p(i) stores the gridcell offset (g-begg+1) set in init_dim_obs_SIF.
    ! fsif_current(:) is the PE-local array of gridcell-average FSIF values
    ! read from the eCLM history file in read_fsif_from_history().
    !
    ! This is analogous to obs_op_SM using state_p(obs_index_p(i)),
    ! except here we read from the history diagnostic array, not the
    ! prognostic state vector — because FSIF is not in the state vector.

      DO i = 1, thisobs%dim_obs_p
          ostate_p(i) = fsif_current(obs_index_p(i))
      END DO

      CALL PDAFomi_gather_obsstate(thisobs, ostate_p, ostate)
      DEALLOCATE(ostate_p)

    END SUBROUTINE obs_op_SIF

    !-------------------------------------------------------------------------------
    !> Initialize local information on the module-type observation
    !!
    !! The routine is called during the loop over all local
    !! analysis domains. It has to initialize the information
    !! about local observations of the module type. It returns
    !! number of local observations of the module type for the
    !! current local analysis domain in DIM_OBS_L and the full
    !! and local offsets of the observation in the overall
    !! observation vector.
    !!
    !! This routine calls the routine PDAFomi_init_dim_obs_l
    !! for each observation type. The call allows to specify a
    !! different localization radius and localization functions
    !! for each observation type and  local analysis domain.
    !!
    SUBROUTINE init_dim_obs_l_SIF(domain_p, step, dim_obs, dim_obs_l)

      ! Include PDAFomi function
      USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l, pi

      ! Include localization radius and local coordinates
      USE mod_assimilation, &
           ONLY: cradius_SIF, locweight, sradius_SIF, screen

      USE enkf_clm_mod, ONLY: state_loc2clm_c_p, clmstatevec_allcol, clmstatevec_only_active

      use shr_kind_mod, only: r8 => shr_kind_r8

      use decompMod , only : get_proc_bounds


#ifdef CLMFIVE
      USE GridcellType, ONLY: grc
      USE ColumnType, ONLY : col
#else
      USE clmtype, ONLY : clm3
#endif
      use clm_varcon, only: spval

      use mod_parallel_pdaf, &
           ONLY: mype_world



      IMPLICIT NONE

  ! *** Arguments ***
      INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
      INTEGER, INTENT(in)  :: step         !< Current time step
      INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
      INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector

      REAL :: coords_l(2)      ! Coordinates of local analysis domain

      real(r8), pointer :: lon(:)
      real(r8), pointer :: lat(:)
      integer, pointer :: mycgridcell(:) !Pointer for CLM3.5/CLM5.0 col->gridcell index arrays

#ifdef CLMFIVE
      ! Obtain CLM lon/lat information
      lon   => grc%londeg
      lat   => grc%latdeg
      ! Obtain CLM column-gridcell information
      mycgridcell => col%gridcell
#else
      lon   => clm3%g%londeg
      lat   => clm3%g%latdeg
      mycgridcell => clm3%g%l%c%gridcell
#endif


  ! **********************************************
  ! *** Initialize local observation dimension ***
  ! **********************************************
    ! count observations within a radius

    if (thisobs%infile==1) then


      if (clmstatevec_allcol==0 .and. clmstatevec_only_active==0) then
        if (lon(state_loc2clm_c_p(domain_p))>180) then
          coords_l(1) = lon(state_loc2clm_c_p(domain_p)) - 360.0
        else
          coords_l(1) = lon(state_loc2clm_c_p(domain_p))
        end if
        coords_l(2) = lat(state_loc2clm_c_p(domain_p))

      else

        ! get coords_l --> coordinates of local analysis domain
        if (lon(mycgridcell(state_loc2clm_c_p(domain_p)))>180) then
        ! if SM should be assimilated. Else, state_loc2clm_c_p is not allocated
          coords_l(1) = lon(mycgridcell(state_loc2clm_c_p(domain_p))) - 360.0
        else
          coords_l(1) = lon(mycgridcell(state_loc2clm_c_p(domain_p)))
        end if
        coords_l(2) = lat(mycgridcell(state_loc2clm_c_p(domain_p)))

      end if

      if (thisobs%disttype==3) then ! if haversine formula in distance calculation, the coordinates have to be converted to radians
        coords_l(1) = coords_l(1) * pi / 180.0
        coords_l(2) = coords_l(2) * pi / 180.0
      end if

    else

      coords_l(1) = spval
      coords_l(2) = spval

    end if

    ! for disttype=3, the cradius and sradius have to passed in meters,
    ! so I multiply by 1000 to be able to put it in km in the input file

    if (thisobs%disttype==3) then
      CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
          locweight, cradius_SIF*1000.0, sradius_SIF*1000.0, dim_obs_l)
    else
      CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
          locweight, cradius_SIF, sradius_SIF, dim_obs_l)
    end if


    END SUBROUTINE init_dim_obs_l_SIF



  !-------------------------------------------------------------------------------
  !> Perform covariance localization for local EnKF on the module-type observation
  !!
  !! The routine is called in the analysis step of the localized
  !! EnKF. It has to apply localization to the two matrices
  !! HP and HPH of the analysis step for the module-type
  !! observation.
  !!
  !! This routine calls the routine PDAFomi_localize_covar
  !! for each observation type. The call allows to specify a
  !! different localization radius and localization functions
  !! for each observation type.
  !!
    SUBROUTINE localize_covar_SIF(dim_p, dim_obs, HP_p, HPH, coords_p)

      ! Include PDAFomi function
      USE PDAFomi, ONLY: PDAFomi_localize_covar

      ! Include localization radius and local coordinates
      USE mod_assimilation, &
           ONLY: cradius_SIF, locweight, sradius_SIF

      use enkf_clm_mod, only: state_pdaf2clm_c_p

      use shr_kind_mod, only: r8 => shr_kind_r8

      USE GridcellType, ONLY: grc
      USE ColumnType, ONLY : col

      IMPLICIT NONE

  ! *** Arguments ***
      INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
      INTEGER, INTENT(in) :: dim_obs               !< Dimension of observation vector
      REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
      REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
      REAL, INTENT(inout) :: coords_p(:,:)         !< Coordinates of state vector elements

      integer :: i

      real(r8), pointer :: lon(:)
      real(r8), pointer :: lat(:)
      integer, pointer :: mycgridcell(:) !Pointer for CLM3.5/CLM5.0 col->gridcell index arrays

#ifdef CLMFIVE
      ! Obtain CLM lon/lat information
      lon   => grc%londeg
      lat   => grc%latdeg
      ! Obtain CLM column-gridcell information
      mycgridcell => col%gridcell
#else
      lon   => clm3%g%londeg
      lat   => clm3%g%latdeg
      mycgridcell => clm3%g%l%c%gridcell
#endif


  ! *************************************
  ! *** Apply covariance localization ***
  ! *************************************

      do i = 1,dim_p
        if (lon(mycgridcell(state_pdaf2clm_c_p(i)))>180) then
          coords_p(1,i) = lon(mycgridcell(state_pdaf2clm_c_p(i))) - 360.0
        else
          coords_p(1,i) = lon(mycgridcell(state_pdaf2clm_c_p(i)))
        end if
        coords_p(2,i) = lat(mycgridcell(state_pdaf2clm_c_p(i)))
      end do



      CALL PDAFomi_localize_covar(thisobs, dim_p, locweight, cradius_SIF, sradius_SIF, &
           coords_p, HP_p, HPH)

    END SUBROUTINE localize_covar_SIF


    SUBROUTINE add_obs_err_SIF(step, dim_obs, C)

        USE mod_parallel_pdaf, &
          ONLY: npes_filter

        use PDAFomi, only: obsdims

        implicit none
        INTEGER, INTENT(in) :: step       ! Current time step
        INTEGER, INTENT(in) :: dim_obs  ! Dimension of observation vector
        REAL, INTENT(inout) :: C(dim_obs,dim_obs) ! Matrix to that
                                        ! observation covariance R is added
        integer :: i, pe, cnt
        INTEGER, ALLOCATABLE :: id_start(:) ! Start index of obs. type in global averall obs. vector
        INTEGER, ALLOCATABLE :: id_end(:)   ! End index of obs. type in global averall obs. vector

        ALLOCATE(id_start(npes_filter), id_end(npes_filter))

        ! Initialize indices --> we only have information about local
        ! obs. dims per PE, so we get the global indices, more
        ! generalizable than using the afrrays initiliazed in
        ! init_dim_obs_SM as we can also consider different
        ! observation types in one observation file. Arrays from
        ! init_dim_obs_pdaf (e.g. obs_nc2pdaf) may not be necessary
        ! anymore, @ Johannes, please have a check here., see also in
        ! PDAFomi_obs_f.F90, there the same code is used
        !
        ! addition: I also use now the obs_pdaf2nc for reordering the
        ! observation covariance matrix to the PDAF internal order
        !
        ! So for an obs type where correlations should be accounted
        ! for, this should not be removed!

        pe = 1
        id_start(1) = 1
        IF (thisobs%obsid>1) id_start(1) = id_start(1) + sum(obsdims(1, 1:thisobs%obsid-1))
        id_end(1)   = id_start(1) + obsdims(1,thisobs%obsid) - 1
        DO pe = 2, npes_filter
          id_start(pe) = id_start(pe-1) + SUM(obsdims(pe-1,thisobs%obsid:))
          IF (thisobs%obsid>1) id_start(pe) = id_start(pe) + sum(obsdims(pe,1:thisobs%obsid-1))
          id_end(pe) = id_start(pe) + obsdims(pe,thisobs%obsid) - 1
        END DO


        cnt = 1
        DO pe = 1, npes_filter
          DO i = id_start(pe), id_end(pe)
            C(i,i) = C(i,i) + 1.0/thisobs%ivar_obs_f(cnt)
            cnt = cnt + 1
          end do
        end do

        DEALLOCATE(id_start, id_end)

    END SUBROUTINE add_obs_err_SIF


    SUBROUTINE init_obscovar_SIF(step, dim_obs, dim_obs_p, covar, m_state_p, isdiag)

        USE mod_parallel_pdaf, &
          ONLY: npes_filter

        use PDAFomi, only: obsdims, map_obs_id

        implicit none
        INTEGER, INTENT(in) :: step                ! Current time step
        INTEGER, INTENT(in) :: dim_obs             ! Dimension of observation vector
        INTEGER, INTENT(in) :: dim_obs_p           ! PE-local dimension of observation vector
        REAL, INTENT(inout) :: covar(dim_obs, dim_obs) ! Observation error covariance matrix
        REAL, INTENT(in)  :: m_state_p(dim_obs_p)  ! PE-local observation vector
        LOGICAL, INTENT(inout) :: isdiag             ! Whether the observation error covar. matrix is diagonal

        integer :: i, pe, cnt
        INTEGER, ALLOCATABLE :: id_start(:) ! Start index of obs. type in global averall obs. vector
        INTEGER, ALLOCATABLE :: id_end(:)   ! End index of obs. type in global averall obs. vector

        ALLOCATE(id_start(npes_filter), id_end(npes_filter))

        ! Initialize indices --> we only have information about local
        ! obs. dims per PE, so we use the same logic as in
        ! add_obs_err_SIF
        pe = 1
        id_start(1) = 1
        IF (thisobs%obsid>1) id_start(1) = id_start(1) + sum(obsdims(1, 1:thisobs%obsid-1))
        id_end(1)   = id_start(1) + obsdims(1,thisobs%obsid) - 1
        DO pe = 2, npes_filter
            id_start(pe) = id_start(pe-1) + SUM(obsdims(pe-1,thisobs%obsid:))
            IF (thisobs%obsid>1) id_start(pe) = id_start(pe) + sum(obsdims(pe,1:thisobs%obsid-1))
            id_end(pe) = id_start(pe) + obsdims(pe,thisobs%obsid) - 1
        END DO

        ! Initialize mapping vector (to be used in
        ! PDAF_enkf_obs_ensemble) --> has to be initialized here, else
        ! there will be errors!
        cnt = 1
        IF (thisobs%obsid-1 > 0) cnt = cnt+ SUM(obsdims(:,1:thisobs%obsid-1))
        DO pe = 1, npes_filter
            DO i = id_start(pe), id_end(pe)
              map_obs_id(i) = cnt
              cnt = cnt + 1
            END DO
        END DO

        cnt = 1
        DO pe = 1, npes_filter
            DO i = id_start(pe), id_end(pe)
              ! the inverse of the observation variance is saved for
              ! each observation, so we do not need any other
              covar(i, i) = covar(i, i) + 1.0/thisobs%ivar_obs_f(cnt)
              ! array here. As we initiliazed the indices for each
              ! process, we also can just take index cnt instead of
              ! complicated mapping between nc and pdaf indices
              cnt = cnt + 1
            ENDDO
        ENDDO

        ! The matrix is diagonal
        ! This setting avoids the computation of the SVD of COVAR
        ! in PDAF_enkf_obs_ensemble
        isdiag = .TRUE.

        DEALLOCATE(id_start, id_end)


    END SUBROUTINE init_obscovar_SIF


    SUBROUTINE prodRinvA_SIF(step, dim_obs_p, rank, obs_p, A_p, C_p)

      INTEGER, INTENT(in) :: step                ! Current time step
      INTEGER, INTENT(in) :: dim_obs_p           ! PE-local dimension of obs. vector
      INTEGER, INTENT(in) :: rank                ! Rank of initial covariance matrix
      REAL, INTENT(in)    :: obs_p(dim_obs_p)    ! PE-local vector of observations
      REAL, INTENT(in)    :: A_p(dim_obs_p,rank) ! Input matrix from analysis routine
      REAL, INTENT(inout)   :: C_p(dim_obs_p,rank) ! Output matrix

      INTEGER :: i, j       ! index of observation component
      INTEGER :: off        ! row offset in A_l and C_l

      off = thisobs%off_obs_f

      do j=1, rank
        do i=1, thisobs%dim_obs_f
          C_p(i+off, j) = thisobs%ivar_obs_f(i) * A_p(i+off, j)
        END DO
      end do

    END SUBROUTINE prodRinvA_SIF


    SUBROUTINE prodRinvA_l_SIF(domain_p, step, dim_obs, rank, obs_l, A_l, C_l)

        use shr_kind_mod, only: r8 => shr_kind_r8
        USE mod_assimilation, &
           ONLY: cradius_SM, locweight, sradius_SM
        use pdafomi, only: PDAFomi_observation_localization_weights

        implicit none

        INTEGER, INTENT(in) :: domain_p             ! Current local analysis domain
        INTEGER, INTENT(in) :: step                 ! Current time step
        INTEGER, INTENT(in) :: dim_obs             ! Dimension of local observation vector, multiple observation types possible,
                                                   ! then we have to access with thisobs_l%dim_obs_l
        INTEGER, INTENT(in) :: rank                 ! Rank of initial covariance matrix
        REAL, INTENT(in)    :: obs_l(dim_obs)     ! Local vector of observations
        REAL, INTENT(inout) :: A_l(dim_obs, rank) ! Input matrix from analysis routine
        REAL, INTENT(out)   :: C_l(dim_obs, rank) ! Output matrix

        INTEGER :: verbose       ! Verbosity flag
        INTEGER :: verbose_w     ! Verbosity flag for weight computation
        INTEGER, SAVE :: domain_save = -1  ! Save previous domain index
        INTEGER :: wtype         ! Type of weight function
        INTEGER :: rtype         ! Type of weight regulation
        REAL, ALLOCATABLE :: weight(:)     ! Localization weights
        REAL, ALLOCATABLE :: A_obs(:,:)    ! Array for a single row of A_l
        REAL    :: var_obs                 ! Variance of observation error

        INTEGER :: i, j

        INTEGER :: off                     ! row offset in A_l and C_l
        INTEGER :: idummy                  ! Dummy to access nobs_all

        real(r8) :: ivariance_obs


        off = thisobs_l%off_obs_l
        idummy = dim_obs

        IF ((domain_p <= domain_save .OR. domain_save < 0) .AND. mype_filter==0) THEN
            verbose = 1
        ELSE
            verbose = 0
        END IF
        domain_save = domain_p

        ! Screen output
        IF (verbose == 1) THEN
            WRITE (*, '(8x, a, f12.3)') &
                '--- Use global rms for observations of ', rms_obs_SIF
            WRITE (*, '(8x, a, 1x)') &
                '--- Domain localization'
            WRITE (*, '(12x, a, 1x, f12.2)') &
                '--- Local influence radius', cradius_SIF

            IF (locweight > 0) THEN
                WRITE (*, '(12x, a)') &
                        '--- Use distance-dependent weight for observation errors'

                IF (locweight == 3) THEN
                    write (*, '(12x, a)') &
                        '--- Use regulated weight with mean error variance'
                ELSE IF (locweight == 4) THEN
                    write (*, '(12x, a)') &
                        '--- Use regulated weight with single-point error variance'
                END IF
            END IF
        ENDIF

        ALLOCATE(weight(thisobs_l%dim_obs_l))
        call PDAFomi_observation_localization_weights(thisobs_l, thisobs, rank, A_l, &
                                         weight, verbose)

        do j=1,rank
            do i=1,thisobs_l%dim_obs_l
                C_l(i+off,j) = thisobs_l%ivar_obs_l(i) * weight(i) * A_l(i+off, j)
            end do
        end do

        deallocate(weight)

    END SUBROUTINE prodRinvA_l_SIF

    SUBROUTINE deallocate_obs_SIF()

        USE PDAFomi, ONLY: PDAFomi_deallocate_obs
        USE PDAFomi_obs_l, ONLY: obs_l_all, firstobs

        implicit none

        if (mype_filter==0) then
            WRITE (*,*) 'Deallocating observations type SM'
        end if
        call PDAFomi_deallocate_obs(thisobs)

        if (allocated(thisobs_l%id_obs_l)) deallocate(thisobs_l%id_obs_l)
        if (allocated(thisobs_l%ivar_obs_l)) deallocate(thisobs_l%ivar_obs_l)
        if (allocated(thisobs_l%distance_l)) deallocate(thisobs_l%distance_l)
        if (allocated(thisobs_l%cradius_l)) deallocate(thisobs_l%cradius_l)
        if (allocated(thisobs_l%sradius_l)) deallocate(thisobs_l%sradius_l)
        if (allocated(thisobs_l%dist_l_v)) deallocate(thisobs_l%dist_l_v)
        if (allocated(thisobs_l%cradius)) deallocate(thisobs_l%cradius)
        if (allocated(thisobs_l%sradius)) deallocate(thisobs_l%sradius)

        if (allocated(obs_l_all)) deallocate(obs_l_all)

        firstobs=0

    END SUBROUTINE deallocate_obs_SIF

END MODULE obs_SIF_pdafomi
#endif
