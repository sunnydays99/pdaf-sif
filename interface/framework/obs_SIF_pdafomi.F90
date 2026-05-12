#ifdef CLMFIVE
MODULE obs_SIF_pdafomi

    USE mod_parallel_pdaf, &
         ONLY: mype_filter    ! Rank of filter process
    USE PDAFomi, &
         ONLY: obs_f, obs_l   ! Declaration of observation data types

    IMPLICIT NONE
    SAVE

    PUBLIC

    ! Variables which are inputs to the module (usually set in init_pdaf)
    LOGICAL :: assim_SIF       !< Whether to assimilate this data type
    REAL    :: rms_obs_SIF     !< Observation error standard deviation (for constant errors)

    ! longitude and latitude of grid cells and observation cells (same pattern as SM)
    INTEGER, ALLOCATABLE :: longxy(:), latixy(:), longxy_obs(:), latixy_obs(:)

    ! -----------------------------------------------------------------------
    ! Module-level SAVE arrays for SIF forward operator.
    !
    ! These are filled once per analysis step in init_dim_obs_SIF via
    ! read_fsif_from_history(), then consumed in obs_op_SIF.
    ! They live here (not in enkf_clm_mod) because they are purely
    ! observation-operator data, not model state data.
    !
    ! Indexed as (1 : endg-begg+1), i.e. local gridcell offset.
    ! -----------------------------------------------------------------------
    REAL, ALLOCATABLE, SAVE :: sifescn_current(:) !< gridcell-avg Nadir Escaping SIF [W/m2/um]
    REAL, ALLOCATABLE, SAVE :: fsno_current(:)    !< gridcell snow fraction [-]

    ! Declare instances of observation data types used here
    TYPE(obs_f), TARGET, PUBLIC :: thisobs      ! full observation
    TYPE(obs_l), TARGET, PUBLIC :: thisobs_l    ! local observation

    !$OMP THREADPRIVATE(thisobs_l)

    !---------------------------------------------------------------------------

CONTAINS

    !! brief Initialize information on the SIF observation type
    !!
    !! Reads gridded satellite SIF observations and eCLM simulated FSIF,
    !! applies physical QC masks (night, snow, no-vegetation), and builds
    !! the PE-local observation arrays for PDAF-OMI.
    !!
    !! Mirrors init_dim_obs_SM exactly at the gridcell level (not column).
    !!
    SUBROUTINE init_dim_obs_SIF(step, dim_obs)

      USE mpi, ONLY: MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM
      USE mod_parallel_pdaf, &
          ONLY: mype_filter, comm_filter, npes_filter, abort_parallel, mype_world
      USE mod_assimilation, &
          ONLY: obs_index_p, obs_filename, &
                obs_pdaf2nc, obs_nc2pdaf, &
                local_dims_obs, local_disp_obs, &
                longxy_obs_floor, latixy_obs_floor, &
                screen, cradius_SIF
      USE PDAFomi,      ONLY: PDAFomi_gather_obs, pi
      USE mod_read_obs, ONLY: multierr, read_obs_nc_type
      USE enkf_clm_mod, ONLY: clmstatevec_allcol, clmstatevec_only_active, &
                               clmstatevec_max_layer, state_clm2pdaf_p, &
                               domain_def_clm
      USE shr_kind_mod, ONLY: r8 => shr_kind_r8
      USE GridcellType, ONLY: grc
      USE clm_varcon,   ONLY: spval, ispval
      USE decompMod,    ONLY: get_proc_bounds, get_proc_global
      USE ColumnType,   ONLY: col
      USE PatchType,    ONLY: patch
      USE mod_tsmp,     ONLY: point_obs, da_print_obs_index

      IMPLICIT NONE

      ! *** Arguments ***
      INTEGER, INTENT(in)    :: step
      INTEGER, INTENT(inout) :: dim_obs

      ! *** Local variables ***
      INTEGER :: i, g, c                  ! counters
      INTEGER :: cnt                      ! counter for filled obs
      INTEGER :: dim_obs_p                ! PE-local obs count
      REAL, ALLOCATABLE :: obs_p(:)       ! PE-local observation values
      REAL, ALLOCATABLE :: obs_g(:)       ! global observation vector (root only)
      REAL, ALLOCATABLE :: ivar_obs_p(:)  ! PE-local inverse obs error variance
      REAL, ALLOCATABLE :: ocoord_p(:,:)  ! PE-local obs coordinates [rad]
      REAL, ALLOCATABLE :: lon_obs(:)     ! obs longitudes
      REAL, ALLOCATABLE :: lat_obs(:)     ! obs latitudes
      REAL, ALLOCATABLE :: dr_obs(:)      ! snapping radius [deg lon, deg lat]
      REAL, ALLOCATABLE :: obserr(:)      ! per-obs error (if multierr==1)
      REAL, ALLOCATABLE :: obscov(:,:)    ! obs covariance (unused for SIF)
      INTEGER, ALLOCATABLE :: layer_obs(:)! obs layer (always 1 for SIF)
      CHARACTER(len=110) :: current_observation_filename
      CHARACTER(len=20)  :: obs_type_name

      INTEGER :: begp, endp, begc, endc, begl, endl, begg, endg
      INTEGER :: numg, numl, numc, nump
      INTEGER :: sum_dim_obs_p, ierror
      REAL(r8), POINTER :: lon(:), lat(:)
      REAL :: deltax, deltay
      LOGICAL :: obs_snapped

      INTEGER :: cg               ! column->gridcell index (used in point_obs block)
      LOGICAL :: newgridcell      ! flag for first column in a gridcell
      LOGICAL :: is_use_dr        ! always .TRUE. for SIF (use deltax/deltay snapping)

      IF (mype_filter==0) &
          WRITE(*,*) 'Assimilate observations - obs type SIF (FSIF at 740 nm)'

      ! Store whether to assimilate this observation type
      IF (assim_SIF) thisobs%doassim = 1

      ! Use geographic (haversine) distance for SIF localisation
      thisobs%disttype = 3
      thisobs%ncoord   = 2

      obs_type_name = 'SIF'
      is_use_dr     = .TRUE.   ! always use deltax/deltay snapping for SIF

      WRITE(current_observation_filename, '(a, i5.5)') TRIM(obs_filename)//'.', step

      ! Root process reads observations from NetCDF
      IF (mype_filter == 0) THEN
          CALL read_obs_nc_type(current_observation_filename, obs_type_name, &
                                dim_obs, obs_g, lon_obs, lat_obs, layer_obs, &
                                dr_obs, obserr, obscov)
      END IF

      CALL mpi_bcast(dim_obs, 1, MPI_INTEGER, 0, comm_filter, ierror)

      ! -----------------------------------------------------------------------
      ! Handle zero-obs case, same as iin obs_SM_pdafomi
      ! -----------------------------------------------------------------------
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

      ! -----------------------------------------------------------------------
      ! Broadcast obs arrays to all filter PEs, same as obs_SM_pdafomi
      ! -----------------------------------------------------------------------
      CALL mpi_bcast(multierr, 1, MPI_INTEGER, 0, comm_filter, ierror)

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

      CALL mpi_bcast(obs_g,     dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
      IF (multierr==1) &
          CALL mpi_bcast(obserr, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
      CALL mpi_bcast(lon_obs,   dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
      CALL mpi_bcast(lat_obs,   dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
      CALL mpi_bcast(dr_obs,    2,       MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
      CALL mpi_bcast(layer_obs, dim_obs, MPI_INTEGER,          0, comm_filter, ierror)

      thisobs%infile = 1
      CALL domain_def_clm(lon_obs, lat_obs, dim_obs, longxy, latixy, longxy_obs, latixy_obs)

      ! Obtain CLM lon/lat arrays
      lon => grc%londeg
      lat => grc%latdeg

      CALL get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
      CALL get_proc_global(numg, numl, numc, nump)

      ! -----------------------------------------------------------------------
      ! Allocate module-level QC arrays and fill them from the live CLM instance.
      ! read_fsif_from_history fills fsif_current, coszen_current, fsno_current
      ! as module-level SAVE arrays indexed (1:endg-begg+1).
      ! It is defined at the bottom of this module.
      ! -----------------------------------------------------------------------
      IF (ALLOCATED(sifescn_current))   DEALLOCATE(sifescn_current)
      IF (ALLOCATED(fsno_current))   DEALLOCATE(fsno_current)
      ALLOCATE(sifescn_current (endg-begg+1))
      ALLOCATE(fsno_current  (endg-begg+1))

      CALL collect_sif_data(begg, endg, begp, endp, begc, endc)
      ! After this call: fsif_current(g-begg+1) + fsno_current(g-begg+1) are populated for all local gridcells.

      ! -----------------------------------------------------------------------
      ! First pass: count PE-local valid observations (dim_obs_p)
      !
      ! SIF maps directly to gridcells, no inner column loop needed.
      ! The loop over g is sufficient; we do not need c tog mapping here.
      ! -----------------------------------------------------------------------
      IF (ALLOCATED(thisobs%id_obs_p)) DEALLOCATE(thisobs%id_obs_p)
      ALLOCATE(thisobs%id_obs_p(1, endg-begg+1))
      thisobs%id_obs_p(1,:) = 0

      dim_obs_p = 0
      cnt       = 1

      DO i = 1, dim_obs
          obs_snapped = .FALSE.
          DO g = begg, endg

              deltax = ABS(lon(g) - lon_obs(i))
              IF (deltax > 180.0) deltax = 360.0 - deltax
              deltay = ABS(lat(g) - lat_obs(i))

              IF ((deltax <= dr_obs(1)) .AND. (deltay <= dr_obs(2))) THEN

                  ! -----------------------------------------------------------
                  ! QC masks — applied before counting this observation
                  ! -----------------------------------------------------------
                  ! 1. Snow-dominated gridcell: satellite SIF unreliable
                  !  mask out high-albedo/snow contamination
                  IF (fsno_current(g-begg+1) > 0.5) CYCLE

                  ! 2. Model SIF invalid (Nighttime / No Veg / Lake)
                  ! sifescn_current <= 0 handles:
                  ! - Nighttime (APAR=0)
                  ! - Bare soil/Urban/Lake
                  ! - Non-growing season
                  IF (sifescn_current(g-begg+1) <= 0.0 .OR. &
                      sifescn_current(g-begg+1) == spval) CYCLE

                  ! 3. Satellite observation negative: physically impossible
                  IF (obs_g(i) <= 0.0) CYCLE

                  ! -----------------------------------------------------------

                  dim_obs_p = dim_obs_p + 1
                  thisobs%id_obs_p(1, g-begg+1) = 1

                  IF (obs_snapped) THEN
                      WRITE(*,*) "TSMP-PDAF mype(w)=", mype_world, &
                          ": ERROR SIF obs snapped at multiple gridcells, i=", i
                      CALL abort_parallel()
                  END IF
                  obs_snapped = .TRUE.
                  cnt = cnt + 1
              END IF
          END DO
      END DO

      ! -----------------------------------------------------------------------
      ! Dimension consistency check and PE-local dim arrays, same as SM
      ! -----------------------------------------------------------------------
      IF (ALLOCATED(ivar_obs_p)) DEALLOCATE(ivar_obs_p)
      ALLOCATE(ivar_obs_p(dim_obs_p))
      IF (ALLOCATED(ocoord_p)) DEALLOCATE(ocoord_p)
      ALLOCATE(ocoord_p(2, dim_obs_p))

      CALL mpi_allreduce(dim_obs_p, sum_dim_obs_p, 1, MPI_INTEGER, MPI_SUM, &
                         comm_filter, ierror)

      IF (sum_dim_obs_p /= dim_obs) THEN
          WRITE(*,*) "TSMP-PDAF mype(w)=", mype_world, &
              ": ERROR SIF: sum of PE-local obs dims does not match dim_obs"
          WRITE(*,*) "sum_dim_obs_p=", sum_dim_obs_p, "dim_obs=", dim_obs
          CALL abort_parallel()
      END IF

      IF (ALLOCATED(local_dims_obs)) DEALLOCATE(local_dims_obs)
      ALLOCATE(local_dims_obs(npes_filter))
      CALL mpi_allgather(dim_obs_p, 1, MPI_INTEGER, local_dims_obs, 1, MPI_INTEGER, &
                         comm_filter, ierror)

      IF (ALLOCATED(local_disp_obs)) DEALLOCATE(local_disp_obs)
      ALLOCATE(local_disp_obs(npes_filter))
      local_disp_obs(1) = 0
      DO i = 2, npes_filter
          local_disp_obs(i) = local_disp_obs(i-1) + local_dims_obs(i-1)
      END DO

      IF (mype_filter==0 .AND. screen > 2) &
          WRITE(*,*) "TSMP-PDAF mype(w)=", mype_world, &
              ": init_dim_obs_SIF: local_disp_obs=", local_disp_obs

      ! -----------------------------------------------------------------------
      ! obs_pdaf2nc / obs_nc2pdaf mapping
      !
      ! (point_obs block): For SIF the mapping is at gridcell level.
      ! We iterate obs to gridcells directly without the inner column loop.
      ! -----------------------------------------------------------------------
      IF (ALLOCATED(obs_pdaf2nc)) DEALLOCATE(obs_pdaf2nc)
      ALLOCATE(obs_pdaf2nc(dim_obs))
      obs_pdaf2nc = 0
      IF (ALLOCATED(obs_nc2pdaf)) DEALLOCATE(obs_nc2pdaf)
      ALLOCATE(obs_nc2pdaf(dim_obs))
      obs_nc2pdaf = 0

      IF (point_obs == 1) THEN

          cnt = 1
          DO i = 1, dim_obs
              obs_snapped = .FALSE.

              DO g = begg, endg

                  deltax = ABS(lon(g) - lon_obs(i))
                  IF (deltax > 180.0) deltax = 360.0 - deltax
                  deltay = ABS(lat(g) - lat_obs(i))

                  IF ((deltax <= dr_obs(1)) .AND. (deltay <= dr_obs(2))) THEN
                      ! Apply same QC as in first pass
                      IF (coszen_current(g-begg+1) < 0.01)   CYCLE
                      IF (fsno_current  (g-begg+1) > 0.5)    CYCLE
                      IF (fsif_current  (g-begg+1) <= 0.0 .OR. &
                          fsif_current  (g-begg+1) == spval) CYCLE
                      IF (obs_g(i) <= 0.0)                   CYCLE

                      ! FIX 4: No ispval/state_clm2pdaf_p check needed for SIF —
                      ! gridcells without vegetation are already masked by
                      ! fsif_current <= 0 above.

                      obs_pdaf2nc(local_disp_obs(mype_filter+1) + cnt) = i
                      obs_nc2pdaf(i) = local_disp_obs(mype_filter+1) + cnt
                      cnt = cnt + 1
                      obs_snapped = .TRUE.
                  END IF
              END DO

              IF (.NOT. obs_snapped) THEN
                  IF (screen > 2) &
                      WRITE(*,*) "TSMP-PDAF mype(w)=", mype_world, &
                          ": SIF obs i=", i, " not snapped on this PE (expected for off-PE obs)"
              END IF
          END DO

      END IF

      ! Collect obs_pdaf2nc / obs_nc2pdaf across PEs (only non-zero on owning PE)
      CALL mpi_allreduce(MPI_IN_PLACE, obs_pdaf2nc, dim_obs, MPI_INTEGER, &
                         MPI_SUM, comm_filter, ierror)
      CALL mpi_allreduce(MPI_IN_PLACE, obs_nc2pdaf, dim_obs, MPI_INTEGER, &
                         MPI_SUM, comm_filter, ierror)

      IF (mype_filter==0 .AND. screen > 2) &
          WRITE(*,*) "TSMP-PDAF mype(w)=", mype_world, &
              ": init_dim_obs_SIF: obs_pdaf2nc=", obs_pdaf2nc

      ! -----------------------------------------------------------------------
      ! Second pass: fill PE-local obs_p, ivar_obs_p, ocoord_p, obs_index_p
      !
      ! obs_index_p(cnt) stores the gridcell offset (g-begg+1).
      ! This is used in obs_op_SIF to index into fsif_current.
      ! -----------------------------------------------------------------------
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
                  IF (fsno_current  (g-begg+1) > 0.5)    CYCLE
                  IF (sifescn_current  (g-begg+1) <= 0.0 .OR. &
                      sifescn_current  (g-begg+1) == spval) CYCLE
                  IF (obs_g(i) <= 0.0)                   CYCLE

                  ! Haversine: coordinates must be in radians
                  ocoord_p(1, cnt) = lon_obs(i) * pi / 180.0
                  ocoord_p(2, cnt) = lat_obs(i) * pi / 180.0

                  ! Store gridcell offset for obs_op_SIF lookup of fsif_current
                  obs_index_p(cnt) = g - begg + 1

                  obs_p(cnt) = obs_g(i)
                  IF (multierr == 1) THEN
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

      ! Clean up
      DEALLOCATE(obs_g, obs_p, ocoord_p, ivar_obs_p)
      IF (ALLOCATED(lon_obs))   DEALLOCATE(lon_obs)
      IF (ALLOCATED(lat_obs))   DEALLOCATE(lat_obs)
      IF (ALLOCATED(dr_obs))    DEALLOCATE(dr_obs)
      IF (ALLOCATED(layer_obs)) DEALLOCATE(layer_obs)
      IF (ALLOCATED(obserr))    DEALLOCATE(obserr)
      IF (ALLOCATED(obscov))    DEALLOCATE(obscov)

    END SUBROUTINE init_dim_obs_SIF


    !---------------------------------------------------------------------------
    !! SIF observation operator
    !!
    !! Maps the model state to SIF observation space by reading from the
    !! module-level fsif_current array (set in init_dim_obs_SIF).
    !!
    !! This is a direct lookup: obs_index_p(i) holds the gridcell offset
    !! (g-begg+1) set during init_dim_obs_SIF. No state_p indexing is needed
    !! because FSIF is a diagnostic variable recomputed each CLM step and is
    !! NOT in the prognostic state vector.
    !!
    SUBROUTINE obs_op_SIF(dim_p, dim_obs, state_p, ostate)

      ! fsif_current is a module-level SAVE array in this module,
      USE mod_assimilation,  ONLY: obs_index_p
      USE PDAFomi_obs_f,     ONLY: PDAFomi_gather_obsstate

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

      ! Direct lookup into the module-level FSIF array
      DO i = 1, thisobs%dim_obs_p
          ostate_p(i) = fsif_current(obs_index_p(i))
      END DO

      CALL PDAFomi_gather_obsstate(thisobs, ostate_p, ostate)
      DEALLOCATE(ostate_p)

    END SUBROUTINE obs_op_SIF


    !---------------------------------------------------------------------------
    !! Initialize local SIF observation dimension 
    !!
    !! Identical structure to init_dim_obs_l_SM; uses cradius_SIF/sradius_SIF.
    !!
    SUBROUTINE init_dim_obs_l_SIF(domain_p, step, dim_obs, dim_obs_l)

      USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l, pi
      USE mod_assimilation, &
          ONLY: cradius_SIF, locweight, sradius_SIF, screen
      USE enkf_clm_mod, ONLY: state_loc2clm_c_p, clmstatevec_allcol, &
                               clmstatevec_only_active
      USE shr_kind_mod, ONLY: r8 => shr_kind_r8
      USE decompMod,    ONLY: get_proc_bounds
      USE GridcellType, ONLY: grc
      USE ColumnType,   ONLY: col
      USE clm_varcon,   ONLY: spval
      USE mod_parallel_pdaf, ONLY: mype_world

      IMPLICIT NONE

      INTEGER, INTENT(in)    :: domain_p
      INTEGER, INTENT(in)    :: step
      INTEGER, INTENT(in)    :: dim_obs
      INTEGER, INTENT(inout) :: dim_obs_l

      REAL :: coords_l(2)
      REAL(r8), POINTER :: lon(:), lat(:)
      INTEGER, POINTER :: mycgridcell(:)

      lon         => grc%londeg
      lat         => grc%latdeg
      mycgridcell => col%gridcell

      IF (thisobs%infile == 1) THEN

          IF (clmstatevec_allcol==0 .AND. clmstatevec_only_active==0) THEN
              IF (lon(state_loc2clm_c_p(domain_p)) > 180.0) THEN
                  coords_l(1) = lon(state_loc2clm_c_p(domain_p)) - 360.0
              ELSE
                  coords_l(1) = lon(state_loc2clm_c_p(domain_p))
              END IF
              coords_l(2) = lat(state_loc2clm_c_p(domain_p))
          ELSE
              IF (lon(mycgridcell(state_loc2clm_c_p(domain_p))) > 180.0) THEN
                  coords_l(1) = lon(mycgridcell(state_loc2clm_c_p(domain_p))) - 360.0
              ELSE
                  coords_l(1) = lon(mycgridcell(state_loc2clm_c_p(domain_p)))
              END IF
              coords_l(2) = lat(mycgridcell(state_loc2clm_c_p(domain_p)))
          END IF

          ! Haversine requires radians
          IF (thisobs%disttype == 3) THEN
              coords_l(1) = coords_l(1) * pi / 180.0
              coords_l(2) = coords_l(2) * pi / 180.0
          END IF

      ELSE
          coords_l(1) = spval
          coords_l(2) = spval
      END IF

      ! cradius/sradius in km → pass in metres (×1000) for haversine
      IF (thisobs%disttype == 3) THEN
          CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
              locweight, cradius_SIF*1000.0, sradius_SIF*1000.0, dim_obs_l)
      ELSE
          CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
              locweight, cradius_SIF, sradius_SIF, dim_obs_l)
      END IF

    END SUBROUTINE init_dim_obs_l_SIF


    !---------------------------------------------------------------------------
    !! Covariance localisation 
    !!
    !! Identical structure to localize_covar_SM; uses cradius_SIF/sradius_SIF.
    !!
    SUBROUTINE localize_covar_SIF(dim_p, dim_obs, HP_p, HPH, coords_p)

      USE PDAFomi, ONLY: PDAFomi_localize_covar
      USE mod_assimilation, &
          ONLY: cradius_SIF, locweight, sradius_SIF
      USE enkf_clm_mod, ONLY: state_pdaf2clm_c_p
      USE shr_kind_mod, ONLY: r8 => shr_kind_r8
      USE GridcellType, ONLY: grc
      USE ColumnType,   ONLY: col

      IMPLICIT NONE

      INTEGER, INTENT(in)    :: dim_p, dim_obs
      REAL, INTENT(inout)    :: HP_p(dim_obs, dim_p)
      REAL, INTENT(inout)    :: HPH(dim_obs, dim_obs)
      REAL, INTENT(inout)    :: coords_p(:,:)

      INTEGER :: i
      REAL(r8), POINTER :: lon(:), lat(:)
      INTEGER, POINTER  :: mycgridcell(:)

      lon         => grc%londeg
      lat         => grc%latdeg
      mycgridcell => col%gridcell

      DO i = 1, dim_p
          IF (lon(mycgridcell(state_pdaf2clm_c_p(i))) > 180.0) THEN
              coords_p(1,i) = lon(mycgridcell(state_pdaf2clm_c_p(i))) - 360.0
          ELSE
              coords_p(1,i) = lon(mycgridcell(state_pdaf2clm_c_p(i)))
          END IF
          coords_p(2,i) = lat(mycgridcell(state_pdaf2clm_c_p(i)))
      END DO

      CALL PDAFomi_localize_covar(thisobs, dim_p, locweight, &
                                   cradius_SIF, sradius_SIF, coords_p, HP_p, HPH)

    END SUBROUTINE localize_covar_SIF


    !---------------------------------------------------------------------------
    ! add_obs_err_SIF, init_obscovar_SIF, prodRinvA_SIF:
    ! Copied from obs_SM_pdafomi: no changes needed as the
    ! logic uses thisobs%obsid / thisobs%ivar_obs_f which are type-generic.
    !---------------------------------------------------------------------------

    SUBROUTINE add_obs_err_SIF(step, dim_obs, C)

      USE mod_parallel_pdaf, ONLY: npes_filter
      USE PDAFomi, ONLY: obsdims

      IMPLICIT NONE
      INTEGER, INTENT(in)    :: step, dim_obs
      REAL,    INTENT(inout) :: C(dim_obs, dim_obs)

      INTEGER :: i, pe, cnt
      INTEGER, ALLOCATABLE :: id_start(:), id_end(:)

      ALLOCATE(id_start(npes_filter), id_end(npes_filter))

      pe = 1
      id_start(1) = 1
      IF (thisobs%obsid > 1) &
          id_start(1) = id_start(1) + SUM(obsdims(1, 1:thisobs%obsid-1))
      id_end(1) = id_start(1) + obsdims(1, thisobs%obsid) - 1
      DO pe = 2, npes_filter
          id_start(pe) = id_start(pe-1) + SUM(obsdims(pe-1, thisobs%obsid:))
          IF (thisobs%obsid > 1) &
              id_start(pe) = id_start(pe) + SUM(obsdims(pe, 1:thisobs%obsid-1))
          id_end(pe) = id_start(pe) + obsdims(pe, thisobs%obsid) - 1
      END DO

      cnt = 1
      DO pe = 1, npes_filter
          DO i = id_start(pe), id_end(pe)
              C(i,i) = C(i,i) + 1.0 / thisobs%ivar_obs_f(cnt)
              cnt = cnt + 1
          END DO
      END DO

      DEALLOCATE(id_start, id_end)

    END SUBROUTINE add_obs_err_SIF


    SUBROUTINE init_obscovar_SIF(step, dim_obs, dim_obs_p, covar, m_state_p, isdiag)

      USE mod_parallel_pdaf, ONLY: npes_filter
      USE PDAFomi, ONLY: obsdims, map_obs_id

      IMPLICIT NONE
      INTEGER, INTENT(in)    :: step, dim_obs, dim_obs_p
      REAL,    INTENT(inout) :: covar(dim_obs, dim_obs)
      REAL,    INTENT(in)    :: m_state_p(dim_obs_p)
      LOGICAL, INTENT(inout) :: isdiag

      INTEGER :: i, pe, cnt
      INTEGER, ALLOCATABLE :: id_start(:), id_end(:)

      ALLOCATE(id_start(npes_filter), id_end(npes_filter))

      pe = 1
      id_start(1) = 1
      IF (thisobs%obsid > 1) &
          id_start(1) = id_start(1) + SUM(obsdims(1, 1:thisobs%obsid-1))
      id_end(1) = id_start(1) + obsdims(1, thisobs%obsid) - 1
      DO pe = 2, npes_filter
          id_start(pe) = id_start(pe-1) + SUM(obsdims(pe-1, thisobs%obsid:))
          IF (thisobs%obsid > 1) &
              id_start(pe) = id_start(pe) + SUM(obsdims(pe, 1:thisobs%obsid-1))
          id_end(pe) = id_start(pe) + obsdims(pe, thisobs%obsid) - 1
      END DO

      cnt = 1
      IF (thisobs%obsid-1 > 0) cnt = cnt + SUM(obsdims(:, 1:thisobs%obsid-1))
      DO pe = 1, npes_filter
          DO i = id_start(pe), id_end(pe)
              map_obs_id(i) = cnt
              cnt = cnt + 1
          END DO
      END DO

      cnt = 1
      DO pe = 1, npes_filter
          DO i = id_start(pe), id_end(pe)
              covar(i,i) = covar(i,i) + 1.0 / thisobs%ivar_obs_f(cnt)
              cnt = cnt + 1
          END DO
      END DO

      isdiag = .TRUE.  ! SIF obs errors are uncorrelated (diagonal R)

      DEALLOCATE(id_start, id_end)

    END SUBROUTINE init_obscovar_SIF


    SUBROUTINE prodRinvA_SIF(step, dim_obs_p, rank, obs_p, A_p, C_p)

      IMPLICIT NONE
      INTEGER, INTENT(in)  :: step, dim_obs_p, rank
      REAL,    INTENT(in)  :: obs_p(dim_obs_p)
      REAL,    INTENT(in)  :: A_p(dim_obs_p, rank)
      REAL,    INTENT(inout) :: C_p(dim_obs_p, rank)

      INTEGER :: i, j, off

      off = thisobs%off_obs_f

      DO j = 1, rank
          DO i = 1, thisobs%dim_obs_f
              C_p(i+off, j) = thisobs%ivar_obs_f(i) * A_p(i+off, j)
          END DO
      END DO

    END SUBROUTINE prodRinvA_SIF


    !---------------------------------------------------------------------------
    !>  Apply R^{-1} A locally for LESTKF
    !!
    !! use cradius_SIF/sradius_SIF.
    !!
    SUBROUTINE prodRinvA_l_SIF(domain_p, step, dim_obs, rank, obs_l, A_l, C_l)

      USE shr_kind_mod, ONLY: r8 => shr_kind_r8
      USE mod_assimilation, ONLY: cradius_SIF, locweight, sradius_SIF
      USE PDAFomi, ONLY: PDAFomi_observation_localization_weights

      IMPLICIT NONE

      INTEGER, INTENT(in)    :: domain_p, step, dim_obs, rank
      REAL,    INTENT(in)    :: obs_l(dim_obs)
      REAL,    INTENT(inout) :: A_l(dim_obs, rank)
      REAL,    INTENT(out)   :: C_l(dim_obs, rank)

      INTEGER, SAVE :: domain_save = -1
      INTEGER :: verbose, i, j, off, idummy
      REAL, ALLOCATABLE :: weight(:)

      off    = thisobs_l%off_obs_l
      idummy = dim_obs

      IF ((domain_p <= domain_save .OR. domain_save < 0) .AND. mype_filter==0) THEN
          verbose = 1
      ELSE
          verbose = 0
      END IF
      domain_save = domain_p

      IF (verbose == 1) THEN
          WRITE(*, '(8x, a, f12.3)') &
              '--- SIF obs rms: ', rms_obs_SIF
          WRITE(*, '(12x, a, f12.2)') &
              '--- SIF local influence radius [km]: ', cradius_SIF
      END IF

      ALLOCATE(weight(thisobs_l%dim_obs_l))
      CALL PDAFomi_observation_localization_weights(thisobs_l, thisobs, rank, A_l, &
                                                    weight, verbose)

      DO j = 1, rank
          DO i = 1, thisobs_l%dim_obs_l
              C_l(i+off, j) = thisobs_l%ivar_obs_l(i) * weight(i) * A_l(i+off, j)
          END DO
      END DO

      DEALLOCATE(weight)

    END SUBROUTINE prodRinvA_l_SIF


    !---------------------------------------------------------------------------
    !! Deallocate all SIF observation arrays
    !!
    !! Also deallocates the module-level SIF diagnostic arrays.
    !!
    SUBROUTINE deallocate_obs_SIF()

      USE PDAFomi,       ONLY: PDAFomi_deallocate_obs
      USE PDAFomi_obs_l, ONLY: obs_l_all, firstobs

      IMPLICIT NONE

      IF (mype_filter==0) &
          WRITE(*,*) 'Deallocating observations type SIF'   

      CALL PDAFomi_deallocate_obs(thisobs)

      IF (ALLOCATED(thisobs_l%id_obs_l))   DEALLOCATE(thisobs_l%id_obs_l)
      IF (ALLOCATED(thisobs_l%ivar_obs_l)) DEALLOCATE(thisobs_l%ivar_obs_l)
      IF (ALLOCATED(thisobs_l%distance_l)) DEALLOCATE(thisobs_l%distance_l)
      IF (ALLOCATED(thisobs_l%cradius_l))  DEALLOCATE(thisobs_l%cradius_l)
      IF (ALLOCATED(thisobs_l%sradius_l))  DEALLOCATE(thisobs_l%sradius_l)
      IF (ALLOCATED(thisobs_l%dist_l_v))   DEALLOCATE(thisobs_l%dist_l_v)
      IF (ALLOCATED(thisobs_l%cradius))    DEALLOCATE(thisobs_l%cradius)
      IF (ALLOCATED(thisobs_l%sradius))    DEALLOCATE(thisobs_l%sradius)

      IF (ALLOCATED(obs_l_all)) DEALLOCATE(obs_l_all)

      ! Deallocate module-level SIF diagnostic arrays (FIX 8)
      IF (ALLOCATED(fsif_current))   DEALLOCATE(fsif_current)
      IF (ALLOCATED(coszen_current)) DEALLOCATE(coszen_current)
      IF (ALLOCATED(fsno_current))   DEALLOCATE(fsno_current)

      firstobs = 0

    END SUBROUTINE deallocate_obs_SIF


    !===========================================================================
    ! read FSIF, COSZEN, FSNO from the live CLM instance
    !===========================================================================

    !> @brief Fill module-level fsif_current, coszen_current, fsno_current
    !!
    !! FSIF is a patch-level variable (ptr_patch in hist_addfld1d).
    !! We compute a patch-area-weighted average over all vegetated patches
    !! within each local gridcell. This matches what a nadir-viewing satellite
    !! observes: the area-weighted top-of-canopy SIF from the gridcell footprint.
    !!
    !! COSZEN and FSNO are already at gridcell / column level and are averaged
    !! over hydrologically active columns.
    !!
    !! All three arrays are indexed as (1 : endg-begg+1) = local offset.
    !! Must be called AFTER ALLOCATE of fsif_current etc. in init_dim_obs_SIF.
    !!
    SUBROUTINE read_fsif_from_history(begg, endg, begp, endp, begc, endc)

      USE shr_kind_mod,  ONLY: r8 => shr_kind_r8
      USE PatchType,     ONLY: patch
      USE ColumnType,    ONLY: col
      USE clm_varcon,    ONLY: spval
      USE clm_instMod,   ONLY: photosyns_inst, waterstate_inst, solarabs_inst
      ! photosyns_inst%fsif_patch        — FSIF [W/m2/um] at 740 nm
      ! waterstate_inst%frac_sno_col     — column snow cover fraction
      ! solarabs_inst%coszen_patch or grc%coszen — cosine solar zenith

      ! If grc%coszen is not available in your eCLM version, use:
      ! USE clm_instMod, ONLY: atm2lnd_inst
      ! atm2lnd_inst%forc_solad_grc(:,1) + forc_solai_grc(:,1) > 0 as proxy
      USE GridcellType,  ONLY: grc

      IMPLICIT NONE

      INTEGER, INTENT(in) :: begg, endg, begp, endp, begc, endc

      INTEGER  :: g, p, c
      REAL(r8) :: wt_sum, fsif_sum, fsno_sum, coszen_val
      INTEGER  :: n_col_active

      DO g = begg, endg

          ! --- FSIF: patch-area-weighted average over vegetated patches ---
          fsif_sum = 0.0_r8
          wt_sum   = 0.0_r8

          DO p = begp, endp
              IF (patch%gridcell(p) == g) THEN
                  ! patch%wtgcell: fractional weight of this patch in the gridcell
                  IF (patch%wtgcell(p) > 0.0_r8 .AND. &
                      photosyns_inst%fsif_patch(p) /= spval .AND. &
                      photosyns_inst%fsif_patch(p) > 0.0_r8) THEN
                      fsif_sum = fsif_sum + photosyns_inst%fsif_patch(p) * patch%wtgcell(p)
                      wt_sum   = wt_sum   + patch%wtgcell(p)
                  END IF
              END IF
          END DO

          IF (wt_sum > 0.0_r8) THEN
              fsif_current(g-begg+1) = REAL(fsif_sum / wt_sum)
          ELSE
              ! Lake, urban, or bare soil gridcell — no SIF possible
              fsif_current(g-begg+1) = 0.0
          END IF

          ! --- COSZEN: gridcell cosine solar zenith angle ---
          ! grc%coszen is populated by CLM at each time step
          coszen_current(g-begg+1) = REAL(grc%coszen(g))

          ! --- FSNO: column-average snow fraction ---
          ! Average over all hydrologically active columns in this gridcell
          fsno_sum    = 0.0_r8
          n_col_active = 0

          DO c = begc, endc
              IF (col%gridcell(c) == g .AND. col%hydrologically_active(c)) THEN
                  fsno_sum     = fsno_sum + waterstate_inst%frac_sno_col(c)
                  n_col_active = n_col_active + 1
              END IF
          END DO

          IF (n_col_active > 0) THEN
              fsno_current(g-begg+1) = REAL(fsno_sum / REAL(n_col_active, r8))
          ELSE
              fsno_current(g-begg+1) = 0.0
          END IF

      END DO

    END SUBROUTINE read_fsif_from_history

END MODULE obs_SIF_pdafomi
#endif
