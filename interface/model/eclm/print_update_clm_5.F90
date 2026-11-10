!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)
!
!This file is part of TSMP-PDAF
!
!TSMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TSMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TSMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!print_update_clm_5.F90: Module for printing updated CLM5 ensemble
!-------------------------------------------------------------------------------------------

#if defined CLMSA
subroutine print_update_clm(ts,ttot) bind(C,name="print_update_clm")

    use, intrinsic :: iso_c_binding, only : c_int
    use shr_kind_mod , only : r8 => shr_kind_r8
    use subgridavemod, only : p2g, c2g
    use domainMod    , only : ldomain
    use clm_varpar   , only : nlevsoi
    use clm_varcon   , only : nameg, spval
    use decompmod    , only : get_proc_global, get_proc_bounds, ldecomp, get_proc_total
    use spmdgathscatmod , only : gather_data_to_master
    use spmdmod      , only : masterproc, npes, mpicom, iam
    use mpi, only : mpi_gatherv, mpi_real8
    use clm_time_manager        , only : get_nstep
    use clm_instMod, only : soilstate_inst, waterstate_inst, &
                         bgc_vegetation_inst
    use netcdf, only : nf90_create, NF90_CLOBBER, nf90_def_dim, nf90_def_var, &
                       NF90_DOUBLE, NF90_FLOAT, nf90_enddef, nf90_open, NF90_WRITE, &
                       nf90_inq_varid, nf90_put_var, nf90_close
    use enkf_clm_mod, only : clmupdate_swc, clmupdate_texture, clmprint_swc, &
                              clmupdate_sif, &              ! SIF addition
                              clm_begp, clm_endp, &         ! SIF: patch bounds
                              sif_idx_of_g                  ! SIF: gridcells in the SIF block
    use PatchType,    only : patch                          ! SIF: patch%gridcell, patch%wtgcell

    implicit none

    integer(c_int), intent(in) :: ts, ttot

    ! *** local variables ***
    integer :: numg, numl, numc, nump
    integer :: begg, endg, begl, endl, begc, endc, begp, endp
    integer :: isec, info, jn, jj, ji, g1, jx
    real(r8), pointer :: swc(:,:)
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    real(r8), pointer :: porgm(:,:)
    real(r8), allocatable :: clmstate_out(:,:,:)
    real(r8), allocatable :: clmstate_out_2d(:,:)     ! SIF: 2-D output array
    integer, dimension(4) :: dimids
    integer, dimension(2) :: dimids_2d                ! SIF: dims without z
    integer, dimension(1) :: il_var_id
    integer :: il_file_id
    integer :: ncvarid(6)                             ! SIF: extended to 6 (was 4)
    integer :: status
    character(len=300) :: update_filename
    integer :: nerror
    integer :: ndlon, ndlat

    ! SIF: local variables for patch -> gridcell averaging and gathering
    integer  :: g, p, ier, pid, numsend
    integer  :: ncells, nlunits, ncols, npfts, ncohorts
    integer  :: numrecvv(0:npes-1), displsv(0:npes-1)
    real(r8), allocatable :: lc_g(:), ln_g(:), wt_g(:)
    real(r8), allocatable :: sif_global(:)


    call get_proc_global(ng=numg, nl=numl, nc=numc, np=nump)
    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ndlon = ldomain%ni
    ndlat = ldomain%nj

    if (masterproc) then
        allocate(clmstate_out(ndlon, ndlat, nlevsoi), stat=nerror)
        if (clmupdate_sif/=0) then
            allocate(clmstate_out_2d(ndlon, ndlat), stat=nerror)
            allocate(sif_global(1:numg), stat=nerror)
        end if
    end if

    if (masterproc) then
        call get_update_filename(update_filename)
        if (ts==1) then
            status = nf90_create(update_filename, NF90_CLOBBER, il_file_id)
            status = nf90_def_dim(il_file_id, "x",  ndlon,  dimids(1))
            status = nf90_def_dim(il_file_id, "y",  ndlat,  dimids(2))
            status = nf90_def_dim(il_file_id, "z",  nlevsoi,dimids(3))
            status = nf90_def_dim(il_file_id, "t",  ttot,   dimids(4))

            ! 2-D dims (no z) for SIF variables
            dimids_2d = [ dimids(1), dimids(2) ]

            if (clmprint_swc==1) then
                status = nf90_def_var(il_file_id, "swc",  NF90_DOUBLE, dimids,   ncvarid(1))
            end if
            if (clmupdate_texture==1) then
                status = nf90_def_var(il_file_id, "sand", NF90_DOUBLE, dimids,   ncvarid(2))
                status = nf90_def_var(il_file_id, "clay", NF90_DOUBLE, dimids,   ncvarid(3))
            end if
            if (clmupdate_texture==2) then
                status = nf90_def_var(il_file_id, "sand", NF90_DOUBLE, dimids,   ncvarid(2))
                status = nf90_def_var(il_file_id, "clay", NF90_DOUBLE, dimids,   ncvarid(3))
                status = nf90_def_var(il_file_id, "orgm", NF90_DOUBLE, dimids,   ncvarid(4))
            end if

            ! ==================================================================
            ! SIF: posterior LEAFC [gC/m2] and LEAFN [gN/m2] as (x, y, t)
            ! ==================================================================
            if (clmupdate_sif/=0) then
                status = nf90_def_var(il_file_id, "LEAFC", NF90_FLOAT, &
                                      [dimids(1), dimids(2), dimids(4)], ncvarid(5))
                status = nf90_def_var(il_file_id, "LEAFN", NF90_FLOAT, &
                                      [dimids(1), dimids(2), dimids(4)], ncvarid(6))
            end if

            status = nf90_enddef(il_file_id)
        else
            status = nf90_open(update_filename, NF90_WRITE, il_file_id)
        end if
    end if

    ! ------------------------------------------------------------------
    ! SWC
    ! ------------------------------------------------------------------
    if (clmprint_swc==1) then
        swc => waterstate_inst%h2osoi_vol_col
        if (masterproc) then
            do jn = 1, nlevsoi
                do g1 = 1, numg
                    ji = mod(ldecomp%gdc2glo(g1)-1, ldomain%ni) + 1
                    jj = (ldecomp%gdc2glo(g1) - 1)/ldomain%ni + 1
                    clmstate_out(ji,jj,jn) = swc(g1, jn)
                end do
            end do
            status = nf90_inq_varid(il_file_id, "swc", ncvarid(1))
            status = nf90_put_var(il_file_id, ncvarid(1), clmstate_out(:,:,:), &
                                  start=[1,1,1,ts], count=[ndlon,ndlat,nlevsoi,1])
        end if
    end if

    ! ------------------------------------------------------------------
    ! Texture: sand, clay, orgm 
    ! ------------------------------------------------------------------
    if ((clmupdate_texture==1) .or. (clmupdate_texture==2)) then
        psand => soilstate_inst%cellsand_col
        pclay => soilstate_inst%cellclay_col
        if (masterproc) then
            do jn = 1, nlevsoi
                do g1 = 1, numg
                    ji = mod(ldecomp%gdc2glo(g1)-1, ldomain%ni) + 1
                    jj = (ldecomp%gdc2glo(g1) - 1)/ldomain%ni + 1
                    clmstate_out(ji,jj,jn) = psand(g1,jn)
                end do
            end do
            status = nf90_inq_varid(il_file_id, "sand", ncvarid(2))
            status = nf90_put_var(il_file_id, ncvarid(2), clmstate_out(:,:,:), &
                                  start=[1,1,1,ts], count=[ndlon,ndlat,nlevsoi,1])
            do jn = 1, nlevsoi
                do g1 = 1, numg
                    ji = mod(ldecomp%gdc2glo(g1)-1, ldomain%ni) + 1
                    jj = (ldecomp%gdc2glo(g1) - 1)/ldomain%ni + 1
                    clmstate_out(ji,jj,jn) = pclay(g1, jn)
                end do
            end do
            status = nf90_inq_varid(il_file_id, "clay", ncvarid(3))
            status = nf90_put_var(il_file_id, ncvarid(3), clmstate_out(:,:,:), &
                                  start=[1,1,1,ts], count=[ndlon,ndlat,nlevsoi,1])
        end if
        if (clmupdate_texture==2) then
            porgm => soilstate_inst%cellorg_col
            if (masterproc) then
                do jn = 1, nlevsoi
                    do g1 = 1, numg
                        ji = mod(ldecomp%gdc2glo(g1)-1, ldomain%ni) + 1
                        jj = (ldecomp%gdc2glo(g1) - 1)/ldomain%ni + 1
                        clmstate_out(ji,jj,jn) = porgm(g1, jn)
                    end do
                end do
                status = nf90_inq_varid(il_file_id, "orgm", ncvarid(4))
                status = nf90_put_var(il_file_id, ncvarid(4), clmstate_out(:,:,:), &
                                      start=[1,1,1,ts], count=[ndlon,ndlat,nlevsoi,1])
            end if
        end if
    end if

    ! ==================================================================
    ! SIF: write updated LEAFC and LEAFN (posterior gridcell averages)
    !
    ! Every PE averages its own patches to its gridcells (begg:endg), weighted
    ! with patch%wtgcell, in one pass; the gridcell fields are then gathered on
    ! masterproc with mpi_gatherv (as in print_inc_clm). Gridcells outside the
    ! SIF block are written as spval.
    ! ==================================================================
    if (clmupdate_sif/=0) then

        ! --- gather layout (same as in print_inc_clm)
        call get_proc_total(iam, ncells, nlunits, ncols, npfts, ncohorts)
        numsend = ncells
        do pid = 0, npes-1
            call get_proc_total(pid, ncells, nlunits, ncols, npfts, ncohorts)
            numrecvv(pid) = ncells
        end do
        displsv(0) = 0
        do pid = 1, npes-1
            displsv(pid) = displsv(pid-1) + numrecvv(pid-1)
        end do

        ! --- PE-local gridcell averages (one pass over the patches)
        allocate(lc_g(begg:endg), ln_g(begg:endg), wt_g(begg:endg))
        lc_g = 0.0_r8; ln_g = 0.0_r8; wt_g = 0.0_r8
        do p = clm_begp, clm_endp
            g = patch%gridcell(p)
            if (sif_idx_of_g(g) == 0) cycle
            if (.not. patch%active(p) .or. patch%wtgcell(p) <= 0.0_r8) cycle
            if (bgc_vegetation_inst%cnveg_carbonstate_inst%leafc_patch(p) /= spval) &
                lc_g(g) = lc_g(g) + bgc_vegetation_inst%cnveg_carbonstate_inst%leafc_patch(p) * patch%wtgcell(p)
            if (bgc_vegetation_inst%cnveg_nitrogenstate_inst%leafn_patch(p) /= spval) &
                ln_g(g) = ln_g(g) + bgc_vegetation_inst%cnveg_nitrogenstate_inst%leafn_patch(p) * patch%wtgcell(p)
            wt_g(g) = wt_g(g) + patch%wtgcell(p)
        end do
        do g = begg, endg
            if (wt_g(g) > 0.0_r8) then
                lc_g(g) = lc_g(g) / wt_g(g)
                ln_g(g) = ln_g(g) / wt_g(g)
            else
                lc_g(g) = spval
                ln_g(g) = spval
            end if
        end do

        call write_sif_field(lc_g, "LEAFC", 5)
        call write_sif_field(ln_g, "LEAFN", 6)

        deallocate(lc_g, ln_g, wt_g)

    end if   ! clmupdate_sif/=0

    ! ------------------------------------------------------------------
    ! Close and deallocate
    ! ------------------------------------------------------------------
    if (masterproc) then
        status = nf90_close(il_file_id)
        deallocate(clmstate_out)
        if (allocated(clmstate_out_2d)) deallocate(clmstate_out_2d)
        if (allocated(sif_global))      deallocate(sif_global)
    end if

contains

    !> gather a PE-local gridcell field (begg:endg) on masterproc and write it at time ts
    subroutine write_sif_field(field_local, varname, ivar)
        real(r8),         intent(in) :: field_local(begg:endg)
        character(len=*), intent(in) :: varname
        integer,          intent(in) :: ivar
        if (masterproc) then
            call mpi_gatherv(field_local(begg), numsend, MPI_REAL8, &
                             sif_global, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
        else
            call mpi_gatherv(field_local(begg), numsend, MPI_REAL8, &
                             0._r8, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
        end if
        if (masterproc) then
            clmstate_out_2d(:,:) = spval
            do g1 = 1, numg
                ji = mod(ldecomp%gdc2glo(g1)-1, ldomain%ni) + 1
                jj = (ldecomp%gdc2glo(g1) - 1)/ldomain%ni + 1
                clmstate_out_2d(ji,jj) = sif_global(g1)
            end do
            status = nf90_inq_varid(il_file_id, varname, ncvarid(ivar))
            status = nf90_put_var(il_file_id, ncvarid(ivar), clmstate_out_2d(:,:), &
                                  start=[1,1,ts], count=[ndlon,ndlat,1])
        end if
    end subroutine write_sif_field

end subroutine print_update_clm


!==============================================================================
subroutine print_inc_clm() bind(C,name="print_inc_clm")
!==============================================================================
! Writes DA increments (posterior - prior) to <caseid>.inc.<inst>.<date>.nc.
! SIF addition: writes LEAFC_INC [gC/m2] and LEAFN_INC [gN/m2], the realised
! gridcell-average increments filled by clm_update_sif in enkf_clm_mod.
!==============================================================================

  use shr_kind_mod , only : r8 => shr_kind_r8
  use domainMod    , only : ldomain
  use clm_varpar   , only : nlevsoi
  use clm_varcon   , only : nameg, spval
  use decompmod    , only : get_proc_global, get_proc_bounds, ldecomp, get_proc_total
  use spmdmod      , only : masterproc, npes, mpicom, iam
  use clm_time_manager        , only : get_nstep
  use clm_instMod, only : soilhydrology_inst, waterstate_inst, atm2lnd_inst
  use netcdf, only : nf90_create, NF90_CLOBBER, nf90_def_dim, nf90_def_var, &
                     NF90_FLOAT, nf90_enddef, nf90_open, NF90_WRITE, &
                     nf90_inq_varid, nf90_put_var, nf90_close
  use ColumnType         , only : col
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use mpi, only: mpi_gatherv, mpi_real8
  ! SIF additions
  use enkf_clm_mod, only : clmupdate_sif, &
                            leafc_inc_g, leafn_inc_g   ! filled by clm_update_sif (begg:endg)

  implicit none

  integer :: numg, numl, numc, nump
  integer :: begg, endg, begl, endl, begc, endc, begp, endp
  integer :: ncells, nlunits, ncols, npfts, ncohorts
  integer :: isec, info, jn, jj, ji, g1, jx, c, l, j, g, index, p, count, count2
  real(r8), pointer :: h2osoi_liq(:,:)
  real(r8), pointer :: h2osoi_ice(:,:)
  real(r8), pointer :: h2osno(:)
  real(r8), pointer :: clmstate_tmp_local(:,:)
  real(r8), pointer :: clmstate_tmp_global(:)
  real(r8), allocatable :: clmstate_out(:,:,:)
  real(r8), allocatable :: clmstate_out_2d(:,:)     ! SIF: 2-D increment output
  integer, dimension(3) :: dimids
  integer, dimension(2) :: dimids_1level
  integer, dimension(1) :: il_var_id
  integer :: il_file_id
  integer :: ncvarid(6)                              ! SIF: extended to 6
  integer :: status
  character(len=300) :: inc_filename
  integer :: nerror, ndlon, ndlat

  integer :: ier, beg
  integer :: numrecvv(0:npes-1)
  integer :: displsv(0:npes-1)
  integer :: numsend, pid
  integer :: count_columns
  real(r8) :: sum_columns
  real(r8), allocatable :: tws_inc(:)

  ! SIF: local buffers for gathering 2-D inc fields
  real(r8), allocatable :: sif_inc_local(:)   ! PE-local (begg:endg) buffer
  real(r8), allocatable :: sif_inc_global(:)  ! global (1:numg) on masterproc

  h2osoi_liq => waterstate_inst%h2osoi_liq_col_inc
  h2osoi_ice => waterstate_inst%h2osoi_ice_col_inc
  h2osno     => waterstate_inst%h2osno_col_inc

  call get_proc_global(ng=numg, nl=numl, nc=numc, np=nump)
  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

  allocate(clmstate_tmp_local(begg:endg, 1:nlevsoi), stat=nerror)
  allocate(tws_inc(begg:endg), stat=nerror)
  tws_inc(begg:endg) = 0._r8

  ndlon = ldomain%ni
  ndlat = ldomain%nj

  if (masterproc) then
      allocate(clmstate_tmp_global(1:numg), stat=nerror)
      allocate(clmstate_out(ndlon, ndlat, nlevsoi), stat=nerror)
      clmstate_out(:,:,:) = nan
      if (clmupdate_sif/=0) then
          allocate(clmstate_out_2d(ndlon, ndlat), stat=nerror)
          clmstate_out_2d(:,:) = 0.0_r8
      end if
  end if

  call get_proc_total(iam, ncells, nlunits, ncols, npfts, ncohorts)
  numsend = ncells
  do pid = 0, npes-1
      call get_proc_total(pid, ncells, nlunits, ncols, npfts, ncohorts)
      numrecvv(pid) = ncells
  end do
  beg = begg
  displsv(0) = 0
  do pid = 1, npes-1
      displsv(pid) = displsv(pid-1) + numrecvv(pid-1)
  end do

  ! ------------------------------------------------------------------
  ! Define NetCDF variables
  ! ------------------------------------------------------------------
  if (masterproc) then
      call get_inc_filename(inc_filename)
      status = nf90_create(inc_filename, NF90_CLOBBER, il_file_id)
      status = nf90_def_dim(il_file_id, "lon", ndlon, dimids(1))
      status = nf90_def_dim(il_file_id, "lat", ndlat, dimids(2))
      status = nf90_def_dim(il_file_id, "z",   nlevsoi, dimids(3))

      dimids_1level = [ dimids(1), dimids(2) ]

      status = nf90_def_var(il_file_id, "SOILLIQ", NF90_FLOAT, dimids,        ncvarid(1))
      status = nf90_def_var(il_file_id, "SOILICE", NF90_FLOAT, dimids,        ncvarid(2))
      status = nf90_def_var(il_file_id, "H2OSNO",  NF90_FLOAT, dimids_1level, ncvarid(3))
      status = nf90_def_var(il_file_id, "TWS",     NF90_FLOAT, dimids_1level, ncvarid(4))

      ! ================================================================
      ! SIF: LEAFC_INC [gC/m2] and LEAFN_INC [gN/m2] as 2-D fields (lon, lat)
      ! ================================================================
      if (clmupdate_sif/=0) then
          status = nf90_def_var(il_file_id, "LEAFC_INC", NF90_FLOAT, dimids_1level, ncvarid(5))
          status = nf90_def_var(il_file_id, "LEAFN_INC", NF90_FLOAT, dimids_1level, ncvarid(6))
      end if

      status = nf90_enddef(il_file_id)
  end if

  ! ------------------------------------------------------------------
  ! SOILLIQ increment
  ! ------------------------------------------------------------------
  clmstate_tmp_local(begg:endg,:) = 0._r8
  do j = 1, nlevsoi
      do g = begg, endg
          count_columns = 0
          sum_columns   = 0
          do c = begc, endc
              if (g==col%gridcell(c) .and. col%hydrologically_active(c) &
                                     .and. j<=col%nbedrock(c)) then
                  sum_columns   = sum_columns + h2osoi_liq(c,j)
                  count_columns = count_columns + 1
              end if
          end do
          if (count_columns > 0) clmstate_tmp_local(g,j) = sum_columns/count_columns
          if (j==1) then
              tws_inc(g) = clmstate_tmp_local(g,j)
          else
              if (clmstate_tmp_local(g,j) /= spval) &
                  tws_inc(g) = tws_inc(g) + clmstate_tmp_local(g,j)
          end if
      end do
  end do
  do jn = 1, nlevsoi
      if (masterproc) then
          call mpi_gatherv(clmstate_tmp_local(beg,jn), numsend, MPI_REAL8, &
                           clmstate_tmp_global, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
      else
          call mpi_gatherv(clmstate_tmp_local(beg,jn), numsend, MPI_REAL8, &
                           0._r8, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
      end if
      if (masterproc) then
          do g1 = 1, numg
              ji = mod(ldecomp%gdc2glo(g1)-1, ldomain%ni) + 1
              jj = (ldecomp%gdc2glo(g1) - 1)/ldomain%ni + 1
              clmstate_out(ji,jj,jn) = clmstate_tmp_global(g1)
          end do
      end if
  end do
  if (masterproc) then
      status = nf90_inq_varid(il_file_id, "SOILLIQ", ncvarid(1))
      status = nf90_put_var(il_file_id, ncvarid(1), clmstate_out(:,:,:), &
                            start=[1,1,1], count=[ndlon,ndlat,nlevsoi])
  end if

  ! ------------------------------------------------------------------
  ! SOILICE increment 
  ! ------------------------------------------------------------------
  clmstate_tmp_local(begg:endg,:) = 0._r8
  do j = 1, nlevsoi
      do g = begg, endg
          count_columns = 0
          sum_columns   = 0
          do c = begc, endc
              if (g==col%gridcell(c) .and. col%hydrologically_active(c) &
                                     .and. j<=col%nbedrock(c)) then
                  sum_columns   = sum_columns + h2osoi_ice(c,j)
                  count_columns = count_columns + 1
              end if
          end do
          if (count_columns > 0) clmstate_tmp_local(g,j) = sum_columns/count_columns
          if (clmstate_tmp_local(g,j) /= spval) &
              tws_inc(g) = tws_inc(g) + clmstate_tmp_local(g,j)
      end do
  end do
  do jn = 1, nlevsoi
      if (masterproc) then
          call mpi_gatherv(clmstate_tmp_local(beg,jn), numsend, MPI_REAL8, &
                           clmstate_tmp_global, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
      else
          call mpi_gatherv(clmstate_tmp_local(beg,jn), numsend, MPI_REAL8, &
                           0._r8, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
      end if
      if (masterproc) then
          do g1 = 1, numg
              ji = mod(ldecomp%gdc2glo(g1)-1, ldomain%ni) + 1
              jj = (ldecomp%gdc2glo(g1) - 1)/ldomain%ni + 1
              clmstate_out(ji,jj,jn) = clmstate_tmp_global(g1)
          end do
      end if
  end do
  if (masterproc) then
      status = nf90_inq_varid(il_file_id, "SOILICE", ncvarid(2))
      status = nf90_put_var(il_file_id, ncvarid(2), clmstate_out(:,:,:), &
                            start=[1,1,1], count=[ndlon,ndlat,nlevsoi])
  end if

  ! ------------------------------------------------------------------
  ! H2OSNO increment 
  ! ------------------------------------------------------------------
  clmstate_tmp_local(begg:endg,:) = 0._r8
  do g = begg, endg
      count_columns = 0
      sum_columns   = 0
      do c = begc, endc
          if (g==col%gridcell(c) .and. col%hydrologically_active(c)) then
              sum_columns   = sum_columns + h2osno(c)
              count_columns = count_columns + 1
          end if
      end do
      if (count_columns > 0) then
          clmstate_tmp_local(g,1) = sum_columns/count_columns
          tws_inc(g) = tws_inc(g) + clmstate_tmp_local(g,1)
      end if
  end do
  if (masterproc) then
      call mpi_gatherv(clmstate_tmp_local(beg,1), numsend, MPI_REAL8, &
                       clmstate_tmp_global, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
  else
      call mpi_gatherv(clmstate_tmp_local(beg,1), numsend, MPI_REAL8, &
                       0._r8, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
  end if
  if (masterproc) then
      do g1 = 1, numg
          ji = mod(ldecomp%gdc2glo(g1)-1, ldomain%ni) + 1
          jj = (ldecomp%gdc2glo(g1) - 1)/ldomain%ni + 1
          clmstate_out(ji,jj,1) = clmstate_tmp_global(g1)
      end do
      status = nf90_inq_varid(il_file_id, "H2OSNO", ncvarid(3))
      status = nf90_put_var(il_file_id, ncvarid(3), clmstate_out(:,:,1), &
                            start=[1,1], count=[ndlon,ndlat])
  end if

  ! ------------------------------------------------------------------
  ! TWS increment 
  ! ------------------------------------------------------------------
  if (masterproc) then
      call mpi_gatherv(tws_inc(beg), numsend, MPI_REAL8, &
                       clmstate_tmp_global, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
  else
      call mpi_gatherv(tws_inc(beg), numsend, MPI_REAL8, &
                       0._r8, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
  end if
  if (masterproc) then
      do g1 = 1, numg
          ji = mod(ldecomp%gdc2glo(g1)-1, ldomain%ni) + 1
          jj = (ldecomp%gdc2glo(g1) - 1)/ldomain%ni + 1
          clmstate_out(ji,jj,1) = clmstate_tmp_global(g1)
      end do
      status = nf90_inq_varid(il_file_id, "TWS", ncvarid(4))
      status = nf90_put_var(il_file_id, ncvarid(4), clmstate_out(:,:,1), &
                            start=[1,1], count=[ndlon,ndlat])
  end if

  ! ==================================================================
  ! SIF: write LEAFC_INC and LEAFN_INC
  !
  ! leafc_inc_g / leafn_inc_g (begg:endg) are reset and filled every DA
  ! cycle in clm_update_sif with the realised gridcell-average increment
  ! (after limiting and clamping); gridcells outside the SIF block are 0.
  ! ==================================================================
  if (clmupdate_sif/=0) then

      allocate(sif_inc_local(begg:endg), stat=nerror)
      if (masterproc) allocate(sif_inc_global(1:numg), stat=nerror)

      sif_inc_local(begg:endg) = leafc_inc_g(begg:endg)
      call gather_and_write_inc("LEAFC_INC", 5)

      sif_inc_local(begg:endg) = leafn_inc_g(begg:endg)
      call gather_and_write_inc("LEAFN_INC", 6)

      deallocate(sif_inc_local)
      if (masterproc) deallocate(sif_inc_global)

  end if   ! clmupdate_sif/=0

  ! ------------------------------------------------------------------
  ! Close and deallocate
  ! ------------------------------------------------------------------
  if (masterproc) then
      status = nf90_close(il_file_id)
      deallocate(clmstate_out)
      deallocate(clmstate_tmp_global)
      if (allocated(clmstate_out_2d)) deallocate(clmstate_out_2d)
  end if
  deallocate(tws_inc)
  deallocate(clmstate_tmp_local)

contains

  !> gather sif_inc_local on masterproc and write it as a 2-D field
  subroutine gather_and_write_inc(varname, ivar)
      character(len=*), intent(in) :: varname
      integer,          intent(in) :: ivar
      if (masterproc) then
          call mpi_gatherv(sif_inc_local(beg), numsend, MPI_REAL8, &
                           sif_inc_global, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
      else
          call mpi_gatherv(sif_inc_local(beg), numsend, MPI_REAL8, &
                           0._r8, numrecvv, displsv, MPI_REAL8, 0, mpicom, ier)
      end if
      if (masterproc) then
          clmstate_out_2d(:,:) = 0.0_r8
          do g1 = 1, numg
              ji = mod(ldecomp%gdc2glo(g1)-1, ldomain%ni) + 1
              jj = (ldecomp%gdc2glo(g1) - 1)/ldomain%ni + 1
              clmstate_out_2d(ji,jj) = sif_inc_global(g1)
          end do
          status = nf90_inq_varid(il_file_id, varname, ncvarid(ivar))
          status = nf90_put_var(il_file_id, ncvarid(ivar), clmstate_out_2d(:,:), &
                                start=[1,1], count=[ndlon,ndlat])
      end if
  end subroutine gather_and_write_inc

end subroutine print_inc_clm
#endif


subroutine get_update_filename(iofile)
    use clm_varctl, only : caseid
    use clm_time_manager, only : get_curr_date, get_prev_date
    implicit none
    character(len=300), intent(inout) :: iofile
    character(len=256) :: cdate
    integer :: day, mon, yr, sec

    call get_prev_date(yr, mon, day, sec)
    write(cdate,'(i4.4,"-",i2.2)') yr, mon
    call get_curr_date(yr, mon, day, sec)
    write(cdate,'(i4.4)') yr
    iofile = trim(caseid)//".update."//trim(cdate)//".nc"
end subroutine get_update_filename


subroutine get_inc_filename(iofile)
    use clm_varctl, only : caseid, inst_suffix
    use clm_time_manager, only : get_curr_date, get_prev_date
    implicit none
    character(len=300), intent(inout) :: iofile
    character(len=256) :: cdate
    integer :: day, mon, yr, sec

    call get_prev_date(yr, mon, day, sec)
    write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr, mon, day, sec
    iofile = trim(caseid)//".inc"//trim(inst_suffix)//"."//trim(cdate)//".nc"
end subroutine get_inc_filename
