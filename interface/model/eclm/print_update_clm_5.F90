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
    use decompmod    , only : get_proc_global, get_proc_bounds, ldecomp
    use spmdgathscatmod , only : gather_data_to_master
    use spmdmod      , only : masterproc
    use clm_time_manager        , only : get_nstep
    use clm_instMod, only : soilstate_inst, waterstate_inst, &
                             cnveg_carbonstate_inst, canopystate_inst
    use netcdf, only : nf90_create, NF90_CLOBBER, nf90_def_dim, nf90_def_var, &
                       NF90_DOUBLE, NF90_FLOAT, nf90_enddef, nf90_open, NF90_WRITE, &
                       nf90_inq_varid, nf90_put_var, nf90_close
    use enkf_clm_mod, only : clmupdate_swc, clmupdate_texture, clmprint_swc, &
                              clmupdate_sif, &           ! SIF addition
                              clm_begp, clm_endp, &      ! SIF: patch bounds
                              num_hactiveg, hactiveg_levels  ! SIF: active gridcells
    use PatchType,    only : patch                        ! SIF: patch%gridcell, patch%wtgcell

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

    ! SIF: local variables for patch→gridcell averaging
    integer  :: count, g, p
    real(r8) :: leafc_sum, tlai_sum, wt_sum


    call get_proc_global(ng=numg, nl=numl, nc=numc, np=nump)
    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ndlon = ldomain%ni
    ndlat = ldomain%nj

    if (masterproc) then
        allocate(clmstate_out(ndlon, ndlat, nlevsoi), stat=nerror)
        if (clmupdate_sif==1) then
            allocate(clmstate_out_2d(ndlon, ndlat), stat=nerror)
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
            ! SIF: define LEAFC and TLAI as 2-D time-series (lon, lat, t)
            ! Units match clm_update_sif bounds: leafc [gC/m2], tlai [m2/m2]
            ! ==================================================================
            if (clmupdate_sif==1) then
                ! Reuse dimids(1:2) + dimids(4) for (x, y, t) 3-D variables
                status = nf90_def_var(il_file_id, "LEAFC", NF90_FLOAT, &
                                      [dimids(1), dimids(2), dimids(4)], ncvarid(5))
                status = nf90_def_var(il_file_id, "TLAI",  NF90_FLOAT, &
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
    ! SIF: write updated LEAFC and TLAI (posterior gridcell averages)
    !
    ! We compute the patch%wtgcell-weighted average over vegetated patches
    ! in each PE-local gridcell, then scatter to masterproc for output.
    ! Same averaging as set_clm_statevec_sif and read_fsif_from_history.
    ! ==================================================================
    if (clmupdate_sif==1) then

        ! --- LEAFC ---
        if (masterproc) clmstate_out_2d(:,:) = real(spval, r8)

        do count = 1, num_hactiveg
            g       = hactiveg_levels(count, 1)
            leafc_sum = 0.0_r8
            wt_sum    = 0.0_r8
            do p = clm_begp, clm_endp
                if (patch%gridcell(p) == g .and. patch%wtgcell(p) > 0.0_r8) then
                    if (cnveg_carbonstate_inst%leafc_patch(p) /= spval) then
                        leafc_sum = leafc_sum + cnveg_carbonstate_inst%leafc_patch(p) &
                                                * patch%wtgcell(p)
                        wt_sum    = wt_sum   + patch%wtgcell(p)
                    end if
                end if
            end do
            if (masterproc) then
                ji = mod(ldecomp%gdc2glo(g-begg+1)-1, ldomain%ni) + 1
                jj = (ldecomp%gdc2glo(g-begg+1) - 1)/ldomain%ni + 1
                if (wt_sum > 0.0_r8) then
                    clmstate_out_2d(ji,jj) = real(leafc_sum / wt_sum)
                else
                    clmstate_out_2d(ji,jj) = 0.0
                end if
            end if
        end do

        if (masterproc) then
            status = nf90_inq_varid(il_file_id, "LEAFC", ncvarid(5))
            status = nf90_put_var(il_file_id, ncvarid(5), clmstate_out_2d(:,:), &
                                  start=[1,1,ts], count=[ndlon,ndlat,1])
        end if

        ! --- TLAI ---
        if (masterproc) clmstate_out_2d(:,:) = real(spval, r8)

        do count = 1, num_hactiveg
            g      = hactiveg_levels(count, 1)
            tlai_sum = 0.0_r8
            wt_sum   = 0.0_r8
            do p = clm_begp, clm_endp
                if (patch%gridcell(p) == g .and. patch%wtgcell(p) > 0.0_r8) then
                    if (canopystate_inst%tlai_patch(p) /= spval) then
                        tlai_sum = tlai_sum + canopystate_inst%tlai_patch(p) &
                                              * patch%wtgcell(p)
                        wt_sum   = wt_sum   + patch%wtgcell(p)
                    end if
                end if
            end do
            if (masterproc) then
                ji = mod(ldecomp%gdc2glo(g-begg+1)-1, ldomain%ni) + 1
                jj = (ldecomp%gdc2glo(g-begg+1) - 1)/ldomain%ni + 1
                if (wt_sum > 0.0_r8) then
                    clmstate_out_2d(ji,jj) = real(tlai_sum / wt_sum)
                else
                    clmstate_out_2d(ji,jj) = 0.0
                end if
            end if
        end do

        if (masterproc) then
            status = nf90_inq_varid(il_file_id, "TLAI", ncvarid(6))
            status = nf90_put_var(il_file_id, ncvarid(6), clmstate_out_2d(:,:), &
                                  start=[1,1,ts], count=[ndlon,ndlat,1])
        end if

    end if   ! clmupdate_sif==1

    ! ------------------------------------------------------------------
    ! Close and deallocate
    ! ------------------------------------------------------------------
    if (masterproc) then
        status = nf90_close(il_file_id)
        deallocate(clmstate_out)
        if (clmupdate_sif==1 .and. allocated(clmstate_out_2d)) &
            deallocate(clmstate_out_2d)
    end if

end subroutine print_update_clm


!==============================================================================
subroutine print_inc_clm() bind(C,name="print_inc_clm")
!==============================================================================
! Writes DA increments (posterior - prior) to <caseid>.inc.<inst>.<date>.nc.
! SIF addition: writes LEAFC_INC and TLAI_INC using the gridcell-level
! increment arrays filled by clm_update_sif in enkf_clm_mod.
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
                            clm_begg, clm_endg, &
                            num_hactiveg, hactiveg_levels, &
                            leafc_inc_g, tlai_inc_g   ! filled by clm_update_sif

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
      if (clmupdate_sif==1) then
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
      ! SIF: define LEAFC_INC and TLAI_INC as 2-D fields (lon, lat)
      ! Units: LEAFC_INC in gC/m2, TLAI_INC in m2/m2
      ! ================================================================
      if (clmupdate_sif==1) then
          status = nf90_def_var(il_file_id, "LEAFC_INC", NF90_FLOAT, dimids_1level, ncvarid(5))
          status = nf90_def_var(il_file_id, "TLAI_INC",  NF90_FLOAT, dimids_1level, ncvarid(6))
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
  ! SIF: write LEAFC_INC and TLAI_INC
  !
  ! leafc_inc_g(g) and tlai_inc_g(g) are filled in clm_update_sif
  ! (enkf_clm_mod) as (new_gridcell_mean - old_gridcell_mean).
  ! They are indexed by CLM global gridcell index g (begg:endg).
  ! We use mpi_gatherv to collect them on masterproc, then write 2-D.
  !
  ! For gridcells NOT in num_hactiveg (lake, urban, bare), the arrays
  ! are initialised to 0 in define_clm_statevec_sif — correct (no update).
  ! ==================================================================
  if (clmupdate_sif==1) then

      allocate(sif_inc_local(begg:endg), stat=nerror)
      if (masterproc) allocate(sif_inc_global(1:numg), stat=nerror)

      ! --- LEAFC_INC ---
      sif_inc_local(begg:endg) = 0.0_r8
      do count = 1, num_hactiveg
          g = hactiveg_levels(count, 1)
          if (g >= begg .and. g <= endg) &
              sif_inc_local(g) = real(leafc_inc_g(g), r8)
      end do

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
              clmstate_out_2d(ji,jj) = real(sif_inc_global(g1))
          end do
          status = nf90_inq_varid(il_file_id, "LEAFC_INC", ncvarid(5))
          status = nf90_put_var(il_file_id, ncvarid(5), clmstate_out_2d(:,:), &
                                start=[1,1], count=[ndlon,ndlat])
      end if

      ! --- TLAI_INC ---
      sif_inc_local(begg:endg) = 0.0_r8
      do count = 1, num_hactiveg
          g = hactiveg_levels(count, 1)
          if (g >= begg .and. g <= endg) &
              sif_inc_local(g) = real(tlai_inc_g(g), r8)
      end do

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
              clmstate_out_2d(ji,jj) = real(sif_inc_global(g1))
          end do
          status = nf90_inq_varid(il_file_id, "TLAI_INC", ncvarid(6))
          status = nf90_put_var(il_file_id, ncvarid(6), clmstate_out_2d(:,:), &
                                start=[1,1], count=[ndlon,ndlat])
      end if

      deallocate(sif_inc_local)
      if (masterproc) deallocate(sif_inc_global)

  end if   ! clmupdate_sif==1

  ! ------------------------------------------------------------------
  ! Close and deallocate
  ! ------------------------------------------------------------------
  if (masterproc) then
      status = nf90_close(il_file_id)
      deallocate(clmstate_out)
      deallocate(clmstate_tmp_global)
      if (clmupdate_sif==1 .and. allocated(clmstate_out_2d)) &
          deallocate(clmstate_out_2d)
  end if
  deallocate(tws_inc)
  deallocate(clmstate_tmp_local)

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
