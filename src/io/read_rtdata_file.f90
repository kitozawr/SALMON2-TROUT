!
!  Copyright 2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!=======================================================================
module read_rtdata_file
  use filesystem, only: open_filehandle
  use inputoutput, only: t_unit_ac,t_unit_time
  implicit none

contains

!===================================================================================================================================

  function count_rows_from_rtdata_file(filename) result(icount)
    ! This function returns the number of data rows
    ! in the `?????_rt.data` formatted file.
    ! Blank lines and comment lines are skipped.
    !
    ! Robustness: a missing / unreadable / empty field file must ABORT with an
    ! actionable message, never fall through with icount=0 (which lets the caller
    ! allocate zero-length arrays and index them out of bounds). The open below
    ! carries an explicit iostat so a file that is not yet visible on this node
    ! -- the classic distributed-filesystem "input read" failure -- is reported
    ! rather than crashing the Fortran runtime.
    implicit none
    character(*), intent(in) :: filename
    character(1024) :: buf
    integer :: icount
    integer ::  flag, fh, ios
    logical :: exists

    inquire(file=trim(filename), exist=exists)
    if (.not. exists) then
        write(*, '(a)') '# ERROR: external-field file not found: '//trim(filename)
        error stop 'count_rows_from_rtdata_file: input field file missing'
    end if

    open(newunit=fh, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*, '(a,i0)') '# ERROR: cannot open external-field file '// &
                           trim(filename)//', iostat=', ios
        error stop 'count_rows_from_rtdata_file: cannot open input field file'
    end if
    icount = 0
    do while (.true.)
        read(fh, '(a)', iostat=flag)  buf
        if (flag == 0) then
            buf = adjustl(buf)
            if (len_trim(buf) < 1) cycle
            if (buf(1:1) == '#' .or. buf(1:1) == '!') cycle
            icount = icount + 1
        else
            exit
        end if
    end do
    close(fh)
    if (icount < 2) then
        write(*, '(a,i0,a)') '# ERROR: external-field file '//trim(filename)// &
                             ' has ', icount, ' data row(s); need >= 2 to interpolate.'
        error stop 'count_rows_from_rtdata_file: input field file too short / empty'
    end if
    return
  end function count_rows_from_rtdata_file

!===================================================================================================================================

  subroutine load_Ac_from_rtdata_file(filename, n_dat, t_dat, Ac_dat)
    ! This subroutine Reads the data of the given number of samples 
    ! from `???_rt.data` file. The empty and comment lines are skipped,
    ! and the values of the first four columns:
    ! time(1) and Ac_ext (2-4) are stored in t_dat, Ac_dat.
    implicit none
    character(*), intent(in) :: filename
    integer, intent(in) :: n_dat
    real(8), intent(out) :: t_dat(n_dat)
    real(8), intent(out) :: Ac_dat(1:3, n_dat)

    integer :: i, fh, ios
    character(1024) :: buf

    open(newunit=fh, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*, '(a,i0)') '# ERROR: cannot open external-field file '// &
                           trim(filename)//', iostat=', ios
        error stop 'load_Ac_from_rtdata_file: cannot open input field file'
    end if
    i = 1
    do while (i <= n_dat)
        read(fh, '(a)', iostat=ios)  buf
        if (ios /= 0) then
            ! Fewer valid rows than count_rows reported: a truncated / partially
            ! staged file (the two passes disagreeing is itself a red flag).
            write(*, '(a,i0,a,i0,a)') '# ERROR: external-field file '//trim(filename)// &
                ' ended after ', i-1, ' of ', n_dat, ' expected rows.'
            error stop 'load_Ac_from_rtdata_file: unexpected end of input field file'
        end if
        buf = adjustl(buf)
        if (len_trim(buf) < 1) cycle
        if (buf(1:1) == '#' .or. buf(1:1) == '!') cycle
        read(buf, *, iostat=ios) t_dat(i), Ac_dat(1, i), Ac_dat(2, i), Ac_dat(3, i)
        if (ios /= 0) then
            write(*, '(a)') '# ERROR: malformed row in external-field file '// &
                            trim(filename)//': '//trim(buf)
            error stop 'load_Ac_from_rtdata_file: malformed data row'
        end if
        i = i + 1
    end do
    close(fh)

    ! Transform unit of Ac field
    Ac_dat(:, :) = Ac_dat(:, :) / t_unit_ac%conv
    t_dat(:) = t_dat(:) / t_unit_time%conv

    return
  end subroutine load_Ac_from_rtdata_file

  ! if you want use this subroutine, please switch subroutine name in
  ! rt/em_filed.f90
  subroutine load_Ac_from_rtdata_file_Ac_tot(filename, n_dat, t_dat, Ac_dat)
    ! This subroutine Reads the data of the given number of samples 
    ! from `???_rt.data` file. The empty and comment lines are skipped,
    ! and the values of the four columns:
    ! time(1) and Ac_tot (8-10) are stored in t_dat, Ac_dat.
    implicit none
    character(*), intent(in) :: filename
    integer, intent(in) :: n_dat
    real(8), intent(out) :: t_dat(n_dat)
    real(8), intent(out) :: Ac_dat(1:3, n_dat)
    real(8) :: tmp(10)

    integer :: i, fh
    character(1024) :: buf

    fh = open_filehandle(trim(filename), status='old')
    i = 1
    do while (i <= n_dat)
        read(fh, '(a)')  buf
        buf = adjustl(buf)
        if (len_trim(buf) < 1) cycle
        if (buf(1:1) == '#' .or. buf(1:1) == '!') cycle
        read(buf, *) t_dat(i), tmp(1:6), Ac_dat(1:3, i)  !Ac_tot
        i = i + 1
    end do
    close(fh)

    ! Transform unit of Ac field
    Ac_dat(:, :) = Ac_dat(:, :) / t_unit_ac%conv
    t_dat(:) = t_dat(:) / t_unit_time%conv
    
    return
  end subroutine load_Ac_from_rtdata_file_Ac_tot

!===================================================================================================================================

end module read_rtdata_file
