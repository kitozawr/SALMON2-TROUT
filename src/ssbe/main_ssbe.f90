subroutine main_ssbe(icomm)
    ! use mpi
    use omp_lib
    use communication
    use multiscale_ssbe
    use realtime_ssbe
    use sfsb_ssbe
    use salmon_global
    implicit none
    integer, intent(in) :: icomm

    select case(trim(theory))
    case ("sbe")
        if (yn_sbe_sfsb == 'y') then
            ! SFSB non-Markovian memory-integral ionization mode [B25]
            call main_sfsb_ssbe(icomm)
        else
            call main_realtime_ssbe(icomm)
        end if
    case ("maxwell_sbe")
        call main_multiscale_ssbe(icomm)
    end select

    call comm_sync_all(icomm)

    return
end subroutine main_ssbe 
