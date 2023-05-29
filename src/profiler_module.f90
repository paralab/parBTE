MODULE profiler
    USE omp_lib
    USE PRECISIONS
    IMPLICIT NONE
    
    !INTEGER, PARAMETER :: int_p  = SELECTED_INT_KIND(8)     !...int_p = kind of integer giving precision to 8 digits
    !INTEGER, PARAMETER :: real_p = SELECTED_REAL_KIND(8)    !...real_p = kind of real giving precision to 8 digits after decimal

    TYPE :: profiler_t
        REAL(real_p) :: seconds
        INTEGER(int_p) :: num_calls
        REAL(real_p) :: pri_seconds
    END TYPE profiler_t

    CONTAINS
    SUBROUTINE profiler_init(p)
        IMPLICIT NONE
        TYPE(profiler_t), INTENT(OUT) :: p

        p%seconds = 0.0
        p%num_calls = 0
        p%pri_seconds = 0.0
    END SUBROUTINE profiler_init

    SUBROUTINE profiler_start(p)
        IMPLICIT NONE
        TYPE(profiler_t), INTENT(OUT) :: p

        p%pri_seconds = omp_get_wtime()
    END SUBROUTINE profiler_start

    SUBROUTINE profiler_stop(p)
        IMPLICIT NONE
        TYPE(profiler_t), INTENT(INOUT) :: p

        !...accumulate to seconds = seconds + (wallclock_at_stop - wallclock_at_start)
        p%seconds = p%seconds - p%pri_seconds
        p%pri_seconds = omp_get_wtime()
        p%seconds = p%seconds + p%pri_seconds

        !...accumulate number of stop calls
        p%num_calls = p%num_calls + 1
    END SUBROUTINE profiler_stop

END MODULE profiler