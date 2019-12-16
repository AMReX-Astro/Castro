module lapack_module

  implicit none

contains

  subroutine dgtsv(N, NRHS, DL, D, DU, B, LDB, INFO)

    implicit none

    integer :: INFO, LDB, N, NRHS

    double precision :: B(LDB,*), D(*), DL(*), DU(*)

    double precision, parameter :: ZERO = 0.0D+0

    integer :: I, J
    double precision :: FACT, TEMP

    INFO = 0
    if (N .lt. 0) then
       INFO = -1
    else if (NRHS .lt. 0) then
       INFO = -2
    else if (LDB .lt. max(1, N)) then
       INFO = -7
    end if

    if (INFO .ne. 0) then
       call xerbla('DGTSV ', -INFO)
       return
    end if

    if (N .eq. 0) return

    if (NRHS .eq. 1) then

       do I = 1, N - 2

          if (abs(D(I)) .ge. abs(DL(I))) then

             if (D(I) .ne. ZERO) then
                FACT = DL(I) / D(I)
                D(I+1) = D(I+1) - FACT*DU(I)
                B(I+1,1) = B(I+1, 1) - FACT*B(I, 1)
             else
                INFO = I
                return
             end if
             DL(I) = ZERO

          else

             FACT = D(I) / DL(I)
             D(I) = DL(I)
             TEMP = D(I+1)
             D(I+1) = DU(I) - FACT*TEMP
             DL(I) = DU(I+1)
             DU(I+1) = -FACT*DL(I)
             DU(I) = TEMP
             TEMP = B(I, 1)
             B(I, 1) = B(I+1, 1)
             B(I+1, 1) = TEMP - FACT*B(I+1, 1)

          end if

       end do

       if (N .GT. 1) then

          I = N - 1

          if (abs(D(I)).ge.abs(DL(I))) then

             if (D(I) .ne. ZERO) then
                FACT = DL(I) / D(I)
                D(I+1) = D(I+1) - FACT*DU(I)
                B(I+1, 1) = B(I+1, 1) - FACT*B(I, 1)
             else
                INFO = I
                return
             end if

          else

             FACT = D(I) / DL(I)
             D(I) = DL(I)
             TEMP = D(I+1)
             D(I+1) = DU(I) - FACT*TEMP
             DU(I) = TEMP
             TEMP = B(I, 1)
             B(I, 1) = B(I+1, 1)
             B(I+1, 1) = TEMP - FACT*B(I+1, 1)

          end if

       end if

       if (D(N) .eq. ZERO) then
          INFO = N
          return
       end if

    else

       do I = 1, N - 2

          if (abs(D(I)) .ge. abs(DL(I))) then

             if (D(I) .ne. ZERO) then
                FACT = DL(I) / D(I)
                D(I+1) = D(I+1) - FACT*DU(I)
                do J = 1, NRHS
                   B(I+1, J) = B(I+1, J) - FACT*B(I, J)
                end do
             else
                INFO = I
                return
             end if

             DL(I) = ZERO

          else

             FACT = D(I) / DL(I)
             D(I) = DL(I)
             TEMP = D(I+1)
             D(I+1) = DU(I) - FACT*TEMP
             DL(I) = DU(I+1)
             DU(I+1) = -FACT*DL(I)
             DU(I) = TEMP
             do J = 1, NRHS
                TEMP = B(I, J)
                B(I, J) = B(I+1, J)
                B(I+1, J) = TEMP - FACT*B(I+1, J)
             end do

          end if

       end do

       if (N .gt. 1) then

          I = N - 1

          if (abs(D(I)) .ge. abs(DL(I))) then

             if (D(I).ne.ZERO) then
                FACT = DL(I) / D(I)
                D(I+1) = D(I+1) - FACT*DU(I)
                do J = 1, NRHS
                   B(I+1, J) = B(I+1, J) - FACT*B(I, J)
                end do
             else
                INFO = I
                return
             end if

          else

             FACT = D(I) / DL(I)
             D(I) = DL(I)
             TEMP = D(I+1)
             D(I+1) = DU(I) - FACT*TEMP
             DU(I) = TEMP
             do J = 1, NRHS
                TEMP = B(I, J)
                B(I, J) = B(I+1, J)
                B(I+1, J) = TEMP - FACT*B(I+1, J)
             end do

          end if

       end if

       if (D(N) .eq. ZERO) then
          INFO = N
          return
       end if

    end if

    if (NRHS .le. 2) then

       J = 1
70     continue
       B(N, J) = B(N, J) / D(N)
       if (N .gt. 1) then
          B(N-1, J) = (B(N-1, J)-DU(N-1)*B(N, J)) / D(N-1)
       end if
       do I = N - 2, 1, -1
          B(I, J) = (B(I, J)-DU(I)*B(I+1, J)-DL(I)* &
                     B(I+2, J)) / D(I)
       end do
       if (J .lt. NRHS) then
          J = J + 1
          go to 70
       end if

    else

       do J = 1, NRHS
          B(N, J) = B(N, J) / D(N)
          if (N .gt. 1) then
             B(N-1, J) = (B(N-1, J)-DU(N-1)*B(N, J)) / &
                          D(N-1)
          end if
          do I = N - 2, 1, -1
             B(I, J) = (B(I, J)-DU(I)*B(I+1, J)-DL(I)* &
                        B(I+2, J)) / D(I)
          end do
       end do

    end if

  end subroutine dgtsv


  subroutine xerbla(SRNAME, INFO)

    implicit none

    character*(*) :: SRNAME
    integer       :: INFO

    write(*, FMT = 999) SRNAME(1:LEN_TRIM(SRNAME)), INFO

    stop

999 format(' ** On entry to ', A, ' parameter number ', I2, ' had ', &
           'an illegal value')

  end subroutine xerbla

end module lapack_module
