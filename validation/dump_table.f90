! Reference table reader: reproduces the read sequence of src/read_sf_file.f90
! verbatim, so we can see what the Fortran actually extracts and whether the
! irregular M_{L-} block desynchronises the reader.
!
! usage: dump_table <table.tbl> <iq2> <iw>
!
! Prints "SF<k> <value>" for the 62 structure functions at grid point
! (iq2, iw), 1-based, exactly as the generator would see them.
program dump_table
    implicit none

    integer, parameter :: nvar1 = 101, nvar2 = 93
    real var1(nvar1), var2(nvar2)
    real sf(62, nvar1, nvar2)
    real dumvar1, dumvar2, dumvar3, dumvar4
    character*80 dummy
    character*2 dummyX
    character*5 dummyY
    real var1tmp1, var1tmp2, var2tmp1, var2tmp2
    integer jvar1, jvar2, iu, iq2, iw, i
    character*512 arg

    iu = 11
    call get_command_argument(1, arg)
    if (len_trim(arg) == 0) then
        write(*,'(a)') 'usage: dump_table <table.tbl> <iq2> <iw>'
        stop 1
    end if
    open(unit=iu, file=trim(arg), status='old')
    call get_command_argument(2, arg)
    read(arg,*) iq2
    call get_command_argument(3, arg)
    read(arg,*) iw

    do jvar1 = 1, nvar1
        do jvar2 = 1, nvar2
            read(iu, FMT = 15, err = 1000) dummyX, var2(jvar2), &
                    dummyY, var1(jvar1)
15          format(A8, f4.2, A7, f7.5)
            if (jvar1 .eq. 1)     var1tmp1 = var1(jvar1)
            if (jvar1 .eq. nvar1) var1tmp2 = var1(jvar1)
            if (jvar2 .eq. 1)     var2tmp1 = var2(jvar2)
            if (jvar2 .eq. nvar2) var2tmp2 = var2(jvar2)

            ! SL+
            read(iu, '(a)', err = 1000) dummy
            read(iu, *, err = 1000) sf(1,jvar1,jvar2), sf(2,jvar1,jvar2), &
                    sf(3,jvar1,jvar2), sf(4,jvar1,jvar2), &
                    sf(5,jvar1,jvar2), sf(6,jvar1,jvar2)
            read(iu, *, err = 1000) sf(7,jvar1,jvar2), sf(8,jvar1,jvar2), &
                    sf(9,jvar1,jvar2), sf(10,jvar1,jvar2), &
                    sf(11,jvar1,jvar2), sf(12,jvar1,jvar2)
            ! SL-
            read(iu, '(a)', err = 1000) dummy
            read(iu, *, err = 1000) dumvar1, dumvar2, &
                    sf(13,jvar1,jvar2), sf(14,jvar1,jvar2), &
                    sf(15,jvar1,jvar2), sf(16,jvar1,jvar2)
            read(iu, *, err = 1000) sf(17,jvar1,jvar2), sf(18,jvar1,jvar2), &
                    sf(19,jvar1,jvar2), sf(20,jvar1,jvar2), &
                    sf(21,jvar1,jvar2), sf(22,jvar1,jvar2)
            ! EL+
            read(iu, '(a)', err = 1000) dummy
            read(iu, *, err = 1000) sf(23,jvar1,jvar2), sf(24,jvar1,jvar2), &
                    sf(25,jvar1,jvar2), sf(26,jvar1,jvar2), &
                    sf(27,jvar1,jvar2), sf(28,jvar1,jvar2)
            read(iu, *, err = 1000) sf(29,jvar1,jvar2), sf(30,jvar1,jvar2), &
                    sf(31,jvar1,jvar2), sf(32,jvar1,jvar2), &
                    sf(33,jvar1,jvar2), sf(34,jvar1,jvar2)
            ! EL-
            read(iu, '(a)', err = 1000) dummy
            read(iu, *, err = 1000) dumvar1, dumvar2, &
                    dumvar3, dumvar4, &
                    sf(35,jvar1,jvar2), sf(36,jvar1,jvar2)
            read(iu, *, err = 1000) sf(37,jvar1,jvar2), sf(38,jvar1,jvar2), &
                    sf(39,jvar1,jvar2), sf(40,jvar1,jvar2), &
                    sf(41,jvar1,jvar2), sf(42,jvar1,jvar2)
            ! ML+
            read(iu, '(a)', err = 1000) dummy
            read(iu, *, err = 1000) dumvar1, dumvar2, &
                    sf(43,jvar1,jvar2), sf(44,jvar1,jvar2), &
                    sf(45,jvar1,jvar2), sf(46,jvar1,jvar2)
            read(iu, *, err = 1000) sf(47,jvar1,jvar2), sf(48,jvar1,jvar2), &
                    sf(49,jvar1,jvar2), sf(50,jvar1,jvar2), &
                    sf(51,jvar1,jvar2), sf(52,jvar1,jvar2)
            ! ML-  (CORRECTED: 2 padding columns, not 4 -- see notes below)
            read(iu, '(a)', err = 1000) dummy
            read(iu, *, err = 1000) dumvar1, dumvar2, &
                    sf(53,jvar1,jvar2), sf(54,jvar1,jvar2), &
                    sf(55,jvar1,jvar2), sf(56,jvar1,jvar2)
            read(iu, *, err = 1000) sf(57,jvar1,jvar2), sf(58,jvar1,jvar2), &
                    sf(59,jvar1,jvar2), sf(60,jvar1,jvar2), &
                    sf(61,jvar1,jvar2), sf(62,jvar1,jvar2)
        enddo
    enddo
    close(iu)

    write(*,'(a)')    'NVAR1_NVAR2'
    write(*,'(2i5)')  nvar1, nvar2
    write(*,'(a,2f9.4)') 'VAR1_MIN_MAX ', var1tmp1, var1tmp2
    write(*,'(a,2f9.4)') 'VAR2_MIN_MAX ', var2tmp1, var2tmp2
    write(*,'(a,2f10.5)') 'AT_Q2_W    ', var1(iq2), var2(iw)
    do i = 1, 62
        write(*,'(a,i3,1x,es24.16)') 'SF', i, sf(i, iq2, iw)
    enddo
    write(*,'(a)') 'OK'
    stop

1000 write(*,'(a,2i6)') 'ERROR_READING jvar1 jvar2', jvar1, jvar2
    stop 2
end
