module REGIONGROW
implicit none

public :: POP, APPEND, REGIONGROW_REAL, REGIONGROW_BIN, REGIONGROW_REAL_LIST, REGIONGROW_BIN_LIST


contains

    !-------------Pop value from array end-------------
    subroutine POP(array,val,cur)
        ! Pop value from array end
        ! <array>: 1d array, int, input.
        ! <val>: int, output, last value in array.
        ! <cur>: int, inout, current location, should point to the last.
        ! NOTE: handle the current location outside this sub.
        implicit none
        integer, intent(inout), dimension(:) :: array
        integer, intent(out) :: val
        integer, intent(inout) :: cur

        if (cur==0) then
            return
        end if

        val=array(cur)
        array(cur)=0

    end subroutine pop


    !------------Append value to array end------------
    subroutine APPEND(array,val,cur)
        ! Append value to array end
        ! NOTE: handle the current location outside this sub.
        implicit none
        integer, intent(inout), dimension(:) :: array
        integer, intent(in) :: val
        integer, intent(inout) :: cur
        integer :: nx

        nx=size(array)
        if (cur==nx) then
            return
        end if

        array(cur)=val

    end subroutine append



    subroutine REGIONGROW_REAL(data,seedy,seedx,eps,res,ny,nx,diagonal)
    ! Region grow search of real data from a given seed and select data in range.

        ! <data>: 2D real, input, input data to search.
        ! <seedy>, <seedx>: int, input, y and x coordinates of seed point where
        !                   search starts.
        ! <eps>: real, input, positive, data in the range of [v - eps, v + eps]
        !        will be selected, where v is the value at seed point.
        ! <res>: 2D int, output, regiongrow search result, same shape as <data>, 
        !        with selected points returned as 1s, and others filled with 0s.
        ! <ny>, <nx>: int, input, shape of <data>.
        ! <diagonal>: int, input, optional, if 1, include diagnoal points as neighbours,
        !             exclude if 0.

        ! Author: guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
        ! Update time: 2017-04-14 09:52:32.

        implicit none
        integer, parameter :: ikind = selected_real_kind(p=15,r=10)

        integer :: ny,nx
        real(kind=ikind), intent(in), dimension(ny,nx) :: data
        integer, intent (in) :: seedy, seedx
        integer, intent(out), dimension(ny,nx) :: res
        integer, intent(inout), optional :: diagonal
        real(kind=ikind), intent(in) :: eps

        integer, dimension(ny*nx) :: queuey
        integer, dimension(ny*nx) :: queuex
        integer, dimension(ny,nx) :: visited
        real(kind=ikind) :: val
        integer :: ii,jj,queuelen,current,cy,cx

        !-------------------Check inputs-------------------
        if (seedy > ny .OR. seedy<=0) then
            write(*,*) '<seedy> out of range.'
            res=0.
            return
        end if

        if (seedx > nx .OR. seedx<=0) then
            write(*,*) '<seedx> out of range.'
            res=0.
            return
        end if

        if (eps<=0) then
            write(*,*) '<eps> needs to be positve.'
            return
        end if

        !----------------------Setup----------------------
        res=0
        visited=0
        queuelen=1
        current=1
        queuey(current)=seedy
        queuex(current)=seedx
        val=data(seedy,seedx)
        visited(seedy,seedx)=1
        res(seedy,seedx)=1


        if (present(diagonal) .eqv. .FALSE.) then
            diagonal=1
        end if

        !----------------Loop through grids----------------
        do while (.TRUE.)
            do ii = -1,1
                do jj = -1,1
                    if (diagonal==0 .AND. ii*jj/=0) then
                        cycle
                    end if
                    cy=queuey(current)+ii
                    cx=queuex(current)+jj

                    if (cy<=0 .OR. cy>ny .OR. cx<=0 .OR. cx>nx) then
                        cycle
                    end if

                    if (visited(cy,cx)==0 .AND. abs(data(cy,cx)-val)<=eps) then
                        res(cy,cx)=1
                        visited(cy,cx)=1
                        queuelen=queuelen+1
                        queuey(queuelen)=cy
                        queuex(queuelen)=cx
                    end if
                end do
            end do

            if (queuelen==current) then
                exit
            end if
            current=current+1
        end do

    end subroutine REGIONGROW_REAL


    subroutine REGIONGROW_BIN(data,seedy,seedx,res,ny,nx,diagonal)
    ! Region grow search of binary data from a given seed.

        ! <data>: 2D int, input, input data to search.
        ! <seedy>, <seedx>: int, input, y and x coordinates of seed point where
        !                   search starts.
        ! <res>: 2D int, output, regiongrow search result, same shape as <data>, 
        !        with selected points returned as 1s, and others filled with 0s.
        ! <ny>, <nx>: int, input, shape of <data>.
        ! <diagonal>: int, input, optional, if 1, include diagnoal points as neighbours,
        !             exclude if 0.

        ! Author: guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
        ! Update time: 2017-04-14 09:52:32.

        implicit none
        integer :: ny,nx
        integer, intent(in), dimension(ny,nx) :: data
        integer, intent (in) :: seedy, seedx
        integer, intent(out), dimension(ny,nx) :: res
        integer, intent(inout), optional :: diagonal

        integer, dimension(ny*nx) :: queuey
        integer, dimension(ny*nx) :: queuex
        integer :: ii,jj,queuelen,current,cy,cx
        integer :: val

        !-------------------Check inputs-------------------
        if (seedy > ny .OR. seedy<=0) then
            write(*,*) '<seedy> out of range.'
            res=0.
            return
        end if

        if (seedx > nx .OR. seedx<=0) then
            write(*,*) '<seedx> out of range.'
            res=0.
            return
        end if

        !----------------Loop through grids----------------
        res=0
        queuelen=1
        current=1
        queuey(current)=seedy
        queuex(current)=seedx
        val=data(seedy,seedx)
        res(seedy,seedx)=1

        if (present(diagonal) .eqv. .FALSE.) then
            diagonal=1
        end if

        do while (.TRUE.)
            do ii = -1,1
                do jj = -1,1
                    if (diagonal==0 .AND. ii*jj/=0) then
                        cycle
                    end if
                    cy=queuey(current)+ii
                    cx=queuex(current)+jj

                    if (cy<=0 .OR. cy>ny .OR. cx<=0 .OR. cx>nx) then
                        cycle
                    end if

                    if (res(cy,cx)==0 .AND. data(cy,cx)==val) then
                        res(cy,cx)=1
                        queuelen=queuelen+1
                        queuey(queuelen)=cy
                        queuex(queuelen)=cx
                    end if
                end do
            end do

            if (queuelen==current) then
                exit
            end if
            current=current+1
        end do

    end subroutine REGIONGROW_BIN





    subroutine REGIONGROW_REAL_LIST(data,seedy,seedx,eps,ry,rx,nr,ny,nx,diagonal)
    ! Region grow search of real data from a given seed and select data in range.
    ! Similar to REGIONGROW_REAL but return result in (y,x) where y and x
    ! are 1d coordinates of selected points.

        ! <data>: 2D real, input, input data to search.
        ! <seedy>, <seedx>: int, input, y and x coordinates of seed point where
        !                   search starts.
        ! <eps>: real, input, positive, data in the range of [v - eps, v + eps]
        !        will be selected, where v is the value at seed point.
        ! <ry>, <rx>: 1d int, y(ny) and x(nx) coordinates of selected points (including seed).
        ! <nr>: int, output, number of valid points in <ry> and <rx>.
        ! <ny>, <nx>: int, input, shape of <data>.
        ! <diagonal>: int, input, optional, if 1, include diagnoal points as neighbours,
        !             exclude if 0.

        ! Author: guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
        ! Update time: 2017-04-14 09:52:32.

        implicit none
        integer, parameter :: ikind = selected_real_kind(p=15,r=10)

        integer :: ny,nx
        real(kind=ikind), intent(in), dimension(ny,nx) :: data
        integer, intent (in) :: seedy, seedx
        integer, intent(out), dimension(ny*nx) :: ry,rx
        integer, intent(out) :: nr
        integer, intent(inout), optional :: diagonal
        real(kind=ikind), intent(in) :: eps

        integer, dimension(ny*nx) :: queuey
        integer, dimension(ny*nx) :: queuex
        integer, dimension(ny,nx) :: visited
        real(kind=ikind) :: val
        integer :: ii,jj,queuelen,cy,cx,neiy,neix

        !-------------------Check inputs-------------------
        if (seedy > ny .OR. seedy<=0) then
            write(*,*) '<seedy> out of range.'
            nr=0
            return
        end if

        if (seedx > nx .OR. seedx<=0) then
            write(*,*) '<seedx> out of range.'
            nr=0
            return
        end if

        if (eps<=0) then
            write(*,*) '<eps> needs to be positve.'
            nr=0
            return
        end if

        !----------------------Setup----------------------
        visited=0
        queuelen=1
        nr=1
        val=data(seedy,seedx)
        visited(seedy,seedx)=1
        call APPEND(queuey,seedy,queuelen)
        call APPEND(queuex,seedx,queuelen)
        call APPEND(ry,seedy,nr)
        call APPEND(rx,seedx,nr)

        if (present(diagonal) .eqv. .FALSE.) then
            diagonal=1
        end if

        !----------------Loop through grids----------------
        do while (.TRUE.)
            if (queuelen==0) then
                exit
            end if
            call POP(queuey,cy,queuelen)
            call POP(queuex,cx,queuelen)
            queuelen=queuelen-1

            !-------------Loop through neighbours-------------
            do ii = -1,1
                do jj = -1,1
                    if (diagonal==0 .AND. ii*jj/=0) then
                        cycle
                    end if

                    neiy=cy+ii
                    neix=cx+jj

                    if (neiy<=0 .OR. neiy>ny .OR. neix<=0 .OR. neix>nx) then
                        cycle
                    end if

                    if (visited(neiy,neix)==0 .AND. abs(data(neiy,neix)-val)<=eps) then
                        visited(neiy,neix)=1
                        queuelen=queuelen+1
                        call APPEND(queuey,neiy,queuelen)
                        call APPEND(queuex,neix,queuelen)
                        nr=nr+1
                        call APPEND(ry,neiy,nr)
                        call APPEND(rx,neix,nr)
                    end if
                end do
            end do
        end do


    end subroutine REGIONGROW_REAL_LIST






    subroutine REGIONGROW_BIN_LIST(data,seedy,seedx,ry,rx,nr,ny,nx,diagonal)
    ! Region grow search of binary data from a given seed.
    ! Similar to REGIONGROW_BIN but return result in (y,x) where y and x
    ! are 1d coordinates of selected points.

        ! <data>: 2D int, input, input data to search.
        ! <seedy>, <seedx>: int, input, y and x coordinates of seed point where
        !                   search starts.
        ! <ry>, <rx>: 1d int, y(ny) and x(nx) coordinates of selected points (including seed).
        ! <nr>: int, output, number of valid points in <ry> and <rx>.
        ! <ny>, <nx>: int, input, shape of <data>.
        ! <diagonal>: int, input, optional, if 1, include diagnoal points as neighbours,
        !             exclude if 0.

        ! Author: guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
        ! Update time: 2017-04-14 09:52:32.


        implicit none
        integer :: ny,nx
        integer, intent(in), dimension(ny,nx) :: data
        integer, intent (in) :: seedy, seedx
        integer, intent(out), dimension(ny*nx) :: ry,rx
        integer, intent(out) :: nr
        integer, intent(inout), optional :: diagonal

        integer, dimension(ny*nx) :: queuey, queuex
        integer, dimension(ny,nx) :: visited
        integer :: val
        integer :: ii,jj,queuelen,cy,cx,neiy,neix

        !-------------------Check inputs-------------------
        if (seedy > ny .OR. seedy<=0) then
            write(*,*) '<seedy> out of range.'
            nr=0
            return
        end if

        if (seedx > nx .OR. seedx<=0) then
            write(*,*) '<seedx> out of range.'
            nr=0
            return
        end if

        !----------------------Setup----------------------
        visited=0
        queuelen=1
        nr=1
        val=data(seedy,seedx)
        visited(seedy,seedx)=1
        call APPEND(queuey,seedy,queuelen)
        call APPEND(queuex,seedx,queuelen)
        call APPEND(ry,seedy,nr)
        call APPEND(rx,seedx,nr)

        if (present(diagonal) .eqv. .FALSE.) then
            diagonal=1
        end if

        !----------------Loop through grids----------------
        do while (.TRUE.)
            if (queuelen==0) then
                exit
            end if
            call POP(queuey,cy,queuelen)
            call POP(queuex,cx,queuelen)
            queuelen=queuelen-1

            !-------------Loop through neighbours-------------
            do ii = -1,1
                do jj = -1,1
                    if (diagonal==0 .AND. ii*jj/=0) then
                        cycle
                    end if

                    neiy=cy+ii
                    neix=cx+jj

                    if (neiy<=0 .OR. neiy>ny .OR. neix<=0 .OR. neix>nx) then
                        cycle
                    end if

                    if (visited(neiy,neix)==0 .AND. data(neiy,neix)==val) then
                        visited(neiy,neix)=1
                        queuelen=queuelen+1
                        call APPEND(queuey,neiy,queuelen)
                        call APPEND(queuex,neix,queuelen)
                        nr=nr+1
                        call APPEND(ry,neiy,nr)
                        call APPEND(rx,neix,nr)
                    end if
                end do
            end do
        end do


    end subroutine REGIONGROW_BIN_LIST

end module
