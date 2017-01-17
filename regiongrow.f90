module REGIONGROW
implicit none

public :: REGIONGROW_REAL, REGIONGROW_BIN


contains

    subroutine REGIONGROW_REAL(data,seedy,seedx,val,eps,res,ny,nx,diagonal)
        ! Region grow from a given seed and return values within range.

        ! <data>: 2d real array. Data to search-select.
        ! <seedy>: int, y coordinate of seed.
        ! <seedx>: int, x coordinate of seed.
        ! <val>: real, target value to match within <data>.
        ! <eps>: real, positive, tolerance range. Values in <data> in the range 
        !        [<val> - <eps>, <val> + <eps>] are regarded matching.
        ! <ny>, <nx>: int, optional, size of <data>. 
        ! <diagonal>: int, optional, if 1, diagonal elements are treated as neighbours, 0 otherwise.
        !
        ! Return <res>: 2d int array, 1s for matching in <data>, 0s otherwise. 


        implicit none
        integer :: ny,nx
        real, intent(in), dimension(ny,nx) :: data
        integer, intent (in) :: seedy, seedx
        real, intent(in) :: val
        integer, intent(out), dimension(ny,nx) :: res
        integer, intent(inout), optional :: diagonal
        real, intent(in) :: eps

        integer, dimension(ny*nx) :: queuey
        integer, dimension(ny*nx) :: queuex
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
        queuelen=1
        current=1
        queuey(current)=seedy
        queuex(current)=seedx

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

                    if (res(cy,cx)==0 .AND. abs(data(cy,cx)-val)<=eps) then
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

    end subroutine REGIONGROW_REAL


    subroutine REGIONGROW_BIN(data,seedy,seedx,val,res,ny,nx,diagonal)
        ! Binary region grow from a given seed.

        ! <data>: 2d int array. Data to search-select.
        ! <seedy>: int, y coordinate of seed.
        ! <seedx>: int, x coordinate of seed.
        ! <val>: int, target value to match within <data>.
        ! <ny>, <nx>: int, optional, size of <data>. 
        ! <diagonal>: int, optional, if 1, diagonal elements are treated as neighbours, 0 otherwise.
        !
        ! Return <res>: 2d int array, with 1s for <data>==<val> and 0s elsewhere. 

        implicit none
        integer :: ny,nx
        integer, intent(in), dimension(ny,nx) :: data
        integer, intent (in) :: seedy, seedx
        integer, intent(in) :: val
        integer, intent(out), dimension(ny,nx) :: res
        integer, intent(inout), optional :: diagonal

        integer, dimension(ny*nx) :: queuey
        integer, dimension(ny*nx) :: queuex
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

        !----------------Loop through grids----------------
        res=0
        queuelen=1
        current=1
        queuey(current)=seedy
        queuex(current)=seedx

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

end module
