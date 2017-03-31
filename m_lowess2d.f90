module m_lowess2d
use m_mrgrnk
use m_qr_solve

! 2D lowess filtering.
!
! EXTERNALS: 
!   - mrgrnk: for array sorting and argsorting
!   - qr_solve: solve general linear least squares using svd.
!
! Compile:
!   - Fortran:
!     gfortran -o m_lowess2d.exe mrgrnk.f90 m_qr_solve.f90 m_lowess2d.f90
!   - f2py:
!     f2py -m m_lowess2d -h m_lowess2d.pyf m_lowess2d.f90
!     f2py -c m_lowess2d.pyf mrgrnk.f90 m_qr_solve.f90 m_lowess2d.f90
!
!
! Author: guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
! Update time: 2017-03-31 16:50:43.

implicit none
private
public LOWESS2D, MEDIAN, POLYFIT2D, REGIONCOPY, GET_ELLIPSE

contains

    subroutine LOWESS2D(z,x,y,nx,ny,span_x,span_y,degree,nsteps,zout,wout)
    ! Compute 2D lowess filtering
    !
    ! <z>: input, real, 2d gridded data to filter.
    ! <x>: input, real, 1d array, x coordinates of <z>.
    ! <y>: input, real, 1d array, y coordinates of <z>.
    ! <nx>, <ny>: input, int, length of <x> and <y>, or number of rows, columns
    !             in <z>, respectively.
    ! <span_x>, <span_y>: input, int, major axes of an ellipse in the x, y
    !                     direction to define the local subset. E.g.
    !                     <span_x>=10, <span_y>=6, an ellipse with HALF axes of
    !                     (5,3) is created as the local subset.
    ! <degree>: input, int, degree of polynomial to fit local subset.
    ! <nsteps>: input, int, number of iterations to robust-fit.
    ! 
    ! <zout>: output, real, 2d filtered data.
    ! <wout>: output, real, 2d local weights.
    !
    ! Author: guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
    ! Update time: 2017-03-31 16:27:04.

        implicit none
        integer, parameter :: ikind = selected_real_kind(p=15,r=10)

        real(kind=ikind), intent(in), dimension(ny,nx) :: z
        real(kind=ikind), intent(in), dimension(nx) :: x
        real(kind=ikind), intent(in), dimension(ny) :: y
        integer, intent(in) :: nx, ny, span_x, span_y, degree, nsteps
        real(kind=ikind), intent(out), dimension(ny,nx) :: zout, wout
        
        !----------------Intermediate vars----------------
        real(kind=ikind), dimension(ny,nx) :: xmesh, ymesh
        real(kind=ikind), dimension(:), allocatable :: xbox, ybox, zbox, &
            & dist, distWeights, zfit, aerr, uu, biWeights, totWeights
        real(kind=ikind) :: xj,yj,mad
        integer, dimension(:), allocatable :: sortidx, bad, badold
        integer, dimension(ny,nx) :: mask
        integer, dimension(:,:), allocatable :: disk
        integer :: a,b,i,j,nbox,p

        !----------------Get mesh and mask----------------
        xmesh=spread(x,1,ny)
        ymesh=spread(y,2,nx)

        a=ceiling(span_x/2.)
        b=ceiling(span_y/2.)
        
        allocate(disk(2*b+1,2*a+1))
        call GET_ELLIPSE(a,b,disk)

        mask=0
        zout=0.
        wout=0.

        !----------------Loop through grids----------------
        do j = 1,ny
            do i = 1,nx
                xj=x(i)
                yj=y(j)

                !-----------------Get local subset-----------------
                call REGIONCOPY(mask,nx,ny,i,j,disk,2*a+1,2*b+1)

                nbox=sum(mask)
                allocate(xbox(nbox))
                allocate(ybox(nbox))
                allocate(zbox(nbox))
                allocate(dist(nbox))
                allocate(sortidx(nbox))
                allocate(distWeights(nbox))
                allocate(zfit(nbox))
                allocate(uu(nbox))
                allocate(biWeights(nbox))
                allocate(bad(nbox))
                allocate(badold(nbox))

                xbox=pack(xmesh,mask==1)
                ybox=pack(ymesh,mask==1)
                zbox=pack(z,mask==1)

                !-----------------Sort by distance-----------------
                dist=sqrt((xbox-xj)**2 + (ybox-yj)**2)
                call mrgrnk(dist,sortidx)
                dist=dist(sortidx)
                xbox=xbox(sortidx)
                ybox=ybox(sortidx)
                zbox=zbox(sortidx)

                !--------Tricube function distance weights--------
                distWeights=(1-(dist/dist(nbox))**3)**3

                !-------------------Initial fit-------------------
                call POLYFIT2D(xbox,ybox,zbox,degree,distWeights,nbox,zfit)

                !--------------------Robust fit--------------------
                do p = 1,nsteps
                    aerr=abs(zfit-zbox)
                    mad=MEDIAN(aerr)
                    uu=(aerr/6./mad)**2
                    where (uu<0)
                        uu=0
                    elsewhere (uu>1)
                        uu=1
                    end where

                    biWeights=(1-uu)**2
                    totWeights=distWeights*biWeights
                    call POLYFIT2D(xbox,ybox,zbox,degree,totWeights,nbox,zfit)

                    if (p==1) then
                        bad=merge(1,0,biWeights<0.34)
                    else
                        badold=bad
                        bad=merge(1,0,biWeights<0.34)
                        if (any(abs(badold-bad)>=0.1) .eqv. .FALSE.) then
                            exit
                        end if
                    end if
                end do

                !----------Save result and clear tmp vars----------
                zout(j,i)=zfit(1)
                wout(j,i)=biWeights(1)
                mask=mask-mask

                deallocate(xbox)
                deallocate(ybox)
                deallocate(zbox)
                deallocate(dist)
                deallocate(sortidx)
                deallocate(distWeights)
                deallocate(zfit)
                deallocate(uu)
                deallocate(biWeights)
                deallocate(bad)
                deallocate(badold)

            end do
        end do

    end subroutine LOWESS2D




    subroutine GET_ELLIPSE(a,b,disk)
    ! Get an ellipse element
    !
    ! <a>, <b>: input, int, half axis length in the x and y directions.
    !
    ! <disk>: output, int, 2d binary array with shape (2*b+1, 2*a+1), filled
    !         with 1s within the ellipse defined by major/minor half axes 
    !         of <a> and <b> and 0s elsewhere. If <a> > <b>, major axis in
    !         the x axis and vice verse.
    !         
        implicit none
        integer, parameter :: ikind = selected_real_kind(p=15,r=10)
        integer, intent(in) :: a, b
        integer, intent(out), dimension(2*b+1,2*a+1) :: disk

        !----------------Intermediate vars----------------
        real(kind=ikind), dimension(2*b+1,2*a+1) :: xmesh, ymesh, ell
        real(kind=ikind), dimension(2*b+1) :: bx
        real(kind=ikind), dimension(2*a+1) :: ax
        real(kind=ikind) :: c
        integer :: i

        ax=[ (i, i=-a, a, 1) ] 
        bx=[ (i, i=-b, b, 1) ] 
        xmesh=spread(ax,1,2*b+1)
        ymesh=spread(bx,2,2*a+1)

        if (a >= b) then
            c=sqrt(float(a**2-b**2))
            ell=sqrt((xmesh-c)**2+ymesh**2) + sqrt((xmesh+c)**2+ymesh**2)
            disk=merge(1,0,ell<=2*a)
        else
            c=sqrt(float(b**2-a**2))
            ell=sqrt(xmesh**2+(ymesh-c)**2) + sqrt(xmesh**2+(ymesh+c)**2)
            disk=merge(1,0,ell<=2*b)
        end if

    end subroutine GET_ELLIPSE




    subroutine REGIONCOPY(slab,nx,ny,xidx,yidx,element,w,h)
    ! Copy a rectangular region to another matrix
    !
    ! <slab>: inout, int, 2d matrix to copy onto.
    ! <nx>, <ny>: input, int, shape of <slab>.
    ! <xidx>, <yidx>: input, int, x and y coordinate in <slab> to paste
    !                 <element>, centered by the center of <element>.
    ! <element>: input, int, 2d array to paste onto <slab>.
    ! <w>, <h>: input, int, shape of <element>.

        implicit none
        integer, intent(inout), dimension(ny,nx) :: slab
        integer, intent(in), dimension(h,w) :: element
        integer, intent(in) :: nx,ny,xidx,yidx,w,h

        !----------------Intermediate vars----------------
        integer :: hh1,yend,hw1,xend,ex1,ex2,ey1,ey2

        if (mod(h,2)==1) then
            hh1=(h-1)/2
            yend=yidx+hh1
        else
            hh1=h/2-1
            yend=yidx+hh1+1
        end if

        if (mod(w,2)==1) then
            hw1=(w-1)/2
            xend=xidx+hw1
        else
            hw1=w/2-1
            xend=xidx+hw1+1
        end if

        !--------Crop element edges if extends out--------
        ex1=1
        ey1=1
        ex2=w
        ey2=h
        if (yidx-hh1<1) then
            if (yend>ny) then
                ey1=hh1-yidx+2
                ey2=ny-yidx+hh1+1
            else
                ey1=hh1-yidx+2
            end if
        else
            if (yend>ny) then
                ey2=ny-yidx+hh1+1
            end if
        end if

        if (xidx-hw1<1) then
            if (xend>nx) then
                ex1=hw1-xidx+2
                ex2=nx-xidx+hw1+1
            else
                ex1=hw1-xidx+2
            end if
        else
            if (xend>nx) then
                ex2=nx-xidx+hw1+1
            end if
        end if

        slab(max(1,yidx-hh1):min(ny,yend), max(1,xidx-hw1):min(nx,xend))=element(ey1:ey2,ex1:ex2)
        
    end subroutine REGIONCOPY




    subroutine POLYFIT2D(x,y,z,degree,weights,n,zfit)
    !------Fit 2d polynomial using weighted least squares------
    !
    ! <x>, <y>: input, real, 1d array of length <n>, coordinates in x and y axes.
    ! <z>: input, real, 1d array of length <n>, function of (x,y).
    ! <degree>: input, int, degree of polynomial to fit.
    ! <weights>: input, real, 1d array of lenght <n>, weights used in the least
    !            squares fit.
    ! <n>: input, int, length of <x>, <y>, <z>, <weights> and <zfit>.
    ! <zfit>: output, real, 1d array of length of <n>, fitted values.
    !
    ! Calls the svd_solve() subroutine for the least squares fitting.

        implicit none
        integer, parameter :: ikind = selected_real_kind(p=15,r=10)
        real(kind=ikind), intent(in), dimension(n) :: x, y, z
        real(kind=ikind), intent(in), dimension(n) :: weights
        integer, intent(in) :: n, degree
        real(kind=ikind), intent(out), dimension(n) :: zfit

        !----------------Intermediate vars----------------
        real(kind=ikind), dimension(n,(degree+1)*(degree+2)/2) :: c, a
        real(kind=ikind), dimension((degree+1)*(degree+2)/2) :: coeff
        real(kind=ikind), dimension(n) :: sw
        integer :: i,j,k,npol

        sw=sqrt(weights)
        zfit=z*sw
        npol=(degree+1)*(degree+2)/2

        !------------------Fill matrix A------------------
        k=1
        do i = 0,degree
            do j = 0,degree-i
                if (i==0 .and. j==0) then
                    c(:,k)=1.
                else
                    c(:,k)=x**j * y**i
                end if
                a(:,k)=c(:,k)*sw
                k=k+1
            end do
        end do

        !--------------------Solve Ax=b--------------------
        call svd_solve(n,npol,a,zfit,coeff)
        !call DGELS('n',n,npol,1,a,n,zcopy,n,work,lwork,ok)
        !!call DGELSY(n,npol,1,a,n,zcopy,n,jpvt,-1.,rank,work,lwork,ok)
        !call DGELSD(n,npol,1,a,n,zcopy,n,s,-1.,rank,work,lwork,iwork,ok)
        !coeff=zcopy(1:npol)

        !-----------------------Fit-----------------------
        zfit=matmul(c,coeff)

    end subroutine POLYFIT2D





    !----------------Find array median----------------
    function MEDIAN(array) result(med)
    ! Find array median
    ! <array>: input, real, 1d array to find median.
    ! <med>: output, real, median of <array>.
        implicit none
        integer, parameter :: ikind = selected_real_kind(p=15,r=10)
        real(kind=ikind), dimension(:) :: array
        real(kind=ikind), allocatable, dimension(:) :: arraycopy
        integer, allocatable, dimension(:) :: sortidx
        real(kind=ikind) :: med
        integer :: n

        n=size(array)
        allocate(sortidx(n))
        allocate(arraycopy(n))
        arraycopy=array
        call mrgrnk(arraycopy,sortidx)
        arraycopy=arraycopy(sortidx)
        if (mod(n,2)==1) then
            med=arraycopy((n+1)/2)
        else
            med=(arraycopy(n/2) + arraycopy(n/2+1))/2
        end if
        deallocate(sortidx)
        deallocate(arraycopy)

    end function MEDIAN

end module m_lowess2d
