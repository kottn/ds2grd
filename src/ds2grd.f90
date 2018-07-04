program main
! _/  _/  _/  _/  _/  _/  _/  _/  _/  _/  _/  _/  _/  _/  _/  _/  _/  _/
implicit none

integer :: ioerr
integer :: i, j, k, ii, jj, kk
integer :: na

character(len= 20), allocatable :: ca(:)
integer, allocatable :: nx(:), ny(:), nm(:)
integer :: nnx, nny, nnm, nst
double precision, allocatable :: diff(:), xulp(:), yulp(:), width(:), height(:)
double precision, allocatable :: x0(:,:), y0(:,:), dp(:,:), fm(:,:), df(:,:)
integer, allocatable :: ir(:,:)
double precision :: xmin, ymin, zmin
double precision :: xmax, ymax, zmax

integer, allocatable :: id(:), fg(:), lu(:)
double precision, allocatable :: xx(:), yy(:), zz(:), nn(:), ft(:), uu(:)

character(len= 99 ) :: dmmy01, dmmy02, charm
character(len= 99 ) :: depthpath, depthfile
character(len= 99 ) :: fmpath, fmfile
character(len= 99 ) :: irpath, irfile
character(len= 99 ) :: deformpath, deformfile, expath
character(len= 99 ) :: datdir, datpath, datfile
character(len= 99 ) :: ascdir, ascpath, ascfile_ID, ascfile_ZZ, ascfile_FT, ascfile_LU
character(len= 99 ) :: wktdir, wktpath, wktfile
character(len=  *), parameter :: infofile = './datainfo.txt'

logical :: datain    = .false., datain_fm = .false.
logical :: datain_ir = .false., datain_df = .false.
logical :: datout    = .false., ascout    = .false., wktout = .false.
! _/  _/  _/  _/  _/  _/  _/  _/  _/  _/  _/  _/  _/  _/  _/  _/  _/  _/

!-- pre.
write(6,*)''
write(6,*)infofile
open(10,file=infofile,status='old')

do
  read(10,'(a1, a)',iostat=ioerr) dmmy01, dmmy02
  if( dmmy01 == '!' )cycle
  charm=trim(adjustl(dmmy01))//trim(adjustl(dmmy02))

  select case (charm)

    case ('<DATAIN>')
      datain = .true.
      read(10,'(a)')depthpath

    case ('<DATAIN_FM>')
      datain_fm = .true.
      read(10,'(a)')fmpath

    case ('<DATAIN_IR>')
      datain_ir = .true.
      read(10,'(a)')irpath

    case ('<DATAIN_DF>')
      datain_df = .true.
      read(10,'(a)')deformpath
      read(10,'(a)')expath

    case ('<DATOUT>')
      datout = .true.
      read(10,'(a)')datdir
      read(10,'(a)')datpath
      datpath = trim(adjustl(datdir))//trim(adjustl(datpath))
      call system( 'mkdir -p "'//trim(datdir)//'" ' )

    case ('<ASCOUT>')
      ascout = .true.
      read(10,'(a)')ascdir
      read(10,'(a)')ascpath
      ascpath = trim(adjustl(ascdir))//trim(adjustl(ascpath))
      call system( 'mkdir -p "'//trim(ascdir)//'" ' )

    case ('<WKTOUT>')
      wktout = .true.
      read(10,'(a)')wktdir
      read(10,'(a)')wktpath
      wktpath = trim(adjustl(wktdir))//trim(adjustl(wktpath))
      call system( 'mkdir -p "'//trim(wktdir)//'" ' )

    case ('<DATALIST>')
      read(10,*)na
      allocate(  ca(na), diff(na), xulp(na), yulp(na), &
                 nx(na), ny(na), width(na), height(na) )

      read(10,*)(ca(ii), diff(ii), xulp(ii), yulp(ii), &
                 nx(ii), ny(ii), width(ii), height(ii), ii = 1, na)

    case ('<END>')
      exit

  end select

end do
close(10)


!-- wkt outp.
if( wktout )then
  write(6,*)''
  write(6,*)''
  write(6,*)'--- WKTOUT ---'
  call makeWKT( wktpath, na, ca, diff, xulp, yulp, width, height )
  write(6,*)''
  write(6,*)''
end if

!-- area loop
allocate( nm(na) )
do i = 1, na

  allocate( x0(nx(i),ny(i)), y0(nx(i),ny(i)), dp(nx(i),ny(i)), &
            fm(nx(i),ny(i)), ir(nx(i),ny(i)), df(nx(i),ny(i)) )


  write(6,*)''
  write(6,*)' Area: ', ca(i)
  write(6,*)'*******************************************'

  !-- adj. upper left point
  !
  ! Chubo XY (node)   ---->  NK xy (cell center)
  !
  !     1 ------> X               1 ----> x
  !   1 @---------@              +---------+      y
  !   | |_|_|_|_|_|            1 |@|_|_|_|@|
  !   | |_|_|_|_|_|            | |_|_|_|_|_|      ^
  !   | |_|_|_|_|_| ny         v |_|_|_|_|_| ny   |
  !   v | | | | | |            y |@| | | |@|      @-- >  x
  !   Y @---------@              +---------+
  !         nx                       nx

  xulp(i) = xulp(i) + diff(i)*0.5
  yulp(i) = yulp(i) - diff(i)*0.5

  write(6,*)'nx, ny: ', nx(i), ny(i)
  write(6,*)''
  write(6,*)'--- DATAIN ---'
  !-- xy gen.
  do jj = 1, ny(i)
    do ii = 1, nx(i)
      x0(ii,jj) = xulp(i) + (ii - 1) * diff(i)
      y0(ii,jj) = yulp(i) - (jj - 1) * diff(i)
    end do
  end do

  !-- inp.
  if( datain )then
    depthfile  = trim(adjustl(depthpath ))//trim(adjustl(ca(i)))//'.dat'
    open(10,file=depthfile,  status='old',action='read')
    write(6,*)depthfile
  end if
  if( datain_fm )then
    fmfile     = trim(adjustl(fmpath    ))//trim(adjustl(ca(i)))//'.dat'
    open(11,file=fmfile,     status='old',action='read')
    write(6,*)fmfile
  end if
  if( datain_ir )then
    irfile     = trim(adjustl(irpath    ))//trim(adjustl(ca(i)))//'.dat'
    open(12,file=irfile,     status='old',action='read')
    write(6,*)irfile
  end if
  if( datain_df )then
    deformfile = trim(adjustl(deformpath))//trim(adjustl(ca(i)))//trim(adjustl(expath))//'.dat'
    open(13,file=deformfile, status='old',action='read')
    write(6,*)deformfile
  end if

  do jj = 1, ny(i)

    !-- Chubo fmt.
    if( datain )then
      read(10,'(10f8.2)')(dp(ii  ,jj), dp(ii+1,jj), dp(ii+2,jj), dp(ii+3,jj), dp(ii+4,jj), &
                          dp(ii+5,jj), dp(ii+6,jj), dp(ii+7,jj), dp(ii+8,jj), dp(ii+9,jj), &
                          ii = 1, nx(i), 10)
    end if
    if( datain_fm )then
      read(11,'(10f8.5)')(fm(ii  ,jj), fm(ii+1,jj), fm(ii+2,jj), fm(ii+3,jj), fm(ii+4,jj), &
                          fm(ii+5,jj), fm(ii+6,jj), fm(ii+7,jj), fm(ii+8,jj), fm(ii+9,jj), &
                          ii = 1, nx(i), 10)
    end if
    if( datain_ir )then
      read(12,'(10i8)'  )(ir(ii  ,jj), ir(ii+1,jj), ir(ii+2,jj), ir(ii+3,jj), ir(ii+4,jj), &
                          ir(ii+5,jj), ir(ii+6,jj), ir(ii+7,jj), ir(ii+8,jj), ir(ii+9,jj), &
                          ii = 1, nx(i), 10)
    end if
    if( datain_df )then
      read(13,'(10f8.2)')(df(ii  ,jj), df(ii+1,jj), df(ii+2,jj), df(ii+3,jj), df(ii+4,jj), &
                          df(ii+5,jj), df(ii+6,jj), df(ii+7,jj), df(ii+8,jj), df(ii+9,jj), &
                          ii = 1, nx(i), 10)
    end if
  end do

  close(10)
  close(11)
  close(12)
  close(13)


  !-- 2 mesh clipping for NK nest mesh
  if( idint(diff(i)) == 1350 )then
    nst = 0
  else if( idint(diff(i)) < 1350 )then
    nst = 2
  end if
  nnx = nx(i) - nst*2
  nny = ny(i) - nst*2
  nnm = nnx*nny
  allocate( id(nnm),                            &
            xx(nnm), yy(nnm), zz(nnm), fg(nnm), &
            ft(nnm), nn(nnm), lu(nnm), uu(nnm) )

  !-- xy sort. & intpo. dataset
  !
  !   Now           ------->     inidom
  !     +---------+                +---------+   ny
  !     |@ → → @|   @-- > nx     |@       @|
  !     |     ↙  |   |            |↑↘   ↑|    ^
  !     |   ↙    |   v            |↑  ↘ ↑|    |
  !     |@  →→ @|  ny            |@       @|    @-- >  nx
  !     +---------+                +---------+

  kk = 0
  do ii = 1+nst, nx(i)-nst
    do jj = ny(i)-nst, 1+nst, -1

      kk = kk + 1
      id(kk) =  kk
      xx(kk) =  x0(ii,jj)
      yy(kk) =  y0(ii,jj)
      zz(kk) =  dp(ii,jj)*(-1)

      nn(kk) = 0.025
      if( datain_fm )then
        nn(kk) = fm(ii,jj)
      end if

      fg(kk) = 0
      ft(kk) = 0.d0
      if( datain_ir )then
             if( 1000 < ir(ii,jj) .and. ir(ii,jj) < 2000 )then
          fg(kk) = 1
          ft(kk) = (ir(ii,jj) - 1000)*0.1
        else if( 2000 < ir(ii,jj) .and. ir(ii,jj) < 3000 )then
          fg(kk) = 2
          ft(kk) = (ir(ii,jj) - 2000)*0.1
        else if( 3000 < ir(ii,jj) .and. ir(ii,jj) < 4000 )then
          fg(kk) = 3
          ft(kk) = (ir(ii,jj) - 3000)*0.1
        end if
      end if

      lu(kk) = 1
      if( zz(kk) < 0.001 )then
        lu(kk) = -1
      end if

      uu(kk) = 0.d0
      if( datain_df )then
        uu(kk) = df(ii,jj)
      end if

    end do
  end do

  !-- outp.
  if( datout )then
    write(6,*)''
    write(6,*)'--- DATOUT ---'
    datfile = trim(adjustl(datpath))//trim(adjustl(ca(i)))//'.dat'
    open(51,file=datfile,status='replace',action='write')
    write(6,*)datfile

    write(51,5000)(id(kk), xx(kk), yy(kk), zz(kk), nn(kk), &
                           fg(kk), ft(kk), lu(kk), uu(kk), kk = 1, nnm)

    close(51)

    xmin = minval( xx(1:nnm) )
    ymin = minval( yy(1:nnm) )
    zmin = minval( zz(1:nnm) )
    xmax = maxval( xx(1:nnm) )
    ymax = maxval( yy(1:nnm) )
    zmax = maxval( zz(1:nnm) )
    write(6,*)'Xmin, Xmax: ',  xmin,  xmax
    write(6,*)'Ymin, Ymax: ',  ymin,  ymax
    write(6,*)'Zmin, Zmax: ',  zmin,  zmax
    write(6,*)'  nx,   ny: ', nx(i), ny(i)
    write(6,*)' nnx,  nny: ',   nnx,   nny
  end if

  if( ascout )then
    ascfile_ZZ = trim(adjustl(ascpath))//trim(adjustl(ca(i)))//'_Elv.asc'
    ascfile_ID = trim(adjustl(ascpath))//trim(adjustl(ca(i)))//'_ID.asc'
    ascfile_FT = trim(adjustl(ascpath))//trim(adjustl(ca(i)))//'_FTEN.asc'
    ascfile_LU = trim(adjustl(ascpath))//trim(adjustl(ca(i)))//'_LandUse.asc'
    write(6,*)''
    write(6,*)'--- ASCOUT (Elv) ---'
    write(6,*)ascfile_ZZ
    call xyz2asc( ascfile_ZZ, xx, yy,       zz, nnm )
    write(6,*)'--- ASCOUT (ID) ---'
    write(6,*)ascfile_ID
    call xyz2asc( ascfile_ID, xx, yy, dble(id), nnm )
    write(6,*)'--- ASCOUT (FTEN) ---'
    write(6,*)ascfile_FT
    call xyz2asc( ascfile_FT, xx, yy,       ft, nnm )
    write(6,*)'--- ASCOUT (LandUse) ---'
    write(6,*)ascfile_LU
    call xyz2asc( ascfile_LU, xx, yy, dble(lu), nnm )
  end if

  deallocate( x0, y0, dp, fm, ir, df )
  deallocate( id, xx, yy, zz, fg, ft, nn, lu, uu )

end do

5000 format(i8,3f15.5, f15.5, i3, f15.5, i3, f15.5)

end program main


!-----------------------------------------------------------------------
subroutine makeWKT( wktfilename, ndata, ca, diff, xulp, yulp, x_dmsize, y_dmsize )
!-----------------------------------------------------------------------
implicit none

integer :: i,j,k
integer ::  nd

double precision, allocatable :: p1(:,:), p2(:,:), p3(:,:), p4(:,:)
double precision :: newdiff
character(len= 50) :: wktfilename2
character(len= 50) :: TMP

character(len=  *), intent(in) :: wktfilename
integer, intent(in) :: ndata
character(len=  *), intent(in) :: ca(ndata)
double precision, intent(in) :: diff(ndata), xulp(ndata), yulp(ndata)
double precision, intent(in) :: x_dmsize(ndata), y_dmsize(ndata)
!-----------------------------------------------------------------------

!-- calc.
allocate(p1(2,ndata),p2(2,ndata),p3(2,ndata),p4(2,ndata))

newdiff = 0.d0
nd = 0
open(50+nd, status='scratch')
do i = 1, ndata

    if( newdiff > diff(i) .or. newdiff < diff(i) )then
      close(50+nd)

      nd = nd +1
      newdiff = diff(i)

      write(TMP,'(a4)') ca(i)

      wktfilename2 = trim(adjustl(wktfilename))//trim(adjustl(TMP))//'.wkt'
      open(50+nd, file=wktfilename2, status='replace')
      write(50+nd, '(a)')'Area,WKT'

      write(6,*) wktfilename2
    end if

    p1(1,i) = xulp(i)
    p2(1,i) = xulp(i)
    p3(1,i) = xulp(i) + x_dmsize(i)
    p4(1,i) = xulp(i) + x_dmsize(i)

    p1(2,i) = yulp(i)
    p2(2,i) = yulp(i) - y_dmsize(i)
    p3(2,i) = yulp(i) - y_dmsize(i)
    p4(2,i) = yulp(i)

    100 format(3a, 2f15.4, a, 2f15.4, a, 2f15.4, a, 2f15.4, a, 2f15.4, a)

    write(50+nd, 100) trim(ca(i)),',', '"POLYGON ((', &
      p1(1,i), p1(2,i), ',', p2(1,i), p2(2,i), ',', &
      p3(1,i), p3(2,i), ',', p4(1,i), p4(2,i), ',', &
      p1(1,i), p1(2,i), '))"'

end do

write(6,*) 'nlist: ', ndata
write(6,*) ' .wkt: ', nd

end subroutine

!-----------------------------------------------------------------------
subroutine xyz2asc(grdfilename, x, y, z, N)
!-----------------------------------------------------------------------
  IMPLICIT NONE
!  CHARACTER(8192) :: csvfilename !kotani
!  CHARACTER(8192) :: grdfilename !kotani
!  CHARACTER(8192) :: hdrfilename !kotani
  CHARACTER(8192) :: line
  REAL(8) :: dummy(100), YY
!  REAL(8),ALLOCATABLE :: x(:), y(:), z(:), zz(:,:) !kotani
  REAL(8),ALLOCATABLE :: zz(:,:)  !kotani
  INTEGER,ALLOCATABLE :: nzz(:,:)
  REAL(8) :: dx, dy
  REAL(8) :: xmax, xmin
  REAL(8) :: ymax, ymin
  REAL(8) :: zmax, zmin
  REAL(8) :: xl, yl
!  INTEGER :: N, nx, ny !kotani
  INTEGER :: nx, ny !kotani
  INTEGER :: iix, iiy
  INTEGER :: i, j, ic, ii, iostatus
  INTEGER :: iX, iY, iZ
  LOGICAL :: ldx = .true., ldy = .true.
  INTEGER :: nonval = -9999
  INTEGER :: nargs
  CHARACTER(12) :: cncols

  !--/ kotani
  CHARACTER(*), intent(in) :: grdfilename
  INTEGER, intent(in) :: N
  REAL(8), intent(in) :: x(N), y(N), z(N)
  !--/
!-----------------------------------------------------------------------

!--/ kotani
!  nargs = iargc()
!  CALL GETARG(1,line); READ(line,*) iX
!  CALL GETARG(2,line); READ(line,*) iY
!  CALL GETARG(3,line); READ(line,*) iZ
!
!  CALL GETARG(4,csvfilename)
!
!  IF(nargs > 4 )THEN
!    CALL GETARG(5,line); READ(line,*) dx
!    ldx = .false.
!  END IF
!  IF(nargs > 5 )THEN
!    CALL GETARG(6,line); READ(line,*) dy
!    ldy = .false.
!  END IF
!
!
!  csvfilename = ADJUSTL(csvfilename)
!
!  ic          = scan(csvfilename,'.',back=.true.)
!  grdfilename = csvfilename(:ic-1)//'.asc'
!
!
!  OPEN(10,file=csvfilename,status='old')
!  !READ(10,'(A)')line
!  N = 0
!  DO
!    READ(10,'(A)',IOSTAT=iostatus)line
!    IF( iostatus/=0 ) EXIT
!    N = N + 1
!  END DO
!
!  WRITE(*,*)'N:',N
!  ALLOCATE(x(N),y(N),z(N))
!
!  CLOSE(10)
!  OPEN(10,file=csvfilename,status='old')
!  !READ(10,'(A)')line
!
!  DO i = 1, N
!    READ(10,'(A)',IOSTAT=iostatus)line
!    READ(line,*) dummy(1:iX-1),x(i)
!    READ(line,*) dummy(1:iY-1),y(i)
!    READ(line,*) dummy(1:iZ-1),z(i)
!  END DO
!--/
  xmin = minval( x(1:N) )
  ymin = minval( y(1:N) )
  zmin = minval( z(1:N) )
  xmax = maxval( x(1:N) )
  ymax = maxval( y(1:N) )
  zmax = maxval( z(1:N) )
  WRITE(*,*)'Xmin,Xmax',xmin,xmax
  WRITE(*,*)'Ymin,Ymax',ymin,ymax
  WRITE(*,*)' Min, Max',zmin,zmax

  IF( ldx )THEN
    DO i = 2, N
      IF( dabs(x(i) - x(1)) > 1.d-6 )THEN
        dx  = dabs( x(i) - x(1) )
        EXIT
      END IF
    END DO
  END IF
  IF( ldy )THEN
    DO i = 2, N
      IF( dabs(y(i) - y(1)) > 1.d-6 )THEN
        dy  = dabs( y(i) - y(1) )
        EXIT
      END IF
    END DO
  END IF

  nx = idnint ( ( xmax - xmin ) / dx + 1 )
  ny = idnint ( ( ymax - ymin ) / dy + 1 )
  WRITE(*,*)'nx, ny: ',nx, ny



  ALLOCATE( ZZ(nx, ny) )
  ALLOCATE(NZZ(nx, ny) )
  DO j = 1, ny
    DO i = 1, nx
      ZZ ( i, j ) = 0.d0
      NZZ( i, j ) = 0
    END DO
  END DO

  DO i = 1, N
    if ( Z(i) == DBLE( nonval ) )cycle
    xl = X(i)
    yl = Y(i)

    iix = idnint ( ( xl - xmin ) / dx + 1 )
    iiy = idnint ( ( yl - ymin ) / dy + 1 )

    !ZZ(iix,iiy) = Z(i)
    ZZ (iix,iiy) = ZZ (iix,iiy) + Z(i)
    NZZ(iix,iiy) = NZZ(iix,iiy) + 1
  END DO

  DO j = 1, ny
    DO i = 1, nx
      IF( NZZ( i, j ) /= 0 )THEN
        ZZ ( i, j ) = ZZ ( i, j ) / NZZ( i, j )
      ELSE
        ZZ ( i, j ) = DBLE( nonval )
      END IF
    END DO
  END DO


  WRITE(*,'(A)')'ASC:'//TRIM(grdfilename) !kotani
  open(60,FILE=grdfilename,STATUS='replace')

  WRITE(line,*)nx
  WRITE                                  (60,'(A)') "ncols "//TRIM(ADJUSTL(line))
  WRITE(line,*)ny
  WRITE                                  (60,'(A)') "nrows "//TRIM(ADJUSTL(line))
  WRITE(line,'(f15.3)') xmin-dx*0.5d0
  WRITE                                  (60,'(A)') "xllcorner "//TRIM(ADJUSTL(line))
  WRITE(line,'(f15.3)') ymin-dy*0.5d0
  WRITE                                  (60,'(A)') "yllcorner "//TRIM(ADJUSTL(line))
  WRITE(line,'(f18.6)') dx
  WRITE                                  (60,'(A)') "cellsize "//TRIM(ADJUSTL(line))
  WRITE                                  (60,'(A)') "nodata_value -9999.000"
  WRITE(line,'(f18.6)') dx
  WRITE                                  (60,'(A)') "cellXsize "//TRIM(ADJUSTL(line))
  WRITE(line,'(f18.6)') dy
  WRITE                                  (60,'(A)') "cellYsize "//TRIM(ADJUSTL(line))

  WRITE(cncols,*) nx
  DO j = ny, 1, -1
    WRITE(60,'('//TRIM(ADJUSTL(cncols))//'e15.7)') ( ZZ(i,j), i=1, nx )
  END DO

  close(60)

  close(60)


END subroutine
