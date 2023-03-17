module nr
   public :: ludcmp, lubksb, tridag
   interface
      subroutine ludcmp(a,n,np,indx,d)
         integer, intent(in) :: n, np
         real, intent(inout), dimension(np,np) :: a
         integer, intent(out), dimension(n) :: indx
         real, intent(out) :: d
      end subroutine ludcmp
      subroutine lubksb(a,n,np,indx,b)
         integer, intent(in) :: n, np
         real, intent(in), dimension(np,np) :: a
         integer, intent(in), dimension(n) :: indx
         real, intent(inout), dimension(n) :: b
      end subroutine lubksb
      subroutine tridag(a,b,c,r,u,n)
         integer, intent(in) :: n
         real, intent(in),  dimension(n) :: a, b, c, r
         real, intent(out), dimension(n) :: u
      end subroutine tridag
   end interface
end module nr

module constants
   implicit none
   ! Note that i runs 1..nx, j runs 0..ny
   integer, parameter, public :: nx=16, ny=16
   real, parameter, public :: dx=3.75e5, dy=6.25e5   ! Gridsize
!   real, parameter, public :: dx=1e6, dy=dx   ! Gridsize
   real, parameter, public :: eps=dx/dy, epsq = eps**2
   real, parameter, public :: heat = 2.0e-3 ! Heating W/kg
   real, parameter, public :: lambdasq = 1.5e-12 ! Stability parameter, m^-2
   real, parameter, public :: rgas = 287.0
   real, parameter, public :: cp = 1004.0
   real, parameter, public :: f0 = 1.0e-4
   real, parameter, public :: beta = 1.6e-11 ! s^{-1} m^{-1}
   real, parameter, public :: p4 = 1.0e5 ! Surface pressure, Pa
   real, parameter, public :: p2 = 0.5*p4
   real, parameter, public :: gamma = lambdasq*dx**2  ! 0.21
   real, parameter, public :: pi = 3.1415926
   real, parameter, public :: k = 4.0e-6 ! s^{-1}  Surface drag
!  Control variables
   real, public, save :: a = 1.0e5  ! m^2/s   Horizonal diffusion
   logical, public, save :: imp = .true.
   real, public, save :: dt1 = 86400.0, dt2 = 3600.0
   real, public, save :: day1 = 131.0, day2 = 164.5
   logical, public, save :: diag_flag = .false.
   real, parameter, public  :: accel = 1.0
end module constants
 

module work
  implicit none
  public :: split, calcvor, calc_zonstream,  relax1, xcalc, diag, diagvv, &
            calc_zvor, relax2
contains

   subroutine split(xt, xz, x)
      use constants
      real, intent(in),  dimension(:,0:) :: xt
      real, intent(out), dimension(0:)   :: xz
      real, intent(out), dimension(:,0:) :: x
      integer :: j
      do j=0,ny
         xz(j) = sum(xt(:,j)) / nx
      end do
      x = xt - spread(xz,dim=1,ncopies=nx)
   end subroutine split

   subroutine calcvor(s1, s3, v1, v3)
      use constants
      real, intent(in),  dimension(:,0:) :: s1, s3
      real, intent(out), dimension(:,0:) :: v1, v3
      integer :: i, j, im, ip, jm, jp

      v1(:,0)  = 0.0
      v1(:,ny) = 0.0
      v3(:,0)  = 0.0
      v3(:,ny) = 0.0
      do j=1,ny-1
         jm = j-1
         jp = j+1
         do i=1,nx
            im = i-1
            if ( im == 0 ) then
               im = nx
            end if
            ip = i+1
            if ( ip == nx+1 ) then
               ip = 1
            end if
            v1(i,j) =       ( s1(ip,j)  + s1(im,j) - 2*s1(i,j) ) +         &
                    epsq  * ( s1(i,jp)  + s1(i,jm) - 2*s1(i,j) ) -         &
                    gamma * ( s1(i,j) - s3(i,j) )
            v3(i,j) =       ( s3(ip,j)  + s3(im,j) - 2*s3(i,j) ) +         &
                    epsq  * ( s3(i,jp)  + s3(i,jm) - 2*s3(i,j) ) + &
                    gamma * ( s1(i,j) - s3(i,j) )
         end do
      end do
   end subroutine calcvor

   subroutine calc_zonstream(v1z, v3z, s1z, s3z)
!     Solve for the zonal mean streamfunction
      use constants
      use nr
      real, intent(in),  dimension(0:) :: v1z, v3z
      real, intent(out), dimension(0:) :: s1z, s3z
      integer, parameter :: nz = 2*ny - 3  ! Number of zonal means to solve for
      real, dimension(nz) :: bmat, junk
      real, dimension(nz,nz) :: amat
      integer :: j
!     For LU solver
      real :: d
      integer, dimension(nz) :: indx

      amat = 0.0
      ! j = 1, level 1
      amat(1,1) = -epsq - gamma
      amat(1,2) = epsq
      ! j = J-1, level 1
      amat(ny-1,ny-2) = epsq
      amat(ny-1,ny-1) = -epsq - gamma
      amat(ny-1,nz) = gamma
      ! j = 2, level 3
      amat(ny,2) = gamma
      amat(ny,ny) = -2.0*epsq - gamma
      amat(ny,ny+1) = epsq
      ! j = J-1, level 3
      amat(nz,ny-1) = gamma
      amat(nz,nz-1) = epsq
      amat(nz,nz) = -epsq - gamma

      !  Level 1
      do j=2,ny-2
         amat(j,j-1) = epsq
         amat(j,j) = -2.0*epsq - gamma
         amat(j,j+1) = epsq
         amat(j,ny-2+j) = gamma
      end do
      !  Level 3
      do j=ny+1,nz-1
         amat(j,j+2-ny) = gamma
         amat(j,j-1) = epsq
         amat(j,j) = -2.0*epsq - gamma
         amat(j,j+1) = epsq
      end do

!!$      print*, "AMAT"
!!$      do j=1,nz
!!$         write(*,"(7f6.2)") amat(j,:)
!!$      end do
!!$
!!$      do j=1,ny-1
!!$         junk(j) = s1z(j)
!!$      end do
!!$      do j=2,ny-1
!!$         junk(ny-2+j) = s3z(j)
!!$      end do
!!$
!!$      junk = matmul(amat,junk)

      do j=1,ny-1
         bmat(j) = v1z(j)
      end do
      do j=2,ny-1
         bmat(ny-2+j) = v3z(j)
      end do

      ! Solve AX=B
      call ludcmp(amat,nz,nz,indx,d)
      call lubksb(amat,nz,nz,indx,bmat)

!!$      do j=1,nz
!!$         write(unit=4,fmt="(i3,2e15.7)") j, bmat(j), junk(j)
!!$      end do

      do j=1,ny-1
         s1z(j) = bmat(j)
!         write(unit=4,fmt="(i3,e15.7)") j, bmat(j)
      end do
      do j=2,ny-1
         s3z(j) = bmat(ny-2+j)
!         write(unit=7,fmt="(i3,e15.7)") j, bmat(ny-2+j)
      end do

!     Apply the BC
      s1z(0) = s1z(1)
      s1z(ny) = s1z(ny-1)
      s3z(0) = 0.0
      s3z(ny) = s3z(ny-1)

   end subroutine calc_zonstream

   subroutine relax1(v1, v3, s1, s3)
      ! Solve for anomaly streamfunction
      use constants
      real, intent(in),  dimension(:,0:) :: v1, v3
      real, intent(inout), dimension(:,0:) :: s1, s3
      integer :: i, j, im, ip, jm, jp
      integer :: iter, iadd, irb
      real :: change, maxdiff, resid

!!$      open(unit=99,file="relax.in",form="unformatted",status="unknown", &
!!$           action="write")
!!$      write(unit=99) v1
!!$      write(unit=99) v3
!!$      write(unit=99) s1
!!$      write(unit=99) s3
!!$      close(unit=99)
      ! Start from the current value of the anomaly streamfunction
      do iter=1,100
         ! Jacobi iteration
         maxdiff = 0.0
         change = 0.0
         do irb=0,1
            do j=1,ny-1
               jm = j-1
               jp = j+1
               do i=1,nx
                  im = i-1
                  if ( im == 0 ) then
                     im = nx
                  end if
                  ip = i+1
                  if ( ip == nx+1 ) then
                     ip = 1
                  end if

                  resid = ( s1(ip,j) + s1(im,j) +     &
                       epsq * ( s1(i,jp) + s1(i,jm) ) - &
                       v1(i,j) + gamma*s3(i,j) ) -      &
                       ( 2.0 + 2.0*epsq + gamma )*s1(i,j)
                  resid = accel*resid / ( 2.0 + 2.0*epsq + gamma )
                  change = change + resid**2 
                  maxdiff = max ( maxdiff, abs(resid) )
                  s1(i,j) = s1(i,j) + resid

                  resid = ( s3(ip,j) + s3(im,j) +     &
                       epsq * ( s3(i,jp) + s3(i,jm) ) - &
                       v3(i,j) + gamma*s1(i,j) ) -      &
                        ( 2.0 + 2.0*epsq + gamma )*s3(i,j)
                  resid = accel*resid / ( 2.0 + 2.0*epsq + gamma )
                  change = change + resid**2
                  maxdiff = max ( maxdiff, abs(resid) )
                  s3(i,j) = s3(i,j) + resid
               end do
            end do
         end do
!            print*, "ITER1", iter, sqrt(change), maxdiff
         ! maxdiff is now only on a single level so halve the convergence 
         ! criterion/
         if ( maxdiff < 0.5*3.75e4 ) then
            ! print*, "ITER1", iter, sqrt(change), maxdiff
            exit
         end if
      end do
      if ( iter >= 100 ) then
         print*, "RELAX1 failed", iter, sqrt(change), maxdiff
         stop
      endif
   end subroutine relax1

   subroutine relax2(x1, x3, dt, v1, v3)
      ! Solve for anomaly vorticity
      use constants
      real, intent(in),  dimension(:,0:) :: x1, x3
      real, intent(in) :: dt
      real, intent(out), dimension(:,0:) :: v1, v3
      real, dimension(nx,0:ny) :: temp1, temp3
      integer :: i, j, im, ip, jm, jp
      integer :: iter
      real :: change1, change3, alpha

      v1 = 0.0
      v3 = 0.0
      alpha = a*dt/dx**2
      do iter=1,100
         ! Jacobi iteration
         do j=1,ny-1
            jm = j-1
            jp = j+1
            do i=1,nx
               im = i-1
               if ( im == 0 ) then
                  im = nx
               end if
               ip = i+1
               if ( ip == nx+1 ) then
                  ip = 1
               end if
               temp1(i,j) = ( alpha*( v1(ip,j) + v1(im,j) + &
                             epsq * ( v1(i,jp) + v1(i,jm) ) ) + &
                             x1(i,j)  ) / &
                            ( 2*alpha*(1.0 + epsq)  + 1.0 )
               temp3(i,j) = ( alpha*( v3(ip,j) + v3(im,j) + &
                             epsq * ( v3(i,jp) + v3(i,jm) ) ) + &
                             x3(i,j)  ) / &
                            ( 2*alpha*(1.0 + epsq)  + 1.0 +1.5*k*dt )
            end do
         end do
         change1 = sum ( ( v1(:,1:ny-1) - temp1(:,1:ny-1) ) **2 )
         v1(:,1:ny-1) = temp1(:,1:ny-1)
         v1(:,0) = v1(:,1)
         v1(:,ny) = v1(:,ny-1)
         change3 = sum ( ( v3(:,1:ny-1) - temp3(:,1:ny-1) ) **2 )
         v3(:,1:ny-1) = temp3(:,1:ny-1)
         v3(:,0) = v3(:,1)
         v3(:,ny) = v3(:,ny-1)
         if ( max(change1, change3) < 1.0 ) then
            ! print*, "ITER2", iter, sqrt(change1), sqrt(change3)
            exit
         end if
      end do
   end subroutine relax2

   subroutine xcalc(v1, v3, v1m, v3m, s1, s3, s1m, s3m, dt, x1, x3)
      use constants
      real, dimension(:,0:), intent(in)  :: v1, v3, v1m, v3m
      real, dimension(:,0:), intent(in)  :: s1, s3, s1m, s3m
      real, intent(in) :: dt
      real, dimension(:,0:), intent(out) :: x1, x3
      integer :: i, j, im, ip, jm, jp
      real :: alpha, b, c, h

      alpha = a*dt/dx**2
      if ( .not. imp ) then
         alpha = 2.0 * alpha
      end if
      b = beta*dx**2*dy
      c = dt/(2.0*dx*dy)
      h = 4*rgas*heat*gamma*dt/(f0*cp)
      x1 = 0.0
      x3 = 0.0
      do j=1,ny-1
         jm = j-1
         jp = j+1
         do i=1,nx
            im = i-1
            if ( im == 0 ) then
               im = nx
            end if
            ip = i+1
            if ( ip == nx+1 ) then
               ip = 1
            end if

            x1(i,j) = v1m(i,j) +                                           &
                      c * ( (v1(ip,j)-v1(im,j))*(s1(i,jp)-s1(i,jm)) -      &
                          (2*b+v1(i,jp)-v1(i,jm))*(s1(ip,j)-s1(im,j)) ) +  &
                      alpha * ( v1m(ip,j)+v1m(im,j)-2*v1m(i,j) +           &
                                epsq*(v1m(i,jp)+v1m(i,jm)-2*v1m(i,j)) ) +  &
                      h*(2*j-ny)/ny

            x3(i,j) = v3m(i,j) +                                          &
                       c * ( (v3(ip,j)-v3(im,j))*(s3(i,jp)-s3(i,jm)) -     &
                           (2*b+v3(i,jp)-v3(i,jm))*(s3(ip,j)-s3(im,j)) ) + &
                      alpha * ( v3m(ip,j)+v3m(im,j)-2*v3m(i,j) +           &
                                epsq*(v3m(i,jp)+v3m(i,jm)-2*v3m(i,j)) ) -  &
                       h*(2*j-ny)/ny
            if ( imp ) then
               x3(i,j) = x3(i,j) -  k*dt*(1.5*v3m(i,j) - v1(i,j) - &
                                          4*gamma*(s1(i,j)-s3(i,j)) )
            else
               x3(i,j) = x3(i,j) -  k*dt*(3.0*v3m(i,j) - v1m(i,j) - &
                                          4*gamma*(s1m(i,j)-s3m(i,j)) )
            end if

         end do
      end do
   end subroutine xcalc

   subroutine calc_zvor(x1z, x3z, dt, v1z, v3z)
!     Solve for the zonal mean vorticity
      use constants
      use nr
      real, intent(in),  dimension(0:) :: x1z, x3z
      real, intent(in) :: dt
      real, intent(out), dimension(0:) :: v1z, v3z
      integer, parameter :: nz = ny-1  ! Number of zonal means to solve for
      real, dimension(nz) :: amat, bmat, cmat, rmat, umat
      integer :: j
      real :: alpha


      alpha = a*dt/dx**2
      ! Level 1
      amat(1) = 0.0
      bmat(1) = -alpha*epsq - 1
      cmat(1) = alpha*epsq
      do j=2,ny-2
         amat(j) = alpha*epsq
         bmat(j) = -2.0*alpha*epsq - 1
         cmat(j) = alpha*epsq
      end do
      amat(nz) = alpha*epsq
      bmat(nz) = -alpha*epsq - 1
      cmat(nz) = 0.0

      rmat = -x1z(1:nz)

      call tridag(amat,bmat,cmat,rmat,umat,nz)

      v1z(1:nz) = umat
      v1z(0) = v1z(1)
      v1z(ny) = v1z(ny-1)

      ! Level 3
      amat(1) = 0.0
      bmat(1) = -alpha*epsq - 1 - 1.5*k*dt
      cmat(1) = alpha*epsq
      do j=2,ny-2
         amat(j) = alpha*epsq
         bmat(j) = -2.0*alpha*epsq - 1 - 1.5*k*dt
         cmat(j) = alpha*epsq
      end do
      amat(nz) = alpha*epsq
      bmat(nz) = -alpha*epsq - 1 - 1.5*k*dt
      cmat(nz) = 0.0

      rmat = -x3z(1:nz)

      call tridag(amat,bmat,cmat,rmat,umat,nz)

      v3z(1:nz) = umat
      v3z(0) = v3z(1)
      v3z(ny) = v3z(ny-1)

   end subroutine calc_zvor
   
   subroutine diag(day, zonal, s1, s3)
      use constants
      real, intent(in) :: day
      logical, intent(in) :: zonal
      real, intent(in), dimension(:,0:) :: s1, s3
      real, dimension(0:ny) :: t2z
      real, dimension(nx,0:ny) :: t2
      real, dimension(1:ny) :: u1z, u3z, v1z, v3z
      real, dimension(nx,1:ny) :: u1t, u3t, u1, u3, v1, v3
      real, dimension(nx,0:ny) :: v1t, v3t
      real :: tke, zke, eke
      integer :: j, jj
      character(len=20) :: fname
      integer, save :: dcount = 1

      if ( zonal ) then
         t2z = f0*(s1(1,:)-s3(1,:))/rgas
         ! Zonal wind, u = -pd{psi}{y}
         u1z = - ( s1(1,1:ny) - s1(1,0:ny-1) ) / dy
         u3z = - ( s3(1,1:ny) - s3(1,0:ny-1) ) / dy

         ! Zonal KE, scaled such that a wind of 1 m/s everywhere gives 10
         zke = 10.0*sum(u1z*u1z + u3z*u3z)/(2*ny)
         ! write(unit=*,fmt="(a,f5.1,4f15.7)") "TZ ", day, maxval(t2z), maxval(u1z), maxval(u3z), zke
      else
         write(unit=3,fmt="(f8.3,2e15.7)") day, s1(1,8), s3(1,8)
         u1t = - ( s1(:,1:ny) - s1(:,0:ny-1) ) / dy
         u3t = - ( s3(:,1:ny) - s3(:,0:ny-1) ) / dy
         v1t = ( s1 - cshift(s1,-1) ) / dx
         v3t = ( s3 - cshift(s3,-1) ) / dx
         ! Average v to get it on the same grid points as u
         do j=ny,1,-1  ! Reverse loop so don't clobber unused values
            v1t(:,j) = 0.25*(v1t(:,j) + cshift(v1t(:,j),1) + &
                 v1t(:,j-1) + cshift(v1t(:,j-1),1) )
         end do
         ! Zonal mean winds
         do j=1,ny
            u1z(j) = sum(u1t(:,j)) / nx
            u3z(j) = sum(u3t(:,j)) / nx
            v1z(j) = sum(v1t(:,j)) / nx
            v3z(j) = sum(v3t(:,j)) / nx
         end do
         tke = 10.0*sum(u1t*u1t + u3t*u3t + v1t(:,1:)*v1t(:,1:) + v3t(:,1:)*v3t(:,1:))/(2*ny*nx)
         u1 = u1t - spread(u1z,dim=1,ncopies=nx)
         u3 = u3t - spread(u3z,dim=1,ncopies=nx)
         v1 = v1t(:,1:ny) - spread(v1z,dim=1,ncopies=nx)
         v3 = v3t(:,1:ny) - spread(v3z,dim=1,ncopies=nx)
         zke = 10.0*sum(u1z*u1z + u3z*u3z)/(2*ny)
         eke = 10.0*sum(u1*u1 + u3*u3 + v1*v1 + v3*v3)/(2*ny*nx)

         t2 = f0*(s1-s3)/rgas
         ! Variables are z1, z2, z4 (surface pressure), T2, u1, u3, v1, v3
         !write(fname,"(a,i3.3)") 'output', dcount
         dcount = dcount + 1
         open(unit=1,file="output",form="unformatted",status="unknown",action="write")
         write(unit=1) day
         write(unit=1) s1
         write(unit=1) s3
         write(unit=1) (1.5*s3(:,ny:0:-1) - 0.5*s1(:,ny:0:-1))*f0
         ! Flip sign so it looks more like temperature when plotted in SH
         write(unit=1) -t2
         write(unit=1) u1t, (/ (0.0,jj=1,nx) /)
         write(unit=1) u3t, (/ (0.0,jj=1,nx) /)
         write(unit=1) v1t(:,1:ny), (/ (0.0,jj=1,nx) /)
         write(unit=1) v3t(:,1:ny), (/ (0.0,jj=1,nx) /)
         close(unit=1)
         ! write(unit=*,fmt="(a,f5.1,6f15.7)") "KE ", day, maxval(t2), maxval(u1t), maxval(u3), zke, eke, tke
      end if

!!$      if ( istep == 130 ) then
!!$         do j=1,ny
!!$            write(unit=7,fmt="(i3,2e15.7)") j, u1z(j), u3z(j)
!!$         end do
!!$      end if

   end subroutine diag

   subroutine diagvv(s1t,s3t, s1tm, s3tm, dt)
      ! Vertical velocity diagnostic
      use constants
      real, intent(in),  dimension(:,0:) :: s1t, s3t, s1tm, s3tm
      real, intent(in) :: dt
      real, dimension(0:ny) :: ds1, ds3, s1, s3, w2, vv
      integer :: j

      ! Just work on zonal means
      ds1 = 0.0 ! s1t(1,:) - s1tm(1,:)
      ds3 = 0.0 !s3t(1,:) - s3tm(1,:)
      s1  = 0.5 * ( s1t(1,:) + s1tm(1,:) )
      s3  = 0.5 * ( s3t(1,:) + s3tm(1,:) )

      w2(0) = 0.0
      w2(ny) = 0.0
      do j=1,ny-1
         ! Jacobian term drops out because using zonal means.
         ! A term also dropped here
         w2(j) = lambdasq*p2/f0 * ( (ds1(j)-ds3(j))/dt + &
              2*rgas*heat*(2*j-ny)/(f0*cp*ny) )
      end do
      vv(0) = 0.0
      do j=1,ny-1
         vv(j) = vv(j-1) - dy/p2 * w2(j)
      end do
      vv(ny) = 0.0
      do j=0,ny
         write(unit=4,fmt="(i3,2e15.7)") j, w2(j), vv(j)
      end do
   end subroutine diagvv

end module work
 
program phillips

!  Implementation of Phillips' model.
!  The time integration is carried forward using the vorticity in a three
!  time level scheme. At each step the streamfunction must be diagnosed from 
!  this. 

!  The streamfunction is denoted by s and vorticity by v.
!  Vector wind components are uc and vc.
!  Total fields are denoted with suffix t and zonal means with suffix z

   use work
   use constants
   use ranlux_m
   implicit none
   real, parameter :: diag_freq = 1./12. ! 0.25
   real, dimension(nx,0:ny) :: s1t, s3t   ! Streamfunction
   real, dimension(nx,0:ny) :: s1tm, s3tm ! Streamfunction at tau-1
   real, dimension(nx,0:ny) :: s1, s3     ! Streamfunction anomalies
   real, dimension(0:ny)    :: s1z, s3z   ! Zonal mean streamfunction
   real, dimension(nx,0:ny) :: v1t, v3t   ! Vorticity
   real, dimension(nx,0:ny) :: v1, v3     ! Vorticity anomalies
   real, dimension(0:ny)    :: v1z, v3z   ! Zonal mean vorticity
   real, dimension(nx,0:ny) :: v1tm, v3tm ! tau - 1 values
!  Temporaries used in timestepping
   real, dimension(nx,0:ny) :: x1t, x3t
   real, dimension(nx,0:ny) :: x1, x3
   real, dimension(0:ny)    :: x1z, x3z
   
   integer :: i, j, im, ip, jm, jp, seed
   real    :: istep, day, swap
   logical :: zonal, addnoise
!   namelist /input/ a, implicit, dt1, dt2, day1, day2, accel

   real :: dt  ! Time step
   real, dimension(nx,ny-1) :: rand
   real, dimension(nx*ny-1) :: rarray

   open(unit=1,file="input",form="formatted",status="old",action="read")
   read(unit=1,fmt=*) seed
!   read(1,nml=input)
   close(unit=1)
!   write(*,nml=input)


!  Initialise at rest
   v1t = 0.0
   v3t = 0.0
   v1tm = v1t
   v3tm = v3t
   s1t = 0.0
   s3t = 0.0
   v1z = 0.
   v3z = 0.

   day = 0.0
   addnoise = .false.

   ! Time loop
   do 

      if ( day < day1 ) then
         zonal = .true. ! 130 day spin up using a one day time step
         dt = dt1
      else
         zonal = .false.
         dt = dt2
         diag_flag = .false.
         if ( .not. addnoise ) then
            call rluxgo(3,seed,0,0)
            call ranlux(rarray,nx*(ny-1))
            rand = reshape(rarray, (/ nx, ny-1 /) )
!            call random_seed()
!            call random_number(rand)
            ! Remove the zonal mean
            do j=1,ny-1
               rand(:,j) = rand(:,j) - sum(rand(:,j))/nx
            end do
            rand = 8.5e6*rand
            ! Define a random streamfunction anomaly
            s1t(:,1:ny-1) = s1t(:,1:ny-1) + rand
            s3t(:,1:ny-1) = s3t(:,1:ny-1) + rand
            ! Convert this to a vorticity anomaly
            call calcvor(s1t, s3t, v1t, v3t)
            addnoise = .true.
            ! Time step has changed so linearly interpolate the previous time
            ! value to ensure a smooth start
            v1tm = v1t - (v1t-v1tm)*dt2/dt1
            v3tm = v3t - (v3t-v3tm)*dt2/dt1
            diag_flag = .true.
         end if
         
      end if

      call split(v1t,v1z,v1)
      call split(v3t,v3z,v3)

!      print*, "V1", v1z(3), v1z(ny-3), v3z(3), v3z(ny-3)
      call calc_zonstream(v1z, v3z, s1z, s3z)

      !  Use relaxation to solve for the anomaly streamfunction
      if ( zonal ) then
         s1 = 0.0
         s3 = 0.0
      else
         call relax1(v1, v3, s1, s3)
      end if
      s1tm = s1t
      s3tm = s3t
      s1t = s1 + spread(s1z,dim=1,ncopies=nx)
      s3t = s3 + spread(s3z,dim=1,ncopies=nx)

      if ( abs ( modulo(day,diag_freq) - nint(modulo(day,diag_freq)) ) < 0.01 ) then
         call diag(day, zonal, s1t, s3t)
         if (.not. zonal ) then
            ! print*, 'Enter to continue'
            read(*,*) 
         end if
     end if
!      call diagvv(s1t, s3t, s1tm, s3tm, dt)

      ! Time stepping
      call xcalc(v1t, v3t, v1tm, v3tm, s1t, s3t, s1tm, s3tm, dt, x1t, x3t)
      call split(x1t,x1z,x1)
      call split(x3t,x3z,x3)

      !  Solve for new zonal mean vorticity from x
      if ( imp ) then
         call calc_zvor(x1z,x3z,dt,v1z,v3z)
      else
         v1z(1:ny-1) = x1z(1:ny-1)
         v1z(0)      = v1z(1)
         v1z(ny)     = v1z(ny-1)
         v3z(1:ny-1) = x3z(1:ny-1)
         v3z(0)      = v3z(1)
         v3z(ny)     = v3z(ny-1)
      end if

      if ( day == 0.0 ) then
         ! Forward step from rest. This simple form for the forward step
         ! assumes that it's starting from rest.
         v1z = 0.5*x1z
         v3z = 0.5*x3z
      else
         ! Update previous value of vorticity
         v1tm = v1t
         v3tm = v3t
      end if

      if ( zonal ) then
         v1 = 0.0
         v3 = 0.0
      else
         if ( imp ) then
            ! Relaxation solver for non-zonal terms
            call relax2(x1, x3, dt, v1, v3)
         else
            v1 = x1
            v3 = x3
         end if
      end if
!      print*, "Vanom", day, maxval(abs(v1))/maxval(abs(v1z))
      v1t = v1 + spread(v1z,dim=1,ncopies=nx)
      v3t = v3 + spread(v3z,dim=1,ncopies=nx)

      day = day + dt/86400.0
      if ( day > day2 ) then
         exit
      end if
   end do

end program phillips
