module constants
   implicit none
   ! Note that i runs 1..nx, j runs 0..ny
   integer, parameter, public :: nx=16, ny=16
   real, parameter, public :: dx=3.75e5, dy=6.25e5   ! Gridsize
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
   real, parameter, public :: k = 4.0e-6 ! s^{-1}  Surface drag
!  Control variables
   real, public, save :: a = 1.0e5  ! m^2/s   Horizonal diffusion
   logical, public, save :: imp = .true.
   real, public, save :: dt1 = 86400.0, dt2 = 3600.0
   real, public, save :: day1 = 131.0, day2 = 200.
   logical, public, save :: diag_flag = .false.
   real, parameter, public  :: accel = 1.0
   ! Tuned to give initial EKE ~ 768
   real, parameter :: noisescale = 1.237e7
end module constants

module work
   implicit none
   public :: split, calcvor, calc_zonstream,  relax1, xcalc, diag, diagvv, &
             calc_zvor, relax2
   integer :: ncid, irec
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
      ! Follow A17 and set end rows to zonal mean of neighbours
      v1(:,0) = sum(v1(:,1))/nx
      v1(:,ny) = sum(v1(:,ny-1))/nx
      v3(:,0) = sum(v3(:,1))/nx
      v3(:,ny) = sum(v3(:,ny-1))/nx
   end subroutine calcvor

   subroutine calc_zonstream(v1z, v3z, s1z, s3z)
!     Solve for the zonal mean streamfunction
      use lapack95
      use constants
      real, intent(in),  dimension(0:) :: v1z, v3z
      real, intent(out), dimension(0:) :: s1z, s3z
      integer, parameter :: nz = 2*ny - 3  ! Number of zonal means to solve for
      real, dimension(nz) :: bmat
      real, dimension(nz,nz) :: amat
      integer :: j
!     For LU solver
      integer, dimension(nz) :: ipiv
      logical, save :: first = .true.

      if (first) then

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
         call getrf(amat, ipiv)

         first = .false.
      end if

      do j=1,ny-1
         bmat(j) = v1z(j)
      end do
      do j=2,ny-1
         bmat(ny-2+j) = v3z(j)
      end do

      ! Solve AX=B
      call getrs(amat, ipiv, bmat)

      do j=1,ny-1
         s1z(j) = bmat(j)
      end do
      do j=2,ny-1
         s3z(j) = bmat(ny-2+j)
      end do

!     Apply the BC
      s1z(0) = s1z(1)
      s1z(ny) = s1z(ny-1)
!     Note the extra restriction on s3(1) which is not solved for.
      s3z(0:1) = 0.0
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
      real :: change1, change3

      ! Start from the current value of the anomaly streamfunction
      do iter=1,100
         ! Jacobi iteration
         maxdiff = 0.0
         change = 0.0
         ! do irb=0,1
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
            ! print*, "ITER1", iter, sqrt(change), maxdiff
         ! end do
         ! maxdiff is now only on a single level so halve the convergence
         ! criterion.
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
         ! Boundary condition A17
         v1(:,0) = sum(v1(:,1))/nx
         v1(:,ny) = sum(v1(:,ny-1))/nx
         change3 = sum ( ( v3(:,1:ny-1) - temp3(:,1:ny-1) ) **2 )
         v3(:,1:ny-1) = temp3(:,1:ny-1)
         v3(:,0) = sum(v3(:,1))/nx
         v3(:,ny) = sum(v3(:,ny-1))/nx
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
      use lapack95
      real, intent(in),  dimension(0:) :: x1z, x3z
      real, intent(in) :: dt
      real, intent(out), dimension(0:) :: v1z, v3z
      integer, parameter :: nz = ny-1  ! Number of zonal means to solve for
      real, dimension(nz) :: amat, bmat, cmat, rmat
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

      call gtsv(amat(2:), bmat, cmat(:nz-1), rmat)

      v1z(1:nz) = rmat
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

      call gtsv(amat(2:), bmat, cmat(:nz-1), rmat)

      v3z(1:nz) = rmat
      v3z(0) = v3z(1)
      v3z(ny) = v3z(ny-1)

   end subroutine calc_zvor

   subroutine zonal_dump(day, v1z, v3z)
      ! Write zonal diagnostics to match Phillips Table 1
      use constants
      real, intent(in) :: day
      real, intent(in), dimension(0:) :: v1z, v3z
      real, dimension(0:ny) :: s1z, s3z, vor1z
      real, dimension(0:ny) :: t2z
      real, dimension(1:ny) :: u1z, u3z
      real ::  u1_t, u2_t, u3_t, u4_t
      real :: tke, zke, eke, zpe, epe
      integer :: j
      character(len=20) :: fname

      ! print*, 'ZONAL V1', v1z(1)
      call calc_zonstream(v1z, v3z, s1z, s3z)
      call zonal_diag(day, s1z, s3z)
      ! Vorticity from eq 22
      ! vor1 = q1 + lambdasq(s1 - s2)
      ! v1 = eta = q dx**2
      vor1z = v1z/dx**2 + lambdasq*(s1z-s3z)

      t2z = f0*(s1z(:)-s3z(:))/rgas
      ! Zonal wind, u = -pd{psi}{y}
      u1z = - ( s1z(1:ny) - s1z(0:ny-1) ) / dy
      u3z = - ( s3z(1:ny) - s3z(0:ny-1) ) / dy
      do j=1,ny-1
         ! Interpolate to the temperature points
         u1_t = 0.5*(u1z(j) + u1z(j+1))
         u3_t = 0.5*(u3z(j) + u3z(j+1))
         u2_t = 0.5*(u1_t + u3_t)
         ! Extrapolation
         u4_t = u3_t + 0.5*(u3_t - u1_t)
         write(unit=*,fmt="(i3,4f10.2,2e14.6)") j, t2z(j), u1_t, u2_t, u4_t, v1z(j), s1z(j)! , 1e4*vor1z(j)
      end do   
   end subroutine zonal_dump

   subroutine zonal_diag(day, s1z, s3z)
      use constants
      real, intent(in) :: day
      real, intent(in), dimension(0:) :: s1z, s3z
      real, dimension(0:ny) :: t2z
      real, dimension(1:ny) :: u1z, u3z

      real :: zke, zpe

      t2z = f0*(s1z(:)-s3z(:))/rgas
      ! Zonal wind, u = -pd{psi}{y}
      u1z = - ( s1z(1:ny) - s1z(0:ny-1) ) / dy
      u3z = - ( s3z(1:ny) - s3z(0:ny-1) ) / dy

      ! Zonal KE, scaled such that a wind of 1 m/s everywhere gives 10
      zke = 10.0*sum(u1z*u1z + u3z*u3z)/(2*ny)

      ! Zonal PE  A21
      ! Sum 1:J-1 / J matches definition of Y operator
      zpe = 5*lambdasq*sum((s1z(1:ny-1) - s3z(1:ny-1))**2) / ny

      write(unit=*,fmt="(a,f6.2,5f15.7)") "TZ ", day, maxval(t2z), maxval(u1z), maxval(u3z), zke, zpe

   end subroutine zonal_diag

   subroutine diag(day, s1, s3)
      use constants
      real, intent(in) :: day
      real, intent(in), dimension(:,0:) :: s1, s3
      real, dimension(0:ny) :: t2z
      real, dimension(nx,0:ny) :: t2
      real, dimension(1:ny) :: u1z, u3z, v1z, v3z
      real, dimension(nx,1:ny) :: u1t, u3t, u1, u3, v1, v3
      real, dimension(nx,0:ny) :: v1t, v3t
      real :: tke, zke, eke, zpe, epe
      integer :: j, jj
      character(len=20) :: fname
      integer, save :: dcount = 1

      u1t = - ( s1(:,1:ny) - s1(:,0:ny-1) ) / dy
      u3t = - ( s3(:,1:ny) - s3(:,0:ny-1) ) / dy
      v1t = ( s1 - cshift(s1,-1) ) / dx
      v3t = ( s3 - cshift(s3,-1) ) / dx

      ! Average v to get it on the same grid points as u
      do j=ny,1,-1  ! Reverse loop so don't clobber unused values
         v1t(:,j) = 0.25*(v1t(:,j) + cshift(v1t(:,j),1) + &
               v1t(:,j-1) + cshift(v1t(:,j-1),1) )
         v3t(:,j) = 0.25*(v3t(:,j) + cshift(v3t(:,j),1) + &
               v3t(:,j-1) + cshift(v3t(:,j-1),1) )
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
      write(unit=*,fmt="(a,f6.2,6f15.7)") "KE ", day, maxval(t2), maxval(u1t), maxval(u3t), zke, eke, tke
      if (eke > 1.e5) then
         stop "EKE too large"
      end if

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

   subroutine create_nc_output()
      use netcdf
      use constants
      integer :: ierr, londim, latdim, levdim, timedim, &
                 lonid, latid, levid, timeid, varid, &
                 vlondim, vlonid, ulatdim, ulatid
      real, dimension(nx) :: glon
      real, dimension(0:ny) :: glat
      integer :: i, j, vid

      ierr = nf90_create ("phillips.nc", NF90_CLOBBER, ncid)
      call check_ncerr(ierr, "Error creating file ")
   
      !  Create dimensions, lon, lat and rec for the output file
      ierr = nf90_def_dim(ncid, 'lon', nx, londim)
      call check_ncerr(ierr,"Error creating londim")
      ierr = nf90_def_dim(ncid, 'lat', ny+1, latdim)
      call check_ncerr(ierr,"Error creating latdim")
      ierr = nf90_def_dim(ncid, 'vlon', nx, vlondim)
      call check_ncerr(ierr,"Error creating vlondim")
      ierr = nf90_def_dim(ncid, 'ulat', ny, ulatdim)
      call check_ncerr(ierr,"Error creating ulatdim")
      ierr = nf90_def_dim(ncid, 'lev', 2, levdim)
      call check_ncerr(ierr,"Error creating levdim")
      ierr = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, timedim)
      call check_ncerr(ierr, 'error creating time dimension')
   
      ierr = nf90_def_var(ncid, 'lon',  NF90_FLOAT, londim, lonid)
      ierr = nf90_put_att(ncid, lonid, "long_name", "longitude")
      ierr = nf90_put_att(ncid, lonid, "units", "degrees_east")

      ierr = nf90_def_var(ncid, 'vlon',  NF90_FLOAT, vlondim, vlonid)
      ierr = nf90_put_att(ncid, vlonid, "long_name", "V longitude")
      ierr = nf90_put_att(ncid, vlonid, "units", "degrees_east")

      ierr = nf90_def_var (ncid, 'lat',  NF90_FLOAT, latdim, latid)
      ierr = nf90_put_att (ncid, latid, "long_name", "latitude" )
      ierr = nf90_put_att (ncid, latid, "units", "degrees_north")

      ierr = nf90_def_var (ncid, 'ulat',  NF90_FLOAT, ulatdim, ulatid)
      ierr = nf90_put_att (ncid, ulatid, "long_name", "U latitude" )
      ierr = nf90_put_att (ncid, ulatid, "units", "degrees_north")

      ierr = nf90_def_var (ncid, 'lev',  NF90_FLOAT, levdim, levid)
      ierr = nf90_put_att (ncid, levid, "long_name", "model level" )

      ierr = nf90_def_var (ncid, 'time',  NF90_FLOAT, timedim, timeid)
      ierr = nf90_put_att (ncid, timeid, "units", "days since 2000-01-01 00:00" )

      ierr = nf90_def_var (ncid, "strm",  NF90_FLOAT, (/londim, latdim, levdim, timedim /), vid)
      ierr = nf90_put_att (ncid, vid, "long_name", "streamfunction")
      ierr = nf90_put_att (ncid, vid, "units", "m2 s-1" )

      ierr = nf90_def_var (ncid, "vor",  NF90_FLOAT, (/londim, latdim, levdim, timedim /), vid)
      ierr = nf90_put_att (ncid, vid, "long_name", "vorticity")
      ierr = nf90_put_att (ncid, vid, "units", "s-1" )

      ierr = nf90_def_var (ncid, "t500",  NF90_FLOAT, (/londim, latdim, levdim, timedim /), vid)
      ierr = nf90_put_att (ncid, vid, "long_name", "air temperature at 500 hPa")
      ierr = nf90_put_att (ncid, vid, "units", "K" )

      ierr = nf90_def_var (ncid, "ps",  NF90_FLOAT, (/londim, latdim, levdim, timedim /), vid)
      ierr = nf90_put_att (ncid, vid, "long_name", "surface pressure")
      ierr = nf90_put_att (ncid, vid, "units", "hPa" )

      ierr = nf90_def_var (ncid, "u",  NF90_FLOAT, (/londim, ulatdim, levdim, timedim /), vid)
      ierr = nf90_put_att (ncid, vid, "long_name", "zonal wind")
      ierr = nf90_put_att (ncid, vid, "units", "m s-1" )

      ierr = nf90_def_var (ncid, "v",  NF90_FLOAT, (/vlondim, latdim, levdim, timedim /), vid)
      ierr = nf90_put_att (ncid, vid, "long_name", "meriodional wind")
      ierr = nf90_put_att (ncid, vid, "units", "m s-1" )

      ierr = nf90_enddef(ncid)

      do i=1,nx
         glon(i) = i
      end do
      ierr = nf90_put_var(ncid, lonid, glon)
      call check_ncerr(ierr)

      do i=1,nx
         glon(i) = i - 0.5
      end do
      ierr = nf90_put_var(ncid, vlonid, glon)
      call check_ncerr(ierr)

      do j=0,ny
         glat(j) = j
      end do
      ierr = nf90_put_var(ncid, latid, glat)
      call check_ncerr(ierr)

      do j=1,ny
         glat(j) = j - 0.5
      end do
      ierr = nf90_put_var(ncid, ulatid, glat(1:))
      call check_ncerr(ierr)

      ierr = nf90_put_var(ncid, levid, [1., 3.])
      call check_ncerr(ierr)

      irec = 0

   end subroutine create_nc_output

   subroutine nc_output(day, v1t, v3t, s1t, s3t)
      use netcdf
      use constants
      real, intent(in) :: day
      real, intent(in), dimension(:,0:) :: v1t, v3t, s1t, s3t
      real, dimension(nx,1:ny) :: u
      real, dimension(nx,0:ny) :: v
      real, dimension(nx,0:ny) :: t2, ps
      integer :: ierr, vid

      irec = irec + 1
      ierr = nf90_inq_varid(ncid, 'vor', vid)
      call check_ncerr(ierr, "Error getting vid for vor")
      ierr = nf90_put_var(ncid, vid, v1t, start=[1, 1, 1, irec], &
                          count=[nx, ny+1, 1, 1])
      ierr = nf90_put_var(ncid, vid, v3t, start=[1, 1, 2, irec], &
                          count=[nx, ny+1, 1, 1])

      ierr = nf90_inq_varid(ncid, 'strm', vid)
      call check_ncerr(ierr, "Error getting vid for strm")
      ierr = nf90_put_var(ncid, vid, s1t, start=[1, 1, 1, irec], &
                          count=[nx, ny+1, 1, 1])
      ierr = nf90_put_var(ncid, vid, s3t, start=[1, 1, 2, irec], &
                          count=[nx, ny+1, 1, 1])

      ierr = nf90_inq_varid(ncid, 'u', vid)
      call check_ncerr(ierr, "Error getting vid for u")
      u = - ( s1t(:,1:ny) - s1t(:,0:ny-1) ) / dy
      ierr = nf90_put_var(ncid, vid, u, start=[1, 1, 1, irec], &
                          count=[nx, ny, 1, 1])
      u = - ( s3t(:,1:ny) - s3t(:,0:ny-1) ) / dy
      ierr = nf90_put_var(ncid, vid, u, start=[1, 1, 2, irec], &
                          count=[nx, ny, 1, 1])

      ierr = nf90_inq_varid(ncid, 'v', vid)
      call check_ncerr(ierr, "Error getting vid for v")
      v = ( s1t - cshift(s1t,-1) ) / dx
      ierr = nf90_put_var(ncid, vid, v, start=[1, 1, 1, irec], &
                          count=[nx, ny+1, 1, 1])
      v = ( s3t - cshift(s3t,-1) ) / dx
      ierr = nf90_put_var(ncid, vid, v, start=[1, 1, 2, irec], &
                          count=[nx, ny+1, 1, 1])

      t2 = f0*(s1t-s3t)/rgas
      ierr = nf90_inq_varid(ncid, 't500', vid)
      call check_ncerr(ierr, "Error getting vid for t500")
      ierr = nf90_put_var(ncid, vid, t2, start=[1, 1, irec], &
                          count=[nx, ny, 1])

      ps = 0.01 * (1.5*s3t - 0.5*s1t)*f0
      ierr = nf90_inq_varid(ncid, 'ps', vid)
      call check_ncerr(ierr, "Error getting vid for ps")
      ierr = nf90_put_var(ncid, vid, ps, start=[1, 1, irec], &
                          count=[nx, ny, 1])

      ierr = nf90_inq_varid(ncid, 'time', vid)
      call check_ncerr(ierr, "Error getting vid for time")
      ierr = nf90_put_var(ncid, vid, day, start=[irec])

      ierr = nf90_sync(ncid)

   end subroutine nc_output

   subroutine check_ncerr(status,mesg)
      use netcdf
      integer, intent(in) :: status
      character(len=*), intent(in), optional :: mesg
      if ( status /= nf90_noerr ) then
         if ( present(mesg) ) then
            print*, mesg
         end if
         print*, trim(nf90_strerror(status))
         stop
      end if
   end subroutine check_ncerr
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
   use msq_rand_mod
   implicit none
   integer, parameter :: diag_freq = 3600
   real, dimension(nx,0:ny) :: s1t, s3t   ! Streamfunction at time tau
   real, dimension(nx,0:ny) :: s1tm, s3tm ! Streamfunction at tau-1
   real, dimension(nx,0:ny) :: s1, s3     ! Streamfunction anomalies
   real, dimension(0:ny)    :: s1z, s3z   ! Zonal mean streamfunction
   real, dimension(nx,0:ny) :: v1t, v3t   ! Vorticity
   real, dimension(nx,0:ny) :: v1tmp, v3tmp   ! Vorticity
   real, dimension(nx,0:ny) :: v1, v3     ! Vorticity anomalies
   real, dimension(0:ny)    :: v1z, v3z   ! Zonal mean vorticity
   real, dimension(nx,0:ny) :: v1tm, v3tm ! tau - 1 values
!  Temporaries used in timestepping
   real, dimension(nx,0:ny) :: x1t, x3t
   real, dimension(nx,0:ny) :: x1, x3
   real, dimension(0:ny)    :: x1z, x3z

   integer :: i, j, im, ip, jm, jp, seed
   integer :: time
   real, parameter :: daylen = 86400.
   real    :: istep, day, swap
   logical :: zonal, addnoise
   integer :: seedsize
   integer, allocatable, dimension(:) :: seed_array
   integer(i8) :: rval
!   namelist /input/ a, implicit, dt1, dt2, day1, day2, accel

   real :: dt  ! Time step
   real, dimension(nx,0:ny) :: rand

!    open(unit=1,file="input",form="formatted",status="old",action="read")
!    read(unit=1,fmt=*) seed
! !   read(1,nml=input)
!    close(unit=1)
! !   write(*,nml=input)

!  Initialise at rest
   v1t = 0.0
   v3t = 0.0
   v1tm = v1t
   v3tm = v3t
   s1t = 0.0
   s3t = 0.0
   v1z = 0.
   v3z = 0.

   time = 0
   day = 0.0
   addnoise = .false.

   call create_nc_output()

   ! Time loop
   do

      if ( day < day1 ) then
         zonal = .true. ! 130 day spin up using a one day time step
         dt = dt1
      else
         ! End of the spin up section
         if (zonal) then
            call zonal_dump(day, v1z, v3z)
         end if
         zonal = .false.
         dt = dt2
         diag_flag = .false.
         if ( .not. addnoise ) then
            ! Time step has changed so linearly interpolate the previous time
            ! value to ensure a smooth start
            v1tm = v1t - (v1t-v1tm)*dt2/dt1
            v3tm = v3t - (v3t-v3tm)*dt2/dt1
            ! call random_seed(size=seedsize)
            ! allocate(seed_array(seedsize))
            ! seed_array = seed
            ! call random_seed(put=seed_array)
            ! call random_number(rand(:,1:ny-1))
            rval = 1111111111_i8
            do j=1,ny-1
               do i=1,nx
                  rval = msq_rand(rval)
                  rand(i,j) = real(rval) / ten10
               end do
            end do
            rand(:,0) = 0.
            rand(:,ny) = 0.
            ! Remove the zonal mean
            do j=1,ny-1
               rand(:,j) = rand(:,j) - sum(rand(:,j))/nx
            end do
            rand = noisescale * rand
            ! Convert this to a vorticity anomaly and add to current 
            ! and previous values
            call calcvor(rand, rand, v1tmp, v3tmp)
            ! print*, 'V1TMP', maxval(abs(v1tmp))
            v1t = v1t + v1tmp
            v3t = v3t + v3tmp
            v1tm = v1tm + v1tmp
            v3tm = v3tm + v3tmp
            addnoise = .true.
            diag_flag = .true.
         end if

      end if

      call split(v1t,v1z,v1)
      call split(v3t,v3z,v3)

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

      if (modulo(time,diag_freq) == 0) then
         if (zonal) then
            call zonal_diag(day, s1z, s3z)
         else
            call diag(day, s1t, s3t)
            call nc_output(day, v1t, v3t, s1t, s3t)
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

      if ( time == 0 ) then
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
      v1t = v1 + spread(v1z,dim=1,ncopies=nx)
      v3t = v3 + spread(v3z,dim=1,ncopies=nx)

      time = time + dt
      day = time/daylen
      if ( day > day2 ) then
         exit
      end if
   end do

end program phillips
