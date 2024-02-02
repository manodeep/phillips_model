# Phillip's model
# Direct translation of fortran code, using fortran arrays starting at one.

import numpy as np
from numpy import linalg
import ranlux

class Grid:
    nx = 16
    ny = 16

class Var():

    # Level 1 and 3 components of 3D variable

    def __init__(self):
        # Total field
        self.l1t = np.zeros( (Grid.nx+1,Grid.ny+1) )
        self.l3t = np.zeros( (Grid.nx+1,Grid.ny+1) )
        # Anomaly (zonal mean removed)
        self.l1  = np.zeros( (Grid.nx+1,Grid.ny+1) )
        self.l3  = np.zeros( (Grid.nx+1,Grid.ny+1) )
        # Zonal means
        self.l1z  = np.zeros( Grid.ny+1 )
        self.l3z  = np.zeros( Grid.ny+1 )

    def dump(self):
        for j in range(0,Grid.ny+1):
            for i in range(1,Grid.nx+1):
		print "%12.4f " % self.l1t[i,j],
            print

    def adump(self):
        for j in range(0,Grid.ny+1):
            for i in range(1,Grid.nx+1):
		print "%12.4f " % self.l1[i,j],
            print

    def settot(self,val):
        if isinstance(val,Var):
            self.l1t[:] = val.l1t
            self.l3t[:] = val.l3t
        else:
            self.l1t[:] = val
            self.l3t[:] = val

    def set(self,val):
        if isinstance(val,Var):
            self.l1[:] = val.l1
            self.l3[:] = val.l3
        else:
            self.l1[:] = val
            self.l3[:] = val

    def calc_zmean(self):
        self.l1z[:] = self.l1t[1:Grid.nx+1,:].mean(axis=0)
        self.l3z[:] = self.l3t[1:Grid.nx+1,:].mean(axis=0)

    def split(self):
        self.calc_zmean()
        self.l1[:] = self.l1t[:] - self.l1z[:]
        self.l3[:] = self.l3t[:] - self.l3z[:]


class Model:
    dx=3.75e5; dy=6.25e5    # Gridsize
    eps=dx/dy; epsq = eps*eps
    heat = 2.0e-3           #  Heating W/kg
    lambdasq = 1.5e-12      # Stability parameter, m^-2
    rgas = 287.0
    cp = 1004.0
    f0 = 1.0e-4
    beta = 1.6e-11          # s^{-1} m^{-1}
    p4 = 1.0e5              # Surface pressure, Pa
    p2 = 0.5*p4
    gamma = lambdasq*dx*dx  # 0.21
    k = 4.0e-6              # s^{-1}  Surface drag

    # Control variables
    a = 1.0e5               # m^2/s   Horizonal diffusion
    imp = True
    diag_flag = False
    accel = 1.0


    def calcvor(self, s, v):
        # This repeats same code for each level. Should the level be another dimension?
        for j in range(1,Grid.ny):
 	    jm = j-1
 	    jp = j+1
            for i in range(1,Grid.nx+1):
                im = i-1
                if im == 0:
                    im = Grid.nx
                ip = i+1;
                if ip == Grid.nx+1:
 		    ip = 1
 		v.l1t[i,j] = ( s.l1t[ip,j]  + s.l1t[im,j] - 2*s.l1t[i,j] ) + \
                    self.epsq * ( s.l1t[i,jp] + s.l1t[i,jm] - 2*s.l1t[i,j] ) - \
                    self.gamma * ( s.l1t[i,j] - s.l3t[i,j] )
 		v.l3t[i,j] = ( s.l3t[ip,j] + s.l3t[im,j] - 2*s.l3t[i,j] ) + \
                    self.epsq * ( s.l3t[i,jp] + s.l3t[i,jm] - 2*s.l3t[i,j] ) + \
                    self.gamma * ( s.l1t[i,j] - s.l3t[i,j] )

    def calc_zonstream(self, v, s):
        # Given vorticity variable as input, solve for the 
        # zonal mean streamfunction

        ny = Grid.ny
        nz = 2*ny - 3  #  Number of zonal means to solve for
        epsq = self.epsq
        gamma = self.gamma

	# Start these arrays from 1 again
        amat = np.zeros((nz+1,nz+1))
        bmat = np.zeros(nz+1)

        # j = 1, level 1
        amat[1,1] = -epsq - gamma
        amat[1,2] = epsq
        # j = J-1, level 1
        amat[ny-1,ny-2] = epsq
        amat[ny-1,ny-1] = -epsq - gamma
        amat[ny-1,nz] = gamma
        # j = 2, level 3
        amat[ny,2] = gamma
        amat[ny,ny] = -2.0*epsq - gamma
        amat[ny,ny+1] = epsq
        # j = J-1, level 3
        amat[nz,ny-1] = gamma
        amat[nz,nz-1] = epsq
        amat[nz,nz] = -epsq - gamma

        #  Level 1
        for j in range(2,ny-1):
            amat[j,j-1] = epsq
            amat[j,j] = -2.0*epsq - gamma
            amat[j,j+1] = epsq
            amat[j,ny-2+j] = gamma
            
        #  Level 3
            for j in range(ny+1, nz):
                amat[j,j+2-ny] = gamma
                amat[j,j-1] = epsq
                amat[j,j] = -2.0*epsq - gamma
                amat[j,j+1] = epsq

#       print*, "AMAT"
#       do j=1,nz
#          write(*,"(7f6.2)") amat(j,:)
#       end do

#       do j=1,ny-1
#          junk(j) = s1z(j)
#       end do
#       do j=2,ny-1
#          junk(ny-2+j) = s3z(j)
#       end do

#      junk = matmul(amat,junk)

        bmat[1:ny] = v.l1z[1:ny]
        for j in range(2,ny):
            bmat[ny-2+j] = v.l3z[j]

        # Solve AX=B
#         print "AMAT"
#         print amat
#         print "BMAT"
#         print bmat
#        print "B shape", bmat.shape
        bmat[1:] = linalg.solve(amat[1:,1:],bmat[1:])

        for j in range(1,ny):
            s.l1z[j] = bmat[j]
        for j in range(2,ny):
            s.l3z[j] = bmat[ny-2+j]

        # Apply the BC
        s.l1z[0] = s.l1z[1]
        s.l1z[ny] = s.l1z[ny-1]
        s.l3z[0] = 0.0
        s.l3z[ny] = s.l3z[ny-1]

    def relax1(self,v, s):
        # Solve for anomaly streamfunction

        nx = Grid.nx
        ny = Grid.ny
        # Start from the current value of the anomaly streamfunction
        for iter in range(100):
            # Jacobi iteration
            maxdiff = 0.0
            change = 0.0
            for irb in range(2):
                for j in range(1,ny):
                    jm = j-1
                    jp = j+1
                    for i in range(1,nx+1):
                        im = i-1
                        if im == 0:
                            im = nx
                        ip = i+1
                        if ip == nx+1:
                            ip = 1

                        resid = ( s.l1[ip,j] + s.l1[im,j] + 
                                  self.epsq*( s.l1[i,jp] + s.l1[i,jm] ) - 
                                  v.l1[i,j] + self.gamma*s.l3[i,j] ) -  \
                                  ( 2.0 + 2.0*self.epsq + self.gamma )*s.l1[i,j]
                        resid = self.accel*resid / ( 2.0 + 2.0*self.epsq + self.gamma )
                        change = change + resid**2 
                        maxdiff = max ( maxdiff, abs(resid) )
                        s.l1[i,j] = s.l1[i,j] + resid

                        resid = ( s.l3[ip,j] + s.l3[im,j] + 
                                  self.epsq*( s.l3[i,jp] + s.l3[i,jm] ) - 
                                  v.l3[i,j] + self.gamma*s.l1[i,j] ) -  \
                                  ( 2.0 + 2.0*self.epsq + self.gamma )*s.l3[i,j]

                        resid = self.accel*resid / ( 2.0 + 2.0*self.epsq + self.gamma )
                        change = change + resid**2
                        maxdiff = max ( maxdiff, abs(resid) )
                        s.l3[i,j] = s.l3[i,j] + resid
            # print "ITER1", iter, np.sqrt(change), maxdiff
            # maxdiff is now only on a single level so halve the convergence 
            # criterion
            if maxdiff < 0.5*3.75e4:
                # print "ITER1 done", iter, np.sqrt(change), maxdiff
                break
        
        if iter >= 100:
            print "RELAX1 failed", iter, np.sqrt(change), maxdiff

    def relax2(self, x, dt, v):
        # Solve for anomaly vorticity
        # real, dimension(nx,0:ny) :: temp1, temp3

        temp = Var()
        nx = Grid.nx
        ny = Grid.ny

        v.l1[:] = 0.0
        v.l3[:] = 0.0
        alpha = self.a*dt/self.dx**2
        for iter in range(100):
            # Jacobi iteration
            for j in range(1,ny):
                jm = j-1
                jp = j+1
                for i in range(1,nx+1):
                    im = i-1
                    if im == 0:
                        im = nx
                    ip = i+1
                    if ip == nx+1:
                        ip = 1
                    temp.l1[i,j] = ( alpha*( v.l1[ip,j] + v.l1[im,j] + 
                                             self.epsq * ( v.l1[i,jp] + v.l1[i,jm] ) ) + 
                                     x.l1[i,j]  ) /    \
                                     ( 2*alpha*(1.0 + self.epsq)  + 1.0 )
                    temp.l3[i,j] = ( alpha*( v.l3[ip,j] + v.l3[im,j] + 
                                             self.epsq * ( v.l3[i,jp] + v.l3[i,jm] ) ) + 
                                     x.l3[i,j]  ) /    \
                                     ( 2*alpha*(1.0 + self.epsq)  + 1.0 +1.5*self.k*dt )
                change1 = np.sum ( ( v.l1[1:,1:] - temp.l1[1:,1:] ) **2 )
                v.l1[:,1:ny] = temp.l1[:,1:ny]
                v.l1[:,0] = v.l1[:,1]
                v.l1[:,ny] = v.l1[:,ny-1]
                change3 = np.sum ( ( v.l3[1:,1:] - temp.l3[1:,1:] ) **2 )
                v.l3[:,1:ny] = temp.l3[:,1:ny]
                v.l3[:,0] = v.l3[:,1]
                v.l3[:,ny] = v.l3[:,ny-1]
                if max(change1, change3) < 1.0:
                    print "ITER2", iter, np.sqrt(change1), np.sqrt(change3)
                    break

    def xcalc(self, v, vm, s, sm, dt, x):

        nx = Grid.nx
        ny = Grid.ny
        alpha = self.a*dt/self.dx**2
        if not self.imp:
            alpha = 2.0 * alpha
            
        b = self.beta*self.dx**2*self.dy
        c = dt/(2.0*self.dx*self.dy)
        h = 4*self.rgas*self.heat*self.gamma*dt/(self.f0*self.cp)
        x.l1t[:] = 0.0
        x.l3t[:] = 0.0
        for j in range(1,ny):
            jm = j-1
            jp = j+1
            for i in range(1,nx+1):
                im = i-1
                if im == 0:
                    im = nx
                ip = i+1
                if ip == nx+1:
                    ip = 1

                x.l1t[i,j] = vm.l1t[i,j] +                                           \
                    c * ( (v.l1t[ip,j]-v.l1t[im,j])*(s.l1t[i,jp]-s.l1t[i,jm]) -             \
                          (2*b+v.l1t[i,jp]-v.l1t[i,jm])*(s.l1t[ip,j]-s.l1t[im,j]) ) +      \
                          alpha * ( vm.l1t[ip,j]+vm.l1t[im,j]-2*vm.l1t[i,j] +           \
                                        self.epsq*(vm.l1t[i,jp]+vm.l1t[i,jm]-2*vm.l1t[i,j]) ) +        \
                                        h*(2*j-ny)/ny

                x.l3t[i,j] = vm.l3t[i,j] +                                             \
                    c * ( (v.l3t[ip,j]-v.l3t[im,j])*(s.l3t[i,jp]-s.l3t[i,jm]) -              \
                           (2*b+v.l3t[i,jp]-v.l3t[i,jm])*(s.l3t[ip,j]-s.l3t[im,j]) ) +       \
                           alpha * ( vm.l3t[ip,j]+vm.l3t[im,j]-2*vm.l3t[i,j] +            \
                                         self.epsq*(vm.l3t[i,jp]+vm.l3t[i,jm]-2*vm.l3t[i,j]) ) -  \
                                         h*(2*j-ny)/ny
                if self.imp:
                    x.l3t[i,j] = x.l3t[i,j] -  self.k*dt*(1.5*vm.l3t[i,j] - v.l1t[i,j] - 
                                               4*self.gamma*(s.l1t[i,j]-s.l3t[i,j]) )
                else:
                    x.l3t[i,j] = x.l3t[i,j] -  self.k*dt*(3.0*vm.l3t[i,j] - vm.l1t[i,j] - 
                                               4*self.gamma*(sm.l1t[i,j]-sm.l3t[i,j]) )

    def tridag(self,A,B,C,R,U,N):
        gam = np.zeros(N+1)
        if B[1] == 0.:
            raise Exception("TRIDAG error")
        bet = B[1]
        U[1] = R[1]/bet
        for j in range(2,N+1):
            gam[j] = C[j-1]/bet
            bet = B[j] - A[j]*gam[j]
            if bet == 0.:
                raise Exception("TRIDAG error")
            U[j] = (R[j]-A[j]*U[j-1])/bet
        for j in range(N-1,0,-1):
            U[j] = U[j]-gam[j+1]*U[j+1]

    def calc_zvor(self, x, dt, v):
        # Solve for the zonal mean vorticity
        ny = Grid.ny
        epsq = self.epsq

        nz = ny-1
        amat = np.zeros(nz+1)
        bmat = np.zeros(nz+1)
        cmat = np.zeros(nz+1)
        rmat = np.zeros(nz+1)
        umat = np.zeros(nz+1)

        alpha = self.a*dt/self.dx**2
        # Level 1
        amat[1] = 0.0
        bmat[1] = -alpha*epsq - 1
        cmat[1] = alpha*epsq
        amat[2:ny-1] = alpha*epsq
        bmat[2:ny-1] = -2.0*alpha*epsq - 1
        cmat[2:ny-1] = alpha*epsq
        amat[nz] = alpha*epsq
        bmat[nz] = -alpha*epsq - 1
        cmat[nz] = 0.0

        rmat[1:] = -x.l1z[1:nz+1]

#         print "CALC_ZVOR A", amat
#         print "CALC_ZVOR B" , bmat
#         print "CALC_ZVOR C" , cmat
#         print "CALC_ZVOR R", rmat

        self.tridag(amat,bmat,cmat,rmat,umat,nz)
#        print "CALC_ZVOR U", umat

        v.l1z[1:nz+1] = umat[1:nz+1]
        v.l1z[0] = v.l1z[1]
        v.l1z[ny] = v.l1z[ny-1]

        # Level 3
        amat[1] = 0.0
        bmat[1] = -alpha*epsq - 1 - 1.5*self.k*dt
        cmat[1] = alpha*epsq
        amat[2:ny-1] = alpha*epsq
        bmat[2:ny-1] = -2.0*alpha*epsq - 1 - 1.5*self.k*dt
        cmat[2:ny-1] = alpha*epsq
        amat[nz] = alpha*epsq
        bmat[nz] = -alpha*epsq - 1 - 1.5*self.k*dt
        cmat[nz] = 0.0

        rmat[1:] = -x.l3z[1:nz+1]

        self.tridag(amat,bmat,cmat,rmat,umat,nz)
        
        v.l3z[1:nz+1] = umat[1:nz+1]
        v.l3z[0] = v.l3z[1]
        v.l3z[ny] = v.l3z[ny-1]

        # end subroutine calc_zvor
   
    def diag(self, day, zonal, s):
#         real, dimension(0:ny) :: t2z
#         real, dimension(nx,0:ny) :: t2
#         real, dimension(1:ny) :: u1z, u3z, v1z, v3z
#         real, dimension(nx,1:ny) :: u1t, u3t, u1, u3, v1, v3
#         real, dimension(nx,0:ny) :: v1t, v3t
#         real :: tke, zke, eke
#         integer :: j, jj
        
        u = Var()
        v = Var()
        nx = Grid.nx; ny = Grid.ny
        dx = self.dx; dy = self.dy

#         print "S1", s.l1t[1,:]
#         print "S1Z", s.l1z[:]

        if zonal:
            # Why not use zonal components here?
            t2z = self.f0*(s.l1t[1,:]-s.l3t[1,:])/self.rgas
            #  Zonal wind, u = -pd{psi}{y}
            u.l1z[1:] = - ( s.l1t[1,1:ny+1] - s.l1t[1,0:ny] ) / dy
            u.l3z[1:] = - ( s.l3t[1,1:ny+1] - s.l3t[1,0:ny] ) / dy

            # Zonal KE, scaled such that a wind of 1 m/s everywhere gives 10
            zke = 10.0*np.sum(u.l1z*u.l1z + u.l3z*u.l3z)/(2*ny)
            print "TZ %5.1f %15.7f %15.7f %15.7f %15.7f" % ( day, t2z.max(), u.l1z.max(), u.l3z.max(), zke)
        else:
            # write(unit=3,fmt="(f8.3,2e15.7)") day, s.l1(1,8), s.l3(1,8)
            if abs(day - round(day,0)) < 0.01 or self.diag_flag:
                u.l1t[:,1:] = - ( s.l1t[:,1:ny+1] - s.l1t[:,0:ny] ) / dy
                u.l3t[:,1:] = - ( s.l3t[:,1:ny+1] - s.l3t[:,0:ny] ) / dy
                for i in range(1,nx+1):
                    im = i-1
                    if im == 0:
                        ix = nx
                    v.l1t[i,:] = ( s.l1t[i,:] - s.l1t[im,:] ) / dx
                    v.l3t[i,:] = ( s.l3t[i,:] - s.l3t[im,:] ) / dx
                #  Average v to get it on the same grid points as u
                for i in range(1,nx+1):
                    ip = i+1
                    if ip > Grid.nx:
                        ip = 1
                        for j in range(Grid.ny, 0, -1):
                            # Reverse loop so don't clobber unused values
                            v.l1t[i][j] = 0.25*(v.l1t[i][j] + v.l1t[ip][j] + 
                                                v.l1t[i][j-1] + v.l1t[ip][j-1])

                # Calculate zonal mean and eddy winds
                u.split()
                v.split()

                tke = 10.0*np.sum(u.l1t[1:,1:]*u.l1t[1:,1:] + u.l3t[1:,1:]*u.l3t[1:,1:] + v.l1t[1:,1:]*v.l1t[1:,1:] + v.l3t[1:,1:]*v.l3t[1:,1:])/(2*ny*nx)

                zke = 10.0*np.sum(u.l1z*u.l1z + u.l3z*u.l3z)/(2*ny)
                eke = 10.0*np.sum(u.l1*u.l1 + u.l3*u.l3 + v.l1*v.l1 + v.l3*v.l3)/(2*ny*nx)

                t2 = self.f0*(s.l1-s.l3)/self.rgas
#             write(unit=1) s1
#             write(unit=1) s3
#             write(unit=1) t2
#             write(unit=1) u1t, (/ (0.0,jj=1,nx) /)
#             write(unit=1) u3t, (/ (0.0,jj=1,nx) /)
#             write(unit=1) v1t(:,1:ny), (/ (0.0,jj=1,nx) /)
#             write(unit=1) v3t(:,1:ny), (/ (0.0,jj=1,nx) /)
                print "KE %5.1f %15.7f %15.7f %15.7f %15.7f %15.7f %15.7f \n" %( day, t2.max(), u.l1t.max(), u.l3t.max(), zke, eke, tke)
#          end if
#       end if
# !!$      if ( istep == 130 ) then
# !!$         do j=1,ny
# !!$            write(unit=7,fmt="(i3,2e15.7)") j, u1z(j), u3z(j)
# !!$         end do
# !!$      end if

                # end subroutine diag

#    subroutine diagvv(s1t,s3t, s1tm, s3tm, dt)
#       ! Vertical velocity diagnostic
#       use constants
#       real, intent(in),  dimension(:,0:) :: s1t, s3t, s1tm, s3tm
#       real, intent(in) :: dt
#       real, dimension(0:ny) :: ds1, ds3, s1, s3, w2, vv
#       integer :: j

#       ! Just work on zonal means
#       ds1 = 0.0 ! s1t(1,:) - s1tm(1,:)
#       ds3 = 0.0 !s3t(1,:) - s3tm(1,:)
#       s1  = 0.5 * ( s1t(1,:) + s1tm(1,:) )
#       s3  = 0.5 * ( s3t(1,:) + s3tm(1,:) )

#       w2(0) = 0.0
#       w2(ny) = 0.0
#       do j=1,ny-1
#          ! Jacobian term drops out because using zonal means.
#          ! A term also dropped here
#          w2(j) = lambdasq*p2/f0 * ( (ds1(j)-ds3(j))/dt + &
#               2*rgas*heat*(2*j-ny)/(f0*cp*ny) )
#       end do
#       vv(0) = 0.0
#       do j=1,ny-1
#          vv(j) = vv(j-1) - dy/p2 * w2(j)
#       end do
#       vv(ny) = 0.0
#       do j=0,ny
#          write(unit=4,fmt="(i3,2e15.7)") j, w2(j), vv(j)
#       end do
#    end subroutine diagvv

    def run(self):

        #  Implementation of Phillips' model.
        # The time integration is carried forward using the vorticity in a three
        # time level scheme. At each step the streamfunction must be diagnosed from 
        #  this. 

        #  The streamfunction is denoted by s and vorticity by v.
        #  Vector wind components are uc and vc.
        #  Total fields are denoted with suffix t and zonal means with suffix z

	day1=131.0; day2 = 160.0
	freq = 0.25
	dt1 = 86400.; dt2 = 3600.
	imp = True

        # All initialised to zero, so model is at rest
        s  = Var() #  ! Streamfunction
        sm = Var() #  ! Streamfunction at tau-1
        v =  Var() #  ! Vorticity
        vm = Var() #  ! Vorticity at tau - 1 values
        # Temporaries used in timestepping
        x = Var()

        seed = 123
        day = 0.0
        addnoise = False

        # Time loop
        while True:
            if day < day1:
                zonal = True # 130 day spin up using a one day time step
                dt = dt1
            else:
                zonal = False
                dt = dt2
                self.diag_flag = False
                if not addnoise:
                    # Noise hasn't been done yet
                    r = ranlux.ranlux()
                    r.rluxgo(3,seed,0,0)
                    rarray = r.ranlux(Grid.nx*(Grid.ny-1))
                    rand = np.zeros((Grid.nx+1,Grid.ny+1))
                    kk = 0
                    for j in range(1,Grid.ny):
                        for i in range(1,Grid.nx+1):
                            rand[i,j] = rarray[kk]
                            kk += 1
                    # Remove the zonal mean
                    rand -= rand[1:,:].mean(axis=0)
                    rand *= 8.5e6
                    #  Define a random streamfunction anomaly
                    s.l1t[:] += rand
                    s.l3t[:] += rand
                    # Convert this to a vorticity anomaly
                    self.calcvor(s, v)
                    addnoise = True
                    # Time step has changed so linearly interpolate the previous time
                    # value to ensure a smooth start
                    vm.l1t[:] = v.l1t - (v.l1t-vm.l1t)*dt2/dt1
                    vm.l3t[:] = v.l3t - (v.l3t-vm.l3t)*dt2/dt1
                    self.diag_flag = True

            v.split()     
            self.calc_zonstream(v, s)
            # print "VS", abs(v.l1z).max(), abs(s.l1z).max()

            #  Use relaxation to solve for the anomaly streamfunction
            if zonal:
                s.l1[:] = 0.0
                s.l3[:] = 0.0
            else:
                self.relax1(v, s)

            sm.settot(s)  # Do a copy
            # Should be a function
            for j in range(Grid.ny+1):
                s.l1t[:,j] = s.l1[:,j] + s.l1z[j]
                s.l3t[:,j] = s.l1[:,j] + s.l3z[j]

            self.diag(day, zonal, s)

            # call diagvv(s1t, s3t, s1tm, s3tm, dt)

            # Time stepping
            # print "SV", v.l3t[5,4], v.l3t[5,5], v.l3t[5,6], s.l3t[5,4], s.l3t[5,5], s.l3t[5,6]
            self.xcalc(v, vm, s, sm, dt, x)
            x.split()
            # print "X", x.l1t[5,5], x.l3t[5,5]
       # print*, "Xanom", day, maxval(abs(x1))/maxval(abs(x1z))

            # Solve for new zonal mean vorticity from x
            if imp:
                self.calc_zvor(x,dt,v)
            else:
                v.l1z[1:ny] = x.l1z[1:ny]
                v.l1z[0]    = v.l1z[1]
                v.l1z[ny]   = v.l1z[ny-1]
                v.l3z[1:ny] = x.l3z[1:ny]
                v.l3z[0]    = v.l3z[1]
                v.l3z[ny]   = v.l3z[ny-1]

            if day == 0.0:
                # Forward step from rest. This simple form for the forward step
                # assumes that it's starting from rest.
                v.l1z = 0.5*x.l1z
                v.l3z = 0.5*x.l3z
            else:
                # Update previous value of vorticity
                vm.settot(v)

            if zonal:
                v.set(0.0)
            else:
                if imp:
                    # Relaxation solver for non-zonal terms
                    self.relax2(x, dt, v)
                else:
                    v.set(x)

            # print*, "Vanom", day, maxval(abs(v1))/maxval(abs(v1z))
            for j in range(Grid.ny+1):
                v.l1t[:,j] = v.l1[:,j] + v.l1z[j]
                v.l3t[:,j] = v.l3[:,j] + v.l3z[j]
            # print "V1Z", abs(v.l1z).max()
#             print "V", abs(v.l1t).max(), abs(v.l3t).max()

            day = day + dt/86400.0
            if day > day2:
                break

if __name__ == '__main__':
    m = Model()
    m.run()
