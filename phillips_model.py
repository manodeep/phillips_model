# Phillip's model
# Direct translation of fortran code, using fortran arrays starting at one.

# The time integration is carried forward using the vorticity in a three
# time level scheme. At each step the streamfunction must be diagnosed from
# this.

# The streamfunction is denoted by s and vorticity by v.
# Vector wind components are uc and vc.
# Total fields are denoted with suffix t and zonal means with suffix z

import numpy as np
from scipy.linalg.lapack import dgtsv, dgetrf, dgetrs
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import netCDF4

def msq_rand(x):
    # Middle 10 digits following Hammer
    return ((x*x) // 10**5) % 10**10

class Grid:
    L = 6.0e6
    W = 5.0e6  # y coord goes from -W to W (Phillips notation])
    nx = 16
    ny = 16
    dx = L / nx   # 3.75e5
    dy = 2*W / ny # 6.25e5    # Gridsize

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
                print(f"{self.l1t[i,j]:12.4f}", end="")
            print()

    def adump(self):
        for j in range(0,Grid.ny+1):
            for i in range(1,Grid.nx+1):
                print(f"{self.l1[i,j]:12.4f}", end="")
            print()

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

    eps = Grid.dx/Grid.dy
    epsq = eps*eps
    heat = 2.0e-3           #  Heating W/kg
    lambdasq = 1.5e-12      # Stability parameter, m^-2
    rgas = 287.0
    cp = 1004.0
    f0 = 1.0e-4
    beta = 1.6e-11          # s^{-1} m^{-1}
    p4 = 1.0e5              # Surface pressure, Pa
    p2 = 0.5*p4
    gamma = lambdasq*Grid.dx**2  # 0.21
    k = 4.0e-6              # s^{-1}  Surface drag

    # Control variables
    a = 1.0e5               # m^2/s   Horizonal diffusion
    diag_flag = False
    accel = 1.0

    noisescale = 7.509e6

    # For netcdf output
    save_netcdf = False
    irec = -1

    day1 = 131.0  # Zonal spin up length
    dt1 = 86400.  # Spin up time step
    day2 = 162.   # Total run length
    dt2 = 7200.   # Time step in regular run
    dt = dt1
    variable_step = True
    diag_freq = 86400

    first_step = True # For solver initialisation

    # All initialised to zero, so model is at rest
    v =  Var() #  ! Vorticity
    vm = Var() #  ! Vorticity at tau - 1 values
    s  = Var() #  ! Streamfunction
    # Temporaries used in timestepping
    x = Var()

    time = 0
    day = 0
    np.seterr(over='raise', invalid='raise', divide='raise')

    def calcvor(self, s, v):
        # This repeats same code for each level. Should the level be another dimension?
        for j in range(1,Grid.ny):
            jm = j-1
            jp = j+1
            for i in range(1,Grid.nx+1):
                im = i-1
                if im == 0:
                    im = Grid.nx
                ip = i+1
                if ip == Grid.nx+1:
                    ip = 1
                v.l1t[i,j] = ( s.l1t[ip,j]  + s.l1t[im,j] - 2*s.l1t[i,j] ) + \
                    self.epsq * ( s.l1t[i,jp] + s.l1t[i,jm] - 2*s.l1t[i,j] ) - \
                    self.gamma * ( s.l1t[i,j] - s.l3t[i,j] )
                v.l3t[i,j] = ( s.l3t[ip,j] + s.l3t[im,j] - 2*s.l3t[i,j] ) + \
                    self.epsq * ( s.l3t[i,jp] + s.l3t[i,jm] - 2*s.l3t[i,j] ) + \
                    self.gamma * ( s.l1t[i,j] - s.l3t[i,j] )
        # Follow A17 and set end rows to zonal mean of neighbours
        v.l1t[:,0] = v.l1t[:,1].mean()
        v.l1t[:,Grid.ny] = v.l1t[:,Grid.ny-1].mean()
        v.l3t[:,0] = v.l3t[:,1].mean()
        v.l3t[:,Grid.ny] = v.l3t[:,Grid.ny-1].mean()

    def calc_zonstream(self, v, s):
        # Given vorticity variable as input, solve for the
        # zonal mean streamfunction

        ny = Grid.ny
        nz = 2*ny - 3  #  Number of zonal means to solve for
        epsq = self.epsq
        gamma = self.gamma

        if self.first_step:
            self.first_step = False
            # Start these arrays from 1 again
            amat = np.zeros((nz+1,nz+1))

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

            self.lu, self.piv, info = dgetrf(amat[1:,1:])

        bmat = np.zeros(nz+1)
        bmat[1:ny] = v.l1z[1:ny]
        for j in range(2,ny):
            bmat[ny-2+j] = v.l3z[j]

        # Solve AX=B
        bmat[1:], info = dgetrs(self.lu, self.piv, bmat[1:])

        for j in range(1,ny):
            s.l1z[j] = bmat[j]
        for j in range(2,ny):
            s.l3z[j] = bmat[ny-2+j]

        # Apply the BC
        s.l1z[0] = s.l1z[1]
        s.l1z[ny] = s.l1z[ny-1]
        s.l3z[0] = 0.0
        # Note the extra restriction on s3(1) which is not solved for.
        s.l3z[1] = 0.0
        s.l3z[ny] = s.l3z[ny-1]

    def relax1(self, v, s):
        # Solve for anomaly streamfunction

        # print("V1 anom", abs(v.l1).max())
        nx = Grid.nx
        ny = Grid.ny
        # Start from the current value of the anomaly streamfunction
        for iter in range(100):
            # Jacobi iteration
            maxdiff = 0.0
            change = 0.0
            # for irb in range(2):
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
            # print("ITER1", iter, np.sqrt(change), maxdiff)
            # maxdiff is now only on a single level so halve the convergence
            # criterion
            if maxdiff < 0.5*3.75e4:
                # print("ITER1 done", iter, np.sqrt(change), maxdiff)
                break
        if iter >= 100:
            raise Exception(f"RELAX1 failed {iter} {np.sqrt(change)} {maxdiff}")

    def relax2(self, x, dt, v):
        # Solve for anomaly vorticity

        temp = Var()
        nx = Grid.nx
        ny = Grid.ny

        v.l1[:] = 0.0
        v.l3[:] = 0.0
        alpha = self.a*dt/Grid.dx**2
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
            # Boundary condition A17
            v.l1[:,0] = v.l1[:,1].mean()
            v.l1[:,ny] = v.l1[:,ny-1].mean()
            change3 = np.sum ( ( v.l3[1:,1:] - temp.l3[1:,1:] ) **2 )
            v.l3[:,1:ny] = temp.l3[:,1:ny]
            v.l3[:,0] = v.l3[:,1].mean()
            v.l3[:,ny] = v.l3[:,ny-1].mean()
            if max(change1, change3) < 1.0:
                # print("ITER2", iter, np.sqrt(change1), np.sqrt(change3))
                break

    def xcalc(self, v, vm, s, dt, x):

        nx = Grid.nx
        ny = Grid.ny
        alpha = self.a*dt/Grid.dx**2
        b = self.beta*Grid.dx**2*Grid.dy
        c = dt/(2.0*Grid.dx*Grid.dy)
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
                x.l3t[i,j] = x.l3t[i,j] -  self.k*dt*(1.5*vm.l3t[i,j] - v.l1t[i,j] -
                                    4*self.gamma*(s.l1t[i,j]-s.l3t[i,j]) )

    def calc_zvor(self, x, dt, v):
        # Solve for the zonal mean vorticity
        ny = Grid.ny
        epsq = self.epsq

        nz = ny-1
        amat = np.zeros(nz-1)
        bmat = np.zeros(nz)
        cmat = np.zeros(nz-1)
        rmat = np.zeros(nz)

        alpha = self.a*dt/Grid.dx**2
        # Level 1
        amat[:] = alpha*epsq
        bmat[0] = bmat[nz-1] = -alpha*epsq - 1
        bmat[1:nz-1] = -2.0*alpha*epsq - 1
        cmat[:] = alpha*epsq
        rmat[:] = -x.l1z[1:nz+1]

        result = dgtsv(amat, bmat, cmat, rmat)
        umat = result[3]
        v.l1z[1:nz+1] = umat[:]
        v.l1z[0] = v.l1z[1]
        v.l1z[ny] = v.l1z[ny-1]

        # Level 3
        # amat[:] = alpha*epsq
        # bmat[0] = bmat[nz-1] = -alpha*epsq - 1 - 1.5*self.k*dt
        # bmat[1:nz-1] = -2.0*alpha*epsq - 1 - 1.5*self.k*dt
        # cmat[:] = alpha*epsq
        bmat -= 1.5*self.k*dt
        rmat[:] = -x.l3z[1:nz+1]

        result = dgtsv(amat, bmat, cmat, rmat)
        umat = result[3]
        v.l3z[1:nz+1] = umat[:]
        v.l3z[0] = v.l3z[1]
        v.l3z[ny] = v.l3z[ny-1]

    def zonal_diag(self, day, s):
        u = Var()
        ny = Grid.ny
        dy = Grid.dy

        t2z = self.f0*(s.l1z[:]-s.l3z[:])/self.rgas
        #  Zonal wind, u = -pd{psi}{y}
        u.l1z[1:] = - ( s.l1z[1:ny+1] - s.l1z[0:ny] ) / dy
        u.l3z[1:] = - ( s.l3z[1:ny+1] - s.l3z[0:ny] ) / dy

        # Zonal KE, scaled such that a wind of 1 m/s everywhere gives 10
        zke = 10.0*np.sum(u.l1z*u.l1z + u.l3z*u.l3z)/(2*ny)

        # Zonal PE  A21
        # Sum 1:J-1 / J matches definition of Y operator
        zpe = 5*self.lambdasq*np.sum((s.l1z[1:ny] - s.l3z[1:ny])**2) / ny

        print("TZ %5.1f %15.7f %15.7f %15.7f %15.7f %15.7f" % ( day, t2z.max(), u.l1z.max(), u.l3z.max(), zke, zpe))

    def calc_T(self):
        return self.f0*(self.s.l1t-self.s.l3t)/self.rgas

    def calc_ps(self):
        return 0.01 * (1.5*self.s.l3t - 0.5*self.s.l1t)*self.f0

    def diag(self, day, s):

        u = Var()
        v = Var()
        vshift = Var()
        nx = Grid.nx; ny = Grid.ny
        dx = Grid.dx; dy = Grid.dy

        u.l1t[:,1:] = - ( s.l1t[:,1:ny+1] - s.l1t[:,0:ny] ) / dy
        u.l3t[:,1:] = - ( s.l3t[:,1:ny+1] - s.l3t[:,0:ny] ) / dy
        for i in range(1,nx+1):
            im = i-1
            if im == 0:
                im = nx
            v.l1t[i,:] = ( s.l1t[i,:] - s.l1t[im,:] ) / dx
            v.l3t[i,:] = ( s.l3t[i,:] - s.l3t[im,:] ) / dx
        #  Average v to get it on the same grid points as u
        for i in range(1,nx+1):
            ip = i+1
            if ip > nx:
                ip = 1
            for j in range(1,ny+1):
                vshift.l1t[i][j] = 0.25*(v.l1t[i][j] + v.l1t[ip][j] +
                                    v.l1t[i][j-1] + v.l1t[ip][j-1])
                vshift.l3t[i][j] = 0.25*(v.l3t[i][j] + v.l3t[ip][j] +
                                    v.l3t[i][j-1] + v.l3t[ip][j-1])
        # print("MAX V", v.l1t.max(), vshift.l1t.max())
        # Calculate zonal mean and eddy winds
        u.split()
        vshift.split()
        v.split()

        # Note factor of 10 here.
        tke = 10.0*np.sum(u.l1t[1:,1:]*u.l1t[1:,1:] + u.l3t[1:,1:]*u.l3t[1:,1:] +
                            vshift.l1t[1:,1:]*vshift.l1t[1:,1:] + vshift.l3t[1:,1:]*vshift.l3t[1:,1:])/(2*ny*nx)

        zke = 10.0*np.sum(u.l1z*u.l1z + u.l3z*u.l3z)/(2*ny)
        # eke = 10.0*np.sum(u.l1*u.l1 + u.l3*u.l3 + vshift.l1*vshift.l1 + vshift.l3*vshift.l3)/(2*ny*nx)
        eke = 10.0*np.sum(u.l1**2 + u.l3**2 + v.l1**2 + v.l3**2)/(2*ny*nx)

        zpe = 5*self.lambdasq*np.sum((s.l1z[1:ny] - s.l3z[1:ny])**2) / ny
        epe = 5*self.lambdasq*np.sum((s.l1[:,1:ny] - s.l3[:,1:ny])**2) / (nx*ny)

        print("KE %6.2f %9.2f %9.2f %9.2f %9.2f" %( day, zke, eke, epe, zpe))
        if eke > 1e5:
            raise Exception("EKE too large")

    def stability_criterion(self, dt, s):
        # Stability criterion (A13)
        smax = 0.
        for i in range(1,Grid.nx+1):
            im = i-1
            if im == 0:
                im = Grid.nx
            ip = i+1
            if ip > Grid.nx:
                ip = 1
            for j in range(1,Grid.ny):
                jm = j-1
                jp = j+1
                smax = max(smax, abs(s.l1t[ip,j]-s.l1t[im,j]) + abs(s.l1t[i,jp]-s.l1t[i,jm]))
                smax = max(smax, abs(s.l3t[ip,j]-s.l3t[im,j]) + abs(s.l3t[i,jp]-s.l3t[i,jm]))
        return 0.5*dt*smax / (Grid.dx*Grid.dy)

    def create_nc_output(self):

        self.ds = netCDF4.Dataset('phillips_model.nc', 'w')
        ds = self.ds
        ds.createDimension('lon', Grid.nx)
        ds.createDimension('lat', 1+Grid.ny)
        ds.createDimension('vlon', Grid.nx)
        ds.createDimension('ulat', Grid.ny)
        ds.createDimension('lev', 2)
        ds.createDimension('time', None)
        v = ds.createVariable('lon', np.float32, ('lon',))
        v.long_name = 'longitude'
        v.units = 'degrees_east'
        v = ds.createVariable('lat', np.float32, ('lat',))
        v.long_name = 'latitude'
        v.units = 'degrees_north'
        v = ds.createVariable('vlon', np.float32, ('vlon',))
        v.long_name = 'V longitude'
        v.units = 'degrees_east'
        v = ds.createVariable('ulat', np.float32, ('ulat',))
        v.long_name = 'U latitude'
        v.units = 'degrees_north'
        v = ds.createVariable('lev', np.float32, ('lev',))
        v.long_name = 'model level'
        v = ds.createVariable('time', np.float32, ('time',))
        v.units = "days since 2000-01-01 00:00"
        v = ds.createVariable('strm', np.float32, ('time', 'lev', 'lat', 'lon'))
        v.long_name = 'streamfunction'
        v.units = 'm2 s-1'
        v = ds.createVariable('vor', np.float32, ('time', 'lev', 'lat', 'lon'))
        v.long_name = 'vorticity'
        v.units = 's-1'
        v = ds.createVariable('t500', np.float32, ('time', 'lat', 'lon'))
        v.long_name = 'air temperature at 500 hPa'
        v.units = 'K'
        v = ds.createVariable('ps', np.float32, ('time', 'lat', 'lon'))
        v.long_name = 'surface pressure'
        v.units = 'hPa'
        v = ds.createVariable('u', np.float32, ('time', 'lev', 'ulat', 'lon'))
        v.long_name = 'zonal wind'
        v.units = 'm s-1'
        v = ds.createVariable('v', np.float32, ('time', 'lev', 'lat', 'vlon'))
        v.long_name = 'meridional wind'
        v.units = 'm s-1'
        v = ds.createVariable('eke', np.float32, ('time',))
        v.long_name = "Eddy kinetic energy (*10)"
        v = ds.createVariable('zke', np.float32, ('time',))
        v.long_name = "Zonal kinetic energy (*10)"
        v = ds.createVariable('epe', np.float32, ('time',))
        v.long_name = "Eddy potential energy (*10)"
        v = ds.createVariable('zpe', np.float32, ('time',))
        v.long_name = "Zonal potential energy (*10)"

        # At latitude 45
        dlon = 360 * Grid.dx / (6.371e6*2*np.pi*np.cos(np.pi/4))
        ds.variables['lon'][:] = np.arange(Grid.nx) * dlon
        lat0 = 45.
        dlat = 360 * Grid.dy / (6.371e6*2*np.pi)
        ds.variables['lat'][:] = 45 + dlat*np.arange(-Grid.ny//2, Grid.ny//2+1)
        ds.variables['vlon'][:] = ds.variables['lon'][:] - 0.5*dlon
        ds.variables['ulat'][:] = ds.variables['lat'][:-1] + 0.5*dlat
        ds.variables['lev'][:] = [1., 3.]

    def nc_output(self, day, v, s):
        # Python variables are (nx, ny) so need to transpose when writing
        # to match netCDF dimensions
        self.irec += 1
        self.ds.variables['vor'][self.irec,0] = v.l1t[1:].T
        self.ds.variables['vor'][self.irec,1] = v.l3t[1:].T
        self.ds.variables['strm'][self.irec,0] = s.l1t[1:].T
        self.ds.variables['strm'][self.irec,1] = s.l3t[1:].T

        nx = Grid.nx; ny=Grid.ny
        u = Var()
        vtmp = Var()
        u.l1t[:,1:] = - ( s.l1t[:,1:ny+1] - s.l1t[:,0:ny] ) / Grid.dy
        u.l3t[:,1:] = - ( s.l3t[:,1:ny+1] - s.l3t[:,0:ny] ) / Grid.dy
        self.ds.variables['u'][self.irec,0] = u.l1t[1:,1:].T
        self.ds.variables['u'][self.irec,1] = u.l3t[1:,1:].T
        for i in range(1,nx+1):
            im = i-1
            if im == 0:
                im = nx
            vtmp.l1t[i,:] = ( s.l1t[i,:] - s.l1t[im,:] ) / Grid.dx
            vtmp.l3t[i,:] = ( s.l3t[i,:] - s.l3t[im,:] ) / Grid.dx
        self.ds.variables['v'][self.irec,0] = vtmp.l1t[1:].T
        self.ds.variables['v'][self.irec,1] = vtmp.l3t[1:].T

        u.split()
        vtmp.split()
        zke = 10.0*np.sum(u.l1z**2 + u.l3z**2)/(2*ny)
        eke = 10.0*np.sum(u.l1**2 + u.l3**2 + vtmp.l1**2 + vtmp.l3**2)/(2*ny*nx)
        zpe = 5*self.lambdasq*np.sum((s.l1z[1:ny] - s.l3z[1:ny])**2) / ny
        epe = 5*self.lambdasq*np.sum((s.l1[:,1:ny] - s.l3[:,1:ny])**2) / (nx*ny)
        self.ds.variables['zke'][self.irec] = zke
        self.ds.variables['eke'][self.irec] = eke
        self.ds.variables['zpe'][self.irec] = zpe
        self.ds.variables['epe'][self.irec] = epe

        t2 = self.f0*(s.l1t-s.l3t)/self.rgas
        self.ds.variables['t500'][self.irec] = t2[1:].T

        ps = 0.01 * (1.5*s.l3t - 0.5*s.l1t)*self.f0
        self.ds.variables['ps'][self.irec] = ps[1:].T

        # Use time since perturbation
        self.ds.variables['time'][self.irec] = day - self.day1

    def spinup(self):

        v = self.v
        vm = self.vm
        s = self.s
        x = self.x

        # Time loop
        while True:
            v.split()
            self.calc_zonstream(v, s)

            s.l1[:] = 0.0
            s.l3[:] = 0.0

            # Should be a function
            for j in range(Grid.ny+1):
                s.l1t[:,j] = s.l1[:,j] + s.l1z[j]
                s.l3t[:,j] = s.l3[:,j] + s.l3z[j]
            if self.time % self.diag_freq == 0:
                self.zonal_diag(self.day, s)

            # Time stepping
            self.xcalc(v, vm, s, self.dt, x)
            x.split()

            # Solve for new zonal mean vorticity from x
            self.calc_zvor(x, self.dt, v)

            if self.time == 0.0:
                # Forward step from rest. This simple form for the forward step
                # assumes that it's starting from rest.
                v.l1z = 0.5*x.l1z
                v.l3z = 0.5*x.l3z
            else:
                # Update previous value of vorticity
                vm.settot(v)

            v.set(0.0)
            for j in range(Grid.ny+1):
                v.l1t[:,j] = v.l1z[j]
                v.l3t[:,j] = v.l3z[j]

            self.time += self.dt
            self.day = self.time/86400.0
            if self.day >= self.day1:
                break

        return vm, v

    def perturb(self):

        vm = self.vm
        v = self.v
        stmp = Var() # For random initialisation
        vtmp = Var() # For random initialisation

        # Time step has changed so linearly interpolate the previous time
        # value to ensure a smooth start
        vm.l1t[:] = v.l1t - (v.l1t-vm.l1t)*self.dt2/self.dt1
        vm.l3t[:] = v.l3t - (v.l3t-vm.l3t)*self.dt2/self.dt1

        rval = 1_111_111_111
        for i in range(1,Grid.nx+1):
            for j in range(1,Grid.ny):
                rval = msq_rand(rval)
                stmp.l1t[i,j] = float(rval) / 10**10
                stmp.l3t[i,j] = stmp.l1t[i,j]
        # Remove the zonal mean and scale
        stmp.split()
        stmp.l1t[:] = self.noisescale*stmp.l1[:]
        stmp.l3t[:] = self.noisescale*stmp.l3[:]
        # Convert this to a vorticity anomaly
        self.calcvor(stmp, vtmp)
        v.l1t[:] += vtmp.l1t
        v.l3t[:] += vtmp.l3t
        vm.l1t[:] += vtmp.l1t
        vm.l3t[:] += vtmp.l3t

    def step(self):

        v = self.v
        vm = self.vm
        s = self.s
        x = self.x

        v.split()
        self.calc_zonstream(v, s)

        #  Use relaxation to solve for the anomaly streamfunction
        self.relax1(v, s)

        # Should be a function
        for j in range(Grid.ny+1):
            s.l1t[:,j] = s.l1[:,j] + s.l1z[j]
            s.l3t[:,j] = s.l3[:,j] + s.l3z[j]
        if self.time % self.diag_freq == 0:
            self.diag(self.day, s)
            if self.save_netcdf:
                self.nc_output(self.day, v, s)

        # Time stepping
        self.xcalc(v, vm, s, self.dt, x)
        x.split()

        # Solve for new zonal mean vorticity from x
        self.calc_zvor(x,self.dt,v)

        if self.time == 0.0:
            # Forward step from rest. This simple form for the forward step
            # assumes that it's starting from rest.
            v.l1z = 0.5*x.l1z
            v.l3z = 0.5*x.l3z
        else:
            # Update previous value of vorticity
            vm.settot(v)

        # Relaxation solver for non-zonal terms
        self.relax2(x, self.dt, v)

        for j in range(Grid.ny+1):
            v.l1t[:,j] = v.l1[:,j] + v.l1z[j]
            v.l3t[:,j] = v.l3[:,j] + v.l3z[j]

        self.time += self.dt
        self.day = self.time/86400.0

        if self.variable_step and self.time % 86400 == 0 and self.dt > 1800:
            stab_crit = self.stability_criterion(self.dt, s)
            # print(f"Stability  {self.day:.2f} {stab_crit:.3f}")
            if stab_crit > 0.9:
                print(f"At {self.day:.2f} {stab_crit:.3f} Adjusting time step to {self.dt-1800}")
                vm.l1t[:] = v.l1t - (v.l1t-vm.l1t)*(self.dt-1800)/self.dt
                vm.l3t[:] = v.l3t - (v.l3t-vm.l3t)*(self.dt-1800)/self.dt
                self.dt -= 1800

def main():
    import time
    m = Model()
    if m.save_netcdf:
        m.create_nc_output()
    m.spinup()
    m.perturb()
    animate = False
    if animate:
        fig, axes = plt.subplots()
        # p = plt.pcolormesh(m.s.l1t.T)
        p = plt.pcolormesh(m.calc_T().T[::-1], vmin=-30, vmax=30, cmap='coolwarm')
        def animate(i):
            m.step()
            # p.set_array(m.s.l1t.T.flatten())
            p.set_array(m.calc_T().T[::-1].flatten())
            axes.set_title(f"Day {m.day:.2f}")
        ani = animation.FuncAnimation(fig, animate, frames=1000,
                                interval=0)

        plt.show()
    else:
        m.dt = m.dt2
        t1 = time.perf_counter()
        while m.day < m.day2:
            m.step()
        t2 = time.perf_counter()
        print("Elapsed time", t2-t1)

# import cProfile
# from pstats import SortKey
# cProfile.run('main()', sort=SortKey.CUMULATIVE)
main()
