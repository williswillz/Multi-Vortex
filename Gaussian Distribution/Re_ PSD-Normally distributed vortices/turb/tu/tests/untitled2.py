# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:21:36 2019

@author: WS1
"""
import numpy as np
from numpy.fft import fftn
from numpy import sqrt, zeros, conj, pi, arange, ones, convolve


# ------------------------------------------------------------------------------

def movingaverage(interval, window_size):
    window = ones(int(window_size)) / float(window_size)
    return convolve(interval, window, 'same')


# ------------------------------------------------------------------------------

def compute_tke_spectrum_1d(u, lx, ly, lz, smooth):
    """
    Given a velocity field u this function computes the kinetic energy
    spectrum of that velocity field in spectral space. This procedure consists of the
    following steps:
    1. Compute the spectral representation of u using a fast Fourier transform.
    This returns uf (the f stands for Fourier)
    2. Compute the point-wise kinetic energy Ef (kx, ky, kz) = 1/2 * (uf)* conjugate(uf)
    3. For every wave number triplet (kx, ky, kz) we have a corresponding spectral kinetic energy
    Ef(kx, ky, kz). To extract a one dimensional spectrum, E(k), we integrate Ef(kx,ky,kz) over
    the surface of a sphere of radius k = sqrt(kx^2 + ky^2 + kz^2). In other words
    E(k) = sum( E(kx,ky,kz), for all (kx,ky,kz) such that k = sqrt(kx^2 + ky^2 + kz^2) ).
    Parameters:
    -----------
    u: 3D array
      The x-velocity component.
    v: 3D array
      The y-velocity component.
    w: 3D array
      The z-velocity component.
    lx: float
      The domain size in the x-direction.
    ly: float
      The domain size in the y-direction.
    lz: float
      The domain size in the z-direction.
    smooth: boolean
      A boolean to smooth the computed spectrum for nice visualization.
    """
    nx = len(u[:, 0, 0])
    ny = len(u[0, :, 0])
    nz = len(u[0, 0, :])

    nt = nx * ny * nz
    n = max(nx, ny, nz)  # int(np.round(np.power(nt,1.0/3.0)))

    uh = fftn(u) / nt

    # tkeh = zeros((nx, ny, nz))
    tkeh = 0.5 * (uh * conj(uh)).real

    length = max(lx, ly, lz)

    knorm = 2.0 * pi / length

    kxmax = nx / 2
    kymax = ny / 2
    kzmax = nz / 2

    wave_numbers = knorm * arange(0, n)
    tke_spectrum = zeros(len(wave_numbers))

    for kx in range(nx):
        rkx = kx
        if kx > kxmax:
            rkx = rkx - nx
        for ky in range(ny):
            rky = ky
            if ky > kymax:
                rky = rky - ny
            for kz in range(nz):
                rkz = kz
                if kz > kzmax:
                    rkz = rkz - nz
                rk = sqrt(rkx * rkx + rky * rky + rkz * rkz)
                k = int(np.round(rk))
                print('k = ', k)
                tke_spectrum[k] = tke_spectrum[k] + tkeh[kx, ky, kz]

    tke_spectrum = tke_spectrum / knorm

    if smooth:
        tkespecsmooth = movingaverage(tke_spectrum, 5)  # smooth the spectrum
        tkespecsmooth[0:4] = tke_spectrum[0:4]  # get the first 4 values from the original data
        tke_spectrum = tkespecsmooth

    knyquist = knorm * min(nx, ny, nz) / 2

    return knyquist, wave_numbers, tke_spectrum


# ------------------------------------------------------------------------------

def compute_tke_spectrum(u, v, w, lx, ly, lz, smooth):
    """
    Given a velocity field u, v, w, this function computes the kinetic energy
    spectrum of that velocity field in spectral space. This procedure consists of the
    following steps:
    1. Compute the spectral representation of u, v, and w using a fast Fourier transform.
    This returns uf, vf, and wf (the f stands for Fourier)
    2. Compute the point-wise kinetic energy Ef (kx, ky, kz) = 1/2 * (uf, vf, wf)* conjugate(uf, vf, wf)
    3. For every wave number triplet (kx, ky, kz) we have a corresponding spectral kinetic energy
    Ef(kx, ky, kz). To extract a one dimensional spectrum, E(k), we integrate Ef(kx,ky,kz) over
    the surface of a sphere of radius k = sqrt(kx^2 + ky^2 + kz^2). In other words
    E(k) = sum( E(kx,ky,kz), for all (kx,ky,kz) such that k = sqrt(kx^2 + ky^2 + kz^2) ).
    Parameters:
    -----------
    u: 3D array
      The x-velocity component.
    v: 3D array
      The y-velocity component.
    w: 3D array
      The z-velocity component.
    lx: float
      The domain size in the x-direction.
    ly: float
      The domain size in the y-direction.
    lz: float
      The domain size in the z-direction.
    smooth: boolean
      A boolean to smooth the computed spectrum for nice visualization.
    """
    nx = len(u[:, 0, 0])
    ny = len(v[0, :, 0])
    nz = len(w[0, 0, :])

    nt = nx * ny * nz
    n = nx  # int(np.round(np.power(nt,1.0/3.0)))

    uh = fftn(u) / nt
    vh = fftn(v) / nt
    wh = fftn(w) / nt

    tkeh = 0.5 * (uh * conj(uh) + vh * conj(vh) + wh * conj(wh)).real

    k0x = 2.0 * pi / lx
    k0y = 2.0 * pi / ly
    k0z = 2.0 * pi / lz

    knorm = (k0x + k0y + k0z) / 3.0
    print('knorm = ', knorm)

    kxmax = nx / 2
    kymax = ny / 2
    kzmax = nz / 2

    # dk = (knorm - kmax)/n
    # wn = knorm + 0.5 * dk + arange(0, nmodes) * dk

    wave_numbers = knorm * arange(0, n)

    tke_spectrum = zeros(len(wave_numbers))

    for kx in range(-nx//2, nx//2-1):
        for ky in range(-ny//2, ny//2-1):
            for kz in range(-nz//2, nz//2-1):
                rk = sqrt(kx**2 + ky**2 + kz**2)
                k = int(np.round(rk))
                tke_spectrum[k] += tkeh[kx, ky, kz]
    # for kx in range(nx):
    #     rkx = kx
    #     if kx > kxmax:
    #         rkx = rkx - nx
    #     for ky in range(ny):
    #         rky = ky
    #         if ky > kymax:
    #             rky = rky - ny
    #         for kz in range(nz):
    #             rkz = kz
    #             if kz > kzmax:
    #                 rkz = rkz - nz
    #             rk = sqrt(rkx * rkx + rky * rky + rkz * rkz)
    #             k = int(np.round(rk))
    #             tke_spectrum[k] = tke_spectrum[k] + tkeh[kx, ky, kz]

    tke_spectrum = tke_spectrum / knorm

    #  tke_spectrum = tke_spectrum[1:]
    #  wave_numbers = wave_numbers[1:]
    if smooth:
        tkespecsmooth = movingaverage(tke_spectrum, 5)  # smooth the spectrum
        tkespecsmooth[0:4] = tke_spectrum[0:4]  # get the first 4 values from the original data
        tke_spectrum = tkespecsmooth

    knyquist = knorm * min(nx, ny, nz) / 2

    return knyquist, wave_numbers, tke_spectrum


# ------------------------------------------------------------------------------

def compute_tke_spectrum_flatarrays(u, v, w, nx, ny, nz, lx, ly, lz, smooth):
    unew = u.reshape([nx, ny, nz])
    vnew = v.reshape([nx, ny, nz])
    wnew = w.reshape([nx, ny, nz])
    k, w, espec = compute_tke_spectrum(unew, vnew, wnew, lx, ly, lz, smooth)
    return k, w, espec

class FileFormats():
  FLAT=1
  XYZ=2
  IJK=3

import numpy as np
from numpy import sin, cos, sqrt, ones, zeros, pi, arange


def generate_isotropic_turbulence(lx, ly, lz, nx, ny, nz, nmodes, wn1, especf):
    """
    Given an energy spectrum, this function computes a discrete, staggered, three
    dimensional velocity field in a box whose energy spectrum corresponds to the input energy
    spectrum up to the Nyquist limit dictated by the grid
    This function returns u, v, w as the axial, transverse, and azimuthal velocities.
    Parameters:
    -----------
    lx: float
      The domain size in the x-direction.
    ly: float
      The domain size in the y-direction.
    lz: float
      The domain size in the z-direction.
    nx: integer
      The number of grid points in the x-direction.
    ny: integer
      The number of grid points in the y-direction.
    nz: integer
      The number of grid points in the z-direction.
    wn1: float
      Smallest wavenumber. Typically dictated by spectrum or domain size.
    espec: functor
      A callback function representing the energy spectrum.
    """

    # generate cell centered x-grid
    dx = lx / nx
    dy = ly / ny
    dz = lz / nz

    # START THE FUN!

    # compute random angles
    phi = 2.0 * pi * np.random.uniform(0.0, 1.0, nmodes)
    nu = np.random.uniform(0.0, 1.0, nmodes)
    theta = np.arccos(2.0 * nu - 1.0)
    psi = np.random.uniform(-pi / 2.0, pi / 2.0, nmodes)

    # highest wave number that can be represented on this grid (nyquist limit)
    wnn = max(np.pi / dx, max(np.pi / dy, np.pi / dz))
    print('I will generate data up to wave number: ', wnn)

    # wavenumber step
    dk = (wnn - wn1) / nmodes

    # wavenumber at cell centers
    wn = wn1 + 0.5 * dk + arange(0, nmodes) * dk

    dkn = ones(nmodes) * dk

    #   wavenumber vector from random angles
    kx = sin(theta) * cos(phi) * wn
    ky = sin(theta) * sin(phi) * wn
    kz = cos(theta) * wn

    # create divergence vector
    ktx = np.sin(kx * dx / 2.0) / dx
    kty = np.sin(ky * dy / 2.0) / dy
    ktz = np.sin(kz * dz / 2.0) / dz

    # Enforce Mass Conservation
    phi1 = 2.0 * pi * np.random.uniform(0.0, 1.0, nmodes)
    nu1 = np.random.uniform(0.0, 1.0, nmodes)
    theta1 = np.arccos(2.0 * nu1 - 1.0)
    zetax = sin(theta1) * cos(phi1)
    zetay = sin(theta1) * sin(phi1)
    zetaz = cos(theta1)
    sxm = zetay * ktz - zetaz * kty
    sym = -(zetax * ktz - zetaz * ktx)
    szm = zetax * kty - zetay * ktx
    smag = sqrt(sxm * sxm + sym * sym + szm * szm)
    sxm = sxm / smag
    sym = sym / smag
    szm = szm / smag

    # verify that the wave vector and sigma are perpendicular
    kk = np.sum(ktx * sxm + kty * sym + ktz * szm)
    print('Orthogonality of k and sigma (divergence in wave space):', kk)

    # get the modes
    km = wn

    espec = especf(km)
    espec = espec.clip(0.0)

    # generate turbulence at cell centers
    um = sqrt(espec * dkn)
    u_ = zeros([nx, ny, nz])
    v_ = zeros([nx, ny, nz])
    w_ = zeros([nx, ny, nz])

    xc = dx / 2.0 + arange(0, nx) * dx
    yc = dy / 2.0 + arange(0, ny) * dy
    zc = dz / 2.0 + arange(0, nz) * dz

    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
                # for every grid point (i,j,k) do the fourier summation
                arg = kx * xc[i] + ky * yc[j] + kz * zc[k] - psi
                bmx = 2.0 * um * cos(arg - kx * dx / 2.0)
                bmy = 2.0 * um * cos(arg - ky * dy / 2.0)
                bmz = 2.0 * um * cos(arg - kz * dz / 2.0)
                u_[i, j, k] = np.sum(bmx * sxm)
                v_[i, j, k] = np.sum(bmy * sym)
                w_[i, j, k] = np.sum(bmz * szm)

    print('done generating turbulence.')
    return u_, v_, w_


def generate_scalar_isotropic_turbulence(lx, ly, lz, nx, ny, nz, nmodes, wn1, especf):
    """
    Given an energy spectrum, this function computes a discrete, staggered, three
    dimensional velocity field in a box whose energy spectrum corresponds to the input energy
    spectrum up to the Nyquist limit dictated by the grid
    This function returns u, v, w as the axial, transverse, and azimuthal velocities.
    Parameters:
    -----------
    lx: float
      The domain size in the x-direction.
    ly: float
      The domain size in the y-direction.
    lz: float
      The domain size in the z-direction.
    nx: integer
      The number of grid points in the x-direction.
    ny: integer
      The number of grid points in the y-direction.
    nz: integer
      The number of grid points in the z-direction.
    wn1: float
      Smallest wavenumber. Typically dictated by spectrum or domain size.
    espec: functor
      A callback function representing the energy spectrum.
    """

    # generate cell centered x-grid
    dx = lx / nx
    dy = ly / ny
    dz = lz / nz

    # START THE FUN!

    # compute random angles
    phi = 2.0 * pi * np.random.uniform(0.0, 1.0, nmodes)
    nu = np.random.uniform(0.0, 1.0, nmodes);
    theta = np.arccos(2.0 * nu - 1.0);
    psi = np.random.uniform(-pi / 2.0, pi / 2.0, nmodes)

    # highest wave number that can be represented on this grid (nyquist limit)
    wnn = max(np.pi / dx, max(np.pi / dy, np.pi / dz))
    print('I will generate data up to wave number: ', wnn)

    # wavenumber step
    dk = (wnn - wn1) / nmodes

    # wavenumber at cell centers
    wn = wn1 + 0.5 * dk + arange(0, nmodes) * dk

    dkn = ones(nmodes) * dk

    #   wavenumber vector from random angles
    kx = sin(theta) * cos(phi) * wn
    ky = sin(theta) * sin(phi) * wn
    kz = cos(theta) * wn

    # get the modes
    km = wn

    espec = especf(km)
    espec = espec.clip(0.0)

    # generate turbulence at cell centers
    um = sqrt(espec * dkn)
    scalar_ = zeros([nx, ny, nz])

    xc = dx / 2.0 + arange(0, nx) * dx
    yc = dy / 2.0 + arange(0, ny) * dy
    zc = dz / 2.0 + arange(0, nz) * dz

    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
                # for every grid point (i,j,k) do the fourier summation
                arg = kx * xc[i] + ky * yc[j] + kz * zc[k] - psi
                bm = 2.0 * um * cos(arg)
                scalar_[i, j, k] = np.sum(bm)

                print('done. I am awesome!')
    return scalar_

import multiprocessing as mp
import numpy as np
import time
from numpy import sin, cos, sqrt, ones, zeros, pi, arange
from numpy import linalg as LA

def compute_turbulence(nthread, dx, dy, dz, psi, um, kx, ky, kz, sxm, sym, szm, nx, ny, nz, nxAll, nyAll, nzAll, ip, jp,
                       kp, q):
    print('Generating turbulence on thread:', nthread)
    t0 = time.time()
    u_ = zeros([nx, ny, nz])
    v_ = zeros([nx, ny, nz])
    w_ = zeros([nx, ny, nz])

    xl = (ip - 1) * nx
    xh = ip * nx
    yl = (jp - 1) * ny
    yh = jp * ny
    zl = (kp - 1) * nz
    zh = kp * nz

    xc = dx / 2.0 + arange(xl, xh) * dx
    yc = dy / 2.0 + arange(yl, yh) * dy
    zc = dz / 2.0 + arange(zl, zh) * dz  # cell centered coordinates

    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
                # for every grid point (i,j,k) do the fourier summation
                arg = kx * xc[i] + ky * yc[j] + kz * zc[k] - psi
                bmx = 2.0 * um * cos(arg - kx * dx / 2.0)
                bmy = 2.0 * um * cos(arg - ky * dy / 2.0)
                bmz = 2.0 * um * cos(arg - kz * dz / 2.0)
                u_[i, j, k] = np.sum(bmx * sxm)
                v_[i, j, k] = np.sum(bmy * sym)
                w_[i, j, k] = np.sum(bmz * szm)

    t1 = time.time()
    print('Thread ', nthread, ' done generating turbulence in ', t1 - t0, 's')
    q.put((ip, jp, kp, u_, v_, w_))
    return ip, jp, kp, u_, v_, w_


def generate_isotropic_turbulence(patches, lx, ly, lz, nx, ny, nz, nmodes, wn1, especf):
    ## grid generation
    # generate cell centered x-grid
    dx = lx / nx
    dy = ly / ny
    dz = lz / nz

    ## START THE FUN!
    # compute random angles
    np.random.seed(0)
    phi = 2.0 * pi * np.random.uniform(0.0, 1.0, nmodes);
    nu = np.random.uniform(0.0, 1.0, nmodes);
    theta = np.arccos(2.0 * nu - 1.0);
    psi = np.random.uniform(-pi / 2.0, pi / 2.0, nmodes);
    alfa = 2.0 * pi * np.random.uniform(0.0, 1.0, nmodes);

    # highest wave number that can be represented on this grid (nyquist limit)
    wnn = max(np.pi / dx, max(np.pi / dy, np.pi / dz));
    print('I will generate data up to wave number: ', wnn)

    # wavenumber step
    dk = (wnn - wn1) / nmodes

    # wavenumber at cell centers
    wn = wn1 + arange(0, nmodes) * dk
    #  wn = wn1 + np.arange(0,nmodes)*dk*np.log(np.arange(0,nmodes) + 1)/np.log(nmodes)
    dkn = ones(nmodes) * dk
    #  dkn = wn[1:nmodes] - wn[0:nmodes-1]
    #  dkn = np.append(dkn,dkn[nmodes-2])

    #   wavenumber vector from random angles
    kx = sin(theta) * cos(phi) * wn
    ky = sin(theta) * sin(phi) * wn
    kz = cos(theta) * wn

    # create divergence vector
    ktx = np.sin(kx * dx / 2.0) / (dx)
    kty = np.sin(ky * dy / 2.0) / (dy)
    ktz = np.sin(kz * dz / 2.0) / (dz)

    #  # Use Davidson's Method to enforce Divergence Free Condition
    #  ktmag = sqrt(ktx*ktx + kty*kty + ktz*ktz)
    #  theta = np.arccos(kzstag/kstagmag)
    #  phi = np.arctan2(kystag,kxstag)
    #  sxm = cos(phi)*cos(theta)*cos(alfa) - sin(phi)*sin(alfa)
    #  sym = sin(phi)*cos(theta)*cos(alfa) + cos(phi)*sin(alfa)
    #  szm = -sin(theta)*cos(alfa)

    # another method to generate sigma = zeta x k_tilde, pick zeta randomly
    #  np.random.seed(3)
    phi1 = 2.0 * pi * np.random.uniform(0.0, 1.0, nmodes);
    nu1 = np.random.uniform(0.0, 1.0, nmodes);
    theta1 = np.arccos(2.0 * nu1 - 1.0);
    zetax = sin(theta1) * cos(phi1)
    zetay = sin(theta1) * sin(phi1)
    zetaz = cos(theta1)
    sxm = zetay * ktz - zetaz * kty
    sym = -(zetax * ktz - zetaz * ktx)
    szm = zetax * kty - zetay * ktx
    smag = sqrt(sxm * sxm + sym * sym + szm * szm)
    sxm = sxm / smag
    sym = sym / smag
    szm = szm / smag

    # verify that the wave vector and sigma are perpendicular
    kk = np.sum(ktx * sxm + kty * sym + ktz * szm)
    print('Orthogonality of k and sigma (divergence in wave space):', kk)

    # get the modes
    km = wn

    # now create an interpolant for the spectrum. this is needed for
    # experimentally-specified spectra
    #  espec = especf(km + dk/2) + especf(km))*0.5
    espec = especf(km)
    espec = espec.clip(0.0)

    # generate turbulence at cell centers
    um = sqrt(espec * dkn)

    #  must use Manager queue here, or will not work
    nxthreads = patches[0];
    nythreads = patches[1];
    nzthreads = patches[2];
    nxt = nx // nxthreads;
    nyt = nx // nythreads;
    nzt = nx // nzthreads;

    manager = mp.Manager()
    mq = manager.Queue()
    pool = mp.Pool(mp.cpu_count())  # assume 2 threads per core

    # fire off workers
    jobs = []
    nthread = 0
    for k in range(1, nzthreads + 1):
        for j in range(1, nythreads + 1):
            for i in range(1, nxthreads + 1):
                nthread = nthread + 1
                job = pool.apply_async(compute_turbulence, (
                nthread, dx, dy, dz, psi, um, kx, ky, kz, sxm, sym, szm, nxt, nyt, nzt, nx, ny, nz, i, j, k, mq))
                jobs.append(job)

    # collect results from the workers through the pool result queue
    print('now collecting results from individual threads...')
    uarrays = []
    varrays = []
    warrays = []
    patches = []
    for job in jobs:
        i, j, k, u, v, w = job.get()
        uarrays.append(u)
        varrays.append(v)
        warrays.append(w)
        patches.append([i, j, k])
    del u, v, w

    pool.terminate()
    pool.close()

    # combine the arrays computed from threads into large arrays
    print('now combining velocity fields generated by the individual threads...')
    uall = zeros([nx, ny, nz])
    vall = zeros([nx, ny, nz])
    wall = zeros([nx, ny, nz])
    nthread = 0
    for k in range(1, nzthreads + 1):
        for j in range(1, nythreads + 1):
            for i in range(1, nxthreads + 1):
                uall[(i - 1) * nxt:i * nxt, (j - 1) * nyt:j * nyt, (k - 1) * nzt:k * nzt] = uarrays[nthread]
                vall[(i - 1) * nxt:i * nxt, (j - 1) * nyt:j * nyt, (k - 1) * nzt:k * nzt] = varrays[nthread]
                wall[(i - 1) * nxt:i * nxt, (j - 1) * nyt:j * nyt, (k - 1) * nzt:k * nzt] = warrays[nthread]
                nthread = nthread + 1

    return uall, vall, wall

#from fileformats import FileFormats
import multiprocessing as mp
import time

class FileFormats():
  FLAT=1
  XYZ=2
  IJK=3
#------------------------------------------------------------------------------

def writefileparallel(u, v, w, dx, dy, dz, fileformat):
  print ('Writing to disk. This may take a while...')
  writeufile = mp.Process(target=writefile, args=('u.txt','x',dx,dy,dz,u, fileformat))
  writeufile.start()
  
  writevfile = mp.Process(target=writefile, args=('v.txt','y',dx,dy,dz,v, fileformat))
  writevfile.start()
  
  writewfile = mp.Process(target=writefile, args=('w.txt','z',dx,dy,dz,w, fileformat))
  writewfile.start()
  
  writeufile.join()
  writevfile.join()
  writewfile.join()
  
#------------------------------------------------------------------------------
  
def writefile(filename, velcomponent, dx, dy, dz, velarray, fileformat):
  t0 = time.time()  

  nx = len(velarray[:,0,0])
  ny = len(velarray[0,:,0])
  nz = len(velarray[0,0,:])  
  
  f = open(filename , 'w')
  zo=[0,0,0]
  if(velcomponent=='x'):
    zo=[0,1,1]
  elif(velcomponent=='y'):
    zo=[1,0,1]
  else:
    zo=[1,1,0]
    
  # loop over the velocity fields generated by each thread
  if (fileformat == FileFormats.XYZ):    
    f.write('%s \n' % 'XYZ')
    xlo = zo[0]*dx/2.0
    ylo = zo[1]*dy/2.0
    zlo = zo[2]*dz/2.0
    for k in range(0,nz):
      for j in range(0,ny):
        for i in range(0,nx):
          x = xlo + i*dx
          y = ylo + j*dy
          z = zlo + k*dz
          u = velarray[i,j,k]              
          f.write('%.16f %.16f %.16f %.16f \n' % (x,y,z,u))        
  elif (fileformat == FileFormats.IJK):
    f.write('%s \n' % 'IJK')    
    for k in range(0,nz):
      for j in range(0,ny):
        for i in range(0,nx):
          u = velarray[i,j,k]              
          f.write('%d %d %d %.16f \n' % (i,j,k,u))
  else:
    f.write('%s \n' % 'FLAT')
    f.write('%d %d %d \n' % (nx, ny, nz))
    for k in range(0,nz):
      for j in range(0,ny):
        for i in range(0,nx):
          u = velarray[i,j,k]              
          f.write('%.16f\n' % u)        
  f.close()
  t1 = time.time()
  print ('Done writing to disk in ', t1 - t0, 's')

from scipy import interpolate
import numpy as np
from numpy import pi, exp
import time
import scipy.io
#from tkespec import compute_tke_spectrum_1d
#import isoturb
#import isoturbo
import matplotlib.pyplot as plt
#from fileformats import FileFormats
#import isoio
plt.interactive(True)

# load an experimental specturm. Alternatively, specify it via a function call
cbcspec = np.loadtxt('cbc_spectrum.txt')
kcbc = cbcspec[:, 0] * 100
ecbc = cbcspec[:, 1] * 1e-6
especf = interpolate.interp1d(kcbc, ecbc, 'cubic')


def cbc_spec(k):
    return especf(k)


def karman_spec(k):
    nu = 1.0e-5
    alpha = 1.452762113
    urms = 0.25
    ke = 40.0
    kappae = np.sqrt(5.0 / 12.0) * ke
    L = 0.746834 / kappae  # integral length scale - sqrt(Pi)*Gamma(5/6)/Gamma(1/3)*1/ke
    #  L = 0.05 # integral length scale
    #  Kappae = 0.746834/L
    epsilon = urms * urms * urms / L
    kappaeta = pow(epsilon, 0.25) * pow(nu, -3.0 / 4.0)
    r1 = k / kappae
    r2 = k / kappaeta
    espec = alpha * urms * urms / kappae * pow(r1, 4) / pow(1.0 + r1 * r1, 17.0 / 6.0) * np.exp(-2.0 * r2 * r2)
    return espec


def power_spec(k):
    Nu = 1 * 1e-3
    L = 0.1
    Li = 1
    ch = 1
    cl = 10
    p0 = 8
    c0 = pow(10, 2)
    Beta = 2
    Eta = Li / 20.0
    ES = Nu * Nu * Nu / (Eta * Eta * Eta * Eta)
    x = k * Eta
    fh = np.exp(-Beta * pow(pow(x, 4) + pow(ch, 4), 0.25) - ch)
    x = k * L
    fl = pow(x / pow(x * x + cl, 0.5), 5.0 / 3.0 + p0)
    espec = c0 * pow(k, -5.0 / 3.0) * pow(ES, 2.0 / 3.0) * fl * fh
    return espec


# ----------------------------------------------------------------------------------------------
# __    __   ______   ________  _______         ______  __    __  _______   __    __  ________
# |  \  |  \ /      \ |        \|       \       |      \|  \  |  \|       \ |  \  |  \|        \
# | $$  | $$|  $$$$$$\| $$$$$$$$| $$$$$$$\       \$$$$$$| $$\ | $$| $$$$$$$\| $$  | $$ \$$$$$$$$
# | $$  | $$| $$___\$$| $$__    | $$__| $$        | $$  | $$$\| $$| $$__/ $$| $$  | $$   | $$
# | $$  | $$ \$$    \ | $$  \   | $$    $$        | $$  | $$$$\ $$| $$    $$| $$  | $$   | $$
# | $$  | $$ _\$$$$$$\| $$$$$   | $$$$$$$\        | $$  | $$\$$ $$| $$$$$$$ | $$  | $$   | $$
# | $$__/ $$|  \__| $$| $$_____ | $$  | $$       _| $$_ | $$ \$$$$| $$      | $$__/ $$   | $$
# \$$    $$ \$$    $$| $$     \| $$  | $$      |   $$ \| $$  \$$$| $$       \$$    $$   | $$
#  \$$$$$$   \$$$$$$  \$$$$$$$$ \$$   \$$       \$$$$$$ \$$   \$$ \$$        \$$$$$$     \$$
# ----------------------------------------------------------------------------------------------

# specify whether you want to use threads or not to generate turbulence
use_parallel = False
patches = [1, 1, 8]
filespec = 'cbc'
whichspec = cbc_spec

# set the number of modes you want to use to represent the velocity.
nmodes = 500
N = 32

# write to file
enableIO = False  # enable writing to file
fileformat = FileFormats.FLAT  # Specify the file format supported formats are: FLAT, IJK, XYZ

# save the velocity field as a matlab matrix (.mat)
savemat = False

# compute the mean of the fluctuations for verification purposes
computeMean = True

# input domain size in the x, y, and z directions. This value is typically
# based on the largest length scale that your data has. For the cbc data,
# the largest length scale corresponds to a wave number of 15, hence, the
# domain size is L = 2pi/15.
lx = 9 * 2.0 * pi / 100.0
ly = 9 * 2.0 * pi / 100.0
lz = 9 * 2.0 * pi / 100.0

# input number of cells (cell centered control volumes). This will
# determine the maximum wave number that can be represented on this grid.
# see wnn below
nx = N  # number of cells in the x direction
ny = N  # number of cells in the y direction
nz = N  # number of cells in the z direction

# enter the smallest wavenumber represented by this spectrum
wn1 = 15  # determined here from cbc spectrum properties

# ------------------------------------------------------------------------------
# END USER INPUT
# ------------------------------------------------------------------------------
t0 = time.time()
phi = isoturb.generate_scalar_isotropic_turbulence(lx, ly, lz, nx, ny, nz, nmodes, wn1, whichspec)
t1 = time.time()
print('it took me ', t1 - t0, ' s to generate the isotropic turbulence.')

dx = lx / nx
dy = ly / ny
dz = lz / nz

#isoio.writefile('u.txt', 'x', dx, dy, dz, u, fileformat)

# if savemat:
#     data = {}  # CREATE empty dictionary
#     data['U'] = u
#     data['V'] = v
#     data['W'] = w
#     scipy.io.savemat('uvw.mat', data)

# compute mean velocities
if computeMean:
    phimean = np.mean(phi)
    print('mean u = ', phimean)

    phifluc = phimean - phi

    # print
    # 'mean u fluct = ', np.mean(ufluc)

    phifrms = np.mean(phifluc * phifluc)

    # print
    # 'u fluc rms = ', np.sqrt(ufrms)
    # print
    # 'v fluc rms = ', np.sqrt(vfrms)
    # print
    # 'w fluc rms = ', np.sqrt(wfrms)

# verify that the generated velocities fit the spectrum
knyquist, wavenumbers, tkespec = compute_tke_spectrum_1d(phi, lx, ly, lz, True)

# compare spectra
# integral comparison:
# find index of nyquist limit
idx = (np.where(wavenumbers == knyquist)[0][0]) - 1

km0 = 2.0 * np.pi / lx
nmodes = 5000
dk0 = (knyquist - km0) / nmodes
exactRange = km0 + np.arange(0, nmodes + 1) * dk0
exactE = np.trapz(karman_spec(exactRange), dx=dk0)
numE = np.trapz(tkespec[0:idx], dx=wavenumbers[0])
# print
# 'diff = ', abs(exactE - numE) / exactE * 100

# analyze how well we fit the input spectrum
# espec = cbc_spec(kcbc) # compute the cbc original spec

# compute the RMS error committed by the generated spectrum
# find index of nyquist limit
idx = (np.where(wavenumbers == knyquist)[0][0]) - 1
exact = whichspec(wavenumbers[4:idx])
num = tkespec[4:idx]
diff = np.abs((exact - num) / exact)
meanE = np.mean(diff)
print('got here ')
# print
# 'Mean Error = ', meanE * 100.0, '%'
# rmsE = np.sqrt(np.mean(diff * diff))
# print
# 'RMS Error = ', rmsE * 100, '%'

# np.savetxt('tkespec_' + filespec + '_' + str(N) + '.txt',np.transpose([wavenumbers,tkespec]))



# fig = plt.figure(figsize=(3.5, 2.6), dpi=100)
# plt.rc("font", size=10, family='serif')
wnn = np.arange(wn1, 2000)
l1, = plt.loglog(kcbc,ecbc, 'k-', label='input')
plt.loglog(wnn, whichspec(wnn), 'k-', label='input')
l2, = plt.loglog(wavenumbers, tkespec, 'bo-', markersize=4, markerfacecolor='w', markevery=1, label='computed')
# plt.axis([8, 10000, 1e-7, 1e-2])
# # plt.xticks(fontsize=12)
# # plt.yticks(fontsize=12)
# plt.axvline(x=knyquist, linestyle='--', color='black')
# plt.xlabel('$\kappa$ (1/m)')
# plt.ylabel('$E(\kappa)$ (m$^3$/s$^2$)')
# plt.grid()
# plt.gcf().tight_layout()
# # plt.title(str(N)+'$^3$')
# # plt.legend(handles=[l1,l2],loc=3)
# # fig.savefig('tkespec_' + filespec + '_' + str(N) + '.pdf')
#
q, ((p1,p2),(p3,p4)) = plt.subplots(2,2)

p1.plot(kcbc, ecbc, 'ob', kcbc, ecbc, '-')
p1.set_title('Interpolated Spectrum')
p1.grid()
p1.set_xlabel('wave number')
p1.set_ylabel('E')

p2.loglog(kcbc, ecbc, '-', wavenumbers, tkespec, 'ro-')
p2.axvline(x=knyquist, linestyle='--', color='black')
p2.set_title('Spectrum of generated turbulence')
p2.grid()

# contour plot
p3.matshow(phi[:,:,nz/2])
p3.set_title('phi')

p4.matshow(phi[:,ny/2,:])
p4.set_title('phi')
# #
plt.show(1)