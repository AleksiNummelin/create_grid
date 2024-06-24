import numpy as np
from scipy import linalg
from scipy.integrate import solve_ivp
import tripolar_funcs as tfuncs
import xarray as xr
#
make_plots=True#True
vectorized=True
outpath='../grid_out/'
#################
# ORIGINAL VALUES
#M =round(720/S);               % Number of grid cells around equator
#N0=round(376/S);               % Number of grid cells towards the south pole
#N1=round(177/S);               % Number of grid cells towards the intersection
#N2=round(216/S);               % Number of grid cells from the intersection towards the north pole
#
#plat= 62;                     % Latitude of the new poles
#plon= -52;
#
# Parameters defining the grid
max_step = 1
rtol = 1E-12
atol = 1E-16
#
S  = 1/4 #0.125
M  = np.round(720/S).astype(int)              # Number of grid cells around equator
#N0 = np.round(360/S).astype(np.int) #np.round(416/S).astype(np.int) #np.round(376/S).astype(np.int)              # Number of grid cells towards the south pole
N0 = np.round(374/S).astype(int) #268
# NOTE THAT IF YOU BRING DOWN THE PLAT (I.E CLOSER TO THE EQUATOR) THEN YOU SHOULD REDUCE N1 AND INCREASE N2
# It is also easy to end up with a situation where you do not get to the north pole. In this case try decreasing 
#N1 = np.round((3)/S).astype(np.int) #N1 = np.round((38)/S).astype(np.int) #N1 = np.round((54)/S).astype(np.int) #np.round(100/S).astype(np.int) #np.round(177/S).astype(np.int)              # Number of grid cells towards the intersection
#N2 = np.round((301)/S).astype(np.int) #N2 = np.round((257)/S).astype(np.int) #N2 = np.round((287)/S).astype(np.int) #np.round(240/S).astype(np.int) #np.round(216/S).astype(np.int)              # Number of grid cells from the intersection towards the north pole
#
# probably best to bring down the o
plat = 45.5 #45 #38 #65 #62                                # Latitude of the new poles
plon = -88 #-90 #90 #-88 #-95 #-52                             # Longitude of the new poles
#
######################################
# tests
# Check that the polar ellipse closes up towards the pole
# moving the pole northward seems to help
# the original values seem to work
#plat = 42 #this works plat=42; N0=428; N1=97; N2=296
#plon = -90
#N0 = np.round(428/S).astype(np.int) #
#
N1 = np.round(40/S).astype(int)
N2 = np.round(232/S).astype(int)
if S<1:
    N2=N2+1
# 235 works
# these almost work
#N1 = np.round(40/S).astype(np.int)  #
#N2 = np.round(230/S).astype(np.int) 
######################################
#
c = 1.08*2*np.tan(np.radians((90-plat)/2))/np.pi # Rate of change of semi-minor axis at the line connecting the new poles
#c = 1.0*2*np.tan(np.radians((90-plat)/2))/np.pi
if S>=1:
   fe = 1/4 #1/4                                         # Meridional divided by zonal grid scale at equator
   ledeg = 16 #15 #30 #22                           # Meridional range of Equatorial grid scale compression
else:
   fe=1
   ledeg = 0

# Set various parameters
#
rad = np.pi/180
deg = 180/np.pi
#
N       = N0+N1+N2
dtheta  = 2*np.pi/M
phi_p   = (90-plat)*rad
theta_p = plon*rad
le      = ledeg*rad

#options = odeset('reltol',1e-12,'abstol',1e-16);
#
# Create angular and radial coordinates where the radial coordinate can
# incorporate grid compression near the equator

theta = np.arange(0,M+1)*dtheta
y     = np.ones((N+1))*np.nan
#
if fe == 1:
  y[range(N0+1)[::-1]] = -np.arange(N0+1)*dtheta
  y[N0:N+1]            = np.arange(N1+N2+1)*dtheta
else:
  # something fishy here
  out1 = solve_ivp(fun=lambda t,y: tfuncs.dydn(t,y,dtheta,fe,le), t_span = (0,-N0), y0=np.array([0]), t_eval=-np.arange(0,N0+1), vectorized=vectorized, max_step=max_step, rtol=rtol, atol=atol) 
  #ode45(@dydn,[0:-1:-N0],0,options);
  out2 = solve_ivp(fun=lambda t,y: tfuncs.dydn(t,y,dtheta,fe,le), t_span = (0,N1+N2), y0=np.array([0]), t_eval=np.arange(0,N1+N2+1), vectorized=vectorized, max_step=max_step, rtol=rtol, atol=atol)
  #ode45(@dydn,[0:N1+N2],0,options);
  y[:N0+1]  = out1.y.squeeze()[::-1]
  y[N0:N+1] = out2.y.squeeze()

r=np.exp(-y)
#
y0=y[0]
y1=y[N0+N1]
y2=y[N]
#
# Compute the coefficients for the polynomial expression of the
# semi-major axis a

D=np.array([[1,       y1,     y1**2,     y1**3,     y1**4],
   [0,        1,     2*y1,    3*y1**2,   4*y1**3],
   [1,       y2,     y2**2,     y2**3,     y2**4],
   [0,        1,     2*y2,    3*y2**2,   4*y2**3],
     [0,        0,        2,       6*y2,  12*y2**2]])

f=np.array([np.exp(-y1), -np.exp(-y1), np.tan(phi_p/2), 0, 0])

#g=(D**(-1)*f);
g = np.dot(linalg.pinv2(D),f)
a0=g[0]
a1=g[1]
a2=g[2]
a3=g[3]
a4=g[4]
a_coeff = [a0, a1, a2, a3, a4]
# Compute the coefficients for the polynomial expression of the
# semi-minor axis b

D=np.array([[1,       y1,     y1**2,    y1**3,     y1**4,      y1**5],
            [0,        1,     2*y1,   3*y1**2,   4*y1**3,    5*y1**4],
            [0,        0,        2,     6*y1,   12*y1**2,   20*y1**3],
            [1,       y2,     y2**2,    y2**3,     y2**4,      y2**5],
            [0,        1,     2*y2,   3*y2**2,   4*y2**3,    5*y2**4],
            [0,        0,        2,     6*y2,   12*y2**2,   20*y2**3]])

f=np.array([np.exp(-y1), -np.exp(-y1), np.exp(-y1), 0, -c, 0])

#g=(D**(-1)*f);
g = np.dot(linalg.pinv2(D),f)
b0=g[0]
b1=g[1]
b2=g[2]
b3=g[3]
b4=g[4]
b5=g[5]
b_coeff = [b0, b1, b2, b3, b4, b5]

print('Latitude of southern boundary of grid: ', np.round((np.pi/2-2*np.arctan(r[0]))*deg,decimals=2))
print('Latitude of grid intersection: ', np.round((np.pi/2-2*np.arctan(r[N0+N1+1]))*deg,decimals=2))
print('Latitude of northern hemisphere poles: ', np.round((np.pi/2-2*np.arctan(tfuncs.a(y2,a_coeff)))*deg,decimals=2))
print('Grid dimension of resulting ocean grid: ' + str(M//2)+'x'+str((N+1)//2))
#
#
#figure(1);set(gcf,'renderer','painters')
#n=N0+N1+1:N+1;
#plot(n,a(y(n)),'r',n,b(y(n)),'b',n,r(n),'k')
#pause

# Create grid south of the intersection

sx  = np.ones((M+1,N+1))*np.nan
sy  = np.ones((M+1,N+1))*np.nan
ang = np.ones((M+1,N+1))*np.nan
sx[:,0:N0+N1+1]  = np.dot(np.cos(theta[:,np.newaxis]),r[0:N0+N1+1][:,np.newaxis].T)
sy[:,0:N0+N1+1]  = np.dot(np.sin(theta[:,np.newaxis]),r[0:N0+N1+1][:,np.newaxis].T)
ang[:,0:N0+N1+1] = 0.

# Create grid north of the intersection

n             = range(N0+N1+1,N+1)
sx[0,n]       = tfuncs.a(y[n],a_coeff)
sy[0,n]       = 0
ang[0,n]      = 0
sx[M//2,n]    = -sx[0,n]
sy[M//2,n]    = 0 
ang[M//2,n] = 0
sx[M,n]     = sx[0,n]
sy[M,n]     = 0
ang[M,n]    = 0

#f_dpdn = partial(dpdn, dtheta=dtheta,fe=fe,le=le,a_coeff=a_coeff,b_coeff=b_coeff)
#
for m in range(1,M//4+1):
    #
    #t,p = solve_ivp(f_dpdn, tspan=np.arange(N1,N1+N2+1), y0=[theta[m] y1], rtol=1e-12, atol=1e-16) #ode45(@dpdn,[N1:N1+N2],[theta(m) y1],options);
    out = solve_ivp(fun=lambda t,y: tfuncs.dpdn(t,y,dtheta,fe,le,a_coeff,b_coeff), t_span = (N1,N1+N2), y0=np.array([theta[m], y1]), t_eval=np.arange(N1,N1+N2+1), vectorized=vectorized, max_step=max_step, rtol=rtol, atol=atol)
    #
    psi = out.y[0,1:]
    yy  = out.y[1,1:]
    #
    a_t = tfuncs.a(yy,a_coeff)
    b_t = tfuncs.b(yy,b_coeff)
    sex = a_t*np.cos(psi) 
    sey = b_t*np.sin(psi)
    a2_t = a_t*a_t
    b2_t = b_t*b_t
    sex2 = sex*sex
    sey2 = sey*sey
    sx[m,N0+N1+1:N+1]  = sex
    sy[m,N0+N1+1:N+1]  = sey
    ang[m,N0+N1+1:N+1] = np.arccos(np.round(np.nanmin([np.ones(N+1-(N0+N1+1)),b2_t*sex2+a2_t*sey2],axis=0)/np.sqrt(((b2_t**2)*sex2+(a2_t**2)*sey2)*(sex2+sey2)),decimals=abs(int(np.log10(rtol)))))
    if np.any(np.isnan(ang[m,N0+N1+1:N+1])):
        print('nans found at:', np.where(np.isnan(ang[m,N0+N1+1:N+1])))
    #
    print('Integrated curve', m-1, np.floor(M/4), out.status)


m = range(0,int(np.ceil(M/4)))
n = range(N0+N1+1,N+1)
ninds,minds = np.meshgrid(n,m)
ninds = ninds.flatten()
minds = minds.flatten()
sx[M//2-minds,ninds]  = -sx[minds,ninds]
sy[M//2-minds,ninds]  =  sy[minds,ninds]
ang[M//2-minds,ninds] = -ang[minds,ninds]
#
m = range(1,M//2)
#
ninds,minds = np.meshgrid(n,m)
ninds =ninds.flatten()
minds =minds.flatten()
#
sx[M-minds,ninds]    =  sx[minds,ninds]
sy[M-minds,ninds]    = -sy[minds,ninds]
ang[M-minds,ninds]   = -ang[minds,ninds]
#
# Extend grid in y-direction
sx = np.concatenate([sx, sx[range(M,-1,-1),:][:,range(N-1,N-3,-1)]],axis=1) #[sx sx(M+1:-1:1,N:-1:N-1)];
sy = np.concatenate([sy, sy[range(M,-1,-1),:][:,range(N-1,N-3,-1)]],axis=1) #[sy sy(M+1:-1:1,N:-1:N-1)];
#
angext       = ang[range(M,-1,-1),:][:,range(N-1,N-3,-1)]-np.pi
ind          = np.where(angext>np.pi)
angext[ind]  = angext[ind]-2*np.pi
ind          = np.where(angext<-np.pi)
angext[ind]  = angext[ind]+2*np.pi
ang          = np.concatenate([ang, angext],axis=1)
#
# Extend grid in x-direction #should be M-1?
sx  = np.concatenate([sx[M-1,:][np.newaxis,], sx],axis=0) #[sx(M,:);sx];
sy  = np.concatenate([sy[M-1,:][np.newaxis,], sy],axis=0) #[sy(M,:);sy];
ang = np.concatenate([ang[M-1,:][np.newaxis,], ang],axis=0) #[ang(M,:);ang];
#
# Convert to longitude/latitude and adjust longitude to get the poles in
# the correct location
lon      = (np.arctan2(sy,sx)+theta_p)*deg
ind      = np.where(lon>180)
lon[ind] = lon[ind]-360
ind      = np.where(lon<=-180)
lon[ind] = lon[ind]+360
lat      = (np.pi/2-2*np.arctan(np.sqrt(sx*sx+sy*sy)))*deg
#
# Test angle
print('Test computation of angles...')
n=int(round(N-N2/4))-1
#
for m in range(M//16,M,M//4):
    dlon=lon[m,n]-lon[m-2,n]
    dlat=lat[m,n]-lat[m-2,n]
    print(np.arctan2(dlat,np.cos(lat[m-1,n]*rad)*dlon)*deg, ang[m-1,n]*deg)

print('lat range:',np.min(lat),np.max(lat))
# Save longitude/latitude/angle data
nj,mj = lon.T.shape
lon_nc = xr.DataArray(lon.T, coords=[np.arange(nj).astype(int), np.arange(mj).astype(int)], dims=['nj', 'mi'],name='lon')
lat_nc = xr.DataArray(lat.T, coords=[np.arange(nj).astype(int), np.arange(mj).astype(int)], dims=['nj', 'mi'],name='lat')
ang_nc = xr.DataArray(ang.T, coords=[np.arange(nj).astype(int), np.arange(mj).astype(int)], dims=['nj', 'mi'],name='ang')
data_out = xr.merge([lon_nc,lat_nc,ang_nc])
data_out.to_netcdf(outpath+'lonlatang_python.nc')


print('Write dimensions to a file')
#
file1 = open("dimensions.h","w+") 
#
file1.write('c --- itdm  = total grid dimension in i direction \n')
file1.write('c --- jtdm  = total grid dimension in j direction \n')
file1.write('      integer    itdm,jtdm \n')
file1.write('      parameter (itdm='+str(M//2)+',jtdm='+str((N+1)//2)+') \n')
file1.write('      integer idm,jdm \n')
file1.write('      parameter (idm=itdm,jdm=jtdm) \n')
file1.write('c --- halo size \n')
file1.write('      integer    nbdy \n')
file1.write('      parameter (nbdy=0) \n')
file1.write('      integer ii,jj,lp \n')
#
file1.close() 
#
#
#
#with open('lonlatang_python.bin', 'wb') as fid:
#    lon.T.tofile(fid)
#    lat.T.tofile(fid)
#    ang.T.tofile(fid)
#
#from array import array
#
#fid = open('lonlatang_python.bin','wb')
#array('d',lon.flatten()).tofile(fid)
#array('d',lat.flatten()).tofile(fid)
#array('d',ang.flatten()).tofile(fid)
#fid.close()
#
#fwrite(fid,lon,'float64');
#fwrite(fid,lat,'float64');
#fwrite(fid,ang,'float64');
#fclose(fid);
#
#import cartopy.crs as ccrs

if make_plots:
    #
    import matplotlib as mpl
    mpl.use('agg')
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    #
    #bathy=xr.open_dataset('orig_bathy.nc')
    #
    #projection1 = ccrs.NearsidePerspective(central_longitude=0, central_latitude=70.0, satellite_height = 1E7)
    projection1=ccrs.Robinson(central_longitude=200)
    #projection1 = ccrs.NorthPolarStereo(central_longitude=0,true_scale_latitude=70)
    figMap,axMap=plt.subplots(nrows=1,ncols=1, sharex=True,sharey=True,figsize=(10,5), subplot_kw={'projection':projection1})
    #axMap.contour(bathy.lon,bathy.lat,bathy.z,levels=[0],colors='gray',transform=ccrs.PlateCarree())
    axMap.coastlines(resolution='50m',color='gray')
    #axMap.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())
    #
    for m in range(0,M//2,int(10/S)):
        tcoords1 = projection1.transform_points(ccrs.PlateCarree(), lon[m,:],lat[m,:])
        clon,clat = tcoords1[:,0], tcoords1[:,1]
        ninds = np.where(np.isfinite(clon))
        clon,clat = clon[ninds],clat[ninds]
        axMap.plot(clon,clat,'k')

    for m in range(M//2+2,M,int(10/S)):
        tcoords1 = projection1.transform_points(ccrs.PlateCarree(), lon[m,:],lat[m,:])
        clon,clat = tcoords1[:,0], tcoords1[:,1]
        ninds = np.where(np.isfinite(clon))
        clon,clat = clon[ninds],clat[ninds]
        axMap.plot(clon,clat,'k')

    #
    for n in range(0,N,int(10/S)):
        jinds=np.argsort(lon[:,n])
        axMap.plot(lon[jinds,n],lat[jinds,n],'k',transform=ccrs.PlateCarree())
    #
    figMap.savefig(outpath+'lonlat_plat'+str(plat)+'_plon'+str(plon)+'_N0'+str(N0)+'_N1'+str(N1)+'_N2'+str(N2)+'_visualized.png',dpi=150)

#if False:
# Plot grid on the sphere
# fig2 = plt.figure() 
#
#m_proj('satellite','lon',0,'lat',50)
##
#m_coast('patch',[.78 .78 .78]);
#
#for m=1:M+1
#  m_plot(lon(m,:),lat(m,:),'k')
#
#for n=1:N+1
#  m_plot(lon(:,n),lat(:,n),'k')#
#
#m_grid

