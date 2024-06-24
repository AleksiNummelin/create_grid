import gsw
import numpy as np
import xarray as xr
from joblib import Parallel, delayed
#from itertools import count
#from shapely.geometry import Point
#from shapely.geometry.polygon import Polygon
#
#point = Point(0.5, 0.5)
#polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])
#print(polygon.contains(point))
#
#mindep = -25.
#
n_jobs  = -1
#
#dpath   = '/cluster/work/users/anu074/create_Hres_grid2/'
dpath = '../grid_out/'
grid    = xr.open_dataset('../grid_out/grid.nc')
nyt,nxt = grid.plat.shape
plat = grid.plat.load()
plon = grid.plon.load()
pclon = grid.pclon.load()
pclat = grid.pclat.load()
#
datain    = xr.open_dataset('../orig_data/pgeo_eur_34.00.nc')
elevation = datain.z.load()
lat       = datain.lat.load()
lon       = datain.lon.load()
#
nys,nxs = elevation.shape
#
odepth = xr.zeros_like(grid.parea).rename('odepth')
ldepth = xr.zeros_like(grid.parea).rename('ldepth')
#
ocount = xr.zeros_like(grid.parea).rename('ocount')
lcount = xr.zeros_like(grid.parea).rename('lcount')
#
def find_depths(odepth_mm,ldepth_mm,ocount_mm,lcount_mm,elevation,lon,lat,plon,plat,pclon,pclat,j,nxt,nys,nxs,c=4):
    '''
    
    '''
    if j%100==0:
        print(j)
    #
    for i in range(nxt):
        dlat  = abs(lat-plat.isel(y=j,x=i))
        dlon  = abs(lon-plon.isel(y=j,x=i))
        jj    = np.where(dlat==np.min(dlat))[0][0]
        ii    = np.where(dlon==np.min(dlon))[0][0]
        jmin  = max(0,jj-c)
        jmax  = min(nys,jj+c+1)
        imin  = ii-c
        imax  = ii+c+1
        if imin<0:
            #elev = xr.concat([elevation[jmin:jmax,imin:nxs],elevation[jmin:jmax,0:imax]],dim='lon')
            #lon1d = xr.concat([lon[imin:nxs],lon[0:imax]],dim='lon')
            elev = elevation.isel(lat=slice(jmin,jmax)).roll(lon=abs(imin)).isel(lon=slice(0,imax))
            lon1d = lon.roll(lon=abs(imin)).isel(lon=slice(0,imax))
        elif imax>nxs:
            #elev = xr.concat([elevation[jmin:jmax,imin:nxs],elevation[jmin:jmax,0:imax-nxs]],dim='lon')
            #lon1d =xr.concat([lon[imin:nxs],lon[0:imax-nxs]],dim='lon')
            dc    = nxs-imax
            elev  = elevation.isel(lat=slice(jmin,jmax)).roll(lon=dc).isel(lon=slice(imin+dc,imax+dc))
            lon1d = lon.roll(lon=dc).isel(lon=slice(imin+dc,imax+dc))
        else:
            elev = elevation[jmin:jmax,imin:imax]
            lon1d = lon[imin:imax]
        #
        lat1d = lat[jmin:jmax]
        lsize=lat1d.size*lon1d.size
        if elev.min()<0:
            #
            colon = np.concatenate([plon.isel(x=i,y=j).values*np.ones((4,1)),pclon.isel(x=i,y=j).values[:,np.newaxis]],axis=-1)
            colat = np.concatenate([plat.isel(x=i,y=j).values*np.ones((4,1)),pclat.isel(x=i,y=j).values[:,np.newaxis]],axis=-1)
            maxdist = gsw.distance(colon,colat,axis=1).max()
            #
            lon2d,lat2d = np.meshgrid(lon1d.values,lat1d.values)
            #if (2*c+1)*(2*c+1) != lon2d.flatten().size:
            #    print((2*c+1)*(2*c+1), lon2d.flatten().size)
            colon = np.concatenate([plon.isel(x=i,y=j).values*np.ones((lsize,1)),lon2d.flatten()[:,np.newaxis]],axis=-1) #*np.ones(((2*c+1)*(2*c+1),1)),lon2d.flatten()[:,np.newaxis]],axis=-1)
            colat = np.concatenate([plat.isel(x=i,y=j).values*np.ones((lsize,1)),lat2d.flatten()[:,np.newaxis]],axis=-1) #np.ones(((2*c+1)*(2*c+1),1)),lat2d.flatten()[:,np.newaxis]],axis=-1)
            dist  = gsw.distance(colon,colat,axis=1).squeeze()
            inds  = np.where(dist<maxdist)[0]
            z     = elev.values.flatten()[inds]
            oinds = np.where(z<0)[0]
            linds = np.where(z>=0)[0]
            odepth_mm[j,i] = z[oinds].sum()
            ldepth_mm[j,i] = z[linds].sum()
            ocount_mm[j,i] = len(oinds)
            lcount_mm[j,i] = len(linds)

print('create memmaps')
odepth_mm = np.memmap(dpath+'odepth.mmap', dtype=float, shape=(grid.plat.shape), mode='w+')
ldepth_mm = np.memmap(dpath+'ldepth.mmap', dtype=float, shape=(grid.plat.shape), mode='w+')
ocount_mm = np.memmap(dpath+'ocount.mmap', dtype=float, shape=(grid.plat.shape), mode='w+')
lcount_mm = np.memmap(dpath+'lcount.mmap', dtype=float, shape=(grid.plat.shape), mode='w+')
#
odepth_mm[:] = np.zeros(grid.plat.shape)
ldepth_mm[:] = np.zeros(grid.plat.shape)
ocount_mm[:] = np.zeros(grid.plat.shape)
lcount_mm[:] = np.zeros(grid.plat.shape)
#
print('loop')
#for j in range(2): #nyt):
#    print(j)
#    #for i in range(nxt):
#        #find_depths(odepth_mm,ldepth_mm,ocount_mm,lcount_mm,elevation,lon,lat,plon,plat,j,i,nys,nxs)
Parallel(n_jobs=n_jobs)(delayed(find_depths)(odepth_mm,ldepth_mm,ocount_mm,lcount_mm,elevation,lon,lat,plon,plat,pclon,pclat,j,nxt,nys,nxs,2) for j in range(nyt))
#
odepth.values[:] = np.array(odepth_mm)
ldepth.values[:] = np.array(ldepth_mm)
ocount.values[:] = np.array(ocount_mm)
lcount.values[:] = np.array(lcount_mm) 
#
#print(odepth.min(),ocount.max(),ldepth.max(),lcount.max())
depth = (odepth/ocount)
depth.values[np.where(ocount<lcount)] = ((odepth+ldepth)/(ocount+lcount)).values[np.where(ocount<lcount)]
#
depth = -1.*(depth.where(depth<=0.).fillna(0.))
depth.to_dataset(name='depth').to_netcdf(dpath+'depth_python.nc')
#
#
# probably fastest to check the distance, but be smart about which points to check so that one doesn't spend too much time looking throuhg all points.
#for j in range(nyt): #range(nyt):
#    print(j)
#    for i in range(nxt):
#         odepth[j,i], ocount[j,i], ldepth[j,i], lcount[j,i] = find_depths(datain.elevation,datain.lon,datain.lat,grid.plon,grid.plat,j,i) 
#         #dlat  = abs(datain.lat-grid.plat.isel(y=j,x=i))
#         #dlon  = abs(datain.lon-grid.plon.isel(y=j,x=i))
#        #jj    = np.where(dlat==np.min(dlat))[0][0]
#         #ii    = np.where(dlon==np.min(dlon))[0][0]
#         #colon = np.concatenate([grid.plon.isel(x=i,y=j).values*np.ones((4,1)),grid.pclon.isel(x=i,y=j).values[:,np.newaxis]],axis=-1)
#         #colat = np.concatenate([grid.plat.isel(x=i,y=j).values*np.ones((4,1)),grid.pclat.isel(x=i,y=j).values[:,np.newaxis]],axis=-1)
#         #maxdist = gsw.distance(colon,colat,axis=1).max()
#         #dist = 0
#         #c    = 1
#         # figure out how large region to consider
#         #while maxdist>dist:
#         #   imax = ii+c
#         #   if imax>=nxs:
#         #       imax = imax-nxs
#         #   #
#         #   colon[0,1] = datain.lon[ii+c]
#         #   colon[1,1] = datain.lon[ii+c]
#         #   colon[2,1] = datain.lon[ii-c]
#         #   colon[3,1] = datain.lon[ii-c]
#         #   colat[0,1] = datain.lat[min(jj+c,nys)]
#         #   colat[1,1] = datain.lat[max(jj-c,0)]
#         #   colat[2,1] = datain.lat[min(jj+c,nys)]
#         #   colat[3,1] = datain.lat[max(jj-c,0)]
#         #   dist       = gsw.distance(colon,colat,axis=1).min()
#         #   c = c+1
#         # load the region and pick up the points that are closer than the bounding box
#         #lon1d = datain.lon[ii-c:ii+c+1].values
#         #lat1d = datain.lat[jj-c:jj+c+1].values
#         #lon2d,lat2d = np.meshgrid(lon1d,lat1d)
#         #colon = np.concatenate([grid.plon.isel(x=i,y=j).values*np.ones(((2*c+1)*(2*c+1),1)),lon2d.flatten()[:,np.newaxis]],axis=-1)
#         #colat = np.concatenate([grid.plat.isel(x=i,y=j).values*np.ones(((2*c+1)*(2*c+1),1)),lat2d.flatten()[:,np.newaxis]],axis=-1)
#         #dist  = gsw.distance(colon,colat,axis=1).squeeze()
#         #inds  = np.where(dist<maxdist)[0]
#         #z     = datain.elevation[jj-c:jj+c+1,ii-c:ii+c+1].values.flatten()[inds]
#         #
#         #odepth[j,i] = odepth[j,i] + z[np.where(z<0)].sum()
#         #ocount[j,i] = ocount[j,i] + len(np.where(z<0)[0])
#         #ldepth[j,i] = ldepth[j,i] + z[np.where(z>=0)].sum()
#         #lcount[j,i] = lcount[j,i] + len(np.where(z>=0)[0])
#         #odepth[jj,ii] = odepth[jj,ii] + max(mindep,datain.elevation[jj-c:jj+c+1,ii-c:ii+c+1].values.flatten()[inds].mean())
#         #
#        poly = Polygon([(grid.pclat[3,j,i],grid.pclon[3,j,i]),(grid.pclat[2,j,i],grid.pclon[2,j,i]),(grid.pclat[1,j,i],grid.pclon[1,j,i]),(grid.pclat[0,j,i],grid.pclon[0,j,i])])
#        po   = Point(grid.plat[j,i],grid.plon[j,i])
#
#mindep = 25.
#jmin = np.where(datain.lat>=grid.pclat.min())[0][0] #do not consider cells south of the plat range
#
#c=0
#for j in range(jmin,nys):
#    for i in range(nxs):
#        c=c+1
#        if c%1E6==0:
#           print(c,c/((nys-jmin)*nxs))
#        #
#        dlat = abs(grid.plat-datain.lat.isel(lat=j))
#        dlon = abs(grid.plon-datain.lon.isel(lon=i))
#        jj,ii = np.where(np.logical_and(dlat==np.min(dlat),dlon==np.min(dlon)))
#        #
#        if datain.elevation.isel(lat=j,lon=i)<0:
#            odepth[jj,ii] = odepth[jj,ii] + datain.elevation.isel(lat=j,lon=i)
#            ocount[jj,ii] = ocount[jj,ii] + 1
#        else:
#            ldepth[jj,ii] = ldepth[jj,ii] + datain.elevation.isel(lat=j,lon=i)
#            lcount[jj,ii] = lcount[jj,ii] + 1
#
#depth = (odepth/ocount)
#depth.values[np.where(ocount<lcount)] = ((odepth+ldepth)/(ocount+lcount)).values[np.where(ocount<lcount)]
#depth.values[np.where(depth>mindep)]  = mindep
#
#depth.to_dataset('depth').to_netcdf('depth_python.nc')
#
