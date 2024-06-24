import xarray as xr
import xesmf as xe
import numpy as np
#
dpath = '../grid_out/'
#
# BLOM GRID
#
grid    = xr.open_dataset(dpath+'grid.nc')
ny,nx = grid.plat.shape
lat_noresm = grid.plat.rename('lat').isel(y=slice(0,ny-1))
lon_noresm = grid.plon.rename('lon').isel(y=slice(0,ny-1))
lon_b_noresm = xr.concat([grid.pclon.isel(nv=0),grid.pclon.isel(nv=1,x=-1)],dim='x')
lon_b_noresm = xr.concat([lon_b_noresm,xr.concat([grid.pclon.isel(nv=3,y=-1),grid.pclon.isel(nv=2,y=-1,x=-1)],dim='x')],dim='y')
lat_b_noresm= xr.concat([grid.pclat.isel(nv=0),grid.pclat.isel(nv=1,x=-1)],dim='x')
lat_b_noresm = xr.concat([lat_b_noresm,xr.concat([grid.pclat.isel(nv=3,y=-1),grid.pclat.isel(nv=2,y=-1,x=-1)],dim='x')],dim='y')
#
lon_b_noresm = lon_b_noresm.isel(y=slice(0,ny)).rename('lon_b').rename({'y':'y_b','x':'x_b'})
lat_b_noresm = lat_b_noresm.isel(y=slice(0,ny)).rename('lat_b').rename({'y':'y_b','x':'x_b'})
#
NorESM_grid = xr.merge([lon_noresm,lat_noresm,lon_b_noresm,lat_b_noresm])
#
# ORIGINAL GRID
#
datain    = xr.open_dataset('../orig_data/paleobathy-topo_34.00Ma-mod-ant.nc') #paleobathy-topo_34.00Ma-gridline2.nc')
elevation = datain.z.load()
lat       = datain.lat.load().rename({'lat':'y'})
lon       = datain.lon.load().rename({'lon':'x'})
#
ny0,nx0 = elevation.shape
dx = 0.1
lon_b   = xr.DataArray(np.arange(-180,180+dx,0.1),dims=('x_b'),name='lon_b')
lat_b   = xr.DataArray(np.arange(-90,90+dx,0.1),dims=('y_b'),name='lat_b')
#
orig_grid = xr.merge([lon,lat,lon_b,lat_b])
#
# create the regridder
regridder = xe.Regridder(orig_grid, NorESM_grid,'conservative',filename=dpath+'orig_to_blom_conservative.nc',reuse_weights=False)
#
ocean_mask    = elevation.where(elevation<0).notnull().astype(float)
elevation_out = regridder(elevation.where(elevation<0,other=0))/regridder(ocean_mask)
elevation_out = -xr.concat([elevation_out,elevation_out.isel(y=slice(-1,ny))[:,::-1]],dim='y')
#
elevation_out.rename('z').to_netcdf(dpath+'depth_python_xesmf.nc')
