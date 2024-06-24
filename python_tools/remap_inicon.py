import xarray as xr
import xesmf as xe
#import sys
#sys.path.append('/cluster/home/anu074/')
#from OMIP_utils import noresm_bounds
#from OMIP_utils import eosben07_tofsig
import numpy as np
from joblib import Parallel, delayed
#
def noresm_bounds(clon,clat):
    '''Return lon,lat bounds based on corner information clon,clat'''
    nv, ny,nx = clon.shape
    lon_b = np.concatenate([clon.values[0,:,:],clon.values[1,:,-1:]],axis=1)
    lon_b = np.concatenate([lon_b,np.concatenate([clon.values[3,-1:,:],clon.values[2,-1:,-1:]],axis=1)],axis=0)
    lat_b = np.concatenate([clat.values[0,:,:],clat.values[1,:,-1:]],axis=1)
    lat_b = np.concatenate([lat_b,np.concatenate([clat.values[3,-1:,:],clat.values[2,-1:,-1:]],axis=1)],axis=0)
    lon_b[np.where(lon_b>180)] = lon_b[np.where(lon_b>180)]-360
    lon_b[np.where(lon_b<-180)] = lon_b[np.where(lon_b<-180)]+360
    #
    lon_b = xr.DataArray(lon_b,dims=('y','x'),name='lon_b')
    lat_b = xr.DataArray(lat_b,dims=('y','x'),name='lat_b')
    #
    return lon_b, lat_b

#def adjust_bottom_depth_loop(j,nx,dz_micom,depth,pdepth,pmask):
#    '''
#    '''
#    iinds=np.where(pmask[j,:]>0)[0]
#    for i in iinds:
#        zi = np.where(depth[:,j,i]<=pdepth[j,i])[0][-1]
#        dz_micom[zi,j,i] = pdepth[j,i]-depth[zi-1,j,i]
#        #
#        dz_micom[zi+1:,j,i] = 0

def adjust_bottom_depth_loop2(j,nx,dz_micom,depth,pdepth,pmask):
    '''
    '''
    iinds=np.where(pmask>0)[0]
    for i in iinds:
        zi = np.where(depth[:,i]<=pdepth[i])[0][-1]
        if zi>0:
            dz_micom[zi,j,i]    = pdepth[i]-depth[zi-1,i]
            dz_micom[zi+1:,j,i] = 0
        else:
            dz_micom[zi+1,j,i]  = pdepth[i]-depth[0,i]
            dz_micom[zi+2:,j,i] = 0


#regrid_mode = 'patch'
regrid_mode = 'nearest_s2d'
inputpath   = '/cluster/shared/noresm/inputdata/ocn/micom/CLIMVOTE_34MA/case1/'
outputpath  = '/cluster/work/users/anu074/create_DOTPaleoGrid_34MA_1/' #'/cluster/shared/noresm/inputdata/ocn/micom/tnx0.125v1/'
dpath       = '/cluster/work/users/anu074/create_DOTPaleoGrid_34MA_1/' #'/cluster/shared/noresm/inputdata/ocn/micom/tnx0.125v1/'
n_jobs=4
#
ini = xr.open_dataset(inputpath+'inicon.nc')
#td  = xr.open_dataset(inputpath+'tidal_dissipation.nc')
#
grid_in = xr.open_dataset(inputpath+'grid.nc')
grid_out= xr.open_dataset(outputpath+'grid.nc')
#
lon_in_b,lat_in_b = noresm_bounds(grid_in.pclon,grid_in.pclat)
lon_out_b,lat_out_b = noresm_bounds(grid_out.pclon,grid_out.pclat)
#{'x':'x_b','y':'y_b'}
grid0=xr.merge([grid_in.plon.rename('lon'),grid_in.plat.rename('lat'),lon_in_b.rename({'x':'x_b','y':'y_b'}),lat_in_b.rename({'x':'x_b','y':'y_b'})])
grid1=xr.merge([grid_out.plon.rename('lon'),grid_out.plat.rename('lat'),lon_out_b.rename({'x':'x_b','y':'y_b'}),lat_out_b.rename({'x':'x_b','y':'y_b'})])
#
# we need to omit the last row because that is simply repetition of the second to last
j=1
yslice = slice(0,len(grid0.y)-j)
yslice1 = slice(0,len(grid1.y)-j)
yslice_b = slice(0,len(grid0.y_b)-j)
yslice1_b = slice(0,len(grid1.y_b)-j)
#
regridder = xe.Regridder(grid0.isel(y=yslice,y_b=yslice_b),grid1.isel(y=yslice1,y_b=yslice1_b),regrid_mode,filename=dpath+'INICON_'+regrid_mode+'_weights.nc',reuse_weights=False)
#
#yslice = slice(0,len(grid0.y)-1)
#yslice1 = slice(0,len(grid1.y)-1)
# regrid and copy the second to last row to the last row
mask1=regridder(grid_in.pmask.astype('float').isel(y=yslice))
for v,var in enumerate(['sigma','dz','saln','temp']):
    print(var)
    if v==0:
        dum    = regridder(ini[var].isel(y=yslice).where(grid_in.pmask.isel(y=yslice)>0).fillna(0.0))/mask1
        dum    = xr.concat([dum,dum.isel(y=slice(len(dum.y)-1,len(dum.y)))[:,::-1,:]],dim='y') 
        iniout = dum.where(grid_out.pmask>0).to_dataset(name=var)
    else:
        dum    = regridder(ini[var].isel(y=yslice).where(grid_in.pmask.isel(y=yslice)>0).fillna(0.0))/mask1
        dum    = xr.concat([dum,dum.isel(y=slice(len(dum.y)-1,len(dum.y)))[:,::-1,:]],dim='y')
        iniout = xr.merge([iniout,dum.where(grid_out.pmask>0).to_dataset(name=var)])

mask2=iniout.sigma.notnull()
nz,ny,nx = iniout.sigma.shape
jinds,iinds = np.where(abs(mask2.isel(z=0)-grid_out.pmask)>0)
##
#iniout['dz'] = iniout.dz.where(iniout.dz>0,other=0.)
# then  loop over those points and fill with nearest neighbor
for j in range(len(jinds)):
    c=1
    if j%500==0:
        print(100*j/len(jinds))
    jj=jinds[j]
    ii=iinds[j]
    yslice=slice(grid_out.jns[jj,ii].values-c-1,grid_out.jnn[jj,ii].values+c)
    xslice=slice(grid_out.inw[jj,ii].values-c-1,grid_out.ine[jj,ii].values+c)
    if xslice.start>xslice.stop or xslice.start<0:
        dum = xr.concat([iniout.saln.isel(y=yslice,x=slice(xslice.start,nx)),iniout.saln.isel(y=yslice,x=slice(0,xslice.stop))],dim='x')
    else:
        dum = iniout.saln.isel(y=yslice,x=xslice)
    #
    while dum.isel(z=0).notnull().max()==0:
        c=c+1
        yslice=slice(grid_out.jns[jj,ii].values-c-1,grid_out.jnn[jj,ii].values+c)
        xslice=slice(grid_out.inw[jj,ii].values-c-1,grid_out.ine[jj,ii].values+c)
        if xslice.start>xslice.stop or xslice.start<0:
            dum = xr.concat([iniout.saln.isel(y=yslice,x=slice(xslice.start,nx)),iniout.saln.isel(y=yslice,x=slice(0,xslice.stop))],dim='x')
        else:
            dum = iniout.saln.isel(y=yslice,x=xslice)
    #
    for v,var in enumerate(['sigma','dz','saln','temp']):
        if xslice.start>xslice.stop or xslice.start<0:
            iniout[var].values[:,jj,ii] = xr.concat([iniout[var].isel(y=yslice,x=slice(xslice.start,nx)),iniout[var].isel(y=yslice,x=slice(0,xslice.stop))],dim='x').mean(dim=('x','y')).values
        else:
            iniout[var].values[:,jj,ii] = iniout[var].isel(y=yslice,x=xslice).mean(dim=('x','y')).values

#iniout['dz'] = iniout.dz.where(iniout.dz>0,other=0.)
#mask3 = iniout.sigma.notnull()
#jinds2,iinds2 = np.where(abs(mask3.isel(z=0)-grid_out.pmask)>0)
#
# FIX BOTTOM DEPTH
print('adjust bottom depth')
dz_micom_mm = np.memmap(dpath+'dz_micom.mmap',  dtype=float, shape=(iniout.dz.shape), mode='w+')
dz_micom_mm[:] = iniout.dz.values.copy()
Parallel(n_jobs=n_jobs)(delayed(adjust_bottom_depth_loop2)(j,nx,dz_micom_mm,iniout.dz.cumsum(dim='z').isel(y=j).values,grid_out.pdepth.isel(y=j).values,grid_out.pmask.isel(y=j).values) for j in range(ny))
#Parallel(n_jobs=n_jobs)(delayed(adjust_bottom_depth_loop)(j,nx,dz_micom_mm,iniout.dz.cumsum(dim='z').values,grid_out.pdepth.values,grid_out.pmask.values) for j in range(ny))
iniout.dz.values = np.array(dz_micom_mm)
print('bottom depth adjusted')
#
#iniout.temp.values = eosben07_tofsig(iniout.sigma.values,iniout.saln.values,pref=2000.e5)
#for var in ['sigma','dz','saln','temp']:
#     iniout['var'] = xr.concat([iniout[var],iniout[var].isel(slice=)])
print('Fix mask')
iniout['sigma'].values = iniout.sigma.where(grid_out.pmask,other=-99.).values
iniout['saln'].values  = iniout.saln.where(grid_out.pmask,other=0.).values
iniout['temp'].values  = iniout.temp.where(grid_out.pmask,other=0.).values
iniout['dz'].values    = iniout.dz.where(grid_out.pmask,other=0.).values
#
print('write file')
iniout = iniout.drop_vars(('lat','lon','z'))
iniout.ffill('z').to_netcdf(outputpath+'inicon_'+regrid_mode+'.nc',format='NETCDF4',encoding={'dz':{'_FillValue':None},'temp':{'_FillValue':None},'saln':{'_FillValue':None},'sigma':{'_FillValue':None}})
# fixing bottom nans with ffill (now implemented above)
# iniout = xr.open_dataset('/projects/NS9869K/noresm/create_grid/DOTPaleo_hres/control_inicon_nearest_s2d.nc')
#outputpath='/projects/NS9869K/noresm/create_grid/DOTPaleo_hres/'
#iniout.ffill('z').to_netcdf(outputpath+'inicon_tnx0.25_DOTPaleo34Ma_control.nc',format='NETCDF3_64BIT',
#                            encoding={'dz':{'_FillValue':None},'temp':{'_FillValue':None},'saln':{'_FillValue':None},
#                                      'sigma':{'_FillValue':None}})
print('done')
# TIDAL MIXING
if False:
    twedon_PD = xr.open_dataset('/cluster/shared/noresm/inputdata/ocn/micom/tnx0.25v4/20170619/tidal_dissipation.nc')
    #
    twedon2 = regridder(twedon_PD.twedon.isel(y=slice(0,1152)).fillna(0))
    twedon2 = xr.concat([twedon2,xr.concat([twedon2.isel(y=slice(ny-1,ny),x=slice(nx//2,nx))[::-1], twedon2.isel(y=slice(ny-1,ny),x=slice(0,nx//2))[::-1]],dim='x')],dim='y')
    twedon2 = xr.DataArray(twedon2.values,dims={'y':ny,'x':nx},name='twedon')
    twedon2.to_dataset(name='twedon').to_netcdf(outputpath+'tidal_dissipation_interpolated.nc',format='NETCDF3_CLASSIC',encoding={'twedon':{'_FillValue':None}})
