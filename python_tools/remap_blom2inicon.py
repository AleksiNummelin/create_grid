import xarray as xr
import xesmf as xe
#import sys
#sys.path.append('/cluster/home/anu074/')
#from OMIP_utils import noresm_bounds
#from OMIP_utils import eosben07_tofsig
import BLOM_utils as butils
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

def adjust_bottom_depth(dz_in,depth,pdepth,pmask):
    '''
    '''
    dz=dz_in.copy()
    if pmask>0:
        zi = np.where(depth<=pdepth)[0][-1]
        if zi>0:
            dz[zi]    = pdepth-depth[zi-1]
            dz[zi+1:] = 0
        else:
            dz[zi+1]  = pdepth-depth[0]
            dz[zi+2:] = 0
    return dz

#regrid_mode = 'patch'
# /nird/projects/NS9874K/noresm/cases/NBF_CLIMVOTE_34MA_1850_case1/ocn/hist/NBF_CLIMVOTE_34MA_1850_case1.micom.hm.2500-01.nc
# /nird/projects/NS9874K/noresm/cases/NBF_CLIMVOTE_34MA_1850_control/ocn/hist/NBF_CLIMVOTE_34MA_1850_control.micom.hm.2500-01.nc
regrid_mode = 'nearest_s2d'
if False: #on Fram
    inputfile   ='/cluster/work/users/anu074/CLIMVOTE_inicon/inputdata/NBF_CLIMVOTE_34MA_1850_control.micom.hm.2500-01.nc'
    inputpath   = '/cluster/shared/noresm/inputdata/ocn/micom/CLIMVOTE_34MA/case1/'
    outputpath  = '/cluster/work/users/anu074/create_DOTPaleoGrid_34MA_1/'
else:
    case='control'
    inputfile   = '/projects/NS9874K/noresm/cases/NBF_CLIMVOTE_34MA_1850_'+case+'/ocn/hist/NBF_CLIMVOTE_34MA_1850_'+case+'.micom.hm.2500-01.nc'
    inputpath   = '/projects/NS9874K/noresm/grids/inicon_DOTPaleo025/'+case+'/'
    outputpath  = '/projects/NS9874K/noresm/grids/inicon_DOTPaleo025/'+case+'/' #'/cluster/shared/noresm/inputdata/ocn/micom/tnx0.125v1/'
    dpath       = '/projects/NS9874K/noresm/grids/inicon_DOTPaleo025/'+case+'/' #'/cluster/shared/noresm/inputdata/ocn/micom/tnx0.125v1/'
n_jobs=8
#
data = xr.open_dataset(inputfile)
P   = data.dp.cumsum('sigma')
#rho = butils.eosben07_rho(P.isel(time=0)/1E4,data.temp.isel(time=0),
#                          data.saln.where(data.saln>0,other=0).isel(time=0)).rename('rho').chunk({'sigma': data.sigma.size,'x':data.x.size,'y':data.y.size}).compute()
# THIS DOESN'T WORK - WE NEED TO REMAP THE DATA TO DENSITY SPACE
rho = (xr.ones_like(data.saln.isel(time=0))*data.sigma).rename('sigma')
#
ini = xr.merge([data.dz.isel(time=0).drop('sigma').drop('time').rename({'sigma':'z'}).to_dataset(name='dz'),
                data.temp.isel(time=0).drop('sigma').drop('time').rename({'sigma':'z'}).to_dataset(name='temp'),
                data.saln.where(data.saln>0,other=0).isel(time=0).drop('sigma').drop('time').rename({'sigma':'z'}).to_dataset(name='saln'),
                rho.drop('sigma').rename({'sigma':'z'}).to_dataset(name='sigma')]) 
#td  = xr.open_dataset(inputpath+'tidal_dissipation.nc')
#
grid_in = xr.open_dataset(inputpath+'grid_old.nc')
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
regridder = xe.Regridder(grid0.isel(y=yslice,y_b=yslice_b),grid1.isel(y=yslice1,y_b=yslice1_b),
                         regrid_mode,filename=dpath+'INICON_'+regrid_mode+'_weights.nc',reuse_weights=True)
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
        if var in ['dz']:
            dum    = regridder(ini[var].cumsum('z').isel(y=yslice).where(grid_in.pmask.isel(y=yslice)>0).fillna(0.0))/mask1
            dum_dz = regridder(ini[var].cumsum('z').isel(y=yslice).where(grid_in.pmask.isel(y=yslice)>0).fillna(0.0))/mask1
            dum = dum.where(dum_dz>0.)
            #dum = xr.concat([xr.zeros_like(dum.isel(z=slice(0,1))),dum],dim='z').diff('z')
        else:
            dum    = regridder(ini[var].isel(y=yslice).where(grid_in.pmask.isel(y=yslice)>0).fillna(0.0))/mask1
        dum    = xr.concat([dum,dum.isel(y=slice(len(dum.y)-1,len(dum.y)))[:,::-1,:]],dim='y')
        if var in ['temp','saln']:
            dum=dum.where(dum!=0.0)
            #dum=dum.ffill('z')
        #elif var in ['dz']:
        #    dum=dum.where(dum>0)
        iniout = xr.merge([iniout,dum.where(grid_out.pmask>0).to_dataset(name=var)])
# make this a loop over all the density classes!
mask2    = iniout.saln.notnull()
mask3    = mask2.isel(z=0).load()
nz,ny,nx = iniout.sigma.shape
jinds,iinds = np.where(abs(mask3.astype(float)-grid_out.pmask.astype(float))>0)
##
#iniout['dz'] = iniout.dz.where(iniout.dz>0,other=0.)
# then  loop over those points and fill with nearest neighbor
# the strategy is to look for ever expanding rectangular area
# here
jns=grid_out.jns.values-1
jnn=grid_out.jnn.values-1
ine=grid_out.ine.values-1
inw=grid_out.inw.values-1
for z in range(nz):
    print(z)
    mask3    = mask2.isel(z=z).load()
    jinds,iinds = np.where(np.logical_and(abs(mask3.astype(float)-grid_out.pmask.astype(float))>0,iniout.dz.isel(z=z)<grid_out.pdepth))
    for j in range(len(jinds)):
        c=1
        if j%500==0:
            print(100*j/len(jinds))
        jj=jinds[j]
        ii=iinds[j]
        yslice=slice(max(0,jns[jj,ii]-c),jnn[jj,ii]+c+1)
        xslice=slice(inw[jj,ii]-c,ine[jj,ii]+c+1)
        if xslice.start>xslice.stop or xslice.start<0:
            dum = xr.concat([mask3.isel(y=yslice,x=slice(xslice.start,nx)),
                             mask3.isel(y=yslice,x=slice(0,xslice.stop))],dim='x')
        else:
            dum = mask3.isel(y=yslice,x=xslice)
        # stop once you find at least one grid cell that is not masked
        while dum.max()==0 and c<200:
            c=c+1
            yslice=slice(max(jns[jj,ii]-c,0),jnn[jj,ii]+c+1)
            xslice=slice(inw[jj,ii]-c,ine[jj,ii]+c+1)
            if xslice.start>xslice.stop or xslice.start<0:
                dum = xr.concat([mask3.isel(y=yslice,x=slice(xslice.start,nx)),
                             mask3.isel(y=yslice,x=slice(0,xslice.stop))],dim='x')
            else:
                dum = mask3.isel(y=yslice,x=xslice)
        #
        #print(c)
        for v,var in enumerate(['dz','saln','temp']):
            if xslice.start>xslice.stop or xslice.start<0:
                 #if var in ['dz']:
                 #    iniout[var].values[:,jj,ii] = xr.concat([iniout[var].cumsum('dz').isel(y=yslice,x=slice(xslice.start,nx)),
                 #                                         iniout[var].cumsum('dz').isel(y=yslice,x=slice(0,xslice.stop))],dim='x').mean(dim=('x','y')).values
                 #else:
                 iniout[var].values[z,jj,ii] = xr.concat([iniout[var].isel(z=z,y=yslice,x=slice(xslice.start,nx)),
                                                     iniout[var].isel(z=z,y=yslice,x=slice(0,xslice.stop))],dim='x').mean(dim=('x','y')).values
            else:
                iniout[var].values[z,jj,ii] = iniout[var].isel(z=z,y=yslice,x=xslice).mean(dim=('x','y')).values
        mask3    = iniout.saln.isel(z=0).notnull().load()

for var in ['temp','saln','dz']:
    iniout[var] = iniout[var].ffill('z')
#iniout['dz'] = iniout.dz.where(iniout.dz>0,other=0.)
#mask3 = iniout.sigma.notnull()
#jinds2,iinds2 = np.where(abs(mask3.isel(z=0)-grid_out.pmask)>0)
#
# FIX BOTTOM DEPTH
print('adjust bottom depth')
#dz_micom_mm = np.memmap(dpath+'dz_micom.mmap',  dtype=float, shape=(iniout.dz.shape), mode='w+')
#dz_micom_mm[:] = iniout.dz.values.copy()
#Parallel(n_jobs=n_jobs)(delayed(adjust_bottom_depth_loop2)(j,nx,dz_micom_mm,iniout.dz.cumsum(dim='z').isel(y=j).values,grid_out.pdepth.isel(y=j).values,grid_out.pmask.isel(y=j).values) for j in range(ny))
#Parallel(n_jobs=n_jobs)(delayed(adjust_bottom_depth_loop)(j,nx,dz_micom_mm,iniout.dz.cumsum(dim='z').values,grid_out.pdepth.values,grid_out.pmask.values) for j in range(ny))
#iniout.dz.values = np.array(dz_micom_mm)
iniout.dz.values[:,:,:] = xr.concat([xr.zeros_like(iniout.dz.isel(z=slice(0,1))),iniout.dz],dim='z').diff('z').values
depth=iniout.dz.cumsum(dim='z')
chunks={'x':nx//n_jobs,'y':ny//n_jobs}
dz = xr.apply_ufunc(adjust_bottom_depth,iniout.dz.chunk(chunks),depth.chunk(chunks),
                    grid_out.pdepth.chunk(chunks),grid_out.pmask.chunk(chunks),
                    input_core_dims=[['z'],['z'],[],[]],
                    output_core_dims=[['z']],
                    output_sizes=({'z':nz}),
                    vectorize=True,
                    dask='parallelized',
                    kwargs={},
                    output_dtypes=[depth.dtype],
).compute()
iniout.dz.values = dz.transpose('z','y','x').values
print('bottom depth adjusted')
#
#iniout.temp.values = eosben07_tofsig(iniout.sigma.values,iniout.saln.values,pref=2000.e5)
#for var in ['sigma','dz','saln','temp']:
#     iniout['var'] = xr.concat([iniout[var],iniout[var].isel(slice=)])
print('Fix mask')
# limit salinty to be 0 at min and recalculate sigma based on t,s
iniout['sigma'].values = iniout.sigma.where(grid_out.pmask,other=-99.).values
iniout['saln'].values  = iniout.saln.where(np.logical_and(grid_out.pmask,iniout.saln>0),other=0.).values
iniout['temp'].values  = iniout.temp.where(grid_out.pmask,other=0.).values
iniout['dz'].values    = iniout.dz.where(grid_out.pmask,other=0.).values
#
print('write file')
iniout = iniout.drop('time')
iniout.to_netcdf(outputpath+'inicon_'+regrid_mode+'.nc',format='NETCDF4',
                 encoding={'dz':{'_FillValue':None},
                           'temp':{'_FillValue':None},
                           'saln':{'_FillValue':None},'sigma':{'_FillValue':None}})
#
print('done')
# TIDAL MIXING
if False:
    twedon_PD = xr.open_dataset('/cluster/shared/noresm/inputdata/ocn/micom/tnx0.25v4/20170619/tidal_dissipation.nc')
    #
    twedon2 = regridder(twedon_PD.twedon.isel(y=slice(0,1152)).fillna(0))
    twedon2 = xr.concat([twedon2,xr.concat([twedon2.isel(y=slice(ny-1,ny),x=slice(nx//2,nx))[::-1], twedon2.isel(y=slice(ny-1,ny),x=slice(0,nx//2))[::-1]],dim='x')],dim='y')
    twedon2 = xr.DataArray(twedon2.values,dims={'y':ny,'x':nx},name='twedon')
    twedon2.to_dataset(name='twedon').to_netcdf(outputpath+'tidal_dissipation_interpolated.nc',format='NETCDF3_CLASSIC',encoding={'twedon':{'_FillValue':None}})
