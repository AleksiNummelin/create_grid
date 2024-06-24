import xarray as xr
import numpy as np
from skimage import measure
#
# force minimum depth - cells that are shallower are made deeper
mindep = 25.
#
fpath ='../grid_out/' # '/cluster/work/users/anu074/create_Hres_grid/'
gridfile=fpath+'grid.nc' #'grid_tnx1v4.nc'
infile = fpath+'depth_python_xesmf.nc'
outfile = fpath+'depth_python_xesmf_regridded.nc'
#icefile = fpath+'../create_Hres_grid2022/N1850frc2NOOBGC_f09_tnx0125v4_test4.cice.h.0036-03.nc'
#icefile2 = fpath+'../create_Hres_grid2022/N1850frc2NOOBGC_f09_tnx0125v4_test4.blom.hm.0036-03.nc'
#ice = xr.open_dataset(icefile)
#ice2 = xr.open_dataset(icefile2)
#icemask = ice.hi.where(ice.hi>0.001).notnull().isel(time=0).compute().astype(int)
#icemask2 = ice2.fice.isel(time=0).where(ice2.fice.isel(time=0)>0,other=0).astype(int)
#icemask2 = icemask2[4:,:]
#
grid=xr.open_dataset(gridfile)
plat=grid.plat
plon=grid.plon
pmask=xr.zeros_like(grid.plat).rename('pmask')
insw=grid.insw-1
jnsw=grid.jnsw-1
ins=grid.ins-1
jns=grid.jns-1
inse=grid.inse-1
jnse=grid.jnse-1
ine=grid.ine-1
jne=grid.jne-1
inne=grid.inne-1
jnne=grid.jnne-1
inn=grid.inn-1
jnn=grid.jnn-1
innw=grid.innw-1
jnnw=grid.jnnw-1
inw=grid.inw-1
jnw=grid.jnw-1
# possible regions where one might have ice
icemask3=pmask.where(np.logical_or(np.logical_or(plat<-55,plat>68),
                     np.logical_or(np.logical_and(plat>50,np.logical_and(plon>90,plon<=180)),
                                   np.logical_and(plat>50,np.logical_and(plon>-120,plon<-18))))).notnull().astype(float)
#
data=xr.open_dataset(infile)
#data=data.rename({'depth':'z'})
#
# FILL IN ALL THE LAKES THAT JUST ONE GRID CELL WIDE/LARGE
def fill_small_lakes(z):
    ''' '''
    mask=z.where(z>0).notnull().astype(float) #ocean = 1, land=0
    jdm,idm = mask.shape
    # north, south, east, west
    narr = xr.concat([mask.roll(y=1,roll_coords=False),mask.roll(y=-1,roll_coords=False),mask.roll(x=1,roll_coords=False),mask.roll(x=-1,roll_coords=False)],dim='dum')
    # northeast, northwest, southeast, southwest
    narr = xr.concat([narr,narr.isel(dum=0).roll(x=1,roll_coords=False),
                      narr.isel(dum=0).roll(x=-1,roll_coords=False),
                      narr.isel(dum=1).roll(x=1,roll_coords=False),
                      narr.isel(dum=1).roll(x=-1,roll_coords=False)],dim='dum')
    #
    narrsum = narr.sum(dim='dum')
    narrmask = 1-narrsum.where(np.logical_and(narrsum<2,mask==1)).notnull().astype(float) # 0 if land or if only 1-gridd cell wide connection to ocean, else 1
    return z*narrmask #

znew = data.z.copy()
#
# we don't actually need this
# because the image processing will take care
# of this as well.
#zold = znew.where(znew>0)
#j=0
#while not znew.equals(zold):
#    print(j)
#    zold=znew.copy()
#    znew=fill_small_lakes(zold)
#    j=j+1
#
if False:
   # CONNECT BLACK SEA - maybe do this based on lat-lon info
   # dardanellit
   #znew.values[1436,929]=60
   #
   znew.values[1431,926]=60
   znew.values[1432,928]=60
   znew.values[1433,929]=60
   znew.values[1434,930]=60
   #
   znew.values[1436,929]=0
   # 
   # bosborus
   znew.values[1442:1448,947]=60
   # CONNECT BALTIC
   znew.values[1632,806]=25
   znew.values[1687:1689,874]=25
   znew.values[1685:1689,873]=25
   # Antarctic 
   znew.values[314,183]=400
   znew.values[328,182]=400
   znew.values[330,184]=400
   znew.values[327,184:186]=400
   znew.values[330,179:181]=400
   #Arctic
   #znew.values[2003,780]=200
   #znew.values[2009:2011,511]=200
   #znew.values[1975,487]=200
   #znew.values[1977,486]=200
   #znew.values[1981,485]=200
   #znew.values[2112:2114,353]=200
   #znew.values[2117,415]=200
   #znew.values[2118,471]=200
   #znew.values[2119,472:475]=200
   #znew.values[2146:2151,563]=200
   #znew.values[2146,559]=200
   #znew.values[2146:2150,555]=200
   #znew.values[2146:2150,554]=200
   #znew.values[2157,512:514]=200
   #znew.values[2040:2042,2467]=200
   #znew.values[2069,2472:2475]=200
   #znew.values[2069:2074,2474]=200
   #znew.values[2069:2074,2475]=200
   #znew.values[2095:2097,2463]=200
   #
   znew.values[2089,2461:2464]=0
   znew.values[2143,554]=0
   znew.values[2113:2116,457]=0
   #
   #znew.values[1633:1636,806]=25
   #znew.values[1635:1636,805]=25
   #znew.values[1645:1649,832]=25
#
# use a image processsing toolkit to connect the ocean and find the lakes
mask = znew.where(znew>0).notnull().fillna(0).astype(float)
ny,nx = mask.shape
# it seems to be better to wrap the Arctic - also padd some extra ocean on the left to connect part of the Gulf of Mexico
# exact choice of this padding will depend on the longitudes of the poles
mask2 = xr.concat([mask.isel(x=slice(0,nx//2)), mask.isel(x=slice(nx//2,nx))[::-1,::-1]],dim='y')
mask3 = xr.concat([mask.isel(x=slice(-nx//4,nx),y=slice(0,ny)), mask.isel(x=slice(0,nx//4),y=slice(0,ny))[::-1,::-1]],dim='y')
mask4 = xr.concat([mask3,mask2],dim='x')
mask4 = xr.concat([mask,mask[::-1,::-1]],dim='y')
mask4 = xr.concat([mask4,mask4],dim='x')
mask_new = xr.ones_like(mask)
#
ny2,nx2=mask4.shape
#
mask4_old = np.zeros(mask4.shape)
#
surroundings=[[jnn,inn],[jnne,inne],[jne,ine],[jnse,inse],[jns,ins],[jnsw,insw],[jnw,inw],[jnnw,innw],[jnn,inn]]
corner_aw=[[jnn,inn],[jne,ine],[jns,ins],[jnw,inw]]
corner_cw=[[jne,ine],[jns,ins],[jnw,inw],[jnn,inn]]
corners=[[jnne,inne],[jnse,inse],[jnsw,insw],[jnnw,innw]]
corners2=[[jnse,inse],[jnsw,insw],[jnnw,innw],[jnne,inne]]
corners3=[[jnsw,insw],[jnnw,innw],[jnne,inne],[jnse,inse]]
corner_aw2=[[jnw,inw],[jnn,inn],[jne,ine],[jns,ins]]
corner_cw2=[[jns,ins],[jnw,inw],[jnn,inn],[jne,ine]]
#
c=0
while not mask4.equals(mask4_old):
    print(c)
    c=c+1
    mask4_old = mask4.copy()
    # get rid of 1-point island
    island_labels = measure.label(1-mask4,connectivity=1)
    island_labels_flat = island_labels.flatten()
    for ilab in range(0,np.max(island_labels_flat)+1):
        i_ind=np.where(island_labels_flat==ilab)[0]
        if i_ind.size<2:
            mask4.values[np.where(island_labels==ilab)]=1
    # fill in lakes - add line of zeros south of southernmost array
    all_labels2 = measure.label(xr.concat([mask4[:1,:]*0.0,mask4,mask4[-1:,:]*0.0],dim='y'),connectivity=1)
    all_labels2 = all_labels2[1:-1,:]
    count1,bins1 = np.histogram(all_labels2,bins=all_labels2.max()+1)
    #count2,bins2 = np.histogram(all_labels2[ny2//2:,:][::-1,::-1],bins=all_labels2[ny2//2:,:].max())
    #all_labels3 = np.concatenate([all_labels2[:ny2//2,:,np.newaxis],all_labels2[ny2//2:,:,np.newaxis][::-1,::-1,:]],axis=-1)
    mask5_1=mask4[:ny2//2,:nx2//2].copy()
    mask5_2=mask4[:ny2//2,nx2//2:].copy()
    mask5_3=mask4[ny2//2:,:nx2//2][::-1,::-1].copy()
    mask5_4=mask4[ny2//2:,nx2//2:][::-1,::-1].copy()
    for j in np.argsort(count1)[:-2]: #this doesn't work, but it is quite close
        mask5_1.values[np.where(all_labels2[:ny2//2,:nx2//2]==j)]=0
        mask5_2.values[np.where(all_labels2[:ny2//2,nx2//2:]==j)]=0
        mask5_3.values[np.where(all_labels2[ny2//2:,:nx2//2][::-1,::-1]==j)]=0
        mask5_4.values[np.where(all_labels2[ny2//2:,nx2//2:][::-1,::-1]==j)]=0
    # reshape
    mask5 = xr.concat([mask5_1,mask5_2,mask5_3,mask5_4],dim='dum').max('dum')
    mask4 = xr.concat([mask5,mask5[::-1,::-1]],dim='y')
    mask4 = xr.concat([mask4,mask4],dim='x')
    #mask_out = np.concatenate([mask4[:ny2//2,nx//4:],mask4[ny2//2:,nx//4:][::-1,::-1]],axis=1)
    #plt.figure(); plt.pcolormesh(mask4);plt.colorbar();
    # make sure all grid cells sufficient connection to each other so that CICE will work on a b-grid.
    mask_out=xr.DataArray(mask5.values.astype(int),dims=('y','x'))
    #
    cmask = xr.concat([mask_out[jnsw,insw], mask_out[jns,ine], mask_out, mask_out[jnw,inw]],dim='dum').min('dum').compute()
    #csum will be <4 if any point surrounding i,j point (or that point itself) will be land - so basically identifies the coastal points
    csum = xr.concat([cmask,cmask[jne,ine],cmask[jnne,inne],cmask[jnn,inn]],dim='dum').sum('dum').compute()
    #nsum will count the ocean points orthogonal to the i,j point (4=all ocean)
    nsum = xr.concat([mask_out[jns,ins],mask_out[jne,ine],mask_out[jnn,inn],mask_out[jnw,inw]],dim='dum').sum('dum').compute()
    #emask1 will be >0 if there will be land surrounding the bottom-left, or top-right corners
    emask1 = xr.concat([cmask,cmask[jnne,inne]],dim='dum').sum('dum').compute()
    #emask1 will be >0 if there will be land surrounding the bottom-right, or top-left corners
    emask2 = xr.concat([cmask[jne,ine],cmask[jnn,inn]],dim='dum').sum('dum').compute()
    #
    mask_out = mask_out.values
    #mask_new.values = mask_out.copy()
    print(csum.dtype)
    # coastal points
    #jinds,iinds=np.where(icemask3.where(np.logical_and(csum.values<4,csum.values>0),other=0))
    jinds,iinds=np.where(csum.where(np.logical_and(csum.values<4,mask_out==1)).notnull().astype(int))
    if not mask4.equals(mask4_old):
        for cc in range(len(jinds)):
            j=jinds[cc]
            i=iinds[cc]
            if cc%10000==0:
                print(cc)
            #check also if one width channel jnw,inw==0 and jne,ine==0
            # or jnn,inn=0 and jns,ins==0, make it 2 depth
            #if (csum[j,i]==0|((csum[j,i]==1&nsum[j,i]==3))| \
            #                 ((csum[j,i]==2)&(emask1[j,i]==2|emask2[j,i]==2))):
            if True:
                    #corners=[[jnne,inne],[jnse,inse],[jnsw,insw],[jnnw,innw]]
                    #corners2=[[jnse,inse],[jnsw,insw],[jnnw,innw],[jnne,inne]]
                    #corners3=[[jnsw,insw],[jnnw,innw],[jnne,inne],[jnse,inse]]
                    #corner_aw=[[jnn,inn],[jne,ine],[jns,ins],[jnw,inw]]
                    #corner_cw=[[jne,ine],[jns,ins],[jnw,inw],[jnn,inn]]
                    #corner_aw2=[[jnw,inw],[jnn,inn],[jne,ine],[jns,ins]]
                    #corner_cw2=[[jns,ins],[jnw,inw],[jnn,inn],[jne,ine]]
                    #
                # fill in 1 width bays
                j2=j; i2=i; dum=True
                while dum:
                    if (mask_out[jnn[j2,i2],inn[j2,i2]]==0 and mask_out[jns[j2,i2],ins[j2,i2]]==0) and (mask_out[jnw[j2,i2],inw[j2,i2]]==0 or mask_out[jne[j2,i2],ine[j2,i2]]==0):
                        mask_out[j2,i2]=0
                        if (mask_out[jnnw[j2,i2],innw[j2,i2]]==0 and mask_out[jnsw[j2,i2],insw[j2,i2]]==0) and mask_out[jnw[j2,i2],inw[j2,i2]]==1:
                            i2=inw[j2,i2]
                        elif (mask_out[jnne[j2,i2],inne[j2,i2]]==0 and mask_out[jnse[j2,i2],inse[j2,i2]]==0) and mask_out[jne[j2,i2],ine[j2,i2]]==1:
                            i2=ine[j2,i2]
                        else:
                            dum=False
                    else:
                        dum=False
                j2=j; i2=i; dum=True
                while dum:
                    if (mask_out[jnw[j2,i2],inw[j2,i2]]==0 and mask_out[jne[j2,i2],ine[j2,i2]]==0) and (mask_out[jnn[j2,i2],inn[j2,i2]]==0 or mask_out[jns[j2,i2],ins[j2,i2]]==0):
                        mask_out[j2,i2]=0
                        if (mask_out[jnnw[j2,i2],innw[j2,i2]]==0 and mask_out[jnne[j2,i2],inne[j2,i2]]==0) and mask_out[jnn[j2,i2],inn[j2,i2]]==1:
                            j2=jnn[j2,i2]
                        elif (mask_out[jnsw[j2,i2],insw[j2,i2]]==0 and mask_out[jnse[j2,i2],inse[j2,i2]]==0) and mask_out[jns[j2,i2],ins[j2,i2]]==1:
                            j2=jns[j2,i2]
                        else:
                            dum=False
                    else:
                        dum=False
                # DEAL WITH CICE B-GRID ISSUES - this part of the code does not converge
                if icemask3[j,i]==1 and True:
                    # open narrow straits for CICE
                    if (mask_out[jnn[j,i],inn[j,i]]==0 and mask_out[jns[j,i],ins[j,i]]==0) and (mask_out[jnw[j,i],inw[j,i]]==1 and mask_out[jne[j,i],ine[j,i]]==1):
                    # open narrow straits for CICE
                    #if icemask3[j,i]==1 and (mask_out[jnw[j,i],inw[j,i]]==1 and mask_out[jne[j,i],ine[j,i]]==1):
                        if mask_out[corners[0][0][j,i],corners[0][1][j,i]]==1 or mask_out[corners[2][0][j,i],corners[2][1][j,i]]==1:
                            mask_out[jnn[j,i],inn[j,i]]=1
                        elif mask_out[corners[1][0][j,i],corners[1][1][j,1]]==1 or mask_out[corners[3][0][j,i],corners[3][1][j,i]]==1:
                            mask_out[jns[j,i],ins[j,i]]=1
                        #else:
                        #    mask_out[jns[j,i],ins[j,i]]=1
                    # fill in 1 width bays
                    #if (mask_out[jnw[j,i],inw[j,i]]==0 or mask_out[jne[j,i],ine[j,i]]==0):
                    #    mask_out[j,i]=0
                    if (mask_out[jnw[j,i],inw[j,i]]==0 and mask_out[jne[j,i],ine[j,i]]==0) and (mask_out[jnn[j,i],inn[j,i]]==1 and mask_out[jns[j,i],ins[j,i]]==1):
                    # open narrow straits for CICE
                    #if icemask3[j,i]==1 and (mask_out[jnn[j,i],inn[j,i]]==1 and mask_out[jns[j,i],ins[j,i]]==1):
                        if mask_out[corners[0][0][j,i],corners[0][1][j,i]]==1 or mask_out[corners[1][0][j,i],corners[1][1][j,i]]==1:
                            mask_out[jne[j,i],ine[j,i]]=1
                        elif mask_out[corners[2][0][j,i],corners[2][1][j,i]]==1 or mask_out[corners[3][0][j,i],corners[3][1][j,i]]==1:
                            mask_out[jnw[j,i],inw[j,i]]=1
                        #else:
                        #    mask_out[jnw[j,i],inw[j,i]]=1
                    # fill in 1 width bays
                    #if (mask_out[jnn[j,i],inn[j,i]]==0 or mask_out[jns[j,i],ins[j,i]]==0):
                    #    mask_out[j,i]=0
                    #Deal with narrow corner points that would cause trouble for CICE - open up if there is ocean around
            #if icemask3[j,i]==1:
                    for pc,corner in enumerate(corners):
                        # diagonal corners just 1-point away where rest of the points are ocean
                        if (mask_out[corner[0][j,i],corner[1][j,i]]==0 and mask_out[corners3[pc][0][j,i],corners3[pc][1][j,i]]==0) and (mask_out[corner_aw[pc][0][j,i],corner_aw[pc][1][j,i]]==1 and mask_out[corner_cw[pc][0][j,i],corner_cw[pc][1][j,i]]==1 and mask_out[corner_aw2[pc][0][j,i],corner_aw2[pc][1][j,i]]==1 and mask_out[corner_cw2[pc][0][j,i],corner_cw2[pc][1][j,i]]==1):
                            mask_out[corner[0][j,i],corner[1][j,i]]=1
                        #
                        if (mask_out[corner[0][j,i],corner[1][j,i]]==0 and mask_out[corner_aw[pc][0][j,i],corner_aw[pc][1][j,i]]==1 and mask_out[corner_cw[pc][0][j,i],corner_cw[pc][1][j,i]]==1) and (mask_out[corner_aw2[pc][0][j,i],corner_aw2[pc][1][j,i]]==0 or mask_out[corner_cw2[pc][0][j,i],corner_cw2[pc][1][j,i]]==0):
                            mask_out[corner[0][j,i],corner[1][j,i]]=1
                    #elif (mask_out[corner_cw[pc][0][j,i],corner_cw[pc][1][j,i]]==1 and mask_out[corner_aw[pc][0][j,i],corner_aw[pc][1][j,i]]==1 and mask_out[corner_cw2[pc][0][j,i],corner_cw2[pc][1][j,i]]==1) and (mask_out[corners[pc][0][j,i],corners[pc][1][j,i]]==0 or mask_out[corners2[pc][0][j,i],corners2[pc][1][j,i]]==0):
                    #    mask_out[corner_cw[pc][0][j,i],corner_cw[pc][1][j,i]]=1
            #
    mask_new.values = mask_out.copy() #np.reshape(mask_out_flat.astype(float),(ny,nx))
    #
    if False:
        mask2 = xr.concat([mask_new.isel(x=slice(0,nx//2)), mask_new.isel(x=slice(nx//2,nx))[::-1,::-1]],dim='y')
        mask3 = xr.concat([mask_new.isel(x=slice(-nx//4,nx),y=slice(0,ny)),mask_new.isel(x=slice(0,nx//4),y=slice(0,ny))[::-1,::-1]],dim='y')
        mask4 = xr.concat([mask3,mask2],dim='x')
    else:
        mask4 = xr.concat([mask_new,mask_new[::-1,::-1]],dim='y')
        mask4 = xr.concat([mask4,mask4],dim='x')
    #plt.figure(); plt.pcolormesh(mask4);plt.colorbar();
    #
#mask_out = np.concatenate([mask4[:ny2//2,nx//4:],mask4[ny2//2:,nx//4:][::-1,::-1]],axis=1)
    #
    #nmask1 = np.logical_and(pmask==1,icemask==1)
    #nmask2 = np.logical_and(csum==1,nsum==3)
    #nmask4 = np.logical_and(csum==2,np.logical_or(emask1==2,emask2==2))
    #emask5 = np.logical_and(emask,np.logicl
    #emask3=((pmask==1)&(icemask==1))&(csum==0|(csum==1&nsum==3)|(csum==2&(emask1==2|emask2==2)))
    #        if (edmask(i,j)==1|edmask(i,j)>3)&icmask(i,j)==1
    #            if  csum==0|
    #                (csum==1&nsum==3)| ...
    #                (csum==2&((cmask(i        ,j        ) ...
    #                +cmask(inne(i,j),jnne(i,j)))==2| ...
    #                (cmask(ine (i,j),jne (i,j)) ...
    #                +cmask(inn (i,j),jnn (i,j)))==2))
    #                switch edmask(i,j)
    #                    case 1
    #                edmask(i,j)=5;
    #                    case 4
    #                edmask(i,j)=6;
    #                otherwise
    #                edmask(i,j)=edmask(i,j)-32;
    # 
    #
#all_labels2 = measure.label(mask4,connectivity=1) # THIS IS A BIT RESTRICTIVE, BUT ONLY ALLOW FOR E-W AND N-S CONNECTIONS
#all_labels3 = np.concatenate([all_labels2[:ny2//2,nx//4:],all_labels2[ny2//2:,nx//4:][::-1,::-1]],axis=1)
#
# land and ocean will have labels 0 and 1. Fill the rest.
#for j in range(2,np.max(all_labels3)+1):
#    mask.values[np.where(all_labels3==j)]=0
#
#icmask(find(plat<-60&pmask==1))=1;
#######################################
# save the  data
depth   = znew.where(mask_out).fillna(0)
# limit minimum depth
depth   = xr.concat([depth,(mindep*xr.ones_like(depth)).where(mask_out).fillna(0)],dim='dum').max(dim='dum')
depth2  = xr.concat([depth.isel(y=slice(0,-1)),depth.isel(y=slice(-2,-1))[:,::-1]],dim='y')
#
dataout = xr.merge([depth2,grid.plon,grid.plat])
# save file
dataout.to_netcdf(outfile,format='NETCDF4')
