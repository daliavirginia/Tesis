import xarray as xr

class GridedDataBases():

    """
    Class for reading and computing anomalies from my nc files. 
    Variables:
        - Zonal wind at 850 hPa (ERA5)
        - Meridional wind at 850 hPa (ERA5)
        - Vertical wind at 700 hPa (ERA5)
        - Vertical integrated moisture divergence (ERA5)
        - Specific humidity at 850 hPa (ERA5)
        - Mean sea level pressure (ERA5)
        - Sea surface temperature (NOAA ERSST)
        - Precipitation (CMAP)
        - Geopotential height at 200, 500 and 850 hPa (ERA5)

    Anomalies are calculated for base period 1981-2010.
    """

    def __init__(self,aniomin,aniomax, region="global",lonmin=None,lonmax=None,latmin=None,latmax=None):
        
        self.aniomin=aniomin
        self.aniomax=aniomax
        self.variables = ["u_850","v_850","w_700","vimd","q_850","msl","sst","pre","z_200","z_500","z_850"]
        self.paths={
            "u_850":"/datos/Datos_Dalia/era5.u.v.w.q.southamerica.1959.2021.nc",
            "v_850":"/datos/Datos_Dalia/era5.u.v.w.q.southamerica.1959.2021.nc",
            "w_700":"/datos/ERA5_updated/mon/ERA5_w700.nc",
            "vimd":"/datos/ERA5_updated/mon/ERA5_vimd_2022.nc",
            "z_200":"/datos/ERA5_updated/mon/ERA5_zg_2022.nc",
            "z_500":"/datos/ERA5_updated/mon/ERA5_zg_2022.nc",
            "z_850":"/datos/ERA5_updated/mon/ERA5_zg_2022.nc",
            "q_850":"/datos/Datos_Dalia/era5.u.v.w.q.southamerica.1959.2021.nc",
            "msl":'/datos/ERA5_updated/mon/ERA5_slp_2022.nc',
            "sst":'/datos/Datos_Dalia/sst.mnmean.nc',
            "precip":'/datos/Datos_Dalia/precip.mon.total.v2018.nc'}

        self.regiones = ["global", "southam","sh","manual"]
        self.region = region
        self.KEYS = ['Verano (DEF)', 'Otoño (MAM)', 'Invierno (JJA)', 'Primavera (SON)']

        if self.region == "global":
            self.lonmin=-180
            self.lonmax=180
            self.latmin=-90
            self.latmax=90
        elif self.region == "southam": 
            self.lonmin=-90
            self.lonmax=-30
            self.latmin=-60
            self.latmax=15           
        elif self.region == "sh":
            self.lonmin=-180
            self.lonmax=180
            self.latmin=-90
            self.latmax=0
        elif self.region == "manual":
            self.lonmin=lonmin
            self.lonmax=lonmax
            self.latmin=latmin
            self.latmax=latmax
        else:
            print("Las opciones de recorte espacial son: %s"%(self.regiones))

    def __str__(self,):

        t="""
        Clase para manipular y obtener facilmente datos reticulados que estan
        en el servidor Vegeta.

        Las variables disponibles son: %s
        Los recortes disponibles son: %s
        El año minimo y máximo se deben declarar en el constructor
        """%(self.variables,self.regiones)

        return t

    def read_nc(self,var):

        """
        Read netCDF file acording to selected variable
        (var). Since databases uses different coords names
        and longitude conventions, the method calls another
        which is special for era5 (read_era5), noaa-ersst(read_ersst)
        or cmap (read_cmap). 
        """

        #Obtengo el path
        path = self.paths[var]
        level=-1

        #Separo la variable donde haya un guion bajo
        splt=var.split("_")

        #Me quedo con la parque me interesa de la variable
        if len(splt)!=1:
            if splt[0]=="z":
                var="z"
                level=splt[1]
            else:
                var=splt[0]

        #Veo de que base de datos es y llamo a la función
        #correspondiente
        print("VARIABLE:%s"%(var))
        era5 = (var=="w" or var=="u" or var=="v" or var=="z" or var=="vimd" or var=="q" or var=="slp")
        cmap = var=="precip"
        noaa = var=="sst"

        if noaa:
            datos = self.read_ersst(var,path,self.aniomin,self.aniomax,level)
        elif cmap:
            datos = self.read_cmap(var,path,self.aniomin,self.aniomax)
        elif era5:
            datos = self.read_era5(var,path,self.aniomin,self.aniomax)
        else:
            print("Ingresaste cualquier variable pa")
            print("Las opciones son: ")
            for v in self.variables:
                print(v)
            datos=None

        return datos
        
    def read_era5(self,path,var,aniomin,aniomax,level=-1):

        timemin="%s-01-01"%(aniomin)
        timemax="%s-12-31"%(aniomax)

        # Abro el archivo usando xarray
        da=xr.open_dataset(path)[var] 

        if var=="z":
            da.sel(level=level)

        #Recorte espacial y temporal
        da= da.sel(
            time=slice(timemin,timemax), 
            latitude=slice(self.latmax, self.latmin),
            longitude=slice(self.lonmin,self.lonmax))


        #Agrupando por meses, le resto la media climatologica mensual a los datos. Asi
        #obtengo anomalías de da. 
        da_est=da.groupby('time.month')
        da_mean=da.sel(time=slice("1981-01-01","2010-12-31")).groupby('time.month').mean('time')
        da_anom=da_est - da_mean

        #Remuestreo de los datos por estación (DJF-MAM-JJA-SON)
        da_anom=da_anom.resample(time='QS-DEC').mean()

        #Guardo los datos separados por estación en un diccionario
        da_anom_est={
            self.KEYS[0]:da_anom.sel(time=da_anom["time.month"]==12),
            self.KEYS[1]:da_anom.sel(time=da_anom["time.month"]==3),
            self.KEYS[2]:da_anom.sel(time=da_anom["time.month"]==6),
            self.KEYS[3]:da_anom.sel(time=da_anom["time.month"]==9)}
        

        #Excluyo el primer y ultimo datos de verano (Tienen meses creados que no estan en los datos)
        da_anom_est[self.KEYS[0]]=da_anom_est[self.KEYS[0]].sel(
            time=slice(da_anom_est[self.KEYS[0]]['time'][0],
            da_anom_est[self.KEYS[0]]['time'][-1]))

        return da_anom_est

    def read_cmap(self,path,var,aniomin,aniomax):

        timemin="%s-01-01"%(aniomin)
        timemax="%s-12-31"%(aniomax)

        if self.lonmin<0 and self.lonmax<0:
            lonmin=self.lonmin+360
            lonmax=self.lonmax+360
        else:
            lonmin=self.lonmin
            lonmax=self.lonmax

        #Abro el archivo usando xarray
        ds_precip=xr.open_dataset(path)[var]
        #Recorte espacial y temporal
        precip_recorte = ds_precip.sel(time=slice(timemin,timemax), lat=slice(self.latmax, self.latmin),
                                    lon=slice(lonmin, lonmax))

        #Uso resample para que la fecha sea la ultima de cada mes, y luego extraigo el día
        precip_recorte = precip_recorte.resample(time='M').mean()
        dias=precip_recorte['time.day'][:]

        precip_recorte = precip_recorte.resample(time='MS').mean()
        #Multiplico la precip mensual por la cantidad de dias (paso de mm/day a mm/month)
        for i in range(len(dias)):
            precip_recorte[i,:,:] = precip_recorte[i,:,:]*dias[i]
        #Agrupando por meses, le resto la media climatologica mensual a los datos. Asi
        #obtengo anomalías de precip. 
        precip=precip_recorte.groupby('time.month')
        precip_mean=precip_recorte.sel(time=slice("1981-01-01","2010-12-31")).groupby('time.month').mean('time')
        precip_anom=precip-precip_mean

        #Remuestreo de los datos por estación (DJF-MAM-JJA-SON)
        precip_anom_est=precip_anom.resample(time='QS-DEC').sum()
        #Guardo los datos separados por estación en un diccionario
        precip_anom_est={self.KEYS[0]:precip_anom.sel(time=precip_anom["time.month"]==12),
                    self.KEYS[1]:precip_anom.sel(time=precip_anom["time.month"]==3),
                    self.KEYS[2]:precip_anom.sel(time=precip_anom["time.month"]==6),
                    self.KEYS[3]:precip_anom.sel(time=precip_anom["time.month"]==9)}
        

        #Excluyo el primer y ultimo datos de verano (Tienen meses creados que no estan en los datos)
        precip_anom_est[self.KEYS[0]]=precip_anom_est[self.KEYS[0]].sel(time=slice(precip_anom_est[self.KEYS[0]]['time'][0],
                                                                    precip_anom_est[self.KEYS[0]]['time'][-1]))
        print(precip_anom_est[self.KEYS[0]])
        del ds_precip, precip_anom, precip_recorte
        
        return precip_anom_est
    
    def read_ersst(self,path,var,aniomin,aniomax,lonmin,lonmax,latmin,latmax):

        timemin="%s-01-01"%(aniomin)
        timemax="%s-12-31"%(aniomax)
        
        if lonmin<0:
            lonmin=lonmin+360
        if lonmax<0:
            lonmax=lonmax+360

        # Abro el archivo usando xarray
        ds_tsm=xr.open_dataset(path)[var]

        #Recorte temporal
        tsm_recorte = ds_tsm.sel(
            time=slice(timemin,timemax),
            lat=slice(latmax,latmin),
            lon=slice())

        #Agrupando por meses, le resto la media climatologica mensual a los datos. Asi
        #obtengo anomalías de tsm. 
        tsm=tsm_recorte.groupby('time.month')
        tsm_media = tsm_recorte.sel(time=slice("1981-01-01","2010-12-31")).groupby('time.month').mean('time')
        tsm_anom = tsm-tsm_media


        #Remuestreo de los datos por estación (DJF-MAM-JJA-SON)
        tsm_anom_est=tsm_anom.resample(time='QS-DEC').mean()

        #Guardo los datos separados por estación en un diccionario
        tsm_anom_est={self.KEYS[0]:tsm_anom.sel(time=tsm_anom["time.month"]==12),
                    self.KEYS[1]:tsm_anom.sel(time=tsm_anom["time.month"]==3),
                    self.KEYS[2]:tsm_anom.sel(time=tsm_anom["time.month"]==6),
                    self.KEYS[3]:tsm_anom.sel(time=tsm_anom["time.month"]==9)}
        

        #Excluyo el primer y ultimo datos de verano (Tienen meses creados que no estan en los datos)
        tsm_anom_est[self.KEYS[0]]=tsm_anom_est[self.KEYS[0]].sel(time=slice(tsm_anom_est[self.KEYS[0]]['time'][0],
                                                                tsm_anom_est[self.KEYS[0]]['time'][-2]))
        #print(tsm_anom_est[KEYS[0]].time)

        del ds_tsm, tsm_anom, tsm_recorte

        return tsm_anom_est

        
    

