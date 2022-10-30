from grideddb import GridedDataBases
import xarray as xr

a = GridedDataBases(1959,2020,region="global")
a.read_nc("w_700") #Dependencies problems??????????????

print(xr.open_dataset("/datos/ERA5_updated/mon/ERA5_w700.nc")) #This works :o 