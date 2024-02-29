import pandas as pd
import numpy as np
# Parses data (i.e Z and LD)

class Data:
    def __init__(self,zfilename,ldfilename,n, approx_bf):
        zfile = readZfile(zfilename)
        self.z = zfile.z
        self.rsid = zfile.rsid
        self.rsid2i = {id:i for i,id in dict(self.rsid).items()} # self.rsid2i[rsXYZ] gives index of rsXYZ in the data
        self.pve = zfile.pve if 'pve' in zfile.columns else None
        ldfile = readLDfile(ldfilename, self.rsid)
        self.ld = ldfile
        print(ldfile)
        print(ldfile.shape)
        self.m = len(self.z)

        if not approx_bf:
            self.z *= 1/np.sqrt(self.z**2 + n - 2) # TODO: find why did we do what we did here??

def readLDfile(ldfilename, rsid = None):
    return pd.read_csv(ldfilename, sep = ' ', header = None).values

def readZfile(zfilename):
    zfile = pd.read_csv(zfilename, sep = ' ', header = None)
    zfile.columns = ['rsid','z','pve'][:len(zfile.columns)]
    return zfile