import numpy as np

def CurveCalc(obj, kt):
    obj.kt = kt
    if obj.DEM.ZFilt is None:
        SZ = obj.DEM.Z
        obj.DEM.ZFilt = SZ
        obj.Filter = [np.nan, np.nan]
    else:
        SZ = obj.DEM.ZFilt
    
    SX, SY = np.meshgrid(obj.DEM.X, obj.DEM.Y)
    du = SX[0, 1] - SX[0, 0]
    dv = SY[1, 0] - SY[0, 0]
    m, n = SZ.shape
    
    SXU = np.ones_like(SX) * obj.DEM.dx
    SXV = np.zeros_like(SX)
    SYU = np.zeros_like(SY)
    SYV = np.ones_like(SY) * obj.DEM.dy
    
    SZU, SZV = np.gradient(SZ, du, dv)
    
    SU = np.zeros((m, n, 3))
    SV = np.zeros((m, n, 3))
    SU[:, :, 0] = SXU
    SU[:, :, 1] = SYU
    SU[:, :, 2] = SZU
    SV[:, :, 0] = SXV
    SV[:, :, 1] = SYV
    SV[:, :, 2] = SZV
    
    E = np.sum(SU * SU, axis=2)
    F = np.sum(SU * SV, axis=2)
    G = np.sum(SV * SV, axis=2)
    
    CUV = np.cross(SU, SV, axis=2)
    AC = np.sqrt(np.sum(CUV ** 2, axis=2))
    NX = CUV[:, :, 0] / AC
    NY = CUV[:, :, 1] / AC
    NZ = CUV[:, :, 2] / AC
    
    NXU, NXV = np.gradient(NX, du, dv)
    NYU, NYV = np.gradient(NY, du, dv)
    NZU, NZV = np.gradient(NZ, du, dv)
    
    NU = np.zeros((m, n, 3))
    NV = np.zeros((m, n, 3))
    NU[:, :, 0] = NXU
    NU[:, :, 1] = NYU
    NU[:, :, 2] = NZU
    NV[:, :, 0] = NXV
    NV[:, :, 1] = NYV
    NV[:, :, 2] = NZV
    
    L = -np.sum(NU * SU, axis=2)
    M = -0.5 * (np.sum(NU * SV, axis=2) + np.sum(NV * SU, axis=2))
    N = -np.sum(NV * SV, axis=2)
    
    K1 = np.zeros_like(SZ)
    K2 = np.zeros_like(SZ)
    K1U = np.zeros_like(SZ)
    K1V = np.zeros_like(SZ)
    K2U = np.zeros_like(SZ)
    K2V = np.zeros_like(SZ)
    
    for i in range(m):
        for j in range(n):
            I = np.array([[E[i, j], F[i, j]], [F[i, j], G[i, j]]])
            II = np.array([[L[i, j], M[i, j]], [M[i, j], N[i, j]]])
            SO = np.linalg.inv(I) @ II
            A = np.max(np.isnan(SO))
            if A == 1:
                KD = np.array([[np.nan, np.nan], [np.nan, np.nan]])
                KM = np.array([[np.nan, np.nan], [np.nan, np.nan]])
            elif A == 0:
                KD, KM = np.linalg.eig(SO)
                if KM[0, 0] < KM[1, 1]:
                    KM = np.rot90(KM, 2)
                    KD = np.fliplr(KD)
            
            K1[i, j] = KM[0, 0]
            K2[i, j] = KM[1, 1]
            K1U[i, j] = KD[0, 0]
            K1V[i, j] = KD[1, 0]
            K2U[i, j] = KD[0, 1]
            K2V[i, j] = KD[1, 1]
    
    K1[np.abs(K1) <= kt] = 0
    K2[np.abs(K2) <= kt] = 0
    
    KG = K1 * K2
    KM = 0.5 * (K1 + K2)
    
    obj.SMAP = np.empty_like(KG) * np.nan
    obj.SDist = [[None] * 9, [None] * 9]
    
    in_ = np.where((KG < 0) & (np.abs(KM) <= 0))
    obj.SMAP[in_] = -4
    obj.SDist[0][0] = ['Perfect Saddles']
    obj.SDist[1][0] = np.size(in_) / np.size(obj.SMAP)
    
    in_ = np.where((KG > 0) & (KM < 0))
    obj.SMAP[in_] = -3
    obj.SDist[0][1] = ['Peaks']
    obj.SDist[1][1] = np.size(in_) / np.size(obj.SMAP)
    
    in_ = np.where((KG < 0) & (KM < 0))
    obj.SMAP[in_] = -2
    obj.SDist[0][2] = ['Antiformal Saddles']
    obj.SDist[1][2] = np.size(in_) / np.size(obj.SMAP)
    
    in_ = np.where((KG == 0) & (KM < 0))
    obj.SMAP[in_] = -1
    obj.SDist[0][3] = ['Antiforms']
    obj.SDist[1][3] = np.size(in_) / np.size(obj.SMAP)
    
    in_ = np.where((KG == 0) & (np.abs(KM) <= 0))
    obj.SMAP[in_] = 0
    obj.SDist[0][4] = ['Planes']
    obj.SDist[1][4] = np.size(in_) / np.size(obj.SMAP)
    
    in_ = np.where((KG == 0) & (KM > 0))
    obj.SMAP[in_] = 1
    obj.SDist[0][5] = ['Synforms']
    obj.SDist[1][5] = np.size(in_) / np.size(obj.SMAP)
    
    in_ = np.where((KG < 0) & (KM > 0))
    obj.SMAP[in_] = 2
    obj.SDist[0][6] = ['Synformal Saddles']
    obj.SDist[1][6] = np.size(in_) / np.size(obj.SMAP)
    
    in_ = np.where((KG > 0) & (KM > 0))
    obj.SMAP[in_] = 3
    obj.SDist[0][7] = ['Basins']
    obj.SDist[1][7] = np.size(in_) / np.size(obj.SMAP)
    
    obj.CMAP = {
        'KG': KG,
        'KM': KM,
        'K1': K1,
        'K2': K2,
        'K1U': K1U,
        'K1V': K1V,
        'K2U': K2U,
        'K2V': K2V
    }
