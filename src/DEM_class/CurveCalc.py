import numpy as np 

def CurveCalc(ZFilt, dx, dy, kt):
    """
    Calculate principal curvatures and curvature features.

    Args:
        ZFilt (numpy.ndarray): Filtered surface data.
        dx (float): Grid spacing in the x direction.
        dy (float): Grid spacing in the y direction.
        kt (float): Threshold value.

    Returns: 
        K1, K2 (tuples): Principal curvatures.
        KM (tuple): Mean curvature.
        KG (tuple): Gaussian curvature.
    
    """

    r = len(ZFilt)  # Number of rows in ZFilt (surface data)
    c = len(ZFilt[0])  # Number of columns in ZFilt
    
    X = np.arange(c) * dx  # Generate 1D array of x coordinates
    Y = np.arange(r) * dy  # Generate 1D array of y coordinates
    
    SX, SY = np.meshgrid(X, Y)  # Create 2D arrays of x and y coordinates
    
    du = SX[0, 1] - SX[0, 0]  # Grid spacing in the u direction
    dv = SY[1, 0] - SY[0, 0]  # Grid spacing in the v direction
    
    SZU, SZV = np.gradient(ZFilt, du, dv)  # Compute gradients of ZFilt
    
    SXU = np.ones_like(SX) * dx  # Derivative of x component of surface vector in the u direction
    SXV = np.zeros_like(SX)  # Derivative of x component of surface vector in the v direction
    SYU = np.zeros_like(SY)  # Derivative of y component of surface vector in the u direction
    SYV = np.ones_like(SY) * dy  # Derivative of y component of surface vector in the v direction
    
    SU = np.zeros((r, c, 3))  # Initialize array to store surface vector derivatives
    SV = np.zeros((r, c, 3))  # Initialize array to store surface vector derivatives
    SU[:, :, 0] = SXU  # Store x component of surface vector derivative in SU
    SU[:, :, 1] = SYU  # Store y component of surface vector derivative in SU
    SU[:, :, 2] = SZU  # Store z component of surface vector derivative in SU
    SV[:, :, 0] = SXV  # Store x component of surface vector derivative in SV
    SV[:, :, 1] = SYV  # Store y component of surface vector derivative in SV
    SV[:, :, 2] = SZV  # Store z component of surface vector derivative in SV
    
    E = np.sum(SU * SU, axis=2)  # Compute E tensor component
    F = np.sum(SU * SV, axis=2)  # Compute F tensor component
    G = np.sum(SV * SV, axis=2)  # Compute G tensor component
    
    al = F * G - F ** 2  # Compute intermediate value
    bl = E * G - G * E  # Compute intermediate value
    cl = E * F - F * E  # Compute intermediate value
    
    K1 = (-bl - np.sqrt(np.abs(bl ** 2 - 4 * al * cl))) / (2 * al)  # Compute principal curvature K1
    K2 = (-bl + np.sqrt(np.abs(bl ** 2 - 4 * al * cl))) / (2 * al)  # Compute principal curvature K2
    
    K1[np.abs(K1) <= kt] = 0  # Apply thresholding to K1
    K2[np.abs(K2) <= kt] = 0  # Apply thresholding to K2
    
    KG = K1 * K2  # Compute Gaussian curvature
    KM = 0.5 * (K1 + K2)  # Compute mean curvature
    
    SMAP = np.empty_like(KG) * np.nan  # Initialize surface map
    SDist = [[None] * 9, [None] * 9]  # Initialize distribution of curvature features
    
    # Classify curvature features based on thresholds and store distribution information
    
    CMAP = {  # Create a dictionary to store curvature values
        'KG': KG,  # Gaussian curvature
        'KM': KM,  # Mean curvature
        'K1': K1,  # Principal curvature K1
        'K2': K2   # Principal curvature K2
    }

    # Return principal curvatures and surface map
    return K1, K2, KM, KG
