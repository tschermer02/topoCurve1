# Reversibility_Test.py
import numpy as np
import matplotlib.pyplot as plt
from topocurve.SpectralFiltering import SpectralFiltering  # adjust import if your path differs

# --- synthetic surface ---
def generate_surface(size=128, scale=10.0):
    x = np.linspace(-np.pi, np.pi, size)
    y = np.linspace(-np.pi, np.pi, size)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(X) * np.cos(Y) * scale
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    return Z, dx, dy

# --- frequency radius grid (unshifted FFT convention) ---
def freq_radius_grid(ny, nx, dy, dx):
    ky = np.fft.fftfreq(ny, d=dy)
    kx = np.fft.fftfreq(nx, d=dx)
    KX, KY = np.meshgrid(kx, ky)
    return np.sqrt(KX**2 + KY**2)

# --- low/high complementary pair using your class’ window/pad/crop (no class edits needed) ---
def reconstruct_via_lp_plus_hp(sf, Z_in, dx, dy, cutoff_max, alpha):
    sf.z_array = Z_in
    sf.dx = [dx, dy]

    padded = sf.padding(alpha)  # detrend -> mirror -> Tukey -> pad
    Ny, Nx = padded.shape

    FZ = np.fft.fft2(padded)
    km = freq_radius_grid(Ny, Nx, dy, dx)

    H_low  = (km <= float(cutoff_max)).astype(float)
    H_high = 1.0 - H_low

    ZL_full = np.fft.ifft2(FZ * H_low ).real
    ZH_full = np.fft.ifft2(FZ * H_high).real

    # crop both to original support
    ZL = ZL_full[(sf.pad_x_max + sf.dim_x): -(sf.pad_x_min + sf.dim_x),
                 (sf.pad_y_max + sf.dim_y): -(sf.pad_y_min + sf.dim_y)]
    ZH = ZH_full[(sf.pad_x_max + sf.dim_x): -(sf.pad_x_min + sf.dim_x),
                 (sf.pad_y_max + sf.dim_y): -(sf.pad_y_min + sf.dim_y)]

    # add plane to each, but only once in the sum
    Z_low  = ZL + sf.plane
    Z_high = ZH + sf.plane
    Z_sum  = Z_low + Z_high - sf.plane
    return Z_sum

def main():
    # original
    Z0, dx, dy = generate_surface(size=128, scale=10.0)

    # create instance without calling __init__(tiff_file)
    sf = SpectralFiltering.__new__(SpectralFiltering)

    # cutoff = 30% of Nyquist; tweak this ratio if you want
    kx_nyq = 0.5 / dx
    ky_nyq = 0.5 / dy
    k_nyq  = min(kx_nyq, ky_nyq)
    cutoff = 0.30 * k_nyq
    alpha  = 0.25  # Tukey window parameter

    # reconstruct (undo) via LP + complement HP
    Z_rec = reconstruct_via_lp_plus_hp(sf, Z0, dx, dy, cutoff, alpha)

    # differences
    diff  = Z_rec - Z0
    mae   = float(np.mean(np.abs(diff)))
    print(f"Reconstruction MAE = {mae:.6e}")

    # plot ONLY: original, reconstructed, |difference|
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)

    im0 = ax[0].imshow(Z0, cmap='viridis');      ax[0].set_title("Original")
    plt.colorbar(im0, ax=ax[0], fraction=0.046, pad=0.04)

    im1 = ax[1].imshow(Z_rec, cmap='viridis');   ax[1].set_title("Reconstructed (LP + HP)")
    plt.colorbar(im1, ax=ax[1], fraction=0.046, pad=0.04)

    im2 = ax[2].imshow(np.abs(diff), cmap='inferno'); ax[2].set_title("|Reconstructed − Original|")
    plt.colorbar(im2, ax=ax[2], fraction=0.046, pad=0.04)

    for a in ax:
        a.set_xticks([]); a.set_yticks([])
    plt.show()

if __name__ == "__main__":
    main()
