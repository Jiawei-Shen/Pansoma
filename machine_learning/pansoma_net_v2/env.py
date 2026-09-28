"""Runtime environment fixes for torch.compile on the cluster's GPU nodes."""
import os
import subprocess
import tempfile


def triton_libcuda():
    """Point Triton at a directory with libcuda.so. Triton builds a small helper with `-lcuda`, which needs the
    unversioned libcuda.so; GPU nodes that ship only the driver's libcuda.so.1 (tequila) fail with "cannot find
    -lcuda" whenever that helper is not already in Triton's cache. Returns the directory used, or None."""
    if os.environ.get("TRITON_LIBCUDA_PATH"):
        return os.environ["TRITON_LIBCUDA_PATH"]
    try:
        listing = subprocess.run(["/sbin/ldconfig", "-p"], capture_output=True, text=True, check=False).stdout
    except OSError:
        return None
    found = [line.split("=>")[-1].strip() for line in listing.splitlines() if "libcuda.so.1 " in line or
             line.rstrip().endswith("libcuda.so.1")]
    if not found or os.path.exists(os.path.join(os.path.dirname(found[0]), "libcuda.so")):
        return None
    d = os.path.join(tempfile.gettempdir(), f"pansoma_libcuda_{os.getuid()}")
    os.makedirs(d, exist_ok=True)
    link = os.path.join(d, "libcuda.so")
    if os.path.islink(link) and os.readlink(link) != found[0]:
        os.unlink(link)
    try:
        os.symlink(found[0], link)
    except FileExistsError:  # another rank made it
        pass
    os.environ["TRITON_LIBCUDA_PATH"] = d
    return d
