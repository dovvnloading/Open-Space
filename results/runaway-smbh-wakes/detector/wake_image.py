"""Read rectangular pixel windows. FITS requires the system CFITSIO library."""

import ctypes as ct
import ctypes.util
from pathlib import Path

import numpy as np


class ImagePlane:
    def __init__(self, path, hdu=0):
        self.path = Path(path)
        self.array = None
        self.handle = None
        if self.path.suffix == ".npy":
            if hdu:
                raise ValueError("NPY images have a single plane at index zero")
            self.array = np.load(self.path, mmap_mode="r", allow_pickle=False)
            if self.array.ndim != 2:
                raise ValueError("Expected a 2D NPY")
            self.height, self.width = self.array.shape
            return
        self.lib = ct.CDLL(ctypes.util.find_library("cfitsio") or "libcfitsio.so")
        self.handle = ct.c_void_p()
        self._call(self.lib.ffopen, ct.byref(self.handle), bytes(self.path), 0)
        try:
            hdutype = ct.c_int()
            self._call(self.lib.ffmahd, self.handle, hdu + 1, ct.byref(hdutype))
            ndim = ct.c_int()
            self._call(self.lib.ffgidm, self.handle, ct.byref(ndim))
            if ndim.value != 2:
                raise ValueError(f"HDU {hdu} does not contain a 2D image")
            dims = (ct.c_long * 2)()
            self._call(self.lib.ffgisz, self.handle, 2, dims)
            self.width, self.height = dims
        except Exception:
            self.close()
            raise

    def _call(self, func, *args):
        status = ct.c_int()
        func(*args, ct.byref(status))
        if status.value:
            message = ct.create_string_buffer(81)
            self.lib.ffgerr(status.value, message)
            raise OSError(f"{self.path}: {func.__name__}: {message.value.decode()}")

    def read(self, x, y, w, h):
        if min(x, y) < 0 or min(w, h) <= 0 or x + w > self.width or y + h > self.height:
            raise ValueError("Window falls outside image")
        if self.array is not None:
            return np.asarray(self.array[y:y+h, x:x+w], dtype=np.float32)
        pair = ct.c_long * 2
        dest = np.empty((h, w), np.float32)
        null = ct.c_float(float("nan"))
        anynull = ct.c_int()
        self._call(self.lib.ffgsv, self.handle, 42, pair(x+1, y+1), pair(x+w, y+h),
                   pair(1, 1), ct.byref(null), ct.c_void_p(dest.ctypes.data), ct.byref(anynull))
        return dest

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    def close(self):
        if self.handle:
            self._call(self.lib.ffclos, self.handle)
            self.handle = None
