"""On-disk encoding; candidate selection always uses unscaled source values."""
import numpy as np

STORAGE_VERSION = 'int8-count-div4-v1'


def int8_quality(value):
    """Preserve the missing-quality sentinel -1 and saturate large qualities."""
    return max(-1, min(127, int(value)))


def encode_count(value):
    """Floor-divide raw counts by four, then saturate to signed-byte range."""
    return min(127, int(value) // 4)


def manifest_dtype(manifest):
    """Read versioned byte encodings and historical int32 v4 shards."""
    dtype = np.dtype(manifest.get('dtype', 'int32'))
    versions = {'int8': STORAGE_VERSION, 'uint8': 'uint8-saturated-v1'}
    if dtype.name in versions:
        if manifest.get('tensor_storage_version') != versions[dtype.name]:
            raise ValueError(f'{dtype} tensors require an explicit storage encoding version')
    elif dtype != np.dtype('int32'):
        raise ValueError('Unsupported v4 tensor dtype: '+str(dtype))
    return dtype
