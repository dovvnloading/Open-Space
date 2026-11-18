"""Continuous nominal angular corridor for an offset-ellipsoid altitude family."""

import numpy as np

from geometry import A, B, intersect, surface


def _direction(vector):
    vector = np.asarray(vector, dtype=float)
    if vector.shape != (3,) or not np.isfinite(vector).all():
        raise ValueError("expected a finite three-vector")
    length = np.linalg.norm(vector)
    if not np.isfinite(length) or length == 0:
        raise ValueError("degenerate direction")
    return vector / length


def _angle(first, second):
    return float(np.arctan2(np.linalg.norm(np.cross(first, second)), first @ second))


def altitude_arc(development_station, receiver_station, longitude, latitude,
                 height_bounds=(40.0, 120.0)):
    """Return unit directions at the ends of the continuous receiver sky arc.

    The first intersections with nested shells move monotonically along the
    development ray. The normalized receiver sightline to the intervening straight
    segment traces the minor great-circle arc between the endpoint directions.
    Offsets model nested ellipsoids, not exact constant-geodetic-height surfaces.
    """
    low, high = height_bounds
    if (not np.isfinite([low, high, longitude, latitude]).all()
            or not 0 <= low <= high
            or not (-180 <= longitude <= 180 and -90 <= latitude <= 90)):
        raise ValueError("invalid altitude or geographic bounds")
    development_station = np.asarray(development_station, dtype=float)
    receiver_station = np.asarray(receiver_station, dtype=float)
    outer_axes = np.array([A + high, A + high, B + high])
    for station in (development_station, receiver_station):
        if station.shape != (3,) or not np.isfinite(station).all():
            raise ValueError("invalid station")
        if np.linalg.norm(station / outer_axes) <= 1:
            raise ValueError("station is inside requested altitude shells")
    mapped_point = surface(longitude, latitude, A + 14, B + 6)
    endpoints = [intersect(development_station, mapped_point, A + height, B + height)
                 for height in (low, high)]
    return tuple(_direction(point - receiver_station) for point in endpoints)


def angular_residual(arc, observed_direction):
    """Return angular distance in radians to a nominal minor arc or point."""
    start, end = map(_direction, arc)
    observed = _direction(observed_direction)
    normal = np.cross(start, end)
    cross_length = np.linalg.norm(normal)
    distance = min(_angle(observed, start), _angle(observed, end))
    if cross_length < 1e-14:
        if start @ end < 0:
            raise ValueError("antipodal endpoints do not define a unique minor arc")
        return distance
    normal /= cross_length
    projection = observed - (normal @ observed) * normal
    if np.linalg.norm(projection) < 1e-14:
        return distance
    projection = _direction(projection)
    if np.cross(start, projection) @ normal >= 0 and np.cross(projection, end) @ normal >= 0:
        return min(distance, _angle(observed, projection))
    return distance
