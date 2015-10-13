
import numpy as np
import math
import random

# See:
# http://www.nayuki.io/res/smallest-enclosing-circle/smallestenclosingcircle.py

# http://stackoverflow.com/q/27673463/1391441

# Returns the smallest circle that encloses all the given points. Runs in
# expected O(n) time, randomized.
# Input: A sequence of pairs of floats or ints, e.g. [(0,5), (3.1,-2.7)].
# Output: A triple of floats representing a circle.
# Note: If 0 points are given, None is returned. If 1 point is given, a circle
# of radius 0 is returned.
#


def _is_in_circle(c, p):
    _EPSILON = 1e-12
    return c is not None and math.hypot(p[0] - c[0], p[1] - c[1]) < c[2] +\
        _EPSILON


def _make_diameter(p0, p1):
    return ((p0[0] + p1[0]) / 2.0, (p0[1] + p1[1]) / 2.0, math.hypot(p0[0] -
            p1[0], p0[1] - p1[1]) / 2.0)


def _cross_product(x0, y0, x1, y1, x2, y2):
    '''
    Returns twice the signed area of the triangle defined by (x0, y0),
    (x1, y1), (x2, y2)
    '''
    return (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)


def _make_circumcircle(p0, p1, p2):
    '''
    Mathematical algorithm from Wikipedia: Circumscribed circle
    '''
    ax = p0[0]
    ay = p0[1]
    bx = p1[0]
    by = p1[1]
    cx = p2[0]
    cy = p2[1]
    d = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2.0
    if d == 0.0:
        return None
    x = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) *
         (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
    y = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) *
         (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
    return (x, y, math.hypot(x - ax, y - ay))


def _make_circle_two_points(points, p, q):
    '''
    Two boundary points known
    '''
    diameter = _make_diameter(p, q)
    if all(_is_in_circle(diameter, r) for r in points):
        return diameter

    left = None
    right = None
    for r in points:
        cross = _cross_product(p[0], p[1], q[0], q[1], r[0], r[1])
        c = _make_circumcircle(p, q, r)
        if c is None:
            continue
        elif cross > 0.0 and (left is None or _cross_product(p[0], p[1], q[0],
                              q[1], c[0], c[1]) > _cross_product(p[0], p[1],
                              q[0], q[1], left[0], left[1])):
            left = c
        elif cross < 0.0 and (right is None or _cross_product(p[0], p[1],
                              q[0], q[1], c[0], c[1]) <
                              _cross_product(p[0], p[1], q[0], q[1],
                                             right[0], right[1])):
            right = c
    return left if (right is None or
                    (left is not None and left[2] <= right[2])) else right


def _make_circle_one_point(points, p):
    '''
    One boundary point known
    '''
    c = (p[0], p[1], 0.0)
    for (i, q) in enumerate(points):
        if not _is_in_circle(c, q):
            if c[2] == 0.0:
                c = _make_diameter(p, q)
            else:
                c = _make_circle_two_points(points[0: i + 1], p, q)
    return c


def make_circle(points):
    # Convert to float and randomize order
    shuffled = [(float(p[0]), float(p[1])) for p in points]
    random.shuffle(shuffled)

    # Progressively add points to circle or recompute circle
    c = None
    for (i, p) in enumerate(shuffled):
        if c is None or not _is_in_circle(c, p):
            c = _make_circle_one_point(shuffled[0: i + 1], p)
    return c


def get_cent_rad(pts_neighbors):
    '''
    Assign a center and a radius for each overdensity/group identified.
    '''
    cent_rad = []
    for g in pts_neighbors:

        pts = np.array(zip(*g))
        mean_pts = np.mean(pts, 0)

        # Translate towards origin
        pts -= mean_pts
        result = make_circle(pts)

        # Move back to correct position.
        c_r = (result[0] + mean_pts[0], result[1] + mean_pts[1], result[2])

        cent_rad.append([c_r[0], c_r[1], c_r[2]])

    return cent_rad
