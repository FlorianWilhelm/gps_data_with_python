# -*- coding: utf-8 -*-

import math
import numpy as np

# RDP algorithm as long as the rdp package is not iterative.
# See https://github.com/fhirschmann/rdp/issues/5
def _DouglasPeucker(points, startIndex, lastIndex, epsilon):
    stk = []
    stk.append([startIndex, lastIndex])
    globalStartIndex = startIndex
    bitArray = np.ones(lastIndex-startIndex+1, dtype=bool)

    while len(stk) > 0:
        startIndex = stk[-1][0]
        lastIndex = stk[-1][1]
        stk.pop()

        dmax = 0.
        index = startIndex

        for i in range(index+1, lastIndex):
            if bitArray[i - globalStartIndex]:
                d = PointLineDistance(points[i], points[startIndex], points[lastIndex])
                if d > dmax:
                    index = i
                    dmax = d
        if dmax > epsilon:
            stk.append([startIndex, index])
            stk.append([index, lastIndex])
        else:
            for i in range(startIndex + 1, lastIndex):
                bitArray[i - globalStartIndex] = False
    return bitArray


def rdp(points, epsilon):
    """
    Ramer-Douglas-Peucker algorithm
    """
    bitArray = _DouglasPeucker(points, 0, len(points)-1, epsilon)
    resList = []
    for i in range(len(points)):
        if bitArray[i]:
            resList.append(points[i])
    return np.array(resList)


def PointLineDistance(point, start, end):
    if np.all(np.equal(start, end)) :
        return np.linalg.norm(point, start)
    n = abs((end[0] - start[0]) * (start[1] - point[1]) - (start[0] - point[0]) * (end[1] - start[1]))
    d = math.sqrt((end[0] - start[0]) * (end[0] - start[0]) + (end[1] - start[1]) * (end[1] - start[1]))
    return n/d


def haversine(coord1, coord2):
    """
    Haversine distance in meters for two (lat, lon) coordinates
    """
    lat1, lon1 = coord1
    lat2, lon2 = coord2
    radius = 6371000 # mean earth radius in meters (GRS 80-Ellipsoid)
    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c
    return d
