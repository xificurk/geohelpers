# -*- coding: utf-8 -*-
"""
Misc geo helpers.

Functions:
    polygon_from_gpx    --- Create Polygon instance from gpx file.

Classes:
    Point               --- Point wrapper (namedtuple).
    Polygon             --- Polygon wrapper.
    Rectangle           --- Rectangle polygon wrapper.
    Mercator            --- Google mercator projection.

"""

__author__ = "Petr Morávek (xificurk@gmail.com)"
__copyright__ = "Copyright (C) 2010 Petr Morávek"
__license__ = "LGPL 3.0"

__version__ = "0.5.0"

from collections import Sequence, namedtuple
import logging
import math
import os.path
from xml.etree import cElementTree as ET

__all__ = ["polygon_from_gpx",
           "Point",
           "Polygon",
           "Rectangle",
           "Mercator"]


_log = logging.getLogger("geohelpers")


############################################################
### Helpers.                                             ###
############################################################

def polygon_from_gpx(file_, precision=0, projector=None):
    """
    Create Polygon instance from gpx file.

    Arguments:
        file_       --- Path to gpx file.

    Keyworded arguments:
        precision   --- Number of digits after decimal point in points coordinates.
        projector   --- Function for reprojecting point coordinates.

    """
    _log.debug("Loading polygon from file {0}.".format(file_))
    points = []
    gpx = ET.parse(file_)
    for element in gpx.findall("{http://www.topografix.com/GPX/1/1}rte/{http://www.topografix.com/GPX/1/1}rtept"):
        x = float(element.get("lon"))
        y = float(element.get("lat"))
        if projector is not None:
            x, y = projector(x, y)
        point = Point(round(x, precision), round(y, precision))
        if len(points) == 0 or point != points[-1]:
            points.append(point)
    return Polygon(points, os.path.basename(file_).replace(".gpx", ""))



############################################################
### Geometry.                                            ###
############################################################

def _orientation(a, b, c):
    return (c.y - a.y) * (b.x - a.x) - (b.y - a.y) * (c.x - a.x)

def _line_segments_intersect(a, b, c, d):
    """ Do non-degenerate line segments ab and cd share a common point? """
    ABC = _orientation(a, b, c)
    if ABC == 0:
        if a.x != b.x:
            if a.x <= c.x <= b.x or b.x <= c.x <= a.x:
                return True
        elif a.y <= c.y <= b.y or b.y <= c.y <= a.y:
            return True
    else:
        ABD = _orientation(a, b, d)
        if ABD == 0:
            if a.x != b.x:
                if a.x <= c.x <= b.x or b.x <= d.x <= a.x:
                    return True
            elif a.y <= c.y <= b.y or b.y <= d.y <= a.y:
                return True
        elif ABC * ABD < 0:
            CDA = _orientation(c, d, a)
            if CDA == 0:
                return True
            elif CDA * _orientation(c, d, b) <= 0:
                return True
    return False


Point = namedtuple("Point", "x y")
"""
Point wrapper.

"""


class Polygon(Sequence):
    """
    Polygon wrapper.

    Essentially a sequence of Point instances representing polygon vertices.

    Attributes:
        name        --- Name of the polygon.
        envelope    --- Rectangle instance containing the whole polygon.

    Methods:
        contains    --- Determine whether the polygon contains point.
        intersects  --- Determine whether the polygon intersects other polygon.

    """

    def __init__(self, points, name="polygon"):
        """
        Arguments:
            points      --- Iterable containing Point instances representing
                            polygon vertices.

        Keyworded arguments:
            name        --- Polygon name.

        """
        points = self._clean_points(points)
        if len(points) < 3:
            raise ValueError("Not enough points for polygon.")
        self.name = str(name)
        self._points = tuple(points)
        min_x = min((p.x for p in points))
        min_y = min((p.y for p in points))
        max_x = max((p.x for p in points))
        max_y = max((p.y for p in points))
        self.envelope = Rectangle(min_x, min_y, max_x, max_y, "{0}.envelope".format(self.name))

    def _clean_points(self, points):
        cleaned_points = []
        for point in points:
            if not isinstance(point, Point):
                raise TypeError("Expected Point instance.")
            # remove repeating vertices
            if len(cleaned_points) > 0 and point == cleaned_points[-1]:
                continue
            cleaned_points.append(point)
        if cleaned_points[0] == cleaned_points[-1]:
            del cleaned_points[-1]
        return cleaned_points

    def __iter__(self):
        return self._points.__iter__()

    def __getitem__(self, index):
        return self._points.__getitem__(index)

    def __len__(self):
        return self._points.__len__()

    def __str__(self):
        return self.name

    def __repr__(self):
        return "{}({}) in ({}, {}, {}, {})".format(self.__class__.__name__, self.name, self.envelope.min_x, self.envelope.min_y, self.envelope.max_x, self.envelope.max_y)

    def contains(self, point):
        """
        Determine whether the polygon contains point.

        Arguments:
            point       --- Point instance.

        """
        if not isinstance(point, Point):
            raise TypeError("Expected Point instance.")
        if point in self._points:
            return True
        if not self.envelope.contains(point):
            return False
        inside = False
        p1 = self._points[-1]
        for p2 in self._points:
            if (p1.y < point.y) ^ (p2.y < point.y):
                x_intersect = (point.y - p2.y) * (p1.x - p2.x) / (p1.y - p2.y) + p2.x
                if point.x < x_intersect:
                    inside = not inside
                elif point.x == x_intersect:
                    return True
            elif (point.y == p1.y == p2.y) and ((p1.x < point.x) ^ (p2.x < point.x)):
                return True
            p1 = p2
        return inside

    def intersects(self, polygon):
        """
        Determine whether the polygon intersects other polygon.

        Arguments:
            polygon     --- Polygon instance.

        """
        if not isinstance(polygon, Polygon):
            raise TypeError("Expected Polygon instance.")
        if not self.envelope.intersects(polygon.envelope):
            return False
        if self.contains(polygon._points[0]) or polygon.contains(self._points[0]):
            return True
        pp1 = polygon._points[-1]
        for pp2 in polygon._points:
            sp1 = self._points[-1]
            for sp2 in self._points:
                if _line_segments_intersect(sp1, sp2, pp1, pp2):
                    return True
                sp1 = sp2
            pp1 = pp2
        return False


class Rectangle(Polygon):
    """
    Rectangle polygon wrapper.

    Attributes:
        min_x       --- Lower x coordinate.
        min_y       --- Lower y coordinate.
        max_x       --- Upper x coordinate.
        max_y       --- Upper y coordinate.
        width       --- Width of the rectangle.
        height      --- Height of the rectangle.
        center      --- Coordinates of center.

    Methods:
        contains    --- Determine whether the rectangle contains point.
        intersects  --- Determine whether the rectangle intersects other polygon.

    """

    def __init__(self, x1, y1, x2, y2, name="rectangle"):
        """
        Arguments:
            x1          --- Left boundary.
            y1          --- Bottom boundary.
            x2          --- Right boundary.
            y2          --- Top boundary.

        Keyworded arguments:
            name        --- Rectangle name.

        """
        if x1 > x2:
            x1, x2 = x2, x1
        elif x1 == x2:
            raise ValueError("Invalid rectangle - left and right coordinates are the same.")
        if y1 > y2:
            y1, y2 = y2, y1
        elif y1 == y2:
            raise ValueError("Invalid rectangle - bottom and top coordinates are the same.")
        self._min_x = x1
        self._max_x = x2
        self._min_y = y1
        self._max_y = y2
        self._points = (Point(x1, y1), Point(x1, y2), Point(x2, y2), Point(x2, y1))
        self.name = str(name)
        self.envelope = self

    @property
    def min_x(self):
        return self._min_x

    @property
    def min_y(self):
        return self._min_y

    @property
    def max_x(self):
        return self._max_x

    @property
    def max_y(self):
        return self._max_y

    @property
    def center(self):
        return Point(0.5 * (self._min_x + self._max_x), 0.5 * (self._min_y + self._max_y))

    @property
    def width(self):
        """ Width of the rectangle. """
        return self._max_x - self._min_x

    @property
    def height(self):
        """ Height of the rectangle. """
        return self._max_y - self._min_y

    def contains(self, point):
        """
        Determine whether the rectangle contains point.

        Arguments:
            point       --- Point instance.

        """
        if not isinstance(point, Point):
            raise TypeError("Expected Point instance.")
        if self._min_x <= point.x <= self._max_x and self._min_y <= point.y <= self._max_y:
            return True
        else:
            return False

    def intersects(self, polygon):
        """
        Determine whether the rectangle intersects other polygon.

        Arguments:
            polygon     --- Polygon instance.

        """
        if isinstance(polygon, Rectangle):
            if self._min_x > polygon.max_x or self._max_x < polygon.min_x or self._min_y > polygon.max_y or self._max_y < polygon.min_y:
                return False
            else:
                return True
        elif isinstance(polygon, Polygon):
            return Polygon.intersects(self, polygon)
        else:
            raise TypeError("Expected Polygon instance.")



############################################################
### Projections.                                         ###
############################################################

class Mercator(object):
    """
    Google mercator projection.

    Note: Mercator x, y pixel coordinates start at upper left corner.

    Attributes:
        R           --- Earth radius.
        zoom        --- Zoom.
        scale       --- Scale in m/px.

    Methods:
        ll2m        --- Convert WGS84 lon, lat to Mercator x, y meters.
        ll2px       --- Convert WGS84 lon, lat to Mercator x, y pixels in current scale.
        m2ll        --- Convert Mercator x, y meters to WGS84 lon, lat.
        px2ll       --- Convert Mercator x, y pixels in current scale to WGS84 lon, lat.

    """

    R = 6378137
    _scale = 1.0
    _zoom = math.log(_scale / (2 * math.pi * R) * 256, 0.5)

    @property
    def zoom(self):
        return self._zoom

    @zoom.setter
    def zoom(self, value):
        self._zoom = float(value)
        self._scale = (2 * math.pi * self.R ) / (256.0 * 2**self._zoom)

    @property
    def scale(self):
        return self._scale

    @scale.setter
    def scale(self, value):
        self._scale = float(value)
        self._zoom = math.log(self._scale / (2 * math.pi * self.R) * 256, 0.5)

    def ll2m(self, lon, lat):
        """
        Convert WGS84 lon, lat to Mercator x, y meters.

        Arguments:
            lon         --- WGS84 longitude.
            lat         --- WGS84 latitude.

        """
        x = lon / 180.0
        y = math.log(math.tan((90.0 + lat) / 360.0 * math.pi )) / math.pi
        if not (-1 <= x <= 1):
            raise ValueError("Longitude out of range.")
        elif not (-1 <= y <= 1):
            raise ValueError("Latitude out of range.")
        x = x * math.pi * self.R
        y = y * math.pi * self.R
        return x, y

    def m2ll(self, x, y):
        """
        Convert Mercator x, y meters to WGS84 lon, lat.

        Arguments:
            x           --- X coordinate.
            y           --- Y coordinate.

        """
        x = x / (math.pi * self.R)
        y = y / (math.pi * self.R)
        lon = x * 180.0
        lat = math.atan(math.exp(y * math.pi)) * 360.0 / math.pi - 90.0
        if not (-180.0 <= lon <= 180.0):
            raise ValueError("X out of range.")
        elif not (-90.0 <= lat <= 90.0):
            raise ValueError("Y out of range.")
        return lon, lat

    def px2m(self, x, y):
        """
        Convert Mercator x, y pixels in current scale to Mercator x, y meters.

        Arguments:
            x           --- X coordinate.
            y           --- Y coordinate.

        """
        x *= self.scale
        y *= self.scale
        x = x - math.pi * self.R
        y = math.pi * self.R - y
        return x, y

    def m2px(self, x, y):
        """
        Convert Mercator x, y meters to Mercator x, y pixels in current scale.

        Arguments:
            x           --- X coordinate.
            y           --- Y coordinate.

        """
        x = math.pi * self.R + x
        y = math.pi * self.R - y
        return x / self._scale, y / self._scale

    def ll2px(self, lon, lat):
        """
        Convert WGS84 lon, lat to Mercator x, y pixels in current scale.

        Arguments:
            lon         --- WGS84 longitude.
            lat         --- WGS84 latitude.

        """
        return self.m2px(*self.ll2m(lon, lat))

    def px2ll(self, x, y):
        """
        Convert Mercator x, y pixels in current scale to WGS84 lon, lat.

        Arguments:
            x           --- X coordinate.
            y           --- Y coordinate.

        """
        return self.m2ll(*self.px2m(x, y))