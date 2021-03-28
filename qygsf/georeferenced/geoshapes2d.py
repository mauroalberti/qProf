from pygsf.geometries.shapes.space2d import Point2D, Segment2D
from pygsf.geometries.shapes.space3d import Line3D
from pygsf.georeferenced.crs import Crs, check_crs
from pygsf.georeferenced.geoshapes3d import GeoLines3D
from pygsf.utils.types import check_type


class GeoPoints2D:
    """
    Collection of points.
    """

    def __init__(self,
                 epsg_code: numbers.Integral,
                 x_array: array,
                 y_array: array
                 ):
        """
        Construct a point list from a set of array values and an EPSG code.

        :param epsg_code: the EPSG code of the points
        :param x_array: the array storing the x values
        :param y_array: the array storing the y values
        """

        check_type(
            var=epsg_code,
            name="EPSG code",
            expected_types=numbers.Integral
        )

        check_type(
            var=x_array,
            name="X array",
            expected_types=array
        )

        check_type(
            var=y_array,
            name="Y array",
            expected_types=array
        )

        array_length = len(x_array)

        if len(y_array) != array_length:
            raise Exception(f"Y array has length {len(y_array)} while X array has length {len(x_array)}")

        self._epsg_code = epsg_code
        self._x_array = x_array
        self._y_array = y_array

    def num_pts(self
                ) -> int:
        """
        Numbers of points.
        """

        return len(self._x_array)

    @classmethod
    def fromPoints(cls,
                   points: List[Point2D],
                   epsg_code: numbers.Integral
                   ):
        """

        :param points: list of points
        :param epsg_code: optional EPSG code
        """

        for ndx, point in enumerate(points):

            check_type(point, "Input point {}".format(ndx), Point2D)

        return GeoPoints2D(
            epsg_code=epsg_code,
            x_array=array('d', [p.x for p in points]),
            y_array=array('d', [p.y for p in points])
        )

    @property
    def xs(self
           ) -> array:
        """
        Returns a copy of the points x values.

        :return: points x values
        """

        return copy(self._x_array)

    @property
    def ys(self
           ) -> array:
        """
        Returns a copy of the points y values.

        :return: points y values
        """

        return copy(self._y_array)

    def pt(self, pt_ndx: numbers.Integral) -> Point2D:
        """
        Extract the point at index pt_ndx.

        :param pt_ndx: point index.
        :return: the extracted Point instance.

        Examples:
        """

        return Point2D(
            x=self._x_array[pt_ndx],
            y=self._y_array[pt_ndx]
        )

    def values_at(self,
        ndx: numbers.Integral
                  ) -> Tuple[float, float]:
        """
        Return the values at given index.

        :param ndx: the index of the point values to extract
        :type ndx: numbers.Integral
        :return: the x and y values
        """

        return (
            self._x_array[ndx],
            self._y_array[ndx]
        )

    def pts(self):

        return [Point2D(*self.values_at(ndx)) for ndx in range(self.num_pts())]

    def __repr__(self) -> str:
        """
        Represents a GeoPoints2D instance as a shortened text.

        :return: a textual shortened representation of a Points instance.
        """

        num_points = self.num_pts()

        if num_points == 0:
            txt = "Empty GeoPoints2D"
        else:
            x1, y1 = self.values_at(0)
            if num_points == 1:
                txt = "GeoPoints2D with unique point: {.4f}.{.4f}".format(x1, y1)
            else:
                x2, y2 = self.values_at(self.num_pts()-1)
                txt = "GeoPoints2D with {} points: ({:.4f}, {:.4f}) ... ({:.4f}, {:.4f})".format(
                    num_points, x1, y1, x2, y2)

        return txt

    def __iter__(self):
        """
        Return each point.
        """

        return (self.pt(ndx) for ndx in range(self.num_pts()))

    @property
    def epsg_code(self) -> numbers.Integral:
        """
        The points EPSG code.

        :return: the points EPSG code
        :rtype: numbers.Integral
        """

        return self._epsg_code

    @property
    def crs(self) -> Crs:
        """
        The points CRS.

        :return: the points CRS
        :rtype: Crs
        """

        return Crs(self.epsg_code)

    def asXyArray(self):
        """
        Convert to a Numpy x-y array
        """

        return np.vstack(
            (
                self.xs,
                self.ys
            )
        ).transpose()

    def add_pt(self, pt) -> None:
        """
        In-place transformation of the original Line instance
        by adding a new point at the end.

        :param pt: the point to add
        :return: nothing
        """

        self._x_array.append(pt.x)
        self._y_array.append(pt.y)

    def add_pts(self,
                pts: List[Point2D]):
        """
        In-place transformation of the original Points instance
        by adding a new set of points at the end.

        :param pts: list of Points.
        """

        check_type(pts, "GeoPoints2D", GeoPoints2D)
        check_crs(self, pts)

        self._x_array.extend(pts.xs)
        self._y_array.extend(pts.ys)

    def x_min(self) -> Optional[numbers.Real]:
        """
        Optional minimum of x values.

        :return: the optional minimum of x values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = GeoPoints2D([[0, 0], [1, 0], [0, 1]])
          >>> l.x_min()
          0.0
          >>> m = GeoPoints2D([])
          >>> m.x_min()
          None
        """

        return np.nanmin(self._x_array) if self.num_pts() > 0 else None

    def x_max(self) -> Optional[numbers.Real]:
        """
        Optional maximum x value.
        """

        return np.nanmax(self._x_array) if self.num_pts() > 0 else None

    def x_mean(self) -> Optional[numbers.Real]:
        """
        Optional mean x value.
        """

        return np.nanmean(self._x_array) if self.num_pts() > 0 else None

    def y_min(self) -> Optional[numbers.Real]:
        """
        Optional minimum y value.
        """

        return np.nanmin(self._y_array) if self.num_pts() > 0 else None

    def y_max(self) -> Optional[numbers.Real]:
        """
        Optional maximum y value.
        """

        return np.nanmax(self._y_array) if self.num_pts() > 0 else None

    def y_mean(self) -> Optional[numbers.Real]:
        """
        Optional mean y value.
        """

        return np.nanmean(self._y_array) if self.num_pts() > 0 else None

    def nanmean_point(self) -> Point2D:
        """
        Returns the nan- excluded mean point of the collection.
        It is the mean point for a collection of point in a x-y frame (i.e., not lat-lon).

        :return: the nan- excluded mean point of the collection.
        :rtype: Point2D
        """

        return Point2D(
            x=np.nanmean(self._x_array),
            y=np.nanmean(self._y_array)
        )

    def segment(self,
        ndx: int
    ) -> Optional[Segment2D]:
        """
        Returns the optional segment starting at index ndx.

        :param ndx: the segment index.
        :return: the optional segment
        """

        if ndx < 0 or ndx >= self.num_pts() - 1:
            return None

        return Segment2D(
            start_pt=self.pt(ndx),
            end_pt=self.pt(ndx + 1)
        )

    def reversed(self) -> 'GeoPoints2D':
        """
        Return a Points instance with reversed point list.

        :return: a new Points instance.
        """

        xs = self._x_array.reverse()
        ys = self._y_array.reverse()

        return GeoPoints2D(
            epsg_code=self.epsg_code,
            x_array=xs,
            y_array=ys
        )


class GeoMPolygon2D:
    """
    A shapely (multi)polygon with EPSG code.

    """

    def __init__(self,
                 shapely_geom: Union[Polygon, MultiPolygon],
                 epsg_code: numbers.Integral
                 ):
        """
        :param shapely_geom: the (multi)polygon
        :type shapely_geom: Union[Polygon, MultiPolygon]
        :param epsg_code: the EPSG code of the two geometries
        :type epsg_code: numbers.Integral
        """

        check_type(
            shapely_geom,
            "Polygon",
            (Polygon, MultiPolygon)
        )

        self._geom = shapely_geom
        self._epsg_code = epsg_code

    @property
    def geom(self):
        return self._geom

    @property
    def epsg_code(self):
        return self._epsg_code

    def intersect_line(self,
                       line: LineString,
                       ) -> GeoLines3D:
        """
        Determine the intersections between a mpolygon and a line.

        :param line: the line
        :type line: shapely.geometry.LineString
        :return: the intersecting lines
        :rtype: GeoLines3D
        """

        lines = GeoLines3D()

        intersections = line.intersection(self.geom)

        if intersections:

            if intersections.geom_type == "LineString":

                inters_ln = line2d_from_shapely(
                    shapely_geom=intersections,
                    epsg_code=self.epsg_code
                )

                lines.append(inters_ln)

            elif intersections.geom_type == "MultiLineString":

                for intersection_line in intersections:

                    inters_ln = line2d_from_shapely(
                        shapely_geom=intersection_line,
                        epsg_code=self.epsg_code
                    )

                    lines.append(inters_ln)

            else:

                pass

        return lines


def line2d_from_shapely(
        shapely_geom: LineString,
        epsg_code: numbers.Integral
) -> GeoPoints2D:
    # Side effects: none
    """
    Create a Line instance from a shapely Linestring instance.

    :param shapely_geom: the shapely input LineString instance
    :type shapely_geom: shapely.geometry.linestring.LineString
    :param epsg_code: the EPSG code of the LineString instance
    :type epsg_code: numbers.Integral
    :return: the converted Line instance
    :rtype: Line3D
    """

    x_array, y_array = shapely_geom.xy

    return Line3D.fromArrays(
        x_array,
        y_array
    )


def line2d_to_shapely(
        src_line: GeoPoints2D
) -> LineString:
    """
    Create a shapely.LineString instance from a Line one.

    :param src_line: the source line to convert to the shapely format
    :type src_line: Line3D
    :return: the shapely LineString instance and the EPSG code
    :rtype: Tuple[LineString, numbers.Integral]
    """

    return LineString(src_line.xy_zipped()), src_line.epsg_code