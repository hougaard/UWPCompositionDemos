using Microsoft.Graphics.Canvas;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Windows.Foundation;
using Windows.UI.Composition;

namespace DeepEarth.Layers
{
    public class Map
    {
        //CompositionVirtualDrawingSurface drawingSurface;

        public ISpatialReference SpatialReference { get; set; }
        public double ActualWidth { get; set; }
        public Point LogicalOrigin { get; set; }
        public Size MapViewPixelSize { get; set; }
        public Size MapViewLogicalSize { get; set; }
        public Size TileSize { get; internal set; }

        public double RotationAngle = 0;
    }

    /// <summary>
    /// <para>
    /// Helper class to transform Logical, Pixel and Geographics coordinates
    /// </para>
    /// <list type="bullet">
    /// <item>Logical - essentially the object is 1 unit, halfway is 0.5</item>
    /// <item>Pixel - Screen pixels, can be relative from the top left of the screen 
    /// or absolute from the top left of the object</item>
    /// <item>Geometry - Longitude, Latitude</item>
    /// </list>
    /// <example>
    /// 
    /// <code title="Transform a Pixel coordinate to a Geographical location"> 
    /// //Within our Map class we have this as a public property:
    /// 
    ///    private CoordTransform _CoordHelper;  
    ///    public CoordTransform CoordHelper
    ///    {
    ///        get
    ///        {
    ///            if (_CoordHelper == null)
    ///            {
    ///                _CoordHelper = new CoordTransform(this);
    ///            }
    ///            return _CoordHelper;
    ///       }
    ///    }
    /// 
    /// //then your code, assuming the map is "_map"
    /// Point geo = _map.CoordHelper.PixelToGeo(new Point(300,400));
    /// </code>
    /// </example>
    /// </summary>
    public class CoordTransform
    {
        private readonly Map _Map;

        /// <summary>
        /// Helper class to transform Logical, Pixel and Geographics coordinates
        /// </summary>
        /// <param name="map">The instance of the current Map object to perform the conversions on</param>
        public CoordTransform(Map map)
        {
            _Map = map;
        }

        #region 2D Point Transformations

        /// <summary>
        /// Convert a Geographical Point (longitude, latitude) to a logic point on the current map.
        /// </summary>
        /// <param name="geographicPoint">Geographical Point (longitude, latitude)</param>
        /// <returns>The logical Point</returns>
        public Point GeoToLogical(Point geographicPoint)
        {
            Point logical = GeographicToLogical(_Map.SpatialReference, geographicPoint);
            return logical;
        }

        /// <summary>
        /// Convert a Geographical Point (longitude, latitude) to a pixel point on the current screen.
        /// </summary>
        /// <param name="geographicPoint">Geographical Point (longitude, latitude)</param>
        /// <returns>Pixel Point</returns>
        public Point GeoToPixel(Point geographicPoint)
        {
            Point logical = GeoToLogical(geographicPoint);
            Point screen = LogicalToPixel(logical);
            return screen;
        }

        /// <summary>
        /// Convert a logic point to a Geographical Point (longitude, latitude) on the current map.
        /// </summary>
        /// <param name="logicalPoint">The logical Point</param>
        /// <returns>Geographical Point (longitude, latitude)</returns>
        public Point LogicalToGeo(Point logicalPoint)
        {
            Point geographic = LogicalToGeographic(_Map.SpatialReference, logicalPoint);
            return geographic;
        }

        /// <summary>
        /// Convert a logic point to a Pixel Point on the current screen at a particular zoom level.
        /// </summary>
        /// <param name="logicalPoint">The logical Point</param>
        /// <returns>Pixel Point</returns>
        public Point LogicalToPixel(Point logicalPoint)
        {
            //Note:  Below is the calculation of avoid the rubberband effect caused by the MSI Animation.
            //The code below will keep the pins snapped to the correct location on MSI.Spring animations.
            //code used to be: return _Msi.LogicalToElementPoint(logicalPoint);

            double scaleFactor = _Map.ActualWidth;
            var offset = new Point(-_Map.LogicalOrigin.X / _Map.MapViewLogicalSize.Width, -_Map.LogicalOrigin.Y / _Map.MapViewLogicalSize.Width);
            double zoomFactor = Math.Log(1 / _Map.MapViewLogicalSize.Width, 2);

            var interimPoint = new Point(offset.X + Math.Pow(2, zoomFactor) * logicalPoint.X, offset.Y + Math.Pow(2, zoomFactor) * logicalPoint.Y);
            var pixel = new Point(interimPoint.X * scaleFactor, interimPoint.Y * scaleFactor);
            return pixel;
        }

        /// <summary>
        /// Convert a pixel point to a Geographical Point (longitude, latitude) on the current screen.
        /// </summary>
        /// <param name="screenPoint">Pixel Point</param>
        /// <returns>Geographical Point (longitude, latitude)</returns>
        public Point PixelToGeo(Point screenPoint)
        {
            Point logical = PixelToLogical(screenPoint);
            Point geographic = LogicalToGeo(logical);
            return geographic;
        }

        /// <summary>
        /// Convert a pixel point to a Logical Point on the current screen.
        /// </summary>
        /// <param name="pixel">Pixel Point</param>
        /// <returns>The logical Point</returns>
        public Point PixelToLogical(Point pixel)
        {
            Point offset = _Map.LogicalOrigin;
            double pixelFactorX = _Map.MapViewLogicalSize.Width / _Map.MapViewPixelSize.Width;
            double pixelFactorY = _Map.MapViewLogicalSize.Height / _Map.MapViewPixelSize.Height;

            Point logical = new Point((pixel.X * pixelFactorX) + offset.X, (pixel.Y * pixelFactorY) + offset.Y);
            return logical;
        }

        internal Point PixelToLogicalIncRotation(Point pixel)
        {
            //consider rotation
            return PixelToLogical(RotatePixelbyMapRotation(pixel));
        }

        /// <summary>
        /// Rotates a Pixel coordinate by the current map rotation
        /// </summary>
        /// <param name="point">Pixel Point</param>
        /// <returns>Pixel Point</returns>
        public Point RotatePixelbyMapRotation(Point point)
        {
            double a = (Math.PI / 180) * _Map.RotationAngle;
            double x = point.X * Math.Cos(a) + point.Y * Math.Sin(a);
            double y = -point.X * Math.Sin(a) + point.Y * Math.Cos(a);
            return new Point(x, y);
        }

        #endregion

        #region 3D Point Transformations

        /// <summary>
        /// Converts a pixel point to a Geographical Point (longitude, latitude) at the supplied Zoom Level.
        /// </summary>
        /// <param name="screenPoint">Pixel Point</param>
        /// <param name="zoomLevel">The zoom level of the map</param>
        /// <returns>Geographical Point (longitude, latitude)</returns>
        public Point PixelToGeoAtZoom(Point screenPoint, int zoomLevel)
        {
            const double projectionOffset = Constants.ProjectionOffset;
            double c = Constants.EarthCircumference / ((1 << zoomLevel) * _Map.TileSize.Width);
            double f = screenPoint.X * c - projectionOffset;
            double g = projectionOffset - screenPoint.Y * c;

            double latitude = RadiansToDegrees(Math.PI / 2 - 2 * Math.Atan(Math.Exp(-g / Constants.EarthRadius)));
            double longitude = RadiansToDegrees(f / Constants.EarthRadius);

            return Clip(new Point(longitude, latitude));
        }

        /// <summary>
        /// Convert a Geographical Point (longitude, latitude) to a pixel point at the supplied Zoom Level.
        /// </summary>
        /// <param name="geographicPoint">Geographical Point (longitude, latitude)</param>
        /// <param name="zoomLevel">The zoom level of the map</param>
        /// <returns>Pixel Point</returns>
        public Point GeoToPixelAtZoom(Point geographicPoint, int zoomLevel)
        {
            const double projectionOffset = Constants.ProjectionOffset;
            double e = Math.Sin(DegreesToRadians(geographicPoint.Y));
            double g = Constants.EarthRadius * DegreesToRadians(geographicPoint.X);
            double h = Constants.EarthRadius / 2 * Math.Log((1 + e) / (1 - e));
            double c = Constants.EarthCircumference / ((1 << zoomLevel) * _Map.TileSize.Width);

            double pixelX = (projectionOffset + g) / c;
            double pixelY = (projectionOffset - h) / c;

            return new Point(pixelX, pixelY);
        }

        #endregion

        #region LogicalView to ZoomLevel Conversions

        /// <summary>
        /// Gets the Logical Viewport for the supplied zoomlevel maintaining the centre point of the map
        /// </summary>
        /// <param name="zoomLevel">The zoom level of the map</param>
        /// <returns>The Logical Viewport</returns>
        public Size ZoomLevelToLogicalView(double zoomLevel)
        {
            //If the zoom is invalid, do nothing. 
            Map map = _Map;
            //TileSource bl = map.BaseLayer.Source;
            double targetWidth = map.MapViewPixelSize.Width / map.TileSize.Width / Math.Pow(2.0, zoomLevel);
            double targetHeight = targetWidth * (map.MapViewPixelSize.Height / map.MapViewPixelSize.Width);
            Size viewPort = new Size(targetWidth, targetHeight);
            return viewPort;
        }

        /// <summary>
        /// Gets the Zoomlevel required to fit the supplied Logical Viewport
        /// </summary>
        /// <param name="viewPort">Logical Viewport</param>
        /// <returns>The zoom level to fit</returns>
        public double LogicalViewToZoomLevel(Size viewPort)
        {
            Map map = _Map;
            double targetZoom = Math.Log(map.MapViewPixelSize.Width / map.TileSize.Width / viewPort.Width, 2);
            return targetZoom;
        }

        #endregion

        #region GetScaleAtZoomLevel

        /// <summary>
        /// Gets the scale of the map at the supplied zoom level at the equator in the unit provided
        /// </summary>
        /// <param name="zoomLevel">The zoom level of the map</param>
        /// <param name="unit">The unit of measure to use</param>
        /// <returns>The scale value</returns>
        public double GetScaleAtZoomLevel(double zoomLevel, Unit unit)
        {
            return GetScaleAtZoomLevel(0, zoomLevel, unit);
        }

        /// <summary>
        /// Gets the scale of the map at the supplied zoom level at the supplied Geographical Point in the unit provided
        /// </summary>
        /// <param name="zoomLevel">The zoom level of the map</param>
        /// <param name="geographicPoint">Geographical Point (longitude, latitude)</param>
        /// <param name="unit">The unit of measure to use</param>
        /// <returns>The scale value</returns>
        public double GetScaleAtZoomLevel(double zoomLevel, Point geographicPoint, Unit unit)
        {
            return GetScaleAtZoomLevel(geographicPoint.Y, zoomLevel, unit);
        }

        /// <summary>
        /// Gets the scale of the map at the supplied zoom level at the supplied Latitude in the unit provided
        /// </summary>
        /// <param name="latitude">The Line of Latitude to calculate at</param>
        /// <param name="zoomLevel">The zoom level of the map</param>
        /// <param name="unit">The unit of measure to use</param>
        /// <returns>The scale value</returns>
        public double GetScaleAtZoomLevel(double latitude, double zoomLevel, Unit unit)
        {
            return GetScaleAtZoomLevel(latitude, zoomLevel, 96, unit);
        }

        /// <summary>
        /// Gets the scale of the map at the supplied zoom level at the supplied Latitude 
        /// at the supplied Sceen DPI in the unit provided
        /// </summary>
        /// <param name="latitude">The Line of Latitude to calculate at</param>
        /// <param name="zoomLevel">The zoom level of the map</param>
        /// <param name="screenDpi">The Dots Per Inch of the current Screen</param>
        /// <param name="unit">The unit of measure to use</param>
        /// <returns>The scale value</returns>
        public double GetScaleAtZoomLevel(double latitude, double zoomLevel, double screenDpi, Unit unit)
        {
            double resolution = GetResolutionAtZoomLevel(latitude, zoomLevel, unit);
            double scale = resolution * screenDpi;

            double unitPerInch = ConvertUnit(1, Unit.Inch, unit);
            return scale / unitPerInch;
        }

        #endregion

        #region GetResolutionAtZoomLevel

        /// <summary>
        /// Gets the resolution of a pixel at the equator at the supplied zoom level in the unit provided
        /// </summary>
        /// <param name="zoomLevel">The zoom level of the map</param>
        /// <param name="unit">The unit of measure to use</param>
        /// <returns>The resolution value</returns>
        public double GetResolutionAtZoomLevel(double zoomLevel, Unit unit)
        {
            return GetResolutionAtZoomLevel(0, zoomLevel, unit);
        }

        /// <summary>
        /// Gets the resolution of a pixel at the supplied Geographical Point (longitude, latitude) 
        /// at the supplied zoom level in the unit provided
        /// </summary>
        /// <param name="zoomLevel">The zoom level of the map</param>
        /// <param name="geographicPoint">Geographical Point (longitude, latitude)</param>
        /// <param name="unit">The unit of measure to use</param>
        /// <returns>The resolution value</returns>
        public double GetResolutionAtZoomLevel(double zoomLevel, Point geographicPoint, Unit unit)
        {
            return GetResolutionAtZoomLevel(geographicPoint.Y, zoomLevel, unit);
        }

        /// <summary>
        /// Gets the resolution of a pixel at the supplied latitude at the supplied zoom level in the unit provided
        /// </summary>
        /// <param name="latitude">The Line of Latitude to calculate at</param>
        /// <param name="zoomLevel">The zoom level of the map</param>
        /// <param name="unit">The unit of measure to use</param>
        /// <returns>The resolution value</returns>
        public double GetResolutionAtZoomLevel(double latitude, double zoomLevel, Unit unit)
        {
            //double mapWidth = _Map.BaseLayer.Source.MapSizeAtTileLevel(Convert.ToInt32(zoomLevel)).Width;
            //double earthRadius = ConvertUnit(Constants.EarthRadius, Unit.Meter, unit);
            //latitude = Clip(latitude, MinLatitude, MaxLatitude);
            //return Math.Cos(latitude * Math.PI / 180) * 2 * Math.PI * earthRadius / mapWidth;
            throw new Exception("Not working yet");
        }

        #endregion

        #region Static Math Calculations

        /// <summary>
        /// Converts a unit of measure to another unit of measure
        /// </summary>
        /// <param name="sourceValue">The value to convert</param>
        /// <param name="sourceUnit">The source unit of measure</param>
        /// <param name="destinationUnit">The unit of measure to convert to</param>
        /// <returns>The converted unit value</returns>
        public static double ConvertUnit(double sourceValue, Unit sourceUnit, Unit destinationUnit)
        {
            if (sourceUnit != destinationUnit)
            {
                switch (sourceUnit)
                {
                    case Unit.Kilometer:
                        sourceValue *= 1000;
                        break;

                    case Unit.Inch:
                        sourceValue *= Constants.INCH_TO_METER;
                        break;

                    case Unit.Mile:
                        sourceValue *= Constants.MILE_TO_METER;
                        break;
                }

                switch (destinationUnit)
                {
                    case Unit.Kilometer:
                        sourceValue /= 1000;
                        break;

                    case Unit.Inch:
                        sourceValue *= Constants.METER_TO_INCH;
                        break;

                    case Unit.Mile:
                        sourceValue *= Constants.METER_TO_MILE;
                        break;
                }
            }

            return sourceValue;
        }

        /// <summary>
        /// Converts Radians to Degrees
        /// </summary>
        /// <param name="r">Value in Radians</param>
        /// <returns>Value in Degrees</returns>
        public static double RadiansToDegrees(double r)
        {
            return (180 / Math.PI) * r;
        }

        /// <summary>
        /// Converts Degrees to Radians
        /// </summary>
        /// <param name="d">Value in Degrees</param>
        /// <returns>Value in Radians</returns>
        public static double DegreesToRadians(double d)
        {
            return (Math.PI / 180) * d;
        }

        /// <summary>
        /// Reduces the number to fit within the supplied bounds
        /// </summary>
        /// <param name="n">The number to clip</param>
        /// <param name="minValue">The minimal allowed value</param>
        /// <param name="maxValue">The maximum allowed value</param>
        /// <returns>The clipped value</returns>
        public static double Clip(double n, double minValue, double maxValue)
        {
            return Math.Min(Math.Max(n, minValue), maxValue);
        }

        /// <summary>
        /// Reduces the Geographical Point to fit within the bounds of the Earth.
        /// </summary>
        /// <param name="geographicPoint">Geographical Point (longitude, latitude)</param>
        /// <returns>Cliped Geographical Point (longitude, latitude)</returns>
        public static Point Clip(Point geographicPoint)
        {
            geographicPoint.X = Clip(geographicPoint.X, Constants.EarthMinLongitude, Constants.EarthMaxLongitude);
            geographicPoint.Y = Clip(geographicPoint.Y, Constants.EarthMinLatitude, Constants.EarthMaxLatitude);
            return geographicPoint;
        }

        #endregion

        #region General Coordinate Transformations
        /*
         * Future work:  Add transformation support from one datum & projection to another.
         * This will become very important with upcoming support for geometric shapes and other,
         * non-WGS WMS data providers.  Currently, all sources are in WGS, so default to Virtual
         * Earth's WGS Mercator Projection
         */

        /// <summary>
        /// Convert a logic point to a Geographical Point (longitude, latitude) for the specified spatial reference.
        /// </summary>
        /// <param name="srs">The spatial reference to use</param>
        /// <param name="logicalPoint">Logical Point</param>
        /// <returns>Geographical Point (longitude, latitude)</returns>
        public Point LogicalToGeographic(ISpatialReference srs, Point logicalPoint)
        {
            return srs.LogicalToGeographic(logicalPoint);
        }

        /// <summary>
        /// Convert a Geographical Point (longitude, latitude) to a logic point for the specified spatial reference.
        /// </summary>
        /// <param name="srs">The spatial reference to use</param>
        /// <param name="geographicPoint">Geographical Point (longitude, latitude)</param>
        /// <returns>Logical Point</returns>
        public Point GeographicToLogical(ISpatialReference srs, Point geographicPoint)
        {
            return srs.GeographicToLogical(geographicPoint);
        }

        #endregion General Coordinate Transformations

    }
    public static class Constants
    {
        //These are Earth specific, not any projection specific
        //See http://en.wikipedia.org/wiki/Geographic_coordinate_system.
        public const double EarthMinLatitude = -85.05112878D;
        public const double EarthMaxLatitude = 85.05112878D;
        public const double EarthMinLongitude = -180D;
        public const double EarthMaxLongitude = 180D;
        public const double EarthCircumference = EarthRadius * 2 * Math.PI;
        public const double HalfEarthCircumference = EarthCircumference / 2;
        public const double EarthRadius = 6378137;

        public const double ProjectionOffset = EarthCircumference * 0.5;

        public const double INCH_TO_METER = 0.0254D;
        public const double METER_TO_INCH = 39.3700787D;
        public const double METER_TO_MILE = 0.000621371192D;
        public const double MILE_TO_METER = 1609.344D;
    }
    public enum Unit
    {
        /// <summary>
        /// 2.54 centimeters
        /// </summary>
        Inch,

        /// <summary>
        /// statute mile, 1,609.344 meters
        /// </summary>
        Mile,

        /// <summary>
        /// distance travelled by light in free space in 1/299,792,458 of a second ;)
        /// </summary>
        Meter,

        /// <summary>
        /// 1000 meters
        /// </summary>
        Kilometer
    }
    public class SpatialReference : ISpatialReference
    {

        public const double HalfPi = 0.159154943091895;
        public const double RadiansToDegrees = 57.2957795130823;

        public string GeoGCS { get; set; }
        public string Datum { get; set; }
        public string DatumAuthority { get; set; }
        public double SpheroidRadius { get; set; }
        public double SpheroidFlattening { get; set; }
        public string SpheroidAuthority { get; set; }
        public double Primem { get; set; }
        public string PrimemAuthority { get; set; }
        public string ProjectionAuthority { get; set; }
        public double FalseEasting { get; set; }
        public double FalseNorthing { get; set; }
        public double CentralMeridian { get; set; }
        public double StandardParallel { get; set; }
        public double LatitudeOfOrigin { get; set; }
        public double AngularUnitOfMeasurement { get; set; }
        public string UnitAuthority { get; set; }
        public string Authority { get; set; }

        /// <summary>
        /// The real world coordinate scale at a given longitude
        /// </summary>
        public double ScaleX { get; set; }

        /// <summary>
        /// The real world coordinate scale at a given latitude
        /// </summary>
        public double ScaleY { get; set; }

        /// <summary>
        /// Logical X offset to centre of earth
        /// </summary>
        public double OffsetX { get; set; }

        /// <summary>
        /// Logical Y offset to centre of earth
        /// </summary>
        public double OffsetY { get; set; }

        /// <summary>
        /// Converts Map UI coordinates to real world (geographic) coordinates for a point
        /// </summary>
        /// <param name="logicalPoint">The logical UI point</param>
        /// <returns>Geographic Point</returns>
        public Point LogicalToGeographic(Point logicalPoint)
        {
            return new Point((logicalPoint.X - OffsetX) * RadiansToDegrees / ScaleX, Math.Atan(Math.Sinh((logicalPoint.Y - OffsetY) / ScaleY)) * RadiansToDegrees);
        }

        /// <summary>
        /// Converts real world (geographic) point coordinates to Map UI coordinates
        /// </summary>
        /// <param name="geographicPoint"></param>
        /// <returns>Logical Point</returns>
        public Point GeographicToLogical(Point geographicPoint)
        {
            double d = Math.Sin(geographicPoint.Y * AngularUnitOfMeasurement);
            return new Point((geographicPoint.X * AngularUnitOfMeasurement * ScaleX) + OffsetX, (0.5 * Math.Log((1.0 + d) / (1.0 - d)) * ScaleY) + OffsetY);
        }
    }
    public interface ISpatialReference
    {

        /// <summary>
        /// A coordinate system based on latitude and longitude. Some geographic coordinate systems are Lat/Lon, 
        /// and some are Lon/Lat. You can find out which this is by examining the axes. 
        /// You should also check the angular units, since not all geographic coordinate systems use degrees.
        /// </summary>
        string GeoGCS { get; set; }

        /// <summary>
        /// This indicates the horizontal datum, which corresponds to the procedure used to measure positions on the surface of the Earth. 
        /// </summary>
        string Datum { get; set; }

        /// <summary>
        /// This indicates the horizontal datum, which corresponds to the procedure used to measure positions on the surface of the Earth.
        /// </summary>
        string DatumAuthority { get; set; }

        /// <summary>
        /// This describes a spheroid, which is an approximation of the Earth's surface as a squashed sphere.
        /// </summary>
        double SpheroidRadius { get; set; }

        /// <summary>
        /// This describes a spheroid, which is an approximation of the Earth's surface as a squashed sphere.
        /// </summary>
        double SpheroidFlattening { get; set; }

        /// <summary>
        /// This describes a spheroid, which is an approximation of the Earth's surface as a squashed sphere.
        /// </summary>
        string SpheroidAuthority { get; set; }

        /// <summary>
        /// This defines the meridian used to take longitude measurements from. 
        /// The units of the longitude must be inferred from the context.
        /// </summary>
        double Primem { get; set; }

        /// <summary>
        /// This defines the meridian used to take longitude measurements from. 
        /// The units of the longitude must be inferred from the context.
        /// </summary>
        string PrimemAuthority { get; set; }

        /// <summary>
        /// This describes a projection from geographic coordinates to projected coordinates.
        /// </summary>
        string ProjectionAuthority { get; set; }

        /// <summary>
        /// The value added to all "x" values in the rectangular coordinate for a map projection. 
        /// This value frequently is assigned to eliminate negative numbers. Expressed in the unit of measure identified in Planar Coordinate Units
        /// </summary>
        double FalseEasting { get; set; }

        /// <summary>
        /// The value added to all "y" values in the rectangular coordinates for a map projection. 
        /// This value frequently is assigned to eliminate negative numbers. 
        /// Expressed in the unit of measure identified in Planar Coordinate Units
        /// </summary>
        double FalseNorthing { get; set; }

        /// <summary>
        /// The line of longitude at the center of a map projection generally used as the basis for constructing the projection
        /// </summary>
        double CentralMeridian { get; set; }

        /// <summary>
        /// The line of constant latitude at which the surface of the Earth and the plane or developable surface intersect
        /// </summary>
        double StandardParallel { get; set; }

        /// <summary>
        /// The latitude chosen as the origin of rectangular coordinate for a map projection
        /// </summary>
        double LatitudeOfOrigin { get; set; }

        /// <summary>
        /// The measurement units used to define the angles of a spheroid or ellipse associated with a specific datum.
        /// For DeepEarth, the datum is usually WGS (World Geodetic System) 1984 and the unit of measurement is a degree
        /// </summary>
        double AngularUnitOfMeasurement { get; set; }

        /// <summary>
        /// The authority body that defines the unit of measurement i.e. European Petroleum Survey Group (EPSG).  For DeepEarth
        /// the unit of measurement is usually degrees and the authority for the datum the map uses, WGS 1984 is EPSG:4326
        /// </summary>
        string UnitAuthority { get; set; }

        /// <summary>
        /// The authority body that defines the standards for the spatial reference parameters.  For DeepEarth,
        /// Spatial Reference is usually WGS 1984 and the authority is EPSG:4326
        /// </summary>
        string Authority { get; set; }

        /// <summary>
        /// Converts a logical Point (0->1) to a geographical coordinate (Longitude, Latitude)
        /// </summary>
        /// <param name="logicalPoint">The logical Point</param>
        /// <returns>The geographical coordinate (Longitude, Latitude)</returns>
        Point LogicalToGeographic(Point logicalPoint);

        /// <summary>
        /// Converts a geographical coordinate (Longitude, Latitude) to a logical Point (0->1)
        /// </summary>
        /// <param name="geographyPoint">The geographical coordinate (Longitude, Latitude)</param>
        /// <returns>The logical Point</returns>
        Point GeographicToLogical(Point geographyPoint);
    }
}