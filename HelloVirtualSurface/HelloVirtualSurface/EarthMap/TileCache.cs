using Microsoft.Graphics.Canvas;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HelloVirtualSurface
{
    public static class TileCache
    {
        public static Dictionary<string, CanvasBitmap> Tiles = new Dictionary<string, CanvasBitmap>();
        public static void AddImage(int ZoomLevel, int tilePositionX, int tilePositionY,CanvasBitmap bitmap)
        {
            try
            {
                Tiles.Add(GetTileKey(ZoomLevel, tilePositionX, tilePositionY), bitmap);
            }
            catch
            {
                // Ignore duplicate adds due to threading...
            }
        }

        public static string GetTileKey(int ZoomLevel, int tilePositionX, int tilePositionY)
        {
            return ZoomLevel.ToString() + "-" +
                   tilePositionX.ToString() + "-" +
                   tilePositionY.ToString();
        }
        private const string TilePath = @"http://map.glomaris.net/tileserver.aspx?Z={Z}&X={X}&Y={Y}";
        public static Uri GetTileUri(int ZoomLevel, int tilePositionX, int tilePositionY)
        {
            string url = string.Empty;
            url = TilePath;
            url = url.Replace("{Z}", ZoomLevel.ToString());
            url = url.Replace("{X}", tilePositionX.ToString());
            url = url.Replace("{Y}", tilePositionY.ToString());
            return new Uri(url);
        }


    }
}
