using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Gnss_Spp
{
     class CommonData
    {
        public const double LIGHT=299792458.0;
        public const double GM = 3.986004418E14;
        public const double pi = 3.1415926535897932;
        public const double OMGe = 7.2921151467E-5;
        public const double RE_WGS84 = 6378137.0;//地球长半轴
        public const double FE_WGS84 = (1.0 / 298.257223563);//地球曲率(a - b) / a
        public const double minElevAngle = 10 * pi / 180.0;//对应的是10度
        public const double FREQ1_CMP = 1.561098E9;         /* BDS B1I     frequency (Hz) */
        public const double FREQ1 = 1.57542E9;          /* L1/E1/B1C  frequency (Hz) */
    }
}
