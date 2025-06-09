using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Remoting.Metadata.W3cXsd2001;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Tasks;



namespace Gnss_Spp
{
     class DataCenter
    {
        public double []Appro_position;// APPROX POSITION XYZ
        public double t_interval;//时间间隔
        public DateTimeOffset t_start;//开始的时间
        public DateTimeOffset t_end;
        public List<Obs_t> Obs_file;//观测文件里的数据
        public Obv_t Obv_file;//导航文件里的数据

    }
}
//gps时，周+相对本周的秒数，起始时间是1980,1,6
struct gpst_t
{
    public int week;
    public  double sec;
}
//在观测数据里的结构
class Obs_t
{
 
    public int num;//记录有多少颗卫星
    public List<Obs> obs;
    public List<int> idx;
}
//一个历元间的一个观测数据
struct Obs
{
    public bool healtflag;//卫星的健康状态
    public int idx;//序号
    public string PRN;//卫星的标识符
    public DateTimeOffset time;
    public double p;//这里简化一下，默认使用的是单频

    //实现构造函数
    public Obs(bool h,int idx1,string prn, DateTimeOffset t,double p1)
    {
        healtflag = h;
        idx = idx1;
        PRN = prn;
        time = t;
        p = p1;
    }
}

struct Orb_par
{
    public string PRN;
    public DateTimeOffset ti;
    public double toe;//周内秒数
    public double cuc, cus, crc, crs, cic,cis;//轨道改正参数
    public double omega0, i0, M0;
    public double detn,deti,detomega;
    public double w, sqrt_as, es;
    public double f0, f1, f2;//时钟修正参数
    public double tgd,sva;//群延迟参数以及定位精度（用户测距精度）
}

struct Obv_t
{

    public List<Orb_par> Ord;
    public double []utc_cmp;  //北斗的时间参数存疑
    public double []ion_cmp;  /* BeiDou iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    
    public Obv_t(int utcSize=8, int ionSize=8)
    {
        Ord = new List<Orb_par>();
        utc_cmp = new double[utcSize];//似乎没用到
        ion_cmp = new double[ionSize];

    }



}