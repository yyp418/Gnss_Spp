using System;
using System.CodeDom.Compiler;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using MathNet.Numerics;
using MathNet.Numerics.Integration;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using static System.Math;


namespace Gnss_Spp
{
     class Calute
    {
        public void calposition(DataCenter data)
        {
            //


        }
        //对每个周期的数据进行计算
        public DenseMatrix cal_epoch(int eph_num, DataCenter data)
        {
            var A = new DenseMatrix(data.Obs_file[eph_num].num, 4);//设置

            for (int i = 0; i < data.Obs_file[eph_num].num; i++)
            {

                Obs temp = data.Obs_file[eph_num].obs[i];//对应的轨道数据
                Orb_par temp_orb = SearchOrb(temp, data.Obv_file);//对应的星历文件
                //对tk的处理
                gpst_t temp_t = ToGpsTime(temp.time);//这里不对！！！
                //gps toolbox
                gpst_t temp_t_orb = ToGpsTime(temp_orb.ti);
                double ttt = temp.p / CommonData.LIGHT;
                double dt0 = temp_t.sec - temp.p / CommonData.LIGHT - temp_t_orb.sec;
                if (dt0 > 302400) dt0 -= 604800;
                if (dt0 < -302400) dt0 += 604800;
                double dtTmp = dt0;
                for (int j = 0; j < 2; j++) dt0 = dtTmp - (temp_orb.f0 + temp_orb.f1 * dt0 + temp_orb.f2 * dt0 * dt0);
                double dt = temp_orb.f0 + temp_orb.f1 * dt0 + temp_orb.f2 * dt0 * dt0;
                double tk = temp_t.sec - temp.p / CommonData.LIGHT - temp_orb.toe - dt;//
                if (tk > 302400) tk -= 604800;
                if (tk < -302400) tk += 604800;
                tk = tk - 14;//转换成北斗时


                double n = temp_orb.detn + Sqrt(CommonData.GM / Pow(temp_orb.sqrt_as, 6));
                double mk = temp_orb.M0 + n * tk;
                if (mk > 2 * CommonData.pi) mk -= 2 * CommonData.pi;
                if (mk < 0) mk += 2 * CommonData.pi;

                double E = 0, ETmp = 10;
                while (Abs(ETmp - E) > 1e-13)
                {
                    ETmp = E;
                    E -= (E - temp_orb.es * Sin(E) - mk) / (1 - temp_orb.es * Cos(E));
                }
                double vk = Acos((Cos(E) - temp_orb.es) / (1 - temp_orb.es * Cos(E)));
                double omegak = vk + temp_orb.w;
                double duk = temp_orb.cus * Sin(2 * omegak) + temp_orb.cuc * Cos(2 * omegak);
                double drk = temp_orb.crs * Sin(2 * omegak) + temp_orb.crc * Cos(2 * omegak);
                double dik = temp_orb.cis * Sin(2 * omegak) + temp_orb.cic * Cos(2 * omegak);

                double uk = omegak + duk;
                double rk = Pow(temp_orb.sqrt_as, 2) * (1 - temp_orb.es * Cos(E)) + drk;
                double ik = temp_orb.i0 + temp_orb.deti * tk + dik;
                double xk = rk * Cos(uk);
                double yk = rk * Sin(uk);

                double omuk = temp_orb.omega0 + (temp_orb.detomega - CommonData.OMGe) * tk - CommonData.OMGe * temp_orb.toe;


                A[i, 0] = xk * Cos(omuk) - yk * Cos(ik) * Sin(omuk);
                A[i, 1] = xk * Sin(omuk) + yk * Cos(ik) * Cos(omuk);
                A[i, 2] = yk * Sin(ik);


                if (temp_orb.PRN == "C01" ||
                    temp_orb.PRN == "C02" ||
                    temp_orb.PRN == "C03" ||
                    temp_orb.PRN == "C04" ||
                    temp_orb.PRN == "C05" ||
                    temp_orb.PRN == "C59" ||
                    temp_orb.PRN == "C60")
                {
                    double f = -5 * CommonData.pi / 180;
                    double p = CommonData.OMGe * tk;
                    var RX = new DenseMatrix(3, 3);
                    var RZ = new DenseMatrix(3, 3);
                    RX[0, 0] = 1; RX[0, 1] = 0; RX[0, 2] = 0;
                    RX[1, 0] = 0; RX[1, 1] = Cos(f); RX[1, 2] = Sin(f);
                    RX[2, 0] = 0; RX[2, 1] = -Sin(f); RX[2, 2] = Cos(f);

                    RZ[0, 0] = Cos(p); RZ[0, 1] = Sin(p); RZ[0, 2] = 0;
                    RZ[1, 0] = -Sin(p); RZ[1, 1] = Cos(p); RZ[1, 2] = 0;
                    RZ[2, 0] = 0; RZ[2, 1] = 0; RZ[2, 2] = 1;
                    var GEOP= new DenseMatrix(3, 1);
                    var P = new DenseMatrix(3, 1);
                    P[0, 0] = A[i, 0]; P[1, 0] = A[i, 1]; P[2, 0] = A[i, 2];
                    GEOP = RZ * RX * P;
                    A[i, 0] = GEOP[0, 0]; A[i, 1]= GEOP[1, 0]; A[i, 2]= GEOP[2, 0];
                }
                //求解钟差
                double F = -2 * Sqrt(CommonData.GM) / CommonData.LIGHT / CommonData.LIGHT;
                double dtr = F * temp_orb.es * temp_orb.sqrt_as * Sin(E);//这里的E是计算得到的偏近点角
                tk= temp.time.ToUnixTimeSeconds() - temp.p / CommonData.LIGHT - dt-temp_orb.ti.ToUnixTimeSeconds();
                //tk = timediff(ts, eph.toc);//这里的ts是观测时间减去距离和插值以及toc
                double dts = temp_orb.f0 + temp_orb.f1 * tk + temp_orb.f2 * tk * tk;
                A[i, 3] = dts * CommonData.LIGHT + dts;//未减去群波延时校正值
            }

            //按最低仰角排除卫星，即过滤掉低仰角卫星
            int error_count = 0;
            for (int i=0;i < data.Obs_file[eph_num].num; i++)
            {

                double[] rr = new double[3];
                double[] rs = new double[3];
                double[] e = new double[3];
                double[] azel=new double [2];//方位角和仰角
                rs[0] = A[i, 0]; rs[1] = A[i, 1]; rs[2] = A[i, 2];
                rr[0] = data.Appro_position[0]; rr[1] = data.Appro_position[1]; rr[2] = data.Appro_position[2];
                double sa;
                double dis = Geodist(rs, rr, out e,out sa);
                Sat2zel( rr, e, out  azel);
                double alph = azel[1];//仰角
                //int error_count = 0;
                //
                if(Abs(alph)<CommonData.minElevAngle || (rs[0] == 0 && rs[1] == 0 && rs[2] == 0))
                {
                    Obs temp = data.Obs_file[eph_num].obs[i];
                    temp.healtflag = false;
                    //data.Obs_file[eph_num].obs[i] = temp;//这里涉及到值类型和引用类型
                    error_count++;
                }
            }
            //data.Obs_file[eph_num].num_true = new int();
            //data.Obs_file[eph_num].num_true = data.Obs_file[eph_num].num - error_count;
            var SatPosClkValid = new DenseMatrix(data.Obs_file[eph_num].num-error_count, 4);

            // data.Obs_file[eph_num].idx = new List<int>();
            //obss[epoIndex].validIndex = new int[obss[epoIndex].n];
            //data.Obs_file[eph_num].idx = new List<int>();
            List<int> idx = new List<int>();
            for (int i=0,j=0;i< data.Obs_file[eph_num].num;i++)
            {
                Obs temp = data.Obs_file[eph_num].obs[i];
                if(temp.healtflag)
                {
                    SatPosClkValid[j, 0] = A[i, 0];
                    SatPosClkValid[j, 1] = A[i, 1];
                    SatPosClkValid[j, 2] = A[i, 2];
                    SatPosClkValid[j, 3] = A[i, 3];
                    j++;
                    idx.Add(i);
                }
            }
            data.Obs_file[eph_num].idx=idx;

            return SatPosClkValid;
        }

        //构建系数阵和常数阵
        public void Fromcoefficient(int eph_num, DenseMatrix SatPosClkValid,DataCenter data,out DenseMatrix B,out DenseMatrix L,out double[] varIon, out double []varTrp)
        {

            var temp = data.Obs_file[eph_num];
            int nValid = SatPosClkValid.RowCount;
            varIon = new double[nValid];
            varTrp = new double[nValid];
            B = new DenseMatrix(nValid, 4);
            L = new DenseMatrix(nValid, 1);

            double []rr=new double [4];
            rr[0] = data.Appro_position[0];
            rr[1] = data.Appro_position[1];
            rr[2] = data.Appro_position[3];
            rr[3] = 0;
            double[] rr_blh=new double[3];
            Ece2blh(rr, out rr_blh);
            for(int i=0;i< nValid;i++)
            {
                int ind = temp.idx[i];//得到原始的序号
                //eph_t eph = seleph(obs.data[ind]);
                Orb_par temp_orb = SearchOrb(temp.obs[ind],  data.Obv_file);//对应的星历文件
                double dts=SatPosClkValid[i, 3]; 
                double[] rs = new double[3];
                rs[0] = SatPosClkValid[i, 0];
                rs[1] = SatPosClkValid[i, 1];
                rs[2] = SatPosClkValid[i, 2];

                double[] vec = new double[3];
                string prn = temp.obs[ind].PRN;
                double SagEffectCor;
                double dis = Geodist(rs, rr, out vec,out SagEffectCor);
                double ro = dis - SagEffectCor;//原始未修正版距离
                B[i, 0] = -vec[0] / ro; B[i, 1] = -vec[1] / ro; B[i, 2] = -vec[2] / ro; B[i, 3] = 1;

                //群波延时校正值TGD
                double P, TGD;
                TGD = temp_orb.tgd;
                P = temp.obs[ind].p - TGD * CommonData.LIGHT;
                //对流层改正
                double[] azel = new double[2];
                Sat2zel(rr, vec, out azel);
                double el = azel[1];//仰角
                double dtrp = tropmodel(rr_blh, azel, 0.7);//0.7是默认的湿度,rr_blh是接收机的blh坐标
                varTrp[i] = dtrp * 0.3;
                gpst_t time = ToGpsTime(temp.obs[ind].time);
                double[] ion = new double[8];
                ion = data.Obv_file.ion_cmp;
                double dion = Ionomodel(rr_blh, azel, ion, time);
                dion *= Pow((CommonData.FREQ1 /CommonData.FREQ1_CMP), 2);//电离层和频率的平方成反比？
                varIon[i] = dion * 0.5;
                double dtr = rr[3];//如果赋值为0那有什么意义？
                double l = P - (dis + dtr - dts + dion + dtrp);
                L[i, 0] = l;
            }
        }


        //完成时间的转换,折磨死人.修正，因为处理的是观测数据，已经是gps时间，无需闰秒
        public gpst_t  ToGpsTime(DateTimeOffset dateTime)
            {
                DateTimeOffset GpsEpoch =
                new DateTimeOffset(1980, 1, 6, 0, 0, 0, TimeSpan.Zero);

                // 先把日期转换为UTC时间
                //DateTimeOffset utcDateTime = dateTime.ToUniversalTime();
              // 计算从GPS起始时间到目标时间的总秒数，同时考虑闰秒
                double totalSeconds = (dateTime - GpsEpoch).TotalSeconds;//？是连续秒

                // 计算周数和周内秒数
                int week = (int)(totalSeconds / (7 * 24 * 60 * 60));
                double secondsOfWeek = totalSeconds % (7 * 24 * 60 * 60);

                gpst_t temp_gpst = new gpst_t();
                temp_gpst.week = week;
                temp_gpst.sec = secondsOfWeek;

                return temp_gpst;
            }
        //输入观测文件，寻找对应的星历文件,粗糙版！
        public Orb_par SearchOrb(Obs OBS, Obv_t ObvFile)
        {
            gpst_t obs_temp_gpst = ToGpsTime(OBS.time);
            Orb_par temp = ObvFile.Ord[0];
            double temp_ti = obs_temp_gpst.sec - temp.toe;
            for (int i=0;i<ObvFile.Ord.Count;i++)
            {
                if(OBS.PRN !=ObvFile.Ord[i].PRN)continue;
                double T= obs_temp_gpst.sec- ObvFile.Ord[i].toe;
                if (Abs(T) >Abs( temp_ti)) continue;
                temp_ti = T;
                temp = ObvFile.Ord[i];
            }
            return temp;
        }

        //计算地理距离
        public double Geodist(double[]rs,double []rr ,out double[]e,out double SagEffectCor)
        {
            double  distance;
            SagEffectCor = CommonData.OMGe * (rs[0] * rr[1] - rs[1] * rr[0]) / CommonData.LIGHT;
            //double r;
            e = new double[3];//看看这里的格式
            if (Sqrt(Pow(rs[0], 2) + Pow(rs[1], 2) + Pow(rs[2], 2)) < CommonData.RE_WGS84) return -1;
            for(int i=0;i<3;i++)
            {
                e[i] = rs[i] - rr[i];
            }
            distance = Sqrt(Pow(e[0], 2) + Pow(e[1], 2) + Pow(e[2], 2));
            for (int i = 0; i < 3; i++) e[i] = e[i] / distance;
            //SagEffectCor = CommonData.OMGe * (rs[0] * rr[1] - rs[1] * rr[0]) / CommonData.LIGHT;
            return distance+ SagEffectCor;
        }

        public bool Ece2blh(double[]rr, out double[]blh)
        {
            blh = new double[3];
            double e2 = CommonData.RE_WGS84 * (2 - CommonData.RE_WGS84);
            double p = Pow(rr[0], 2) + Pow(rr[1], 2);
            blh[1] = p*p> 1E-12 ? Atan(rr[1]/rr[0]):0.0;
            double phi = 0; double temp=0;
            double N,h;
            while (Abs(temp-phi) > 1E-4)
            {
                temp = phi;
                N = CommonData.RE_WGS84 / Sqrt(1-e2*Sin(phi));
                h = p / Cos(phi) - N;
                phi= Atan(rr[2] / p * (1 - e2 * (N/(N+h))));
            }
            blh[0] = phi;
            blh[2] = p / Cos(blh[0]) - CommonData.RE_WGS84 / Sqrt(1 - e2 * Sin(phi));
            return true;
        }

        //ecf转东北天坐标系
        public void Sat2zel(double[] rr, double[] v, out double[] ael)
        {

            ael = new double[2];
            bool sta = Ece2blh(rr, out double[] blh);
            double B = blh[0], L = blh[1], H = blh[2];
            var S = new DenseMatrix(3, 3);
            S[0, 0] = -Sin(L); S[0, 1] = Cos(L); S[0, 2] = 0.0;
            S[1, 0] = -Sin(B) * Cos(L); S[1, 1] = -Sin(B) * Sin(L); S[1, 2] = Cos(B);
            S[2, 0] = Cos(B) * Cos(L); S[2, 1] = Cos(B) * Sin(L); S[2, 2] = Sin(B);

            var enu = new DenseMatrix(3, 1);
            var mat_v= new DenseMatrix(3, 1);
            mat_v[0, 0] = v[0]; mat_v[1, 0] = v[1]; mat_v[2, 0] = v[2];
            enu = S * mat_v;

            double norm = Sqrt(Pow(enu[0,0],2)+ Pow(enu[1, 0], 2)+ Pow(enu[2, 0], 2));
            double alph = (Pow(enu[0, 0], 2) + Pow(enu[1, 0], 2)) < 1E-12 ? 0.0 : Atan((enu[0,0]/ norm)/ (enu[1,0] / norm));
            ael[0] = alph;
            ael[1] = Asin(enu[2,0] / norm);
        }

        //对流层改正模型,humi是相对湿度
        public double tropmodel(double []pos,double[]azel,double humi)
        {
            double temp0 = 15;
            double hgt,pres,temp,e,z,trph,trpw;
            if (pos[2] < -100 || pos[2] > 10000 || azel[1]<=0) return 0.0;
            hgt = pos[2] < 0.0 ? 0.0 : pos[2];
            pres = 1013.25 * Pow(1 - 2.2557E-5 * hgt, 5.2568);
            temp = temp0 - 6.5 * E - 3 * hgt + 273.16;
            e = 6.108 * Exp((17.15 * temp - 4684.0) / (temp-38.45))*humi;

            //z为天顶角
             z = CommonData.pi / 2-azel[2];
            trph = 0.0022768 * pres / (1.0 - 0.00266 * Cos(2.0 * pos[0]) - 0.00028 * hgt / 1E3) / Cos(z);
            trpw = 0.002277 * (1255.0 / temp + 0.05) * e / Cos(z);
            return trph + trpw;
        }
        
        //电离层改正模型，克罗布歇，
        public double Ionomodel(double[]pos,double[]azel,double[] ion,gpst_t t)
        {
            double[] ion_default = new double[8];
            ion_default = new double[]{ 0.1118E-07, -0.7451E-08, -0.5961E-07, 0.1192E-06,
        0.1167E+06, -0.2294E+06, -0.1311E+06, 0.1049E+07 };

            double tt, f, psi, phi, lam, amp, per, x;
            int week;

            double norm = 0;
            for (int i = 0; i < 8; i++) norm += Pow(ion[i], 2);

            if (pos[2] < -1E3 || azel[1] <= 0) return 0.0;// 如果没有电离层参数，默认参数
            if (norm <= 0.0) ion = ion_default;

            psi = 0.0137 / (azel[1] / PI + 0.11) - 0.022;//计算地心角
            phi = pos[0] / PI + psi * Cos(azel[0]);//计算穿刺点地理维度
            if (phi > 0.416) phi = 0.416;
            else if (phi < -0.416) phi = -0.416;
            lam = pos[1] / PI + psi * Sin(azel[0]) / Cos(phi * PI);//经度度
            phi += 0.064 * Cos((lam - 1.617) * PI);//计算穿刺点地磁纬度
            //计算穿刺点地方时
            tt = 43200.0 * lam + t.sec;
            tt -= Floor(tt / 86400.0) * 86400.0; /* 0<=tt<86400 */
            f = 1.0 + 16.0 * Pow(0.53 - azel[1] / PI, 3.0);//投影系数
            /* ionospheric delay */
            amp = ion[0] + phi * (ion[1] + phi * (ion[2] + phi * ion[3]));
            per = ion[4] + phi * (ion[5] + phi * (ion[6] + phi * ion[7]));
            amp = amp < 0.0 ? 0.0 : amp;
            per = per < 72000.0 ? 72000.0 : per;
            x = 2.0 * PI * (tt - 50400.0) / per;
            return CommonData.LIGHT * f * (Abs(x) < 1.57 ? 5E-9 + amp * (1.0 + x * x * (-0.5 + x * x / 24.0)) : 5E-9);
        }
    }
}
