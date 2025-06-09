using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Windows.Forms;
using static System.Windows.Forms.VisualStyles.VisualStyleElement.TaskbarClock;

namespace Gnss_Spp
{
     class FileHelper
    {
        /*读取观测文件*/
        public void OpenObsfile(DataCenter data)
        {
            try
            {
                OpenFileDialog opf = new OpenFileDialog();
                opf.Filter = "文本文件（*.22o）|*.22o";
                opf.Title = "请输入要导入观测值的文件";
                if (opf.ShowDialog() == DialogResult.OK)
                {
                    var reader = new StreamReader(opf.FileName); 

                    /* APPROX POSITION XYZ测站的近似坐标*/
                    string line;
                    int linecount1 = 9;
                    while (linecount1>0)//c#无法将int直接转换成bool类型
                    {
                        line = reader.ReadLine();
                        linecount1--;
                    }
                    line = reader.ReadLine();
                    var lines = line?.Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                    data.Appro_position = new double[3];
                    data.Appro_position[0] = double.Parse(lines[0]);
                    data.Appro_position[1] = double.Parse(lines[1]);
                    data.Appro_position[2] = double.Parse(lines[2]);

                    /* 读取时间间隔*/
                    int linecount2 = 13;
                    while (linecount2 > 0)//c#无法将int直接转换成bool类型
                    {
                        line = reader.ReadLine();
                        linecount2--;
                    }
                    line = reader.ReadLine();
                    var lines1= line?.Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries); ;
                    data.t_interval = double.Parse(lines1[0]);

                    /* read 'TIME OF FIRST OBS',观测开始的时间*/
                    line = reader.ReadLine();
                    var lines2= line?.Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                    data.t_start = str2DateTimeOffset(lines2);

                    /* read 'TIME OF LAst OBS',观测结束的时间*/
                    line = reader.ReadLine();
                    var lines3 = line?.Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                    data.t_end = str2DateTimeOffset(lines3);

                    /* 读取每个周期的观测数据*/
                    while ((line = reader.ReadLine()) != null)
                    {
                        if (line.Trim().Equals("END OF HEADER"))
                        {
                            break; // 找到结束标记，跳出循环
                        }
                    }

                    data.Obs_file = new List<Obs_t>();

                    TimeSpan timeDiff = data.t_end - data.t_start;
                    // 将时间差转换为总秒数（double类型）
                    double totalSeconds = timeDiff.TotalSeconds;
                    int epoch = (int)(totalSeconds / data.t_interval)+1;//这里可能出问题
                    for(int i=0;i< epoch;i++)
                    {
                        line = reader.ReadLine();
                        var lines4 = line?.Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Skip(1).ToArray();

                        DateTimeOffset timeEpoch = str2DateTimeOffset(lines4);
                        int num_oneepoch = int.Parse(lines4[7]);//一个历元间的数据条数
                        List<Obs> Temobs = new List<Obs>();//对于一个周期的数据进行暂存


                        int BDS_count = 0;
                        for(int j=0;j< num_oneepoch;j++)
                        {
                            line = reader.ReadLine();
                            lines4= line?.Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                            Obs temp_obs = new Obs(true, j, lines4[0], timeEpoch, double.Parse(lines4[1]));
                           
                            Temobs.Add(temp_obs);
                            //对于北斗数据进行处理
                            if (lines4[0].Substring(0, 1)=="C")
                            {
                                BDS_count += 1;
                                //Temidex.Add(BDS_count - 1);
                            }
                        }
                        List<Obs> BDS_obs = Temobs.Take(BDS_count).ToList();
                        
                        Obs_t temp_obst = new Obs_t();
                        temp_obst.num = BDS_obs.Count;
                        temp_obst.obs = BDS_obs;
                        //temp_obst.idx = Temidex;
                        data.Obs_file.Add(temp_obst);
                    }
                }
                else
                {
                    MessageBox.Show("数据导入失败！请重试");
                }
            }
            catch (Exception)
            {
                MessageBox.Show("文件导入失败，请重新导入！");
                throw;
            }

        }

        /* 读取导航文件*/
        public void OpenNavfile(DataCenter data)
        {
            try
            {
                OpenFileDialog opf = new OpenFileDialog();
                opf.Filter = "文本文件（*.22p）|*.22p";
                opf.Title = "请输入要导入的导航文件";
                if (opf.ShowDialog() == DialogResult.OK)
                {
                    var reader = new StreamReader(opf.FileName);
                    string line;
                    int countCLines = 0;
                    data.Obv_file.Ord = new List<Orb_par>();
                    while ((line = reader.ReadLine()) != null)
                    {
                        if (line.StartsWith("C"))
                        {
                            countCLines++;
                        }
                    }
                    //重置文件指针到开头
                    reader.BaseStream.Seek(0, SeekOrigin.Begin);
                    /*读取从第六行开始的电离层8个参数*/
                    int linecount1 = 5;
                    while (linecount1 > 0)
                    {
                        line = reader.ReadLine();
                        linecount1--;
                    }
                    //GPSA
                    line = reader.ReadLine();
                    var lines = line?.Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                    data.Obv_file.ion_cmp = new double[8];
                    data.Obv_file.ion_cmp[0] = double.Parse(lines[1]);
                    data.Obv_file.ion_cmp[1] = double.Parse(lines[2]);
                    data.Obv_file.ion_cmp[2] = double.Parse(lines[3]);
                    data.Obv_file.ion_cmp[3] = double.Parse(lines[4]);
                    //GPSB
                    line = reader.ReadLine();
                    lines = line?.Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                    data.Obv_file.ion_cmp[4] = double.Parse(lines[1]);
                    data.Obv_file.ion_cmp[5] = double.Parse(lines[2]);
                    data.Obv_file.ion_cmp[6] = double.Parse(lines[3]);
                    data.Obv_file.ion_cmp[7] = double.Parse(lines[4]);

                    while ((line = reader.ReadLine()) != null)
                    {
                        if (line.StartsWith("C"))
                        {

                            line = Regex.Replace(line, @"(?<=\d)-(?=\d)", " -");
                            lines = line?.Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                            Orb_par temp_orb = new Orb_par();
                            temp_orb.PRN = lines[0];
                            lines = line?.Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Skip(1).ToArray();
                            //DateTimeOffset timeEpoch = str2DateTimeOffset(lines4);
                            temp_orb.ti = str2DateTimeOffset(lines);
                            //卫星钟误差修正参数
                            temp_orb.f0 = double.Parse(lines[6]);
                            temp_orb.f1 = double.Parse(lines[7]);
                            temp_orb.f2 = double.Parse(lines[8]);
                            //C_rs (meters)
                            //-Delta n(radians / sec)
                            //- M0(radians)
                            line = reader.ReadLine();
                            line = Regex.Replace(line, @"(?<=\d)-(?=\d)", " -");
                            lines = line?.Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                            temp_orb.crs = double.Parse(lines[1]);
                            temp_orb.detn= double.Parse(lines[2]);
                            temp_orb.M0= double.Parse(lines[3]);

                            line = reader.ReadLine();
                            line = Regex.Replace(line, @"(?<=\d)-(?=\d)", " -");
                            lines = line?.Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                            /*- C_uc (radians)
                                - e Eccentricity
                                - C_us (radians)
                                - sqrt(A) (sqrt(m))*/
                            temp_orb.cuc= double.Parse(lines[0]);
                            temp_orb.es= double.Parse(lines[1]);
                            temp_orb.cus= double.Parse(lines[2]);
                            temp_orb.sqrt_as= double.Parse(lines[3]);

                            line = reader.ReadLine();
                            line = Regex.Replace(line, @"(?<=\d)-(?=\d)", " -");
                            lines = line?.Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                            //  -T_oe Time of Ephemeris(sec of GPS wk)周内秒数
                            //- C_ic(radians)
                            //- OMEGA0(radians)
                            //- C_is(radians)
                            //tm的又再整时间，真的讨厌，去死吧
                            temp_orb.toe= double.Parse(lines[0]);
                            temp_orb.cic= double.Parse(lines[1]);
                            temp_orb.omega0= double.Parse(lines[2]);
                            temp_orb.cis= double.Parse(lines[3]);
                            // -i0(radians)
                            //- C_rc(meters)
                            //- omega(radians)
                            //- OMEGA DOT(radians / sec)
                            line = reader.ReadLine();
                            line = Regex.Replace(line, @"(?<=\d)-(?=\d)", " -");
                            lines = line?.Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                            temp_orb.i0= double.Parse(lines[0]);
                            temp_orb.crc= double.Parse(lines[1]);
                            temp_orb.w= double.Parse(lines[2]);
                            temp_orb.detomega= double.Parse(lines[3]);

                            line = reader.ReadLine();
                            line = Regex.Replace(line, @"(?<=\d)-(?=\d)", " -");
                            lines = line?.Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                            temp_orb.deti = double.Parse(lines[0]);
                            //存在一个tgd，详情见4.3.1
                            line = reader.ReadLine();
                            line = Regex.Replace(line, @"(?<=\d)-(?=\d)", " -");
                            lines = line?.Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                            temp_orb.sva= double.Parse(lines[0]);
                            temp_orb.tgd= double.Parse(lines[2]);

                            line = reader.ReadLine();
                            data.Obv_file.Ord.Add(temp_orb);
                        }

                    }


                }
                else
                {
                    MessageBox.Show("数据导入失败！请重试");
                }
            }
            catch (Exception)
            {
                MessageBox.Show("文件导入失败，请重新导入！");
                throw;
            }


        }
        public DateTimeOffset str2DateTimeOffset(string[] lines2)
        {
            int year, mon, day, hour, min, sec;
            double d_sec;
            year = int.Parse(lines2[0]);
            mon = int.Parse(lines2[1]);
            day = int.Parse(lines2[2]);
            hour = int.Parse(lines2[3]);
            min = int.Parse(lines2[4]);
            d_sec = double.Parse(lines2[5]);
            sec = (int)d_sec;
            int fraction = (int)((d_sec - sec) * 10000000);
            DateTimeOffset t = new DateTimeOffset(year, mon, day, hour, min, sec, 0,TimeSpan.Zero).AddTicks(fraction);
            return t;
        }
    }
}
