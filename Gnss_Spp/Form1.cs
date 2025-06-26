using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Diagnostics;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Gnss_Spp
{
    public partial class Form1 : Form
    {
        DataCenter data=new DataCenter();
        public Form1()
        {
            InitializeComponent();
        }

        private void 观测文件ToolStripMenuItem_Click(object sender, EventArgs e)
        {
            Stopwatch stopwatch = Stopwatch.StartNew();
            FileHelper fhelper=new FileHelper();
            fhelper.OpenObsfile(data);
            stopwatch.Stop();
            Console.WriteLine($"代码执行时间: {stopwatch.ElapsedMilliseconds} 毫秒");
        }

        private void 导航文件ToolStripMenuItem_Click(object sender, EventArgs e)
        {
            Stopwatch stopwatch = Stopwatch.StartNew();
            FileHelper fhelper = new FileHelper();
            fhelper.OpenNavfile(data);
            stopwatch.Stop();
            Console.WriteLine($"代码执行时间: {stopwatch.ElapsedMilliseconds} 毫秒");
        }

        private void sPPToolStripMenuItem_Click(object sender, EventArgs e)
        {
            Calute cal = new Calute();
            DenseMatrix SatPosClkValid = cal.cal_epoch(0, data);
            DenseMatrix B;
            DenseMatrix L;
            cal.Fromcoefficient(0,  SatPosClkValid, data, out  B, out L, out double[] varIon, out double[] varTrp);
        }
    }
     

}
