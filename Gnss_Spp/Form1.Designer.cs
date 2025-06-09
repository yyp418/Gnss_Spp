namespace Gnss_Spp
{
    partial class Form1
    {
        /// <summary>
        /// 必需的设计器变量。
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// 清理所有正在使用的资源。
        /// </summary>
        /// <param name="disposing">如果应释放托管资源，为 true；否则为 false。</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows 窗体设计器生成的代码

        /// <summary>
        /// 设计器支持所需的方法 - 不要修改
        /// 使用代码编辑器修改此方法的内容。
        /// </summary>
        private void InitializeComponent()
        {
            this.components = new System.ComponentModel.Container();
            this.menuStrip1 = new System.Windows.Forms.MenuStrip();
            this.文件ToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.观测文件ToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.导航文件ToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.contextMenuStrip1 = new System.Windows.Forms.ContextMenuStrip(this.components);
            this.计算ToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.sPPToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.menuStrip1.SuspendLayout();
            this.SuspendLayout();
            // 
            // menuStrip1
            // 
            this.menuStrip1.ImageScalingSize = new System.Drawing.Size(20, 20);
            this.menuStrip1.Items.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.文件ToolStripMenuItem,
            this.计算ToolStripMenuItem});
            this.menuStrip1.Location = new System.Drawing.Point(0, 0);
            this.menuStrip1.Name = "menuStrip1";
            this.menuStrip1.Size = new System.Drawing.Size(800, 28);
            this.menuStrip1.TabIndex = 0;
            this.menuStrip1.Text = "menuStrip1";
            // 
            // 文件ToolStripMenuItem
            // 
            this.文件ToolStripMenuItem.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.观测文件ToolStripMenuItem,
            this.导航文件ToolStripMenuItem});
            this.文件ToolStripMenuItem.Name = "文件ToolStripMenuItem";
            this.文件ToolStripMenuItem.Size = new System.Drawing.Size(53, 24);
            this.文件ToolStripMenuItem.Text = "文件";
            // 
            // 观测文件ToolStripMenuItem
            // 
            this.观测文件ToolStripMenuItem.Name = "观测文件ToolStripMenuItem";
            this.观测文件ToolStripMenuItem.Size = new System.Drawing.Size(224, 26);
            this.观测文件ToolStripMenuItem.Text = "观测文件";
            this.观测文件ToolStripMenuItem.Click += new System.EventHandler(this.观测文件ToolStripMenuItem_Click);
            // 
            // 导航文件ToolStripMenuItem
            // 
            this.导航文件ToolStripMenuItem.Name = "导航文件ToolStripMenuItem";
            this.导航文件ToolStripMenuItem.Size = new System.Drawing.Size(224, 26);
            this.导航文件ToolStripMenuItem.Text = "导航文件";
            this.导航文件ToolStripMenuItem.Click += new System.EventHandler(this.导航文件ToolStripMenuItem_Click);
            // 
            // contextMenuStrip1
            // 
            this.contextMenuStrip1.ImageScalingSize = new System.Drawing.Size(20, 20);
            this.contextMenuStrip1.Name = "contextMenuStrip1";
            this.contextMenuStrip1.Size = new System.Drawing.Size(61, 4);
            // 
            // 计算ToolStripMenuItem
            // 
            this.计算ToolStripMenuItem.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.sPPToolStripMenuItem});
            this.计算ToolStripMenuItem.Name = "计算ToolStripMenuItem";
            this.计算ToolStripMenuItem.Size = new System.Drawing.Size(53, 24);
            this.计算ToolStripMenuItem.Text = "计算";
            // 
            // sPPToolStripMenuItem
            // 
            this.sPPToolStripMenuItem.Name = "sPPToolStripMenuItem";
            this.sPPToolStripMenuItem.Size = new System.Drawing.Size(224, 26);
            this.sPPToolStripMenuItem.Text = "SPP";
            this.sPPToolStripMenuItem.Click += new System.EventHandler(this.sPPToolStripMenuItem_Click);
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(8F, 15F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(800, 450);
            this.Controls.Add(this.menuStrip1);
            this.MainMenuStrip = this.menuStrip1;
            this.Name = "Form1";
            this.Text = "Form1";
            this.menuStrip1.ResumeLayout(false);
            this.menuStrip1.PerformLayout();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        private System.Windows.Forms.MenuStrip menuStrip1;
        private System.Windows.Forms.ToolStripMenuItem 文件ToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem 观测文件ToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem 导航文件ToolStripMenuItem;
        private System.Windows.Forms.ContextMenuStrip contextMenuStrip1;
        private System.Windows.Forms.ToolStripMenuItem 计算ToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem sPPToolStripMenuItem;

        #endregion

        //private System.Windows.Forms.ToolStripContainer toolStripContainer1;
    }
}

