using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using MersenneTwister;

namespace Demo
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }
        struct BenchData { public float time; public float speed; public string type; public string mode; }
        private readonly DataTable dataTable = new DataTable("Default");
        private byte[] testFile;
        private int[] tTest;
        private int kTest;
        private CancellationTokenSource tokenSource;
        private Action finishBench = null; 
        private delegate void SafeCallDelegate(BenchData data);
        private delegate void SafeCallDelegateText(string text);
        private Chart chart;
        private readonly List<BenchData> benchDatas = new List<BenchData>();
        private async void OpenFileBtn_Click(object sender, EventArgs e)
        {
            openFileBtn.Enabled = false;
            using (OpenFileDialog openFile = new OpenFileDialog())
            {
                if (openFile.ShowDialog() == DialogResult.OK)
                {
                    using (FileStream fstream = File.OpenRead(openFile.FileName))
                    {
                        testFile = new byte[fstream.Length + 256];
                        await fstream.ReadAsync(testFile, 0, (int)fstream.Length);
                        label1.Text = string.Format("Размер файла: {0} Кб", fstream.Length >> 10);
                        Array.Clear(testFile, (int)fstream.Length, 256);
                        fstream.Close();
                    }
                    blockSzCB.Enabled = true;
                    blockSzCB.SelectedIndex = 0;
                    testBtn.Enabled = true;
                }
            }
            openFileBtn.Enabled = true;
        }

        private void BlockSzCB_SelectedIndexChanged(object sender, EventArgs e)
        {
            if (testFile == null || testFile.Length == 0) return;
            float blockCount = 1 + testFile.Length / int.Parse(blockSzCB.SelectedItem.ToString());
            char modifier = ' ';
            if (blockCount > 1e+6f)
            {
                blockCount /= 1e+6f;
                modifier = 'M';
            }
            if (blockCount > 1e+3f)
            {
                blockCount /= 1e+3f;
                modifier = 'k';
            }
            textBox1.Text = string.Format("Число блоков в файле: {0:F1}{1}{2}", 
                blockCount, modifier, Environment.NewLine);
        }

        private async void TestBtn_Click(object sender, EventArgs e)
        {
            if (testFile == null || testFile.Length == 0) return;
            if (finishBench != null)
            {
                finishBench.Invoke();
                return;
            }
            testBtn.Text = "Отмена";
            chart.Series
                .ToList()
                .ForEach(x => x.Points.Clear());
            benchDatas.Clear();
            kTest = int.Parse(blockSzCB.SelectedItem.ToString());
            switch (kTest)
            {
                case 32: //k=32
                    tTest = new int[] { 8, 16, 32, 48 };
                    break;
                case 64: //k=64
                    tTest = new int[] { 8, 16, 28, 40 };
                    break;
                case 128: //k=128
                case 191: //k=191
                    tTest = new int[] { 8, 16, 24, 32 };
                    break;
            }
            dataTable.Rows.Clear();
            tTest.ToList()
                .ForEach(x =>
                {
                    int n = x * 2 + kTest;
                    dataTable.Rows.Add(string.Format("({0},{1}){2}кодир.", n, kTest, Environment.NewLine), 0, 0, 0);
                    dataTable.Rows.Add(string.Format("({0},{1}){2}чист.дек.", n, kTest, Environment.NewLine), 0, 0, 0);
                    dataTable.Rows.Add(string.Format("({0},{1}){2}декод.", n, kTest, Environment.NewLine), 0, 0, 0);
                });
            chart.DataBind();

            tokenSource = new CancellationTokenSource();
            finishBench = () => { tokenSource.Cancel(); };
            try
            {
                foreach (var x in tTest)
                    await RunBenchmark(x, tokenSource.Token);
            }
            catch (OperationCanceledException)
            {
                MessageBox.Show("Отменено пользователем");
            }
            finishBench = null;
            testBtn.Text = "Тест";
        }
        private async Task RunBenchmark(int t, CancellationToken token)
        {
            int fileSz = testFile.Length - 256;
            int blockCount = 1 + fileSz / kTest, n = kTest + t * 2;
            byte[] tempFile = new byte[blockCount * n];
            byte[] origFile = new byte[testFile.Length];
            RSCS[] rs = new RSCS[] { 
                new RSCS((byte)n, (byte)kTest, RSCS.Coders.ALU),
                null, 
                null };
            try
            {
                rs[1] = new RSCS((byte)n, (byte)kTest, RSCS.Coders.SSSE3);
            }
            catch (PlatformNotSupportedException)
            {
                rs[1] = null;
            }
            try
            {
                rs[2] = new RSCS((byte)n, (byte)kTest, RSCS.Coders.AVX2);
            }
            catch (PlatformNotSupportedException)
            {
                rs[2] = null;
            }

            await Task.Run(() => {
                WriteTextSafe(string.Format("\tТестирование RS({0}, {1}){2}", n, kTest, Environment.NewLine));
                Stopwatch stopwatch = new Stopwatch();
                foreach (var coder in rs)
                {
                    if (coder == null) continue;
                    Array.Clear(testFile, 0, testFile.Length);
                    stopwatch.Restart();
                    for (int j = 0; j < blockCount; j++)
                    {
                        Buffer.BlockCopy(testFile, j * kTest, tempFile, j * n, kTest);
                        coder.EncodeBuffer(ref tempFile[j * n]);
                        if (token.IsCancellationRequested) break;
                    }
                    if (token.IsCancellationRequested) throw new OperationCanceledException("Canceled by user");
                    stopwatch.Stop();
                    WriteTextSafe(new BenchData()
                    {
                        time = stopwatch.ElapsedMilliseconds / 1000f,
                        speed = fileSz / 1024f / stopwatch.ElapsedMilliseconds,
                        mode = string.Format("({0},{1}){2}кодир.", n, kTest, Environment.NewLine),
                        type = coder.Version
                    });
                    Array.Clear(origFile, 0, origFile.Length);
                    stopwatch.Restart();
                    for (int j = 0; j < blockCount; j++)
                    {
                        coder.DecodeBuffer(ref tempFile[j * n]);
                        Buffer.BlockCopy(tempFile, j * n, origFile, j * kTest, kTest);
                        if (token.IsCancellationRequested) break;
                    }
                    if (token.IsCancellationRequested) throw new OperationCanceledException("Canceled by user");
                    stopwatch.Stop();
                    WriteTextSafe(new BenchData() 
                    {
                        time = stopwatch.ElapsedMilliseconds / 1000f,
                        speed = fileSz / 1024f / stopwatch.ElapsedMilliseconds,
                        mode = string.Format("({0},{1}){2}чист.дек.", n, kTest, Environment.NewLine),
                        type = coder.Version 
                    });
                    if (!origFile.SequenceEqual(testFile))
                        WriteTextSafe(string.Format("{0} чист.декод. неверно{1}", coder.Version, Environment.NewLine));

                    int errorSum = 0;
                    //Introduce some errors
                    for (int k = 0; k < blockCount; k++)
                    {
                        int errors = Randoms.Next(0, t + 1);
                        for (int m = 0; m < errors; m++)
                        {
                            ushort randn = (ushort)Randoms.Next(ushort.MinValue, ushort.MaxValue);
                            int idx = k * n + randn % n; //Random position within a block
                            randn >>= 8; //Upper byte as a random value
                            tempFile[idx] ^= (byte)(randn != 0 ? randn : 1);
                        }
                        errorSum += errors;
                    }
                    Array.Clear(origFile, 0, origFile.Length);
                    stopwatch.Restart();
                    for (int j = 0; j < blockCount; j++)
                    {
                        coder.DecodeBuffer(ref tempFile[j * n]);
                        Buffer.BlockCopy(tempFile, j * n, origFile, j * kTest, kTest);
                        if (token.IsCancellationRequested) break;
                    }
                    if (token.IsCancellationRequested) throw new OperationCanceledException("Canceled by user");
                    stopwatch.Stop();
                    WriteTextSafe(new BenchData()
                    {
                        time = stopwatch.ElapsedMilliseconds / 1000f,
                        speed = fileSz / 1024f / stopwatch.ElapsedMilliseconds,
                        mode = string.Format("({0},{1}){2}декод.", n, kTest, Environment.NewLine),
                        type = coder.Version
                    });
                    Array.Clear(origFile, fileSz, 256);
                    if (!origFile.SequenceEqual(testFile))
                        WriteTextSafe(string.Format("{0} полн.декод. неверно{1}", coder.Version, Environment.NewLine));
                }
            });
        }
        private void WriteTextSafe(BenchData data)
        {
            if (textBox1.InvokeRequired)
            {
                var d = new SafeCallDelegate(WriteTextSafe);
                textBox1.Invoke(d, new object[] { data });
            }
            else
            {
                benchDatas.Add(data);
                textBox1.AppendText(string.Format("{0} {1}: {2}Мб/с{3}", 
                    data.type, data.mode.Replace(Environment.NewLine, " "), data.speed, Environment.NewLine));
                if (radioButton1.Checked)
                {
                    dataTable.Rows
                        .Find(data.mode)[data.type] = data.speed;
                }
                else
                {
                    float maxSpeed = benchDatas
                        .FindAll(x => x.mode.Equals(data.mode))
                        .Max(y => y.speed);
                    benchDatas
                        .Where(x => x.mode.Equals(data.mode))
                        .ToList()
                        .ForEach(y =>
                            dataTable.Rows
                                .Find(y.mode)[y.type] = y.speed * 100 / maxSpeed
                        );
                }
                chart.DataBind();
            }
        }
        private void WriteTextSafe(string text)
        {
            if (textBox1.InvokeRequired)
            {
                var d = new SafeCallDelegateText(WriteTextSafe);
                textBox1.Invoke(d, new object[] { text });
            }
            else
            {

                textBox1.AppendText(text);
            }
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            chart = new Chart
            {
                Dock = DockStyle.Fill
            };
            panel1.Controls.Add(chart);

            ChartArea chartArea = new ChartArea("Default");
            chartArea.AxisX.LabelAutoFitMaxFontSize = 8;
            chartArea.AxisX.IsMarksNextToAxis = true;
            chartArea.AxisX.LabelStyle.Interval = 1;
            chartArea.AxisX.MajorGrid.LineDashStyle = ChartDashStyle.Dash;
            chartArea.AxisY.LabelAutoFitMaxFontSize = 8;
            chartArea.AxisY.Title = "Мб/с";
            chart.ChartAreas.Add(chartArea);
            Legend legend = new Legend("Default")
            {
                BackColor = Color.Transparent,
                Docking = Docking.Bottom
            };
            chart.Legends.Add(legend);
            dataTable.Columns.Add("Name", typeof(string));
            dataTable.PrimaryKey = new DataColumn[] { dataTable.Columns["Name"] };
            dataTable.Columns.Add("ALU", typeof(float));
            dataTable.Columns.Add("SSSE3", typeof(float));
            dataTable.Columns.Add("AVX2", typeof(float));
            chart.Series.Add(new Series("ALU")
            {
                XValueMember = "Name",
                XValueType = ChartValueType.String,
                YValueMembers = "ALU",
                YValueType = ChartValueType.Single,
                ChartType = SeriesChartType.Column
            });
            chart.Series.Add(new Series("SSSE3")
            {
                XValueMember = "Name",
                XValueType = ChartValueType.String,
                YValueMembers = "SSSE3",
                YValueType = ChartValueType.Single,
                ChartType = SeriesChartType.Column
            });
            chart.Series.Add(new Series("AVX2")
            {
                XValueMember = "Name",
                XValueType = ChartValueType.String,
                YValueMembers = "AVX2",
                YValueType = ChartValueType.Single,
                ChartType = SeriesChartType.Column
            });
            chart.DataSource = dataTable;
        }

        private void Form1_ResizeEnd(object sender, EventArgs e)
        {
            
        }

        private void RadioButton1_CheckedChanged(object sender, EventArgs e)
        {
            if (radioButton1.Checked)
            {
                chart.ChartAreas.First().AxisY.Title = "Мб/с";
                benchDatas.ForEach(x => dataTable.Rows
                    .Find(x.mode)[x.type] = x.speed);
            }
            else
            {
                chart.ChartAreas.First().AxisY.Title = "%";
                benchDatas
                    .GroupBy(x => x.mode)
                    .Select(y => new KeyValuePair<string, float>(y.Key, y.Max(z => z.speed)))
                    .ToList()
                    .ForEach(z => benchDatas
                        .Where(u => u.mode.Equals(z.Key))
                        .ToList()
                        .ForEach(v => dataTable.Rows
                            .Find(v.mode)[v.type] = v.speed * 100 / z.Value)
                    );
            }
            chart.DataBind();
        }
    }
}
