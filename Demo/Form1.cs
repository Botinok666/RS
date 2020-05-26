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
        private string testFileName;
        private long testFileLength;
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
                        testFileLength = fstream.Length;
                        testFile = new byte[fstream.Length + 256];
                        await fstream.ReadAsync(testFile, 0, (int)fstream.Length);
                        label1.Text = string.Format("Размер файла: {0} Кб", fstream.Length >> 10);
                        Array.Clear(testFile, (int)fstream.Length, 256);
                        fstream.Close();
                    }
                    blockSzCB.Enabled = true;
                    blockSzCB.SelectedIndex = 0;
                    testBtn.Enabled = true;
                    testFileName = openFile.FileName;
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
            List<RSCS> rs = new List<RSCS>();
            Enum.GetValues(typeof(RSCS.Coders))
                .OfType<RSCS.Coders>()
                .ToList()
                .ForEach(x => {
                    try
                    {
                        rs.Add(new RSCS((byte)n, (byte)kTest, x));
                    }
                    catch (PlatformNotSupportedException)
                    { }
                });

            await Task.Run(() => {
                WriteTextSafe(string.Format("\tТестирование RS({0}, {1}){2}", n, kTest, Environment.NewLine));
                Stopwatch stopwatch = new Stopwatch();
                foreach (var coder in rs)
                {
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
            Enum.GetNames(typeof(RSCS.Coders))
                .ToList()
                .ForEach(x => {
                    dataTable.Columns.Add(x, typeof(float));
                    chart.Series.Add(new Series(x)
                    {
                        XValueMember = "Name",
                        XValueType = ChartValueType.String,
                        YValueMembers = x,
                        YValueType = ChartValueType.Single,
                        ChartType = SeriesChartType.Column
                    });
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

        private void EncodeBtn_Click(object sender, EventArgs e)
        {
            if (testFile == null || testFile.Length < 1024)
                return;
            int fileSz = testFile.Length - 256;
            int t = 32, k = int.Parse(blockSzCB.SelectedItem.ToString());
            int blockCount = 1 + fileSz / k, n = k + t * 2;
            byte[] tempFile = new byte[blockCount * n];
            RSCS rs = new RSCS((byte)n, (byte)k);

            for (int j = 0; j < blockCount; j++)
            {
                Buffer.BlockCopy(testFile, j * k, tempFile, j * n, k);
                rs.EncodeBuffer(ref tempFile[j * n]);
            }

            using SaveFileDialog saveFile = new SaveFileDialog
            {
                FileName = Path.ChangeExtension(testFileName, "bin")
            };
            if (saveFile.ShowDialog() != DialogResult.OK)
                return;
            using (BinaryWriter writer = new BinaryWriter(File.Create(saveFile.FileName)))
            {
                writer.Write(n);
                writer.Write(k);
                writer.Write(Path.GetExtension(testFileName));
                writer.Write(testFileLength);
                writer.Write(tempFile);
                writer.Flush();
            }
            WriteTextSafe("Файл " + Path.GetFileName(saveFile.FileName) + " сохранён" + Environment.NewLine);
        }

        private void CorruptBtn_Click(object sender, EventArgs e)
        {
            using OpenFileDialog openFile = new OpenFileDialog
            {
                Filter = "RS files (*.bin)|*.bin"
            };
            if (openFile.ShowDialog() != DialogResult.OK)
                return;
            if (!File.Exists(openFile.FileName))
                return;
            using BinaryReader reader = new BinaryReader(File.OpenRead(openFile.FileName));
            int n = reader.ReadInt32();
            int k = reader.ReadInt32();
            string ext = reader.ReadString();
            long length = reader.ReadInt64();
            byte[] tempFile = reader.ReadBytes((int)(reader.BaseStream.Length - reader.BaseStream.Position));
            //Introduce some errors
            int t = (n - k) >> 1;
            for (int j = 0; j < tempFile.Length / n; j++)
            {
                int errors = Randoms.Next(0, t + 1);
                for (int m = 0; m < errors; m++)
                {
                    ushort randn = (ushort)Randoms.Next(ushort.MinValue, ushort.MaxValue);
                    int idx = j * n + randn % n; //Random position within a block
                    randn >>= 8; //Upper byte as a random value
                    tempFile[idx] ^= (byte)(randn != 0 ? randn : 1);
                }
            }
            string corruptName = Path.Combine(Path.GetDirectoryName(openFile.FileName), 
                Path.GetFileNameWithoutExtension(openFile.FileName) + "-corrupt.bin");
            using BinaryWriter writer = new BinaryWriter(File.Create(corruptName));
            writer.Write(n);
            writer.Write(k);
            writer.Write(ext);
            writer.Write(length);
            writer.Write(tempFile);
            writer.Flush();
            WriteTextSafe("Файл " + Path.GetFileName(corruptName) + " сохранён" + Environment.NewLine);
        }

        private void DecodeBtn_Click(object sender, EventArgs e)
        {
            using OpenFileDialog openFile = new OpenFileDialog
            {
                Filter = "RS files (*.bin)|*.bin"
            };
            if (openFile.ShowDialog() != DialogResult.OK)
                return;
            if (!File.Exists(openFile.FileName))
                return;
            using BinaryReader reader = new BinaryReader(File.OpenRead(openFile.FileName));
            int n = reader.ReadInt32();
            int k = reader.ReadInt32();
            string ext = reader.ReadString();
            long length = reader.ReadInt64();
            byte[] tempFile = reader.ReadBytes((int)(reader.BaseStream.Length - reader.BaseStream.Position));
            byte[] origFile = new byte[length];

            RSCS rs = new RSCS((byte)n, (byte)k);
            int blocks = tempFile.Length / n;
            int errors = 0;
            for (int j = 0; j < blocks - 1; j++)
            {
                errors += rs.DecodeBuffer(ref tempFile[j * n]);
                Buffer.BlockCopy(tempFile, j * n, origFile, j * k, k);
            }
            errors += rs.DecodeBuffer(ref tempFile[(blocks - 1) * n]);
            Buffer.BlockCopy(tempFile, (blocks - 1) * n, origFile, (blocks - 1) * k, origFile.Length - (blocks - 1) * k);
            WriteTextSafe("Исправлено " + errors.ToString() + " ошибок" + Environment.NewLine);

            string restoredName = Path.Combine(Path.GetDirectoryName(openFile.FileName),
                Path.GetFileNameWithoutExtension(openFile.FileName) + "-restored" + ext);
            using BinaryWriter writer = new BinaryWriter(File.Create(restoredName));
            writer.Write(origFile);
            writer.Flush();
            WriteTextSafe("Файл " + Path.GetFileName(restoredName) + " сохранён" + Environment.NewLine);
        }
    }
}
