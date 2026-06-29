

A deep learning-based system for detecting anomalous transposable element (TE) distribution patterns in genomic scaffolds. This tool combines three complementary methods: CNN-based autoencoder, pattern analysis, and outlier detection to identify scaffolds with unusual TE distributions.

## 🌟 Features
+ **Multi-method Detection**: Combines CNN autoencoder, pattern analysis, and statistical outlier detection
+ **Automated Training**: Trains on normal TE distributions and learns to identify anomalies
+ **Flexible Modes**: Supports separate training/testing or complete pipeline
+ **Comprehensive Output**: Generates visualizations, detailed reports, and CSV results
+ **Batch Processing**: Handles multiple input files and large datasets efficiently



## 🔧 Installation
### Prerequisites
+ Python 3.8 or higher
+ CUDA-capable GPU (optional)
+ 8GB+ RAM (16GB+ recommended)

### Step 1: Clone the Repository
```bash
git clone https://github.com/yourusername/te-anomaly-detection.git
cd te-anomaly-detection
```

### Step 2: Create Virtual Environment
```bash
# Using conda (recommended)
conda create -n te_detection python=3.9
conda activate te_detection

# Or using venv
python -m venv te_detection_env
source te_detection_env/bin/activate  # On Windows: te_detection_env\Scripts\activate
```

### Step 3: Install Dependencies
```bash
# For GPU support (CUDA 11.7+)
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu117

# For CPU only
pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu

# Install other dependencies
pip install -r requirements.txt
```

**requirements.txt:**

```plain
numpy>=1.21.0
pandas>=1.3.0
matplotlib>=3.4.0
seaborn>=0.11.0
scikit-learn>=1.0.0
torch>=2.0.0
```

### Step 4: Verify Installation
```bash
python te_detection.py --help
```

---

## 🚀 Quick Start
### Example 1: Complete Pipeline (Train + Test)
```bash
python te_detection.py \
  --mode both \
  --train-dir ./data/train/ \
  --test-dir ./data/test/ \
  --epochs 100 \
  --batch-size 16
```

### Example 2: Training Only
```bash
python te_detection.py \
  --mode train \
  --train-dir ./data/train/ \
  --epochs 100 \
  --lr 0.0001
```

### Example 3: Testing with Pre-trained Model
```bash
python te_detection.py \
  --mode test \
  --test-dir ./data/test/ \
  --z-threshold 5.0
```

---

## 📊 Input Data Format
### File Structure
Input files should be **tab-separated** text files (`.txt`) with the following columns:

| Column | Description | Example |
| --- | --- | --- |
| `Chr` | Chromosome/scaffold name | scaffold_1 |
| `From` | Window start position (bp) | 0 |
| `To` | Window end position (bp) | 100000 |
| `LTR` | LTR retrotransposon density | 0.123 |
| `LINE` | LINE element density | 0.045 |
| `SINE` | SINE element density | 0.012 |
| `DNA` | DNA transposon density | 0.067 |
| `unknown` | Unknown TE density | 0.003 |
| `TE` | Total TE density | 0.250 |


### Example File Format
```plain
Chr	From	To	LTR	LINE	SINE	DNA	unknown	TE
scaffold_1	0	100000	0.123	0.045	0.012	0.067	0.003	0.250
scaffold_1	100000	200000	0.134	0.039	0.015	0.071	0.002	0.261
scaffold_1	200000	300000	0.128	0.042	0.013	0.069	0.004	0.256
...
scaffold_2	0	100000	0.089	0.056	0.018	0.045	0.001	0.209
...
```

### Requirements
+ **100 windows per scaffold**: Each scaffold should be divided into exactly 100 equal-length windows
+ **Minimum scaffold length**: 10 Mb (configurable with `--min-length`)
+ **Minimum windows**: 15 windows (configurable internally)
+ **Tab-separated**: Columns must be separated by tabs
+ **No header row for **`**whole_genome**`: Exclude genome-wide summary rows

### Directory Structure
```plain
data/
├── train/
│   ├── species1_TE_distribution.txt
│   ├── species2_TE_distribution.txt
│   └── species3_TE_distribution.txt
└── test/
    ├── sample1_TE_distribution.txt
    └── sample2_TE_distribution.txt
```

---

## 📖 Usage
### Command-Line Interface
```bash
python te_detection.py [OPTIONS]
```

### Interactive Mode
If you run without `--mode`, the program will prompt you:

```bash
python te_detection.py

# Output:
# 🧬 TE Transposon Anomaly Detection System
# ============================================================
# 
# 📋 Select running mode:
#   [1] Training - Train new model and save
#   [2] Testing - Detect using trained model
#   [3] Complete - Train + Test
# ------------------------------------------------------------
# 
# 👉 Enter choice (1/2/3):
```

---

## ⚙️ Parameters
### Required Parameters
| Parameter | Description | Example |
| --- | --- | --- |
| `--mode` | Running mode: `train`<br/>, `test`<br/>, or `both` | `--mode test` |


### File Paths
| Parameter | Description | Default |
| --- | --- | --- |
| `--train-dir` | Training data directory | None (will prompt) |
| `--test-dir` | Testing data directory | None (will prompt) |


### Training Parameters
| Parameter | Description | Default |
| --- | --- | --- |
| `--epochs` | Number of training epochs | 100 |
| `--batch-size` | Batch size for training | 16 |
| `--lr` | Learning rate | 0.0001 |
| `--min-length` | Minimum scaffold length (bp) | 10000000 |


### Detection Thresholds
| Parameter | Description | Default | Range |
| --- | --- | --- | --- |
| `--z-threshold` | Z-score threshold for outlier detection | 4.5 | > 0 |
| `--cnn-doubt-percentile` | CNN doubt threshold percentile | 85 | 0-100 |
| `--cnn-anomaly-percentile` | CNN anomaly threshold percentile | 90 | 0-100 |


### Examples
```bash
# High sensitivity (more anomalies detected)
python te_detection.py --mode test --test-dir ./data/ \
  --z-threshold 3.0 \
  --cnn-doubt-percentile 80 \
  --cnn-anomaly-percentile 85

# Low sensitivity (fewer anomalies)
python te_detection.py --mode test --test-dir ./data/ \
  --z-threshold 5.0 \
  --cnn-doubt-percentile 90 \
  --cnn-anomaly-percentile 95

# Fast training (for testing)
python te_detection.py --mode train --train-dir ./data/ \
  --epochs 50 \
  --batch-size 32

# GPU training with custom learning rate
python te_detection.py --mode train --train-dir ./data/ \
  --epochs 200 \
  --lr 0.0005
```

---

## 📁 Output Files
### Generated Files
After running the detection, you'll find:

```plain
results/
├── te_anomaly_results.csv        # Complete results for all scaffolds
├── anomalous_scaffolds.csv       # Only anomalous scaffolds
└── training_curve.png            # Training/validation loss curve (train mode)

best_anomaly_model.pth            # Trained model file
```

### Result Columns
**te_anomaly_results.csv:**

| Column | Description | Values |
| --- | --- | --- |
| `Scaffold` | Scaffold name | scaffold_1 |
| `CNN_Score` | Reconstruction error from autoencoder | 0.0 - 1.0+ |
| `Pattern_Score` | Distribution pattern anomaly score | 0.0 - 1.0 |
| `Pattern_Type` | Detected pattern type | valley_pattern_anomaly, normal_peak_pattern, etc. |
| `Pattern_Status` | Pattern-based classification | normal, doubt, anomaly |
| `CNN_Status` | CNN-based classification | normal, doubt, anomaly |
| `Has_Outliers` | Whether scaffold has outlier windows | True/False |
| `Final_Status` | Combined final classification | **normal**, **doubt**, **anomaly** |
| `Reason` | Explanation for classification | Contains 3 delta-TE outliers |
| `Source_File` | Original input file | species1_TE_distribution.txt |


### Visualization
The system generates plots showing:

+ TE composition (LTR, LINE, SINE, DNA, Unknown)
+ Total TE distribution along scaffold
+ Anomaly classification and scores

**Only non-normal scaffolds (anomalous or doubtful) are plotted** (max 22 plots).

---

## 📚 Examples
### Example 1: Training on Multiple Species
```bash
# Prepare training data from multiple species
mkdir -p data/train
cp species1_TE.txt species2_TE.txt species3_TE.txt data/train/

# Train model
python te_detection.py \
  --mode train \
  --train-dir data/train/ \
  --epochs 150 \
  --batch-size 16 \
  --lr 0.0001

# Output:
# 📂 Scanning directory: data/train/
#    Found: 3 files
# ...
# ✅ Model training completed
# 📁 Best model saved to: best_anomaly_model.pth
```

### Example 2: Detecting Anomalies in New Samples
```bash
# Test on new samples
python te_detection.py \
  --mode test \
  --test-dir data/test/ \
  --z-threshold 4.5

# Output:
# Total test scaffolds: 45
# Detected anomalies: 7
# Anomaly rate: 15.56%
# 
# 📋 Anomalous scaffold details:
# ================================================================================
#  1. scaffold_12 (source: sample1_TE_distribution.txt)
#     CNN score: 0.2341
#     Pattern score: 0.8500
#     Pattern type: valley_pattern_anomaly
#     Reason: Pattern detected anomaly
#     Has outliers: No
# ...
```

### Example 3: Adjusting Sensitivity
```bash
# More sensitive detection
python te_detection.py \
  --mode test \
  --test-dir data/test/ \
  --z-threshold 3.5 \
  --cnn-doubt-percentile 75 \
  --cnn-anomaly-percentile 85

# Less sensitive (stricter)
python te_detection.py \
  --mode test \
  --test-dir data/test/ \
  --z-threshold 5.5 \
  --cnn-doubt-percentile 90 \
  --cnn-anomaly-percentile 95
```

---

## 🔬 How It Works
### Three-Method Detection System
The system combines three complementary anomaly detection approaches:

#### 1. **CNN Autoencoder** (Reconstruction-based)
+ Trains a convolutional autoencoder on normal TE distributions
+ Learns to compress and reconstruct typical patterns
+ High reconstruction error = anomalous pattern
+ **Best for**: Global distribution abnormalities

#### 2. **Pattern Analysis** (Rule-based)
+ Analyzes TE distribution shape (peak vs valley patterns)
+ Normal pattern: High TE in middle, low at ends (chromosome arms)
+ Anomalous pattern: High TE at ends, low in middle (telomeric clusters?)
+ **Best for**: Specific biological patterns

#### 3. **Outlier Detection** (Statistical)
+ Computes ΔTE (change in TE density between adjacent windows)
+ Uses Z-score to identify sudden spikes or drops
+ Threshold: |Z| > 4.5 (configurable)
+ **Best for**: Localized TE hotspots or insertion events

### Classification Logic
```plain
if has_outliers:
    Final Status = ANOMALY
    Reason = "Contains N delta-TE outliers"
else:
    Pattern Status = f(pattern_score)  # 0-0.3: normal, 0.3-0.7: doubt, 0.7-1.0: anomaly
    CNN Status = f(cnn_score, thresholds)  # Based on trained percentiles
    
    if Pattern=ANOMALY or CNN=ANOMALY:
        Final Status = ANOMALY
    elif Pattern=DOUBT or CNN=DOUBT:
        Final Status = DOUBT
    else:
        Final Status = NORMAL
```

### Training Process
1. **Data splitting**: 60% train, 40% validation
2. **Normalization**: StandardScaler on TE densities
3. **Model**: CNN encoder (3 conv layers) + MLP decoder
4. **Loss**: MSE reconstruction loss
5. **Early stopping**: Saves best model based on validation loss
6. **Threshold computation**: Percentiles computed on validation set reconstruction errors

### Model Architecture
```plain
Input: [batch, 1, 100] (100 windows per scaffold)
  ↓
CNN Encoder:
  Conv1d(1→64, k=8) + ReLU
  Conv1d(64→128, k=12, s=2) + ReLU
  Conv1d(128→256, k=8, s=2) + ReLU
  AdaptiveAvgPool1d(8)
  Flatten → [batch, 2048]
  ↓
MLP Decoder:
  Linear(2048→1024) + ReLU
  Linear(1024→512) + ReLU
  Linear(512→100)
  ↓
Output: [batch, 100] (reconstructed TE distribution)
```

---

## 🐛 Troubleshooting
### Common Issues
#### 1. **"No matching files found"**
```bash
# Error: No matching files found: /path/to/data/*.txt

# Solution: Check that files exist and have .txt extension
ls /path/to/data/
# Ensure files are named like: species_TE_distribution.txt
```

#### 2. **"CUDA out of memory"**
```bash
# Error: RuntimeError: CUDA out of memory

# Solution: Reduce batch size
python te_detection.py --mode train --train-dir ./data/ --batch-size 8
```

#### 3. **"Model file not found"**
```bash
# Error: Model file best_anomaly_model.pth not found

# Solution: Train the model first
python te_detection.py --mode train --train-dir ./data/train/
```

#### 4. **"All encodings failed"**
```bash
# Error: All encodings failed, skipping file

# Solution: Check file encoding and format
file your_data.txt
# Should be ASCII or UTF-8 text

# Convert if needed:
iconv -f GBK -t UTF-8 input.txt > output.txt
```

#### 5. **Scaffolds are filtered out**
```bash
# Warning: Filtered out: 50 scaffolds

# Reason: Scaffolds too short or too few windows
# Solution: Adjust --min-length parameter
python te_detection.py --mode test --test-dir ./data/ --min-length 5000000
```

---

## 💡 Tips and Best Practices
### For Training
1. **Use diverse normal samples**: Include multiple species/populations with typical TE distributions
2. **Training size**: Recommended 50+ scaffolds for robust model
3. **Validation split**: Default 40% is good; adjust for small datasets
4. **Epochs**: 100-200 epochs usually sufficient; monitor validation loss
5. **GPU usage**: Training on GPU is ~10-50x faster

### For Testing
1. **Consistent preprocessing**: Test data must have same format as training data
2. **Threshold tuning**: Adjust based on false positive/negative rates
3. **Multiple runs**: Test with different thresholds to find optimal sensitivity
4. **Visual inspection**: Always check plotted anomalies for false positives

### Performance Optimization
```bash
# Fast training (testing purposes)
python te_detection.py --mode train --train-dir ./data/ \
  --epochs 50 --batch-size 32

# High-quality training (production)
python te_detection.py --mode train --train-dir ./data/ \
  --epochs 200 --batch-size 16 --lr 0.00005

# Batch testing (many samples)
python te_detection.py --mode test --test-dir ./large_dataset/ \
  --batch-size 64  # Faster inference
```

---

## 📊 Expected Runtime
| Task | Dataset Size | GPU | CPU | Time |
| --- | --- | --- | --- | --- |
| Training | 100 scaffolds, 100 epochs | RTX 3090 | i9-12900K | 5 min / 45 min |
| Training | 500 scaffolds, 100 epochs | RTX 3090 | i9-12900K | 20 min / 3 hrs |
| Testing | 100 scaffolds | RTX 3090 | i9-12900K | 10 sec / 1 min |
| Testing | 1000 scaffolds | RTX 3090 | i9-12900K | 1 min / 10 min |


---

## 📝 Citation
If you use this tool in your research, please cite:

```plain
@software{te_anomaly_detection,
  author = {Your Name},
  title = {TE Anomaly Detection System},
  year = {2025},
  url = {https://github.com/yourusername/te-anomaly-detection}
}
```

---

## 🤝 Contributing
Contributions are welcome! Please feel free to submit a Pull Request.

### Development Setup
```bash
git clone https://github.com/yourusername/te-anomaly-detection.git
cd te-anomaly-detection
pip install -r requirements-dev.txt
```

### Running Tests
```bash
pytest tests/
```

---


## 🙏 Acknowledgments
+ PyTorch team for the deep learning framework
+ scikit-learn for preprocessing tools
+ matplotlib/seaborn for visualization

---

---

**Last updated**: 2025-01-XX

