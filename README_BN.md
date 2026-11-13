# TriDiseaseRank-X

**একাধিক Gene Expression Dataset বিশ্লেষণ ও exploratory biomarker prioritization-এর জন্য R Shiny ভিত্তিক স্থানীয় Web Application।**

এই প্রজেক্টে মূল R Research Script-গুলোকে একটি সহজ Browser Interface-এর সঙ্গে যুক্ত করা হয়েছে। একসঙ্গে একাধিক Counts CSV ও Metadata CSV আপলোড করে DEG, WGCNA, Cross-Dataset Comparison, Expression-Correlation ভিত্তিক GNN, Machine Learning, Pathway Enrichment এবং Five-Layer Evidence Ranking চালানো যায়।

> **গুরুত্বপূর্ণ:** এই সফটওয়্যারের ফলাফল গবেষণামূলক ও exploratory। তিনটি bundled GEO Dataset-এর biological comparison একই নয়। কোনো ranked gene-কে সরাসরি clinically validated biomarker বলা যাবে না।

## অন্য Windows কম্পিউটারে চালানোর নিয়ম

### যা লাগবে

- Windows 10/11 (64-bit)
- R for Windows (ডেভেলপমেন্টে R 4.6.x ব্যবহার করা হয়েছে)
- প্রথমবার Package ইনস্টল করার জন্য Internet
- বড় Analysis-এর জন্য ৮ GB RAM ন্যূনতম; ১৬ GB বা বেশি সুপারিশ করা হয়
- Source Package compile করতে হলে Rtools45 লাগতে পারে
- RStudio বা VS Code বাধ্যতামূলক নয়

### সহজ Installation

১. GitHub Repository থেকে **Code → Download ZIP** নির্বাচন করে Download ও Extract করো।

২. R ইনস্টল করো: <https://cran.r-project.org/bin/windows/base/>

৩. Extract করা `TriDiseaseRank-X` Folder খুলে `setup.bat`-এ Double Click করো। এটি প্রয়োজনীয় CRAN ও Bioconductor Package ইনস্টল করবে এবং Torch runtime প্রস্তুত করার চেষ্টা করবে।

৪. Setup শেষ হলে `run_app.bat`-এ Double Click করো।

৫. Browser-এ TriDiseaseRank-X খুলবে। ব্যবহার করার সময় Command Window বন্ধ করবে না।

### Manual পদ্ধতি

Project Folder-এ PowerShell বা Command Prompt খুলে চালাও:

```powershell
Rscript scripts\setup_windows.R
Rscript scripts\check_environment.R
Rscript scripts\launch_app.R
```

Torch Runtime-এ সমস্যা হলে R Console-এ চালাও:

```r
torch::install_torch()
```

তারপর R Restart করে Environment Check আবার চালাও।

## CSV Format

### Counts CSV

- প্রতিটি Row একটি Gene।
- `gene_id` এবং `gene_symbol` কলাম থাকতে হবে।
- পরের কলামগুলো Sample ID এবং Raw non-negative integer count হবে।

### Metadata CSV

- `Samples` ও `Condition` কলাম থাকতে হবে।
- `Condition`-এ `Control` ও `Case` ব্যবহার করতে হবে।
- Metadata-এর Sample ID অবশ্যই Counts CSV-এর Sample Column-এর সঙ্গে মিলতে হবে।
- Pairing filename prefix দিয়ে হয়। যেমন:
  - `GSE120852_counts_data.csv`
  - `GSE120852_metadata.csv`

`data/example/` Folder-এ Development-এ ব্যবহৃত তিনটি GEO Dataset-এর উদাহরণ রাখা হয়েছে। চাইলে নিজের compatible CSV আপলোড করতে পারবে।

## Analysis চালানোর ধাপ

১. Counts ও Metadata CSV Upload এবং Validate।
২. সব Dataset-এর DEG চালাও।
৩. WGCNA চালাও।
৪. Cross-Dataset Comparison চালাও।
৫. GNN চালাও (প্রথমে Quick Preview)।
৬. ML চালাও (প্রথমে Quick Test; Full 113-model mode বেশি সময় ও RAM নেয়)।
৭. Pathway Analysis ও Five-Layer Ranking চালাও।
৮. Dashboard থেকে CSV ও Excel Report Download করো।

## Folder Structure

```text
TriDiseaseRank-X/
├── app.R
├── R/                  # Web analysis modules
├── original/           # Original research R scripts
├── data/example/       # Example GEO CSV files
├── www/                # CSS ও web assets
├── scripts/            # Setup, check, launch
├── setup.bat
├── run_app.bat
└── README.md
```

## Output

Analysis-এর Output স্থানীয়ভাবে `TriDiseaseRank_results/` Folder-এ তৈরি হয়। এই Folder GitHub-এ রাখা হয় না, কারণ এখানে বড় ও ব্যবহারকারীর নিজস্ব Analysis Output থাকতে পারে।

## Troubleshooting

- **Rscript পাওয়া যাচ্ছে না:** R ইনস্টল করে Terminal আবার খুলে `where Rscript` পরীক্ষা করো।
- **Package Install Error:** নির্দিষ্ট Package-এর Error পড়ো। Source Build করতে Rtools45 লাগতে পারে। সমস্যা সমাধানের পর `setup.bat` আবার চালাও।
- **Torch Error:** R Console-এ `torch::install_torch()` চালিয়ে R Restart করো।
- **`R/` বা `original/` ফাইল পাওয়া যাচ্ছে না:** Project Root থেকে `run_app.bat` চালাও।
- **Pathway Result ফাঁকা:** Significant Enrichment না পাওয়া বৈধ ফলাফল হতে পারে। Optional Database Package না থাকলেও Core Analysis চলতে পারে।
- **Analysis ধীর:** Quick Preview/Quick Test ও কম Candidate Gene দিয়ে শুরু করো। Full GNN ও 113-model ML-এ বেশি RAM ও সময় লাগে।

## বৈজ্ঞানিক সীমাবদ্ধতা

- বর্তমান GNN Gene Expression Correlation থেকে Network তৈরি করে; এটি নিজে থেকে curated STRING/BioGRID/IntAct PPI database ব্যবহার করে না।
- তিনটি GEO Dataset-এর Study Design ও Biological Contrast ভিন্ন।
- Internal Cross-validation স্বাধীন External Validation-এর বিকল্প নয়। Cross-validation-এর আগে feature selection হলে performance optimistic হতে পারে।
- Five-Layer Score exploratory prioritization score; এটি clinical probability বা diagnostic score নয়।

মূল Research Script-গুলো `original/` Folder-এ সংরক্ষিত আছে। Web Integration-এর কোড `R/` Folder-এ রয়েছে।
