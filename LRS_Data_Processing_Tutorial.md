# Long-Read Sequencing (LRS) Data Processing Tutorial

Welcome to the Human RNome Project! This tutorial is designed as a step-by-step guide for anyone—from an undergraduate student to a Principal Investigator—to process Long-Read Sequencing (LRS) data. 

In this first part of the tutorial, we will cover how to download the raw data from our portal and run the automated WDL pipelines to perform basecalling and initial processing.

## Step 1: Downloading the Raw Data

In this step, we will download the raw data files generated directly by the sequencing machines. Before we can analyze the RNA, we need this raw electrical signal data (stored in `.pod5` files) to extract the actual nucleotide sequences and map chemical modifications.

1. Open your web browser and navigate to the project's data repository: 
   [https://doi.org/10.25585/DOE-HRP/3377574](https://doi.org/10.25585/DOE-HRP/3377574)
2. Follow the instructions on the portal to download the raw Long-Read Sequencing files (typically `.pod5` format).
3. Save these files into a dedicated directory on your local machine or high-performance computing (HPC) cluster. For this tutorial, we will assume you saved them in `/path/to/pod5_dir`.

> [!TIP]
> Depending on the size of the data, you may need to use tools like Globus or `wget`/`curl` if downloading directly to a remote server or cluster.

## Step 2: Getting the Code (Cloning the Repository)

Next, you will download a complete copy of the project's code to your computer. We need these custom scripts and pipelines developed by the consortium to analyze the data you just downloaded.

Open your terminal (Command Prompt/PowerShell on Windows, Terminal on Mac/Linux) and run the following command to download the code:

```bash
git clone https://github.com/Human-RNome-Project/HRP-benchmarking-project.git
```

**Understanding the command:** 
- `git clone` is a command from the Git version control software that tells your computer to reach out to the internet, copy an entire code repository, and place it exactly as it appears into a new folder on your computer.

This will create a new folder called `HRP-benchmarking-project` in your current directory containing all of the scripts and pipelines.

## Step 3: Navigating to the LRS Pipelines

Now, we need to move our terminal session into the specific folder where the Long-Read Sequencing pipelines live. The terminal needs to be "inside" the correct folder to run the scripts located there, otherwise it won't know where to look.

We use the `cd` (change directory) command to move between folders in the terminal.

1. First, move into the main project folder you just downloaded:
   ```bash
   cd HRP-benchmarking-project
   ```
2. Next, navigate deep into the data processing folders where the LRS WDL pipelines live:
   ```bash
   cd data_processing_pipelines/LRS_Reanalysis/wdl-pipelines-main/LRS/
   ```

You will see three folders here: `rRNA_pipeline`, `polyA_pipeline`, and `tRNA_pipeline`. For this tutorial, we will use the **rRNA pipeline** as an example, but the steps are nearly identical for the others. Move into that folder:

```bash
cd rRNA_pipeline
```

**Understanding the command:**
- `cd` stands for "change directory". It tells the computer to open the folder you specify and make it your current working location in the terminal.

## Step 4: Generating the Input JSON

Before running the pipeline, we must create a configuration file that tells it where your raw data is located. WDL (Workflow Description Language) pipelines act like automated assembly lines, so they need a standardized set of instructions (formatted as a JSON file) to know exactly which files to process and what to name the outputs.

We have provided a helpful bash script to generate this file automatically. To demonstrate this on a real sample from our dataset, we will use **`HRP_A_017_native_rRNA_001`** (you can find this sample ID and its biological details in the `HRP_MetaData_A_LRS.tsv` metadata sheet).

Run the following command, replacing `/path/to/pod5_dir` with the actual folder where you downloaded this specific sample's data:

```bash
./generate_rRNA_inputs.sh /path/to/pod5_dir HRP_A_017_native_rRNA_001
```

**Understanding the command:**
- `./generate_rRNA_inputs.sh` executes a small program (a bash script) that we wrote for you.
- `/path/to/pod5_dir` is the first input to the script. It tells the script where you saved your downloaded data in Step 1.
- `HRP_A_017_native_rRNA_001` is the second input. This is the exact sample ID from the metadata sheet, which the script will use to appropriately name all of your final output files.

This script will automatically write a new configuration file located at `rRNA/inputs_HRP_A_017_native_rRNA_001.json`.

> [!NOTE]
> The `generate_rRNA_inputs.sh` script assumes standard paths for reference genomes and Dorado basecalling models. If your environment is set up differently, you may need to open the script in a text editor and adjust those paths before running it.

## Step 5: Preparing Docker Images (Local Runs Only)

We also need to download the required software containers before starting the pipeline. To ensure the pipeline works seamlessly on different computers, it uses "Docker containers," which act like mini virtual computers with the exact software pre-installed. For a local run, your computer needs to download these images from the internet first so the pipeline engine (like `miniwdl`) can access them.

If you are running this locally (not on a supercomputer like Perlmutter), ensure you have Docker installed and running on your machine. Run the following commands in your terminal to pull the necessary images:

```bash
# Dorado basecaller
docker pull ontresearch/dorado@sha256:c8f356489fa8b44b31beba841b84d2879de2088e

# Seqtagger demultiplexing
docker pull lpryszcz/seqtagger:latest

# Minimap2 alignment
docker pull nanozoo/minimap2:2.28--9e3bd01

# NanoComp QC
docker pull luxendr13/nanocomp:0.6.0

# Modkit modification calling
docker pull ontresearch/modkit@sha256:489d708a48c66368e5d1e118538e5dca68203a64
```

**Understanding the command:**
- `docker pull` tells Docker to reach out to the internet (Docker Hub) and download the specific software image so it is ready and waiting on your hard drive for the WDL pipeline to use.

## Step 6: Running the Pipeline (Basecalling)

Finally, we will launch the automated software that reads the raw signals and converts them into RNA sequences. Because the raw `.pod5` files are simply electrical signals, we use a machine learning tool called **Dorado** (the basecaller) to translate those electrical squiggles into A, C, G, and U nucleotides while detecting chemical modifications on them.

You are now ready to run the primary pipeline! The pipeline (`pipeline_SCATTER_jaws_rRNA.wdl`) will take your raw `.pod5` files and automatically perform demultiplexing and Dorado basecalling (as well as subsequent alignment).

You can run this WDL pipeline using different workflow execution engines depending on your local setup. These engines act as the "managers" that read the WDL file and execute the steps inside it.

### Option A: Local Execution with miniwdl (Recommended for local machines)
If you are running this on a local machine or standard server, `miniwdl` is highly recommended.

```bash
miniwdl run pipeline_SCATTER_jaws_rRNA.wdl --input rRNA/inputs_HRP_A_017_native_rRNA_001.json
```

**Understanding the command:**
- `miniwdl run` tells the `miniwdl` engine to start an automated job.
- `pipeline_SCATTER_jaws_rRNA.wdl` is the master blueprint detailing all the software steps required.
- `--input rRNA/inputs_HRP_A_017_native_rRNA_001.json` provides the configuration file you generated in Step 4, telling the blueprint exactly which data to process.

### Option B: Local Execution with Cromwell
Alternatively, you can use Cromwell:

```bash
java -jar cromwell.jar run pipeline_SCATTER_jaws_rRNA.wdl -i rRNA/inputs_HRP_A_017_native_rRNA_001.json
```

### Option C: HPC Execution with JAWS (NERSC/Perlmutter)
If you are working on the NERSC Perlmutter supercomputer, use JAWS:

```bash
jaws submit --no-cache pipeline_SCATTER_jaws_rRNA.wdl rRNA/inputs_HRP_A_017_native_rRNA_001.json perlmutter --tag "HRP_A_017_native_rRNA_001_run"
```

> [!IMPORTANT]
> The pipeline utilizes Docker containers for all tasks. Because you pulled the images in Step 5, your local execution engine (`miniwdl` or `cromwell`) will automatically detect and use these pre-downloaded images to execute the software steps.

## Step 7: Generating Input JSON for BAM Merging (Replicates)

Now, we will prepare a configuration file that tells the secondary pipeline which replicates to combine. Combining replicates into a single dataset increases our statistical power and confidence when calling modifications.

Once your initial basecalling pipelines have finished successfully, you will have aligned `.bam` files. The next step is to merge multiple biological or technical replicates into a single comprehensive `.bam` file and run `modkit` on it.

Unlike Step 4, there is no automatic bash script for this. You must hand-write an input JSON file using a text editor.

Create a new file called `inputs_merge_modkit_grid_rRNA.json` and paste the following template into it:

```json
{
  "ont_rRNA_merge_modkit_grid.bamfiles": [
    "/path/to/replicate1.merged.transcriptome.aligned.sorted.bam",
    "/path/to/replicate2.merged.transcriptome.aligned.sorted.bam"
  ],
  "ont_rRNA_merge_modkit_grid.bamindices": [
    "/path/to/replicate1.merged.transcriptome.aligned.sorted.bam.bai",
    "/path/to/replicate2.merged.transcriptome.aligned.sorted.bam.bai"
  ],
  "ont_rRNA_merge_modkit_grid.sample_id": "my_merged_rRNA_sample",
  "ont_rRNA_merge_modkit_grid.reference": "/path/to/hs_rRNAs_NR_046235.fa",
  "ont_rRNA_merge_modkit_grid.cpus": 12,
  "ont_rRNA_merge_modkit_grid.mod_thresholds": [
    0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99
  ]
}
```

**Understanding the file:**
- You must replace the paths in `bamfiles` and `bamindices` with the actual file paths produced by Step 6.
- The `sample_id` determines the prefix for your merged output.
- The `mod_thresholds` array specifies the grid of modification probability thresholds `modkit` will evaluate.

## Step 8: Running the BAM Merge and Modkit Pipeline

In this final step, we will launch the automated software that merges your alignments and performs the final modification pileup. Modkit evaluates every nucleotide position across the merged data to determine if a modification exists based on the probability thresholds you set.

With your JSON prepared, you can now launch the second WDL pipeline (`pipeline_rRNA_bam_merge_modkit_SCATTER.wdl`) using the same execution engine you chose in Step 6.

### Option A: Local Execution with miniwdl

```bash
miniwdl run pipeline_rRNA_bam_merge_modkit_SCATTER.wdl --input inputs_merge_modkit_grid_rRNA.json
```

**Understanding the command:**
- `miniwdl run` tells the `miniwdl` engine to start the job.
- `pipeline_rRNA_bam_merge_modkit_SCATTER.wdl` is the secondary blueprint that instructs the system to combine the BAM files and run `modkit`.
- `--input inputs_merge_modkit_grid_rRNA.json` provides the configuration you just wrote.

### Option B: Local Execution with Cromwell

```bash
java -jar cromwell.jar run pipeline_rRNA_bam_merge_modkit_SCATTER.wdl -i inputs_merge_modkit_grid_rRNA.json
```

### Option C: HPC Execution with JAWS (NERSC/Perlmutter)

```bash
jaws submit --no-cache pipeline_rRNA_bam_merge_modkit_SCATTER.wdl inputs_merge_modkit_grid_rRNA.json perlmutter --tag "rRNA_merge_grid_run"
```

---
**Congratulations!** You have successfully downloaded the raw LRS data, basecalled it, and merged your replicates to perform high-quality modification calling using Modkit. Your final pileup files are now ready to be formatted and filtered for downstream figure generation!
