<div style="border: 2px dashed #6c757d; padding: 10px; border-radius: 10px; background-color: #f8f9fa; text-align: left; margin-bottom: 10px;">
  <p style="font-size: 18px; color: #343a40; font-family: 'Courier New', Courier, monospace;">
    <strong>@author:</strong> James V. Talwar<br>
    <strong>@adapted by:</strong> TJ Sears for 2023 Class<br>
    <strong>@adapted by:</strong> Adam Klie for 2024 Class
  </p>
</div>


# <div align="center"><b>Installations</b></div>

In the practical section of this bootcamp, we will need access to several software programs. Unfortunately, these programs don't come preinstalled on every machine and we need to install them ourselves.

This can be a non-trivial task in bioinformatics, as many programs have intracate and nasty dependencies. Below, we walk you through a common way to install bioinformatic software using a package manager called `miniconda`.

# 1) Log into TSCC:

Follow the instructions as outlined in step 2: [2_TSCC](2_TSCC.md) to log in to TSCC.

# 2) Start an Interactive Session

If you successfully logged into TSCC, you may notice that before your cursor, you may see something like the following:
```bash
[etrain82@login2 ~]$
```
`login2` actually means that that you are currently accessing something called the "login node." It is the machine that every other person who logs in to TSCC is also using. 

<div style="border: 2px solid #ff211d; padding: 15px; border-radius: 10px; background-color: #ffffff;">
  <p style="color: #000000; font-family: Arial, sans-serif;">
    Important Note: It is critical that you **DO NOT** run any computationally intensive code on the login node. Otherwise you might make even simple commands like navigating around files painfully slow for fellow users and you most likely draw the ire of the TSCC support team.
  </p>
</div>

Instead, you will request a different machine and set of resources just for you (aren't you special). This can be done in a few ways, but for the purposes of installing software, we will use an interactive session.

To request resources in something called an "interactive session" (will be covered later), **just for the the purposes of installation**, copy the following command to your terminal and hit enter:

```bash
srun -N 1 -n 2 --mem 2G -t 1:00:00 -p hotel -q hotel -A htl191 --pty bash
```

Let's break down this command
- `srun` - command that is used to ask the cluster for some resources to do *something*. That *something* will be defined at the end of the command.
- `-N 1` - give me 1 node. A node is a computer that has multiple processors and memory. I don't think I've ever asked for more than 1 node at a time.
- `-n 2` - give me 2 processors. Each node has multiple processors and you can ask for more than one if you need them. This is useful if something can be parallelized.
- `--mem 2G` - give me 2 gigabytes of memory. TOTAL amount of memory that you are asking for.
- `-t 1:00:00` - give me 1 hour of time to run. This is in the format of hours:minutes:seconds
- `-p hotel` - give me resources from the `hotel` *partition*. A partition is a group of nodes that have similar properties which follow a pay as you go model. We only have access to the `hotel` partition.
- `-q hotel` - give me *quality of service* from the `hotel` *queue*. We only have access to the `hotel` QOS so don't worry about this for now.
- `-A htl191` - give me resources from this *allocation*. We have a set of amount of resources that are *allocated* to us via this identifier.
- `--pty bash` - give me a bash shell to work in. This is the *something* that we are asking the cluster to do. Give me a playground to run commands in.

If successful you should see something like the following:

```bash
srun: job 1974598 queued and waiting for resources
srun: job 1974598 has been allocated resources
```
with the job number changed accordingly.

You might also note that the text before the prompt changed too. You should see something like:
```bash
bash-5.1$
```
This indicates that you are now on a separate machine with resources dedicated to you and you alone (the 5-1 is actually indicating that we are using bash version 5.1)

# 3) Install Miniconda

Installing packages that are compatible can be a painful process for even seasoned bioinformaticians. Luckily for us, there are so called ["package managers"](https://en.wikipedia.org/wiki/Package_manager) that make it a whole lot easier. [Many package managers exist](https://en.wikipedia.org/wiki/List_of_software_package_management_systems), but we will be using [Miniconda](https://docs.conda.io/en/latest/miniconda.html) for this bootcamp, as its very flexible and lightweight.

Miniconda is pretty easy to install itself. Start by copying the the following file to your home directory:

```bash
scp /tscc/nfs/home/hkcarter/Miniconda3-latest-Linux-x86_64.sh ~/.
```

This file is a bash script (set of code instructions) that will install Miniconda on your account. To run the script, type the following command:
```bash
cd ~
bash Miniconda3-latest-Linux-x86_64.sh
```

Press enter when prompted

<div style="border: 2px solid #4CAF50; padding: 15px; border-radius: 10px; background-color: #e8f5e9;">
  <h2 style="color: #388E3C; font-family: Arial, sans-serif;">
    &#128218; Pro Tip: Organizing external programs
  </h2>
  <p style="color: #1B5E20; font-family: Arial, sans-serif;">
    Organizing code, data, and projects is a critical skill for any bioinformatician. There are numerous ways to organize your files, and you'll stumble across many different structures. The key is to be consistent and ensure that your organization is logical, easy to understand, and well-documented.
    <br><br>
    I prefer to download any external programs or software into a folder called <code>opt</code> in my home directory. You can read more where the <code>opt</code> name came from <a href="https://www.baeldung.com/linux/opt-directory#:~:text=The%20FHS%20defines%20%2Fopt%20as,external%20or%20third%2Dparty%20software">here</a>. When the Miniconda installer prompts you to select an installation directory, specify a path other than your home directory, like <code>~/opt/miniconda3</code>. This will keep your home directory clean and organized.
  </p>
</div>

After hitting enter, you will then be presented with the license/terms and conditions. If you want to skip to the end, hit `q` and then accept the license terms by typing _yes_. 

Miniconda will then present you with an installation location. It should be presented as `tscc/nfs/home/etrain##/miniconda3`. Press enter to confirm the location (or specify your own if you want to do a bit of organization.

Miniconda is now installing! This may take a bit so don't get frustrated. Leave your terminal open and let this run.

![image](https://github.com/user-attachments/assets/4808b7c3-5310-4f4d-b709-b926196e9329)

Time to check if this worked. Type:
    
```bash
conda --version
```

You should see the following output:

```bash
conda 24.5.0
```

# 4. Setting up your "base" environment

Miniconda works by putting downloaded software into containers known as [environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#managing-environments).

This allows you to create different containers/environments that have different purposes.

When you first install Miniconda, you are given a default environment called `base`. This is the environment that you are in when you first login to TSCC.

We will need something called Jupyter notebooks to do analyses for this course (don't worry if you don't know what those are for now). To install Jupyter, run the following command:

```bash
conda install -c conda-forge jupyter jupyterlab
```

<div style="border: 2px solid #4CAF50; padding: 15px; border-radius: 10px; background-color: #e8f5e9;">
  <h2 style="color: #388E3C; font-family: Arial, sans-serif;">
    &#128218; Pro Tip
  </h2>
  <p style="color: #1B5E20; font-family: Arial, sans-serif;">
    Up until recently, conda performance was painfully slow. This has since been remedied to some extent, but I'm not sure if it's been fixed entirely. A much faster alternative to conda is a package manager called <a href="https://mamba.readthedocs.io/">mamba</a> that is a drop-in replacement for conda.
    <br><br>
    <strong>To install mamba on an existing installation of conda, use the following:</strong>
    <br>
    <code>conda install -n base --override-channels -c conda-forge mamba 'python_abi=*=*cp*'</code>
  </p>
</div>


# 5. Creating an environment for running an `rna-seq` analysis

During bootcamp, we will be honing our bioinformatic skills using an RNA-seq analysis. We will get more into the details of this in the first couple days, but for now, we need to install some software that we will need later.

As a rule of thumb I install very little software in my base environment. Not doing so can lead to a bloated `base` environment that can cause performance issues. The one exception to that is Jupyter which will allow us to run notebooks on TSCC.

Instead, we will create a new environment specifically for this bootcamp. Run the following command
```bash
conda create -n 2024-mstp-bootcamp python=3.11 r-base=4.3.1 -y
```

Let's break this down
- `conda create -n 2024-mstp-bootcamp` - create a new environment called `2024-mstp-bootcamp`
- `python=3.11` - install python version 3.11 in this environment
- `r-base=4.3.1` - install R version 4.3.1 in this environment
- `-y` - automatically say yes to any prompts

We can now 'activate' (enter) the environment we just created:
```bash
conda activate 2024-mstp-bootcamp
```

Once in the bootcamp environment, we can install *most* all the necessary packages for RNA-seq analysis:

```bash
conda install -c conda-forge -c bioconda numpy pandas matplotlib seaborn STAR fastqc samtools bzip2 subread -y
```

Let's break this down
- `conda install -c conda-forge -c bioconda` - install packages from the conda-forge and bioconda channels
- `numpy pandas matplotlib seaborn` - install the python packages numpy, pandas, matplotlib, and seaborn
- `STAR fastqc samtools bzip2 subread` - install the programs STAR, fastqc, samtools, bzip2, and subread

Some packages are not available via conda and instead can be installed via the Python package manager [`pip`](https://pip.pypa.io/en/stable/).

Lucky for us, `pip` comes default when a new Python environment is created in conda and conda and pip are very compatible. To install the packages we want, all we have to do is:
```bash
pip install decoupler pydeseq2 scanpy sanbomics gseapy PyWGCNA
```

Great! Hopefully these ran successfully for you. We will talk more about the packages and what they are used for in the actual bootcamp.

There is one last thing we need to do. Jupyter notebooks have no way of knowing where these programs are unless we tell it. First, we need to install something called ipykernel:
```bash
conda install -c anaconda ipykernel -y
```

This let's us create a "kernel" (Jupyter jargon) that knows where the software we just installed lives:
    
```bash
python -m ipykernel install --user --name 2024-mstp-bootcamp --display-name "Python 3.11 R 4.3.1 2024-mstp-bootcamp"
```

One last time for this notebook, let's break this down:
- `python -m ipykernel install` - run the command to install a new kernel
- `--user` - install the kernel for the current user only, as opposed to system-wide
- `--name 2024-mstp-bootcamp` - name the kernel `2024-mstp-bootcamp`, this should match the conda environment name
- `--display-name "Python 3.11 R 4.3.1 2024-mstp-bootcamp"` - display the kernel as "Python 3.11 R 4.3.1 2024-mstp-bootcamp" in Jupyter

# DONE!
Congratulations! You have successfully installed Miniconda, Jupyter, and most of the software you need for the bootcamp. You should be all set to go for Day 1. Feel free to email me with any questions!

---
