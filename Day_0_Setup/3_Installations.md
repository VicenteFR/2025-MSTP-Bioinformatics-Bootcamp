<div style="border: 2px dashed #6c757d; padding: 10px; border-radius: 10px; background-color: #f8f9fa; text-align: left; margin-bottom: 10px;">
  <p style="font-size: 18px; color: #343a40; font-family: 'Courier New', Courier, monospace;">
    <strong>@author:</strong> James V. Talwar<br>
    <strong>@adapted by:</strong> TJ Sears for 2023 Class<br>
    <strong>@adapted by:</strong> Adam Klie for 2024 Class<br>
    <strong>@adapted by:</strong> Vicente Fajardo Rosas for 2025 Class
  </p>
</div>


# <div align="center"><b>Installations</b></div>

In the practical section of this bootcamp, we will need access to several software programs. Unfortunately, these programs don't come preinstalled on every machine and we need to install them ourselves.

This can be a non-trivial task in bioinformatics, as many programs have intracate and nasty dependencies. Below, we walk you through a common way to install bioinformatic software using a package manager called Miniconda.

# 1) Log in to TSCC:

Follow the instructions as outlined in step 2: [2_TSCC](2_TSCC.md) to log in to TSCC.

# 2) Start an Interactive Session

When you successfully log in to TSCC, you may notice something like the following before your cursor:
```bash
[etrain82@login2 ~]$
```
`login2` means that that you are currently accessing something called the "login node." It is the machine that every other person who logs in to TSCC is also using. 

<div style="border: 2px solid #ff211d; padding: 15px; border-radius: 10px; background-color: #ffffff;">
  <p style="color: #000000; font-family: Arial, sans-serif;">
    Important Note: It is critical that you **DO NOT** run any computationally intensive code on the login node. Otherwise you might make even simple commands (like navigating around files) painfully slow for fellow users and you most likely draw the ire of the TSCC support team.
  </p>
</div>

Instead, you will request a different machine and set of resources just for you (aren't you special). This can be done in a few ways, but for the purposes of installing software, we will use an "interactive session" (to be covered later).

To request resources for an "interactive session", **just for the the purposes of installation**, copy the following command to your terminal and hit enter:

```bash
srun -N 1 -n 2 --mem 2G -t 1:00:00 -p hotel -q hotel -A htl191 --pty bash
```

If successful you should see something like the following:

```bash
srun: job 1974598 queued and waiting for resources
srun: job 1974598 has been allocated resources
```
with the job number changed accordingly.

You might also note that the text before the prompt changed too. You should see something like:
```bash
[etrain82@tscc-11-2 ~]$
```
This indicates that you are now on a separate machine with resources dedicated to you and you alone.

Let's break down this command
- `srun` - command that is used to ask the cluster for some resources to do *something*. That *something* will be defined at the end of the command.
- `-N 1` - give me 1 node. A node is a computer that has multiple processors and memory. I don't think I've ever asked for more than 1 node at a time.
- `-n 2` - give me 2 processors. Each node has multiple processors and you can ask for more than one if you need them. This is useful if something can be parallelized.
- `--mem 2G` - give me 2 gigabytes of memory. TOTAL amount of memory that you are asking for.
- `-t 2:00:00` - give me 2 hours of time to run. This is in the format of hours:minutes:seconds
- `-p hotel` - give me resources from the `hotel` *partition*. A partition is a group of nodes that have similar properties which follow a pay as you go model. We only have access to the `hotel` partition.
- `-q hotel` - give me *quality of service* from the `hotel` *queue*. We only have access to the `hotel` QOS so don't worry about this for now.
- `-A htl191` - give me resources from this *allocation*. We have a set of amount of resources that are *allocated* to us via this identifier.
- `--pty bash` - give me a bash shell to work in. This is the *something* that we are asking the cluster to do. Give me a playground to run commands in.

Don't worry if you didn't follow a lot of that. We will cover much of this the first day and some of it you don't ever really need to understand.

You are now in your own safe space to perform all the installations your heart desires.

# 3) Install Miniconda

Installing multiple packages that are all compatible with each other can be a painful process for even seasoned bioinformaticians. Luckily for us, there exist package installation aides called ["package managers"](https://en.wikipedia.org/wiki/Package_manager) that make our lives a whole lot easier. [Many package managers exist](https://en.wikipedia.org/wiki/List_of_software_package_management_systems), but we will be using [Miniconda](https://docs.conda.io/en/latest/miniconda.html) for this bootcamp as its very flexible and lightweight.

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

After hitting enter, you will then be presented with the license/terms and conditions. If you want to skip to the end, hit `q` and then accept the license terms by typing _yes_. 

Miniconda will then present you with an installation location as `/tscc/nfs/home/etrain##/miniconda3`. Press enter to confirm the location (or specify your own if you want to do a bit of organization -- see the pro tip at the end of this doc).

Miniconda is now installing! This may take a bit so don't get frustrated. Leave your terminal open and let this run.

![image](https://github.com/user-attachments/assets/4808b7c3-5310-4f4d-b709-b926196e9329)

Finally, Miniconda will ask if you want to run `conda init` to configure your account to automatically use conda on each login. Type __yes__ and hit enter.

Time to check if this worked. Type:
    
```bash
source ~/.bashrc
```

Followed by:

```bash
conda --version
```

You should see the following output:

```bash
conda 24.5.0
```

You should also see your prompt change to something like:
```bash
(base) [etrain82@tscc-11-2 ~]$
```

# 4. Setting up your "base" environment

Miniconda works by putting downloaded software into containers known as [environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#managing-environments). This allows you to create different containers/environments that have different purposes.

When you first install Miniconda, you are put in a default environment called `base`. This is the environment that you are in now and will be whenever you log in to TSCC.

In base, we will need something called Jupyter notebooks for later in the bootcamp (don't worry if you don't know what those are for now). To install Jupyter, run the following command:

```bash
conda install -c conda-forge jupyter jupyterlab -y
```

# 5. Creating an environment for running an `rna-seq` analysis

During bootcamp, we will be honing our bioinformatic skills using an RNA-seq analysis. We will get more into the details in the first couple days, but for now, we need to install some software.

As a rule of thumb I install very little in my base environment (Jupyter being the exception). This is to avoid a bloated `base` environment that can cause performance issues.

Instead, we will create a new environment specifically for this bootcamp. Run the following command
```bash
conda create -n 2025-mstp-bootcamp python=3.11 r-base=4.3.1 -y
```

Let's break this down
- `conda create -n 2025-mstp-bootcamp` - create a new environment called `2025-mstp-bootcamp`
- `python=3.11` - install python version 3.11 in this environment
- `r-base=4.3.1` - install R version 4.3.1 in this environment
- `-y` - automatically say yes to any prompts

We can now 'activate' (enter) the environment we just created:
```bash
conda activate 2025-mstp-bootcamp
```

You should see your prompt change to something like:
```bash
(2025-mstp-bootcamp) [etrain82@tscc-11-2 ~]$
```

This indicates that we are in the bootcamp environment, we can now install *most* of the necessary packages for RNA-seq analysis:

```bash
conda install -c conda-forge -c bioconda numpy pandas matplotlib seaborn STAR fastqc samtools bzip2 subread scanpy gseapy -y
```

Let's break this down
- `conda install -c conda-forge -c bioconda` - install packages from the conda-forge and bioconda channels
- `numpy pandas matplotlib seaborn` - install the python packages numpy, pandas, matplotlib, and seaborn
- `STAR fastqc samtools bzip2 subread` - install the programs STAR, fastqc, samtools, bzip2, and subread

Some packages are not available via conda and instead can be installed via the Python package manager [`pip`](https://pip.pypa.io/en/stable/).

Lucky for us, `pip` comes default when a new Python environment is created in conda, and conda and pip are very compatible. To install the packages we want, all we have to do is:
```bash
pip install decoupler pydeseq2 sanbomics PyWGCNA
```

Great! Hopefully these ran successfully for you. We will talk more about the packages and what they are used for in the actual bootcamp.

There is one last thing we need to do. Jupyter notebooks have no way of knowing where these programs are unless we tell them. We need to install something called ipykernel:
```bash
conda install -c anaconda ipykernel -y
```

and then create a "kernel" (Jupyter jargon) that knows where the software we just installed lives:
    
```bash
python -m ipykernel install --user --name 2025-mstp-bootcamp --display-name "Python 3.11 R 4.3.1 2025-mstp-bootcamp"
```

One last time for this notebook, let's break this down:
- `python -m ipykernel install` - run the command to install a new kernel
- `--user` - install the kernel for the current user only, as opposed to system-wide
- `--name 2025-mstp-bootcamp` - name the kernel `2025-mstp-bootcamp`, this should match the conda environment name
- `--display-name "Python 3.11 R 4.3.1 2025-mstp-bootcamp"` - display the kernel as "Python 3.11 R 4.3.1 2025-mstp-bootcamp" in Jupyter

You can now exit the interactive session by typing `exit` or `CTRL-D`.

# DONE!
Congratulations! You have successfully installed Miniconda, Jupyter, and most of the software you need for the bootcamp. You should be all set to go for Day 1. Feel free to email me with any questions!

---

<div style="border: 2px solid #4CAF50; padding: 15px; border-radius: 10px; background-color: #e8f5e9;">
  <h2 style="color: #388E3C; font-family: Arial, sans-serif;">
    &#128218; Pro Tip 1: Using a screen for your installations (or any other task).
  </h2>
  <p style="color: #1B5E20; font-family: Arial, sans-serif;">
    While installing software, you may face network issues that interrupt your ssh connection and ultimately ruin your installation. Though a bit unlikely, this could happen to you while running any script and not only this installation script.<br>
    To overcome this type of an issue, bash screens come in handy. A screen is like a detachable session window inside your session. It lets you run programs or scripts (like the one we're using for our installations here) that keep running even if you close your laptop or lose your connection.<br>
    To play with screens, we use the command <code>screen</code>. You can check your active screens with <code>screen -ls</code>. When running that for the first time in your session, you will get the following message <code>No Sockets found in /run/screen/S-USER.</code>, which rightfully indicates you don't have any active screens yet. To create a new one, just type <code>screen -S SCREEN.NAME</code> replacing <code>SCREEN.NAME</code> with whatever brief description you'd prefer according to your task (e.g., <code>install</code>). This will take you directly inside the screen (a new, non-connection-dependent session/terminal). In there, you can proceed to run the bash script to install <code>miniconda</code>.<br>
    To exit the screen without killing it and its active processes (i.e., "detaching the screen"), just type <code>ctrl</code>+<code>A</code> and <code>D</code> (D stands for "detach"). You can go back to your screen ("reattach it") by typing <code>screen -r install</code>, replacing <code>install</code> with your screen name (which you can check again, if forgotten, with <code>screen -ls</code>).<br>
    Finally, when you're done with the screen, you can kill it by typing <code>ctrl</code>+<code>D</code>.
  </p>
</div>
<br>
<div style="border: 2px solid #4CAF50; padding: 15px; border-radius: 10px; background-color: #e8f5e9;">
  <h2 style="color: #388E3C; font-family: Arial, sans-serif;">
    &#128218; Pro Tip 2: Organizing external programs
  </h2>
  <p style="color: #1B5E20; font-family: Arial, sans-serif;">
    Organizing code, data, and projects is a critical skill for any bioinformatician. There are numerous ways to organize your files, and you'll stumble across many different structures. The key is to be consistent and ensure that your organization is logical, easy to understand, and well-documented. I prefer to download any external programs or software into a folder called <code>opt</code> in my home directory. You can read more where the <code>opt</code> name came from <a href="https://www.baeldung.com/linux/opt-directory#:~:text=The%20FHS%20defines%20%2Fopt%20as,external%20or%20third%2Dparty%20software">here</a>. When the Miniconda installer prompts you to select an installation directory, specify a path other than your home directory, like <code>/tscc/nfs/home/your_username/opt/miniconda3</code> (be sure to swap out <code>your_username</code> for your actual username). Note, if you already installed miniconda in your home directory, **don't install it again**.
  </p>
</div>
