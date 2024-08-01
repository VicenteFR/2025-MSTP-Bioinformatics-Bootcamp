# <div align="center"><b>Opening Jupyter Notebooks on TSCC</b></div>


![image.png](../static/Day_1/lights.png)

The command line is a bioinformatician's work horse, but doesn't offer the most interactive or visual interface. [Jupyter](https://jupyter.org/) notebooks are a commonly used alternative which allows you to load in data, manipulate it, make pretty figures, export your final data and figures, and annotate your workflow all in one place!

# Running Jupyter notebooks

Normally, we could use the following command to run a Jupyter notebook:

```bash
jupyter notebook
```

This would generate a Jupyter notebook session that we could access via a simple URL.

However, since we are on a remote server (TSCC), our local machines have no way of knowing where the Jupyter session lives.

Lucky for us, the nice folks over at TSCC have set up a program for us to run Jupyter notebooks on the server and access them on our local machines.

# Using `galyleo` to connect to Jupyter notebooks

`galyleo` is a utility command to help you launch Jupyter notebooks on high-performance computing (HPC) systems. To use `galyleo`, you need to first load the `galyleo` module on TSCC:

```bash
module load galyleo
```

Note that we don't have to install `galyleo` because it is already comes pre-installed on TSCC.

It can then be run in much the same way as a interactive job on TSCC. Here is an example command:

```bash
galyleo launch --account htl191 --cpus 1 --memory 4 --time-limit 1:00:00 --partition hotel --qos hotel
```

Let's break that down:
- `galyleo launch` is the command to start a Jupyter notebook session
- `--account htl191` specifies the allocation to use. Remember this is the allocation we used for running the interactive job for installing software.
- `--cpus 1` specifies the number of CPUs to use. This is the number of CPUs that will be allocated to the Jupyter notebook session. For the purposes of this bootcamp, we will only need 1 CPU.
- `--memory 4` specifies the amount of memory to use. This is the amount of memory that will be allocated to the Jupyter notebook session. For the purposes of this bootcamp, we will only need 4 GB of memory.
- `--time-limit 1:00:00` specifies the time limit for the Jupyter notebook session. 
- `--partition hotel` specifies the partition to use. This is the same partition we used for the interactive job for installing software.
- `--qos hotel` specifies the quality of service to use. This is the same quality of service we used for the interactive job for installing software.

For more details and a full list of possible options: https://github.com/mkandes/galyleo

If successful, you should see something like this:

```bash
Submitted Jupyter launch script to Slurm. Your SLURM_JOB_ID is 1994809.
Success! Token linked to jobid.
Please copy and paste the HTTPS URL provided below into your web browser.
Do not share this URL with others. It is the password to your Jupyter notebook session.
Your Jupyter notebook session will begin once compute resources are allocated to your job by the scheduler.
https://shopping-strife-reexamine.tscc-user-content.sdsc.edu?token=d3f6d542acbfa6dcbccfcb0b9c9bd779
```

As the instructions indicate, copy and paste the URL into your web browser. You should see a Jupyter notebook session open up in your browser.

## Starting a notebook

Galyleo actually loads something called [JupyterLab](https://jupyter.org/try-jupyter/lab/). JupyterLab is a flexible and powerful user interface for programming called an Integrated Development Environment or IDE. You can access new Jupyter notebooks via JupyterLab, open existing ones, and perform many other tasks. Once opened, you should see a GUI interface like this:

![image.png](../static/Day_1/jupyterlab.png)

Most of the time, you will want to start a new notebook. This can be done anytime by clicking on the `+` sign in the top left corner and selecting any of the options under the `Notebook` section.

Here is where our previous installations come in handy. You can now use the Jupyter notebook to run Python code, and you can use the Python packages we installed to analyze your data.

Play around with Jupyter, notebooks, and the environment as a whole. We will work through these together initially.

# Ending a Jupyter notebook session

When you are done with your Jupyter notebook session, you can close the tab in your browser. Note, this only closes your connection to the session and will not end the session. You can always return to this URL to resume the session, as the session will continue to run on the server until the time limit is reached.

# DONE!

---

# PRO  TIP: Uuse vscode
