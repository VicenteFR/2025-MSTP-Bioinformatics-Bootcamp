<div style="border: 2px dashed #6c757d; padding: 10px; border-radius: 10px; background-color: #f8f9fa; text-align: left; margin-bottom: 20px;">
  <p style="font-size: 18px; color: #343a40; font-family: 'Courier New', Courier, monospace;">
    <strong>@author:</strong> Adam Klie
  </p>
</div>

# <div align="center"><b>Logging in to the Triton Shared Computing Cluster (TSCC)</b></div>

Much of bioinformatics involves working with data or algorithms that require more computational resources than available on a personal computer (PC). In this scenario, it is common to use a high-performance computing (HPC) cluster. 

The [Triton Shared Computing Cluster (TSCC)](https://www.sdsc.edu/support/user_guides/tscc.html) is an HPC cluster on campus that is used by many research groups. 


We have set up training accounts for you to use on TSCC. This three step guide will walk you through the steps to log in to TSCC and access your training account.

# 1) Open a terminal

- **For Mac users**, you can access the Terminal using Finder: `Applications - Utilities - Launch Terminal`
- **For Linux users**: `Applications - System tools - Terminal`
- **For Windows users** who followed the instructions in the previous notebook: `Start - Applications - Ubuntu`

# 2) Log in

In order to log in to TSCC, you will need to use two-factor authentication via DUO. See [this guide](https://support.ucsd.edu/services?id=kb_article_view&sys_kb_id=dba41d798776d11c947a0fa8cebb3527&sysparm_article=KB0020168) to get this set-up on one of your devices.

Once DUO is set-up, use `ssh <your_username>@login.tscc.sdsc.edu` to log in to TSCC. Replace `<your_username>` with the username provided to you in the introductory email. For example:

```bash
ssh etrain82@login.tscc.sdsc.edu
```

would log in to the account `etrain82`.

You will be prompted to enter your password. This is your UCSD active directory (AD) password that you use to log in to most things UCSD.

Two-factor authentication will then prompt you to confirm your login. Follow the instructions to confirm.

If you successfully log in, you will see something like:
```bash
[etrain82@login2 ~]$
```
before your prompt. This means you are now logged in to TSCC.

To exit TSCC, type `exit` and press enter. You can alos use `CTRL + D`.

# 3) Granting admin access to your account

It may be useful for TAs to have admin access to your account for troubleshooting purposes. Run the following command to grant that access:

```bash
chmod -R o+rx /tscc/nfs/home/<your_username>
```

Again replacing `<your_username>` with your provided username.
# DONE!
You can now move onto the next step: [3_Installations](3_Installations.md)

---

# Pro tips

Most of what will be taught in this bootcamp will be the basic versions of tools, commands, and concepts. If you have some familiarity with bioinformatics, or even programming in other disciplines, you may be interested on expanding on these basics. I've added a few sections to many of the guides throughout the course labeled "Pro tips" meant to tackle best practices, efficiency hacks and more.

<div style="border: 2px solid #ff211d; padding: 15px; border-radius: 10px; background-color: #ffffff;">
  <p style="color: #000000; font-family: Arial, sans-serif;">
    Important Note: These pro tips are completely optional and NOT NCESSARY for completing the bootcamp tasks. If you are someone who is just getting started with bioinformatics, I would recommend skipping these sections for now and coming back to them later.
  </p>
</div>

Pro tips will appear in outlined boxes that look like this:

<div style="border: 2px solid #4CAF50; padding: 15px; border-radius: 10px; background-color: #e8f5e9;">
  <h2 style="color: #388E3C; font-family: Arial, sans-serif;">
    &#128218; Pro Tip: Download iTerm2 (for macOS users)
  </h2>
  <p style="color: #1B5E20; font-family: Arial, sans-serif;">
    I'm not a huge fan of Mac's built in terminal. I prefer iTerm2. It's not necessary for this course but it's a much better terminal emulator, in both appearance and functionality.<br><br>
    Install iTerm2 here: <a href="https://iterm2.com/" style="color: #1E88E5;">https://iterm2.com/</a><br><br>
    There are many other terminal emulators out there, some of which are may be better than iTerm2. Feel free to explore and find one that you like.
  </p>
</div>


<div style="border: 2px solid #4CAF50; padding: 15px; border-radius: 10px; background-color: #e8f5e9;">
  <h2 style="color: #388E3C; font-family: Arial, sans-serif;">
    &#128218; Pro Tip: creating an SSH configuration for TSCC
  </h2>
  <p style="color: #1B5E20; font-family: Arial, sans-serif;">
    Create a configuration file for logging in to TSCC: <a href="https://www.sdsc.edu/support/user_guides/tscc.html" style="color: #337ab7;">https://www.sdsc.edu/support/user_guides/tscc.html</a>. See the section *Set up multiplexing for TSCC Host*<br>
    This will prevent you from having to DUO authenticate every time you log via the same terminal session and also allow you to use VSCode (see later pro tip).
  </p>
</div>

### In your local pc open or create this file: `~/.ssh/config`

Type the following command exactly as it appears below, then press Enter: `vi ~/.ssh/config`

Copy and paste the following into the terminal. Make sure that the entire thing gets copied and pasted.
```bash
#TSCC Account  

Host tscc

HostName login.tscc.sdsc.edu  

User etrain82

ControlPath ~/.ssh/%r@%h:%p  
        ControlMaster auto  
        ControlPersist 10  
```
Make sure to swap out `etrain82` with your username that I provided to you.

Type `:wq` to save the file and exit.
