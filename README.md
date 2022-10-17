# Team Tec-Monterrey 2022 Software Tool 🧬 

## Description  📝

# rnatrix
Find optimal sRNA binding sequence for a given sRNA-mRNA pair and predict its regulatory behavior.
Binding sequence based on Vazquez-Anderson et al.[1] and Kumar et al. [2].
Prediction achieved through a neural network model trained on our dataset [3]. 

https://2022.igem.wiki/tec-monterrey/software

## Installation ⚙️ 

### Mac OS X. Apple Silicon macs (M1, M2…) 🍎  


**I. INSTALL CONDA AND CONFIGURE BIOCONDA** 📓

All of the used python libraries are installed through Bioconda; a package manager with focus in bioinformatics and 
computational biology.

1. **Install Miniconda**

Miniconda is a lightweight version of conda, it contains only conda, python, and a small number of essential packages. 
Go to: https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.pkg

Once downloaded, double click on the .pkg file and proceed according to the instructions
It is important to note that you shouldn’t use the Intel x86 version of miniconda, since we had trouble installing 
the needed packages.

2. **Configure Bioconda**

Open your terminal and execute, line by line the following commands

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config —-add channels conda-forge

❗️**NOTE**: If you go to https://bioconda.github.io/, you’ll notice that you should execute a fourth command:

    conda config --set channel_priority strict

However, we recommend you to not execute it, since it caused problems when we were trying to install all the packages

3. **Create a conda environment**

Open your terminal, deactivate the base conda environment 

    conda deactivate

Create a new environment with Python 3.9, but specifying to conda that you wish to use the x86 architecture rather 
than the native ARM64 one

    CONDA_SUBDIR=osx-64 conda create -n rnatrix python=3.9

❗️**NOTE**: You probably noticed that the miniconda version we’re using is ARM64, but we need to specify that 
the architecture of our environment should be x86 since some packages, like ViennaRNA are not available for 
ARM64 architectures. Also, we highly recommend you to install Python 3.9, since Python 3.10 is not compatible 
with our software.

Next, activate the created environment and specify the architecture

    conda activate rnatrix
    conda config --env --set subdir osx-64

**II. INSTALL ALL PACKAGES** 📚

1. Update pip. This is necessary to perform a successful installation of Nupack package, as specified below

        python3 -m pip install -U pip

2. Install Jupyter lab and matplotlib

        conda install --update-all pip matplotlib jupyterlab

3. Install the RNA packages

        conda install -c bioconda viennarna
        conda install -c bioconda barriers

4. Download nupack files from: https://github.com/beliveau-lab/NUPACK

❗️**NOTE**: If you go to: https://docs.nupack.org/start/#maclinux-installation, 
you’ll see that according to the official manual, the standard way to install Nupack is registering on 
the official website and then downloading the 4.0.0.28 (and further) versions, however since version 4.0.0.27 
Nupack is made to run natively on ARM64 architecture, and this won’t work because we’ve specified previously 
to conda that our preferred architecture will be x86

5. Depending upon the language of our operative system, we select the downloads folder in the terminal

        cd Downloads

❗️**NOTE**: If you’re already on the Downloads folder, please skip step 5

6. Install Nupack using pip

        python3 -m pip install -U nupack -f ~/<Downloads>/nupack-<VERSION>/package

❗️**NOTE**: Replace <Downloads> by the real name of our downloads folder 
(or the direction in which Nupack install files are located, and <VERSION> by the numerical version 
of our chosen Nupack package (by default is 4.0.0.23)

7. Verify that the Nupack installation was successful

        python3 -m pip install -U pytest
        python3 -m pytest -v –pyargs nupack

In our therminal should appear a similar message:

        ================== 41 passed in 14.07s ==================


### INTEL MACS (x86 ARCHITECTURE)🍎


**I. INSTALL CONDA AND CONFIGURE BIOCONDA** 📓

All of the used python libraries are installed through Bioconda; a package manager 
with focus in bioinformatics and computational biology

1. **Install Miniconda**

    Miniconda is a lightweight version of conda, it contains only conda, python, and a small number of essential 
    packages. Go to: https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.pkg

    Once downloaded, double click on the .pkg file and proceed according to the instructions
    It is important to note that you **must** use the Intel x86 version of miniconda

2. **Configure Bioconda**

    Open your terminal and execute, line by line the following commands
    
        conda config --add channels defaults
        conda config --add channels bioconda
        conda config —-add channels conda-forge

    **NOTE**: If you go to https://bioconda.github.io/, you’ll notice that you should execute a fourth command:
    
        conda config --set channel_priority strict

    However, we recommend you to not execute it, since it caused problems when we were trying to install all 
    the packages


3. **Create a conda environment**

    Open your terminal, deactivate the base conda environment 
    
        conda deactivate

    Create a new environment with Python 3.9, but specifying to conda that you wish to use 
    the x86 architecture rather than the ARM64 one
    
        CONDA_SUBDIR=osx-64 conda create -n rnatrix python=3.9

    ❗️**NOTE**: You probably noticed that the miniconda version we’re using is x86, but we need to specify 
    that the architecture of our environment should be x86 since some packages, like Nupack can cause 
    trouble if the environment architecture is not previously specified. Also, we highly recommend you 
    to install Python 3.9, since Python 3.10 is not compatible with our software

    Next, activate the created environment and specify the architecture
    conda activate rnatrix
    
        conda config --env --set subdir osx-64

**II. INSTALL ALL PACKAGES** 📚

1.  Update pip. This is necessary to perform a successful installation of Nupack package, as specified below

        python3 -m pip install -U pip

2.  Install Jupyter lab and matplotlib

        conda install --update-all pip matplotlib jupyterlab

3.  Install the RNA packages

        conda install -c bioconda viennarna
        conda install -c bioconda barriers

4.  Download nupack files from: https://github.com/beliveau-lab/NUPACK

    ❗️**NOTE**: If you go to: https://docs.nupack.org/start/#maclinux-installation, you’ll see that according 
    to the official manual, the standard way to install nupack is registering on the official website and 
    then downloading the 4.0.0.28 (and further) versions, however since version 4.0.0.27 nupack is made to run natively on ARM64 architecture, and this won’t work because we’ve specified previously to conda that our preferred architecture will be x86

5.  Depending upon the language of our operative system, we select the downloads folder in the terminal

        cd Downloads

    ❗️**NOTE**: If you’re already on the Downloads folder, please skip step 5 

6.  Install Nupack using pip

        python3 -m pip install -U nupack -f ~/<Downloads>/nupack-<VERSION>/package

    ❗️**NOTE**: Replace <Downloads> by the real name of our downloads folder (or the direction in which Nupack 
    install files are located, and <VERSION> by the numerical version of our chosen Nupack package 
    (by default is 4.0.0.23)

7.  Verify that the Nupack installation was successful

        python3 -m pip install -U pytest
        python3 -m pytest -v –pyargs nupack
   
    In our therminal should appear a similar message:
    
        ================== 41 passed in 14.07s ==================



### UBUNTU 20.04+ 🐧

This one is the simplest method, since bioconda packages can run natively 
on ubuntu with x86 architecture.


**I. INSTALL CONDA AND CONFIGURE BIOCONDA**📓

1.  **Install Miniconda**

    Miniconda is a lightweight version of conda, it contains only conda, python, 
    and a small number of essential packages. 
    Go to: https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
    and download the installer

2. **Open your terminal and go to Downloads folder**

       cd Downloads

3.  **Run the installcation command**

         bash Miniconda3-latest-Linux-x86_64.shdf

    Follow the prompts on the in staller screens

4.  **Configure Bioconda**

    Open your terminal and execute, line by line the following commands
        
        conda config --add channels defaults
        conda config --add channels bioconda
        conda config —-add channels conda-forge

    ❗️**NOTE**: If you go to https://bioconda.github.io/, 
    you’ll notice that you should execute a fourth command:
    
        conda config --set channel_priority strict

    However, we recommend you to not execute it, 
    since it caused problems when we were trying to install all the packages

5.  **Create a conda environment**

    Open your terminal, deactivate the base conda environment 
    
        conda deactivate

    Create a new environment with Python 3.9
    
        conda create -n rnatrix python=3.9

    ❗️**NOTE**: We highly recommend you to install Python 3.9, 
    since Python 3.10 is not compatible with our software

    Next, activate the created environment and specify the architecture
    conda activate rnatrix


**II. INSTALL ALL PACKAGES**📚

1.  Update pip. This is necessary to perform a successful installation of Nupack package, 
    as specified below

        python3 -m pip install -U pip

2.  Install Jupyter lab and matplotlib

        conda install --update-all pip matplotlib jupyterlab

3.  Install the RNA packages

        conda install -c bioconda viennarna
        conda install -c bioconda barriers

4.  Register on Nupack website http://www.nupack.org/downloads/register. 
    You’ll receive by email an username and a password for further Nupack downloading

5.  Download Nupack installation files from: http://www.nupack.org/downloads/source.
    You’ll be asked to provide your account information

    ❗️**NOTE**: Download the latest version of Nupack (4+)

6.  Depending upon the language of our operative system, we select the downloads folder in the terminal

        cd Downloads

    ❗️**NOTE**: If you’re already on the Downloads folder, please skip step 6

7.  Install Nupack using pip

        python3 -m pip install -U nupack -f ~/<Downloads>/nupack-<VERSION>/package

    ❗️**NOTE**: Replace <Downloads> by the real name of our downloads folder 
    (or the direction in which Nupack install files are located), and <VERSION> 
    by the numerical version of our chosen Nupack package (by default is 4.0.0.28)

8.  Verify that the Nupack installation was successful

        python3 -m pip install -U pytest
        python3 -m pytest -v –pyargs nupack

    In our therminal should appear a similar message:
    
         ================== 41 passed in 14.07s ==================



### WINDOWS 10 AND 11 🪟 

Bioconda only supports 64-bit Linux and MacOS operating systems; however in Windows 10 and Windows 11, 
we are able to run Linux packages and scripts from command line using Windows Subsystem for Linux 2

**I. Install Windows Subsystem for Linux**📕

1.  **Check Windows version**
    You must be running Windows 10 version 2004 (build 19041 and higher) or Windows 11

2.  **Install Ubuntu 20.04 from windows store**
    Go to https://apps.microsoft.com/store/detail/ubuntu-20045-lts/9MTTCL66CPXJ and install Ubuntu 20.04.5 
    This is not an operating system, but is a terminal able to run Linux commands.

    ❗️**NOTE**: Probably are newer versions of Ubuntu in the windows store, however we have not tried them and 
    we recommend the 20.04.5 version, because we have already tested it successfully.

3.  Install Windows Subsystem for Linux

    Run Windows Powershell as administrator and execute
    
         wsl --install

    ❗️**NOTE**: All of the following steps must be completed using the installed Ubuntu terminal


**II. Install Miniconda and Configure Bioconda**📓

All of the used python libraries are installed through Bioconda; 
a package manager with focus in bioinformatics and computational biology

1. Install Miniconda

    Miniconda is a lightweight version of conda, it contains only conda, python, and a small number of essential packages. 
    Go to: https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
    and download the installer

    ❗️**NOTE**: Don’t download the Windows installer, remember we’re using WSL for executing Linux packages in Windows

2.  Open your Ubuntu terminal and excecute

        explorer.exe .

    This command will open up the system folder in which your Linux files will be stored

    ❗️**NOTE**: Don´t forget the dot “.” after explorer.exe

3.  Open your standard Windows Downloads folder. Move the downloaded Miniconda installer to your Linux 

        \\wsl.localhost\Ubuntu<version>\home\<user>

    You should move the installer to this directory. Replace the <version> with your ubuntu version and <user> with your user name

4.  Run the installation command

        bash Miniconda3-latest-Linux-x86_64.sh

    Follow the prompts on the installer screens
    ❗️**NOTE**: Run for file name

5.  Configure Bioconda
    Open your terminal and execute, line by line the following commands

        conda config --add channels defaults
        conda config --add channels bioconda
        conda config --add channels conda-forge

    ❗️**NOTE**: If you go to https://bioconda.github.io/, you’ll notice that you should execute a fourth command:
    conda config --set channel_priority strict

    However, we recommend you to not execute it, since it caused problems when we were trying to install all the packages

6.  Create a conda environment

    Open your terminal, deactivate the base conda environment 
    
        conda deactivate

    Create a new environment with Python 3.9
    
        conda create -n rnatrix python=3.9

    ❗️**NOTE**: We highly recommend you to install Python 3.9, since Python 3.10 is not compatible with our software

    Next, activate the created environment and specify the architecture
    conda activate rnatrix

**III. Install all packages**📚

7.  Update pip. This is necessary to perform a successful installation of Nupack package, as specified below

        python3 -m pip install -U pip

8.  Install Jupyter lab and matplotlib

        conda install --update-all pip matplotlib jupyterlab

9.  Install the RNA packages

        conda install -c bioconda viennarna
        conda install -c bioconda barriers

10. Register on Nupack website http://www.nupack.org/downloads/register. 
    You’ll receive by email an username and a password for further Nupack downloading

11. Download Nupack installation files from: http://www.nupack.org/downloads/source. 
    You’ll be asked to provide your account information

    ❗️**NOTE**: Download the latest version of Nupack (4+)

12. Depending upon the language of our operative system, we select the downloads folder in the terminal
    Move the nupack-<version> file to your user path. !Be sure to extract the zipfile and dont create another file!

        \\wsl.localhost\Ubuntu<version>\home\<user>

13. Install Nupack using pip

        python3 -m pip install -U nupack -f ~/<Downloads>/nupack-<VERSION>/package

    ❗️**NOTE**: Replace <Downloads> by the real name of our downloads folder 
    (or the direction in which Nupack install files are located), 
    and <VERSION> by the numerical version of our chosen Nupack package (by default is 4.0.0.28)

14. Verify that the Nupack installation was successful
    
        python3 -m pip install -U pytest
        python3 -m pytest -v –pyargs nupack

    In our therminal should appear a similar message:
    
        ================== 41 passed in 14.07s ==================

15. Install required packages and dependencies

        python3 -m pip install -r requirements.txt

## Download

Use conda to clone the rnatrix repository 

        git clone https://gitlab.igem.org/2022/software-tools/tec-monterrey.git

This will create a folder in your home/<user> named tec-monterrey

## Usage 📃
There are a few ways to use the program: 

    - sRNA design from mRNA 
    - Prediction of regulatory behavior from designed sRNA
    - Regulatory prediction from sRNA-mRNA pair

1. Open Anaconda and activate rnatrix environment
        
        conda activate rnatrix

2.  Change your working directory to 'tec-monterrey' folder using 'cd' command

    It should look something like this:
    
        (rnatrix)<user>@<PC>: cd ~/tec-monterrey

3. Inside 'tec-monterrey/inputs' you will find two subfolders 'sRNA' and 'mRNA'.
    The program will look for either a single text file inside 'mRNA' or a single text file inside both 'sRNA' and 'mRNA'.
    1. In case you only have a file in your 'mRNA' subfolder, the program will both design and predict its behavior.
    2. If you were to also have a file inside 'sRNA' subfolder, the program would then only predict the behavior of their interaction.
    
4. The program expects as a first argument the length of the desired sRNA to be designed, and as a second argument an energy threshold (kcal/mol).

**!!! Make sure to remove .gitkeep from both 'input' subfolders !!!**
Example:

    python rnatrix.py 24 5

5. Several flags are available to tweak the programs behavior

    python rnatrix.py --help

## Authors and acknowledgment 📨
**Code:** José Manuel Hernández de Labra, Hugo Fidel Campos Espinoza, Juan Alberto Alfar Almaguer

**Research:** Juan Alberto Alfar Almaguer, José Manuel Hernández de Labra, Hugo Fidel Campos Espinoza

**Data collection:** José Manuel Hernández de Labra, Hugo Fidel Campos Espinoza, Juan Alberto Alfar Almaguer, Juan Manuel Domínguez Larrieta, Andrea Zarely Benítez Obregón, Alejandra Velázquez Mañas, Samyr Humberto Nacif López, Sofia Elena Goitia Favela, Emilio Fabian Ortiz, Ana Paola Morales Mendoza, Ana Paula Rodríguez Cavazos, Andrea González Hidalgo
