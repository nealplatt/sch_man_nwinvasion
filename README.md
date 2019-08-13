## ___Schistosoma mansoni___ in the New World

pulication information...

authors...

affiliations...

abstract...

 

---

### Reproducibility
The analyses in the manuscript can be reproduced within a Jupyter Notebook. To begin, complete the following steps in order to replicate the software environment and directory structure.

**Clone the ```git``` repository.**
```
git clone <repository>
cd <repository>
```

**Install Anaconda3**
```
wget https://repo.anaconda.com/archive/Anaconda3-2019.07-Linux-x86_64.sh
bash Anaconda3-2019.07-Linux-x86_64.sh
```

**Install GATK 4.1.2.0**
```
wget -P bin/ https://github.com/broadinstitute/gatk/releases/download/4.1.2.0/gatk-4.1.2.0.zip
unzip bin/gatk-4.1.2.0.zip -d bin/
```
Then modify the included gatk4 yml file as follows:
```
cp bin/gatk-4.1.2.0/gatkcondaenv.yml config/sch_man_nwinvasion-gatk4-env.yml
```
Change the name of the env from ```name: gatk``` to ```name: sch_man_nwinvasion-gatk4```
Update package location from   ```- gatkPythonPackageArchive.zip``` to ```/full/path/to/dir/gatkPythonPackageArchive.zip```

**Install remaining software via ```conda```, build envs, and activate**
Create envs.

```
conda env create -f config/sch_man_nwinvasion-pipeline-env.yml
conda env create -f config/sch_man_nwinvasion-gatk4-env.yml
conda env create -f config/sch_man_nwinvasion-nbanalyses-env.yml
```

Now activate the main analyses environment.
```
conda activate sch_man_nwinvasion-nbanalyses
```


**Start the Jupyter Notebook server**
```
jupyter notebook 
```

Load the notebook.  From there analyses, should be completed serially.

---

In the notebook make sure all of the relevant info is present:
- dowloand data from SRA, exomes, and Sanger SMV7
 

