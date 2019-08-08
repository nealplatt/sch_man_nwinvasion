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

**Install the remaining software via ```conda``` and activate the environment**
```
conda env create -f config/env.yml
conda activate sch_man_nwinvasion
```

**Start the Jupyter Notebook server**
```
jupyter notebook 
```

Then load the notebook.  From there analyses, should be completed serially.


