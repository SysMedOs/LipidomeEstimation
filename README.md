# LipidomeEstimation

This repository contains two python scripts for the estimation of lipidome in the unox/ox form.

Tested in Python 2.7, 3.6 and 3.7 with `pandas` and `scipy` on Windows an Linux

### How to run the prediction

Navigate to the program folder

+ For the prediction of **unox lipidome**
```
$ python unoxlipidome_estimation.py
```

+ For the prediction of **ox lipidome**
```
$ python oxlipidome_estimation.py
```

The output will be displayed in the terminal.


### Default values

The default FA list is in the `data/FA_list.xlsx`


The default lipid classes considered is defined as below:

```python
usr_lipid_classes = {
    # list of lipid classes with 1 FA
    'x1': ['FA', 'CholesterolEster', 'LPA', 'LPC', 'LPE', 'LPG', 'LPI', 'LPS',
           'Monoacylglycerol', 'Ceramide', 'Sphingolipid'],
    # list of lipid classes with 2 FA
    'x2': ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'Diacylglycerol'],
    'x3': ['Triacylglycerol'],  # list of lipid classes with 3 FA
    'x4': ['Cardiolipin'],  # list of lipid classes with 4 FA
}
``` 

The default lipid classes considered is defined as below:
```python
usr_mod_dct = {
    'm_ocp': ['Aldehyde', 'CarboxylicAcid'],  # list of cleavage terminal
    'm_oap': ['OH', 'OOH', 'KETO', 'EPOXY'],  # list of Oxygen addition modifications
    'm_p': ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'],  # list of Prostane Rings
    'm_o': ['D-IsoK', 'E-IsoK', 'TXA', 'TXB'],  # list of Other types
}
``` 

The `usr_lipid_classes` and `usr_mod_dct` can be changed in the `if __name__ == '__main__':` section.


### License

+ This package is Dual-licensed

    * For academic and non-commercial use: [GPLv2 License](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html):
    
    * For commercial use: please contact the develop team by email.

+ Please cite our publication in an appropriate form.
    * Publication submitted

### Fundings
We acknowledge all projects that supports the development of LipidHunter:

+ BMBF - Federal Ministry of Education and Research Germany:

    https://www.bmbf.de/en/

+ e:Med Systems Medicine Network:

    http://www.sys-med.de/en/

+ SysMedOS Project :

    https://home.uni-leipzig.de/fedorova/sysmedos/
