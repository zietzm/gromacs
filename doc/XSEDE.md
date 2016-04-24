# Using XSEDE

* xsede.org is the central hub for all supercomputers
* You need an 'allocation' for any job. In terms of SU ('Service units') - cpu * hour
* Can use WinSCP/nautilus
* Need to use batch jobs -> BASH scripts

Example below:
```
#!/bin/bash
#PBC -l ncpus = 16



```
* Request cpu must be multiple of 16
* Will get output log for each script run.
