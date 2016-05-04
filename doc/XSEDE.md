# Using XSEDE

* xsede.org is the central hub for all supercomputers
* You need an 'allocation' for any job. In terms of SU ('Service units') - cpu * hour
* Can use WinSCP/nautilus
* Need to use batch jobs -> BASH scripts

We have to use Globus for large file transfers.

cd Downloads/globusconnectpersonal-2.2.1/
./globusconnect
(press connect)
go to Globus online
Search endpoint "XSEDE Gordon"
Log in using XSEDE
Press 'Start a transfer from XSEDE Gordon'


Example below:
```
#!/bin/bash
#PBC -l ncpus = 16



```
* Request cpu must be multiple of 16
* Will get output log for each script run.


Environment variables on Gordon

```
[mzietz@gordon-ln1 ~]$ printenv
PDSHROOT=/opt/pdsh
MANPATH=/opt/intel/composer_xe_2013_sp1.2.144/man:/opt/intel/composer_xe_2013.1.117/man:/opt/torque/man:
HOSTNAME=gordon-ln1.sdsc.edu
IPPROOT=/opt/intel/composer_xe_2013.1.117/ipp
INTEL_LICENSE_FILE=/opt/flexlm/license/license-intel.txt
TERM=xterm
SHELL=/bin/bash
HISTSIZE=1000
SSH_CLIENT=165.123.152.247 41958 22
LIBRARY_PATH=/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64:/opt/intel/composer_xe_2013.1.117/compiler/lib/intel64:/opt/intel/composer_xe_2013.1.117/ipp/../compiler/lib/intel64:/opt/intel/composer_xe_2013.1.117/ipp/lib/intel64:/opt/intel/composer_xe_2013.1.117/compiler/lib/intel64:/opt/intel/composer_xe_2013.1.117/mkl/lib/intel64:/opt/intel/composer_xe_2013.1.117/tbb/lib/intel64//cc4.1.0_libc2.4_kernel2.6.16.21:/opt/gridengine/lib/lx26-amd64:/opt/intel/composer_xe_2013.1.117/debugger/lib/intel64:/opt/intel/composer_xe_2013.1.117/mpirt/lib/intel64
FPATH=/opt/intel/composer_xe_2013_sp1.2.144/mkl/include
GMXMAN=/opt/gromacs/share/man
QTDIR=/usr/lib64/qt-3.3
QTINC=/usr/lib64/qt-3.3/include
ROCKSROOT=/opt/rocks/share/devel
SSH_TTY=/dev/pts/10
ANT_HOME=/opt/rocks
USER=mzietz
LD_LIBRARY_PATH=/opt/gromacs/lib:/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64:/opt/gnu/gmp/lib:/opt/gnu/mpc/lib:/opt/gnu/gcc/lib64:/opt/gnu/mpfr/lib:/opt/gnu/lib:/opt/gnu/lib64:/opt/mvapich2/intel/ib/lib:/opt/intel/composer_xe_2013.1.117/compiler/lib/intel64:/opt/intel/composer_xe_2013.1.117/ipp/../compiler/lib/intel64:/opt/intel/composer_xe_2013.1.117/ipp/lib/intel64:/opt/intel/composer_xe_2013.1.117/compiler/lib/intel64:/opt/intel/composer_xe_2013.1.117/mkl/lib/intel64:/opt/intel/composer_xe_2013.1.117/tbb/lib/intel64//cc4.1.0_libc2.4_kernel2.6.16.21:/opt/gridengine/lib/lx26-amd64:/opt/intel/composer_xe_2013.1.117/debugger/lib/intel64:/opt/intel/composer_xe_2013.1.117/mpirt/lib/intel64
LS_COLORS=rs=0:di=01;34:ln=01;36:mh=00:pi=40;33:so=01;35:do=01;35:bd=40;33;01:cd=40;33;01:or=40;31;01:mi=01;05;37;41:su=37;41:sg=30;43:ca=30;41:tw=30;42:ow=34;42:st=37;44:ex=01;32:*.tar=01;31:*.tgz=01;31:*.arj=01;31:*.taz=01;31:*.lzh=01;31:*.lzma=01;31:*.tlz=01;31:*.txz=01;31:*.zip=01;31:*.z=01;31:*.Z=01;31:*.dz=01;31:*.gz=01;31:*.lz=01;31:*.xz=01;31:*.bz2=01;31:*.tbz=01;31:*.tbz2=01;31:*.bz=01;31:*.tz=01;31:*.deb=01;31:*.rpm=01;31:*.jar=01;31:*.rar=01;31:*.ace=01;31:*.zoo=01;31:*.cpio=01;31:*.7z=01;31:*.rz=01;31:*.jpg=01;35:*.jpeg=01;35:*.gif=01;35:*.bmp=01;35:*.pbm=01;35:*.pgm=01;35:*.ppm=01;35:*.tga=01;35:*.xbm=01;35:*.xpm=01;35:*.tif=01;35:*.tiff=01;35:*.png=01;35:*.svg=01;35:*.svgz=01;35:*.mng=01;35:*.pcx=01;35:*.mov=01;35:*.mpg=01;35:*.mpeg=01;35:*.m2v=01;35:*.mkv=01;35:*.ogm=01;35:*.mp4=01;35:*.m4v=01;35:*.mp4v=01;35:*.vob=01;35:*.qt=01;35:*.nuv=01;35:*.wmv=01;35:*.asf=01;35:*.rm=01;35:*.rmvb=01;35:*.flc=01;35:*.avi=01;35:*.fli=01;35:*.flv=01;35:*.gl=01;35:*.dl=01;35:*.xcf=01;35:*.xwd=01;35:*.yuv=01;35:*.cgm=01;35:*.emf=01;35:*.axv=01;35:*.anx=01;35:*.ogv=01;35:*.ogx=01;35:*.aac=01;36:*.au=01;36:*.flac=01;36:*.mid=01;36:*.midi=01;36:*.mka=01;36:*.mp3=01;36:*.mpc=01;36:*.ogg=01;36:*.ra=01;36:*.wav=01;36:*.axa=01;36:*.oga=01;36:*.spx=01;36:*.xspf=01;36:
GROMACSHOME=/opt/gromacs
MKLHOME=/opt/intel/composer_xe_2013_sp1.2.144/mkl
ROCKS_ROOT=/opt/rocks
CPATH=/opt/intel/composer_xe_2013_sp1.2.144/mkl/include
GMXDATA=/opt/gromacs/share
NLSPATH=/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64/locale/%l_%t/%N
MAIL=/var/spool/mail/mzietz
PATH=/opt/gromacs/bin:/opt/gnu/bin:/opt/gnu/gcc/bin:/opt/mvapich2/intel/ib/bin:/opt/intel/composer_xe_2013.1.117/bin/intel64:/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/usr/java/latest/bin:/opt/sdsc/bin:/opt/maui/bin:/opt/torque/bin:/opt/torque/sbin:/opt/torque/bin:/opt/torque/sbin:/state/partition1/catalina/bin:/opt/pdsh/bin:/opt/rocks/bin:/opt/rocks/sbin:/home/mzietz/bin
TBBROOT=/opt/intel/composer_xe_2013.1.117/tbb
PWD=/home/mzietz
_LMFILES_=/opt/modulefiles/compilers/intel/2013.1.117:/opt/modulefiles/mpi/.intel/mvapich2_ib/1.9:/opt/modulefiles/applications/gnubase/1.0:/opt/modulefiles/applications/mkl/11.1:/opt/modulefiles/applications/gromacs/5.0.2
JAVA_HOME=/usr/java/latest
LANG=en_US.UTF-8
MODULEPATH=/opt/modulefiles/mpi/.intel:/opt/modulefiles/applications/.intel:/opt/modulefiles/mpi:/opt/modulefiles/compilers:/opt/modulefiles/applications:/usr/share/Modules/modulefiles:/etc/modulefiles
LOADEDMODULES=intel/2013.1.117:mvapich2_ib/1.9:gnubase/1.0:mkl/11.1:gromacs/5.0.2
LM_LICENSE_FILE=/opt/flexlm/license/license.dat
HISTCONTROL=ignoredups
KRB5CCNAME=FILE:/tmp/krb5cc_513715_CnAY5B
SHLVL=1
HOME=/home/mzietz
ROLLSROOT=/opt/rocks/share/devel/src/roll
MPIHOME=/opt/mvapich2/intel/ib
GMXBIN=/opt/gromacs/bin
LOGNAME=mzietz
QTLIB=/usr/lib64/qt-3.3/lib
CVS_RSH=ssh
GMXLDLIB=/opt/gromacs/lib
SSH_CONNECTION=165.123.152.247 41958 198.202.104.118 22
MODULESHOME=/usr/share/Modules
MKL_ROOT=/opt/intel/composer_xe_2013_sp1.2.144/mkl
LESSOPEN=|/usr/bin/lesspipe.sh %s
INCLUDE=/opt/intel/composer_xe_2013_sp1.2.144/mkl/include
INTELHOME=/opt/intel/composer_xe_2013.1.117
G_BROKEN_FILENAMES=1
BASH_FUNC_module()=() {  eval `/usr/bin/modulecmd bash $*`
}
_=/usr/bin/printenv
```