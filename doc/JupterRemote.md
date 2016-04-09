# Running Jupyter Notebook remotely
One of the main advantages of using IPython/Jupyter/Rodeo is that you can
pull it up remotely, even on Windows machines. This is possible because
these programs run in the browser. 

On ssh, simply type 
```
jupyter notebook --no-browser --port=7000
``` 
obviously the port number is arbitrary, and, hopefully, so too should be the
program we choose to be available by ssh.

Once that command has been issued by ssh, exit that ssh and locally, type
```
ssh -N -f -L localhost:6000:localhost:7000 USER@REMOTELOCATION
``` 
In the above, 6000 is the local and 7000 is the remote. 

Now to access ipython/jupyter/rodeo, go to a browser and type
```
http://localhost:6000
```
This runs the notebook remotely, and we can now have the notebook on our 
local computer, even though the notebook is actually open remotely.
