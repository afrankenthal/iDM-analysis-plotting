# iDM-analysis-plotting
iDM analysis code (plotting and cutflows)

### Run jupyter on the LPC

It is possible to run available versions of jupyter on the LPC. This is in lieu of having to install jupyter on your local machine. The idea is to run the jupyter notebook server on the LPC and connect to it via your local computer browser. For that to work you need to SSH-tunnel the notebook server's port.

To set this up:

1) SSH tunnel. SSH to LPC machines and add port tunneling:

    ```shell
    $ ssh -L 8888:localhost:8888 user@cmslpc-sl6.fnal.gov
    ```

   8888 is the default port for jupyter notebooks.
2) Jupyter notebook server. Once inside an LPC machine, there are two options for enabling jupyter notebooks:

  - Python2: this is the easiest way and comes bundled with CMSSW, but it only offers Python2 support. After running `cmsenv` inside a    CMSSW release, type: `jupyter notebook --no-browser`
  
  - Python3: Python3 has a lot of nifty features that are worth using, but it doesn't come with CMSSW except for the very latest releases (10.1.X I believe). To enable it in the LPC (note this is outside CMSSW):
  
    ```shell
    source /cvmfs/sft.cern.ch/lcg/views/LCG_92python3/x86_64-slc6-gcc62-opt/setup.sh
    export PYTHONPATH=/cvmfs/sft.cern.ch/lcg/views/LCG_92python3/x86_64-slc6-gcc62-opt/lib/python3.6/site-packages:$PYTHONPATH
    ```
  
     Release LCG_92 already comes with jupyter too, so after sourcing it you can just type `jupyter notebook --no-browser` to run the server.
     
     **NOTE:** If you choose Python3 and then afterwards set up a CMSSW environment, it will mess with your jupyter configuration. If you need to use both at the same time, make sure CMSSW is set up _before_ the LCG_92 release. 
     
  
  
3) Access the notebook server on your browser. After the notebook server is set-up, it will give you a link to open (in the form `http://localhost:8888...`). Copy that link and paste it on your browser and you'll enter the jupyter notebook environment, and you're ready to go.


### Scripts
BeamHalo: plots of track quality information and table of beam halo summary
GenKinematicsNew: Studies of the Gen information. 
PlotSignalBkgs: Plots the cutflow for the signal and background from data files containing the summary from SROptimizationAnalysisFull
QCDCorrelationStudies:
SROptimizationAnalysisFull: Loads all the backgrounds scripts and dumps the information into data files to be read by PlotSignalBkgs
