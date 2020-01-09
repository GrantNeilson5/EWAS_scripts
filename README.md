How to make a lienar model for EWAS run in parallel on ISCA 

1) write your script and set number of processors you want to parelalize ie cl<-makeCluster(12)
2) when you submitt your job to run make sure you have called the same number of processors in the .sh 
PBS -l procs=12 # Number of processors

3) USe the EWAS submission scrip tto submit the script from the shell command line

4) make sure you submit as a high memory node as otherwise it wont run PBS -l feature=highmem
5) follow the script to run it
