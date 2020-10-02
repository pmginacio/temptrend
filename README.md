# temptrend

Utility to compute temperature anomaly trends from the GHCN dataset.

## install with git

python2, with numpy and matplotlib are required.
I could install those easily with,

    pip numpy matplotlib

assuming you have pip installed.

Create a directory of your choice to clone the git repo,

    mkdir ~/programs/temptrend
    cd ~/programs/temptrend
    git clone --recurse-submodules git@github.com:pmginacio/temptrend.git .

Then you can simply run,

    ./temptrend.py -h 
    
for a description of how to use the utility.

## run with docker

To simply run it use,
    
    docker -it pmginacio/temptrend -h

to see a description of the program.

Docker containers are a great way to make sure an application runs the same in all sorts of environments.
This comes at the cost of some complexity when running some applications inside the container.

Temptrend outputs plots of the temperature anomaly trend as PDF files.
However all data in docker containers disappears after the containers stops.
In order to access the plots after the the container exits, we need to give the container a local folder to store the results.
I like to create a local directory to run temptrend, e.g., ~/programs/temptrend

    mkdir ~/programs/temptrend
    cd ~/programs/temptrend

To look at the plots after the container is done, a local directory has to be mounted into the container.
Then the container app will write the output plots there.

    mkdir results
    docker run -it -v ~/programs/temptrend/results:/usr/src/temptrend/results pmginacio/temptrend -h

The first time a given year is requested from the GHCN dataset, temptrend will download it and preprocess it.
The download time is quite long, therefore temptrend keeps a local cache of the GHNC dataset.
The next time this year is requested, the local data is already available.
In order to persist the local cache between calls another mount point is required,

    mkdir data
    docker run -it -v ~/programs/temptrend/data:/usr/src/temptrend/data -v ~/programs/temptrend/results:/usr/src/temptrend/results pmginacio/temptrend -h

Finally, by default docker runs as the root user. 
This means that all files and directories created by the container are owned by root.
This is annoying when we need to delete something. 
Ideally, files and directories should be created as the same user which calls the docker program.
That can be done as,

    docker run -it -v ~/programs/temptrend/data:/usr/src/temptrend/data -v ~/programs/temptrend/results:/usr/src/temptrend/results --user "$(id -u):$(id -g)" pmginacio/temptrend -h

All of this leads to a quite long command. 
One idea is to save it in a executable file and forget about it,

    echo '#!/bin/bash -eu
    [[ ! -d data ]] && mkdir data
    [[ ! -d results ]] && mkdir results
    docker run -it -v $(pwd)/data:/usr/src/temptrend/data -v $(pwd)/results:/usr/src/temptrend/results --user "$(id -u):$(id -g)" pmginacio/temptrend $@' > temptrend
    chmod u+x temptrend

Now run it can be run much shorter as,

    temptrend -t 50 60 2010:2019

## be smart and get to work

Temptrend is simple utility which computes temperature anomaly trends on the basis of the GHCN dataset.
The program has basic functionality and I use it to illustrate some of the design principles of processing large datasets.

The full dataset is a list of temperatures over a geographic grid for a period of ~70 years and counting.
The data is available in text form from their FTP server.
Tthe whole dataset is about 700MB compressed, not a big dataset by any standard.
Nonetheless, it is big enough to present some opportunities for optimization.

Temptrend implements a local cache for the GHCN dataset.
For any given year, execution will be slow only for the first time it is requested, since the data has to be downloaded from the server.
Afterwards, the local cache is used as long as it is present. 
If needed, the local cache can be partially or totally deleted and temptrend will simple download the missing data again.
Local caching minimizes the download speed bottleneck.

The next speed bottleneck is reading and parsing text files, which is very slow.
Once a yearly dataset is downloaded, temptrend loads the daily grids and stores them as compressed binary data.
Compressed binary data can be read extremely fast which removes another speed bottleneck.
It also has small footprint; the whole dataset should take up about 80MB, a factor of 10x smaller than the original compressed text file dataset.

## scale, scale, scale

Quite some thought has been put into temptrend to make it scalable to distributed processing.
Any program written with this in mind is ready to tackle large datasets with big computational loads.
Temptrend makes extensive use map() and reduce().
With this structure temptrend is readily to parallelize computations across multiple worker threads.
These functional operators (map, filter and reduce) also allow for the resource use to be kept to a minimum.
There is no need to keep intermediate results.
For example, one can load 60 years of daily grids as

    lgrds = map(ghcn.load, ldates_1950_to_2010)

However, doing so would likely consume all available RAM in my laptop.
Instead I extend the pipeline to read, consume and discard the daily grids without needed to keep them in memory.

The core of temptrend is the single instruction

    sumA = reduce(np.add, 
        map(grd_to_trend_ls, 
          map(ghcn.load, ldates), itertools.repeat(ldates[len(ldates)/2],len(ldates))))

where the daily grids are read, a LS accumulation matrix is computed and the resulting matrices are accumulated in sumA.
This is done quickly and the amount of memory is low and constant with respect to the length of the dataset considered.

## least squares?

Temptrend has a relatively simple task to fit a linear trend y=mx+b to each of the input grid nodes.
Naive implementations of this would require the whole dataset to be loaded in memory.
Then loop each grid node to compute a least-squares trend estimate at each one.
This would be slow, not scalable and would quickly consume all available memory.
Instead the trend estimation problem has to be reformulated in a distributed manner.
A normal matrix and rhs vector can be built simultaneously for all nodes of a daily grid.
Then the normal matrix and rhs vectors can be accumulated, with a call reduce(), in order to combine all the information at once.
With this approach, temptrend is ready to distribute arbitrary chunks of the dataset for computation to multiple workers possibly in different machines.
Temptrend is ready to accept datasets of any size, as opposed to a naive implementation of LS trend estimation.

## numpy arrays for the win!

Wherever possible, temptrend makes use of vectorized operations instead of looping arrays.
A good example of this is the grids.Grid.pack() function which returns a list of x,y,z values from the 2D grid.
A naive implementation would loop the grid and return x,y,z values.
Instead, temptrend builds and reshapes X and Y coordinate arrays which are then flattened together with the 2D gridded values. 
The final step is to keep only the finite values of the grid which is done with logical indexing,

    X = np.tile(self.x.reshape(self.nx,1),(1,self.ny))
    Y = np.tile(self.y,(self.nx, 1))
    xyz = np.hstack((X.reshape(n,1),Y.reshape(n,1),self.z.reshape(n,1)))
    xyz = xyz[np.isfinite(xyz[:,2]),:]

This approach is one order of magnitude faster than a naive loop implementation.
