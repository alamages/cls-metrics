# Clustering evaluation metrics

The RI, NMI and conductance metrics are implemented using Cython. Note that conductance is implemented for unweighted and undirected graph. There are example graph and community files under the data/ directory.

Dependencies
--------------------

 * [Cython](http://cython.org/)
 * [numpy](http://www.numpy.org/)

Compile
------------

The code works for both python2 and python3

`git clone https://github.com/alamages/cls-metrics.git`  
`cd cls-metrics`  
`make` 

Data
---------------

The code supports reads undirected unweighted graphs. The input graph must be provided in this format:

[node\_id] [neighbor\_id1] [neighbor\_id2] ....

The node ids must in an increasing sequence starting from index 1 up to n (where n is the number of nodes). Same with the neighbor ids.

Same for the community/cluster files the format is:  

[node\_id] [community\_id]


Examples of input files can be found under the data/ folder.

Examples
------
 
 `# to print the scripts' arguments use -h/--help argument:`  
 `./conductance.py -h`
 
 `# conductance:`  
 `./conductance.py -g data/graph -c data/clusters`

`# ri and nmi metrics:`  
`./ri_nmi.py -g data/community -c data/clusters`