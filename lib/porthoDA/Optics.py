# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
The optic algorithm perform a clustering based on the density of the data points
"""

import numpy as np

class Optics:
    """ Optic class to hold the algorithm
    """
    def __init__(self, min_samples, epsilon=float('inf') ):
        """ initialise argument
        
        Parameters
        ----------
        min_samples : int
            the minimal number of data point required in a cluster to make it
            a core distance cluster
        epsilon : float
            default set to infinite but for practical purpose use a floating
            distance
        """
        # neighboring radius
        self.epsilon = epsilon
        # minimum number of point in a neighborhood
        self.min_samples = min_samples

    def run(self, X):
        """ the main function
        
        Parameters
        ----------
        X : numpy array
            the data point to cluster
        """
        self.n = X.shape[0]       
        self.cd = np.ones( self.n  )  * float('inf')       
        self.processed = np.zeros( self.n, dtype=bool )
        self.rd = np.ones( self.n  )  * float('inf')        
        self.unprocessed = range(self.n)
        self.ordered = []
        # while there is unprocessed data point ...
        while self.unprocessed:
            pos = self.unprocessed[0]
            # marking currend as processed and look for neighbors
            
            self.processed[pos] = True
            self.unprocessed.remove(pos)
            self.ordered.append(pos)
            
            point_neighbors = np.where(X[pos] <= self.epsilon)[0]
            # check if current point is a core_distance point :
            # number of neighbor of p = min_samples - 1

            # the distance from a point to its nth neighbor (n = min_cluser_size)
            
            if self.cd[pos] < float('inf') : 
                core_dist = self.cd[pos]
            elif len(point_neighbors) >= self.min_samples - 1:
                sorted_neighbors = sorted([X[pos,m] for m in point_neighbors])
                self.cd[pos] = sorted_neighbors[self.min_samples - 2]
                core_dist = self.cd[pos]
            else :
                core_dist = float('inf')

            if core_dist < float('inf') :
                
                # update reachability_distance for each unprocessed neighbor
                
                seeds = []
                for p in [m for m in point_neighbors if not self.processed[m] ] :
                    # find new reachability distance new_rd
                    # if rd is null, keep new_rd and add n to the seed list
                    # otherwise if new_rd < old rd, update rd
                    new_rd = max( self.cd[pos], X[pos,p] )
                    if self.rd[p] == float('inf') :
                        self.rd[p] = new_rd
                        seeds.append(p)
                    elif new_rd < self.rd[p]:
                        self.rd[p] = new_rd
                # for every unprocessed neighbors
                while(seeds):
                    # find the neighbor n with smallest reachability distance
                    seeds.sort(key=lambda n: self.rd[n])
                    n = seeds.pop(0)
                    if not self.processed[n] :
                        # mark n as processed and find neighbors
                        self.processed[n] = True
                        self.unprocessed.remove(n)
                        self.ordered.append(n)
                        
                        n_neighbors = np.where( X[n] <= self.epsilon )[0]
                        # current neighbor is a core_distance?

                        if self.cd[n] < float('inf') : 
                            core_dist_n = self.cd[n]
                        elif len(n_neighbors) >= self.min_samples - 1:
                            sorted_neighbors = sorted([X[n,m] for m in n_neighbors])
                            self.cd[n] = sorted_neighbors[self.min_samples - 2]
                            core_dist_n = self.cd[n]
                        else :
                            core_dist_n = float('inf')
                            
                        if core_dist_n < float('inf') :
                            # update reachability_distance for each of n's neighbors
                            for p in [m for m in n_neighbors if not self.processed[m] ] :
                                # find new reachability distance new_rd
                                # if rd is null, keep new_rd and add n to the seed list
                                # otherwise if new_rd < old rd, update rd
                                new_rd = max( self.cd[n], X[n,p] )
                                if self.rd[p] == float('inf') :
                                    self.rd[p] = new_rd
                                    seeds.append(p)
                                elif new_rd < self.rd[p]:
                                    self.rd[p] = new_rd
        # when all points have been processed
        # return the ordered list, the reachability distance and the core distance
        return self.ordered, self.rd, self.cd
        
    def cluster(self, epsilon_prime):
        """ after ordering by neighbor, cluster the data point according to 
        a new epsilon distance
        
        Parameter
        ---------
        epsilon_prime : float
            the second epsilon distance , sould be inferior to the first one
            
            
        Return
        ------
        labels : numpy array
            the labeling list containing the cluster label of each data point
        """
        clusterid = 0
        labels = np.ones( len(self.ordered) ) * -1.0        
        separators = []
        i = 0
        while i < len(self.ordered) - 1:
            j = i + 1
            obi = self.ordered[i]
            obj = self.ordered[j]
            rdi = self.rd[obi]
            rdj = self.rd[obj]
            
            # use an upper limit to separate the clusters
            if rdi > epsilon_prime:
                separators.append( obi )
            elif rdj > epsilon_prime:
                separators.append(  obj )
                i += 1
            i += 1
        # don't forget the last one
        if separators[-1] != len(self.ordered) :
            separators.append( len(self.ordered) ) 
        #for start, end in separators :
        clusterid = 0 
        for i in range(len(separators) - 1):
            start = separators[i]
            end = separators[i + 1]
            if end - start > self.min_samples:
                labels[self.ordered[start:end]] = clusterid
                clusterid += 1 
        return labels

    #def extractDBSCAN( self, epsilon_p ) :
        #clusterid = 0
        #labels = np.ones( len(self.ordered) ) * -1.0
        #for i in range(len(self.ordered)) :
            #obi = self.ordered[i]
            #if self.rd[ obi ] > epsilon_p :
                #if self.cd[ obi ] <= epsilon_p :
                    #clusterid += 1 
                    #labels[ obi ] = clusterid 
            #else :
                #labels[ obi ] = clusterid 
        #return  labels
