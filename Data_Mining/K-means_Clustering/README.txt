*********** HIGH LEVEL VIEW OF CODE ************

1. In main function we are generating random data points in data array.
2. In main function we are setting first two centroid in "centroids" array. First centroid is random and second is farthest from first.
3. Then we are passing data points, dimension,and other required values in kmeans function. 
3. In kkmeans function we are forming the clusters and their centroids
4. Now, we are passing this data to bipartition function which is going through each data points in cluster.
5. While loop is running until it finds same clusters and centroids.
6. Then in main function we are passing query points which is to be found as testing.
8. After that, search function is trying to find the nearst clusters



********** STEPS TO TEST THE PROGRAM *************

1. Open the main file in c IDE
2. According to the requirement, change the values "dim", "ndata", "total_clusters" in the main function.
3. The definition of "dim", "ndata", "total_clusters" is mentioned by the professor. ("total_clusters" is same as "k")
4. The output is displaying each cluster its starting point and the size of it, the nearest point.
