*********** HIGH LEVEL VIEW OF CODE ************

1. In main function we are generating random data points in data array.
2. Then we are passing data points, dimension,and other required values in kdtree function. 
3. In kdtree we are generating mean and variance. Also, we choose choosen_dim.
4. Now, we are passing this data to bipartition function which is going through each data points in cluster.
5. bipartition function is Spliting data in values to choosen_dim.
6. Then, We again comes to kdtree function and calculate centroid of cluster.
7. Then in main function we are passing query_pt which is to be found as testing.
8. After that, search_kdtree function is trying to find the value.
9. If search_kdtree does not found point, we are trying to find minimum distance and nearest point in cluster



********** STEPS TO TEST THE PROGRAM *************

1. Open the program in c IDE
2. According to the requirement, change the values "dim", "ndata", "kk" in the main function.
3. The definition of "dim", "ndata", "kk" is mentioned by the professor.
4. The output is displaying each cluster its starting point and the size of it, distance between the nearest point in the cluster 
   and given point and the nearest point,max variance for each cluster and its dimension.
