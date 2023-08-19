*********** HIGH LEVEL VIEW OF CODE ************

1. In main function we are generating random data points in data array.
2. Then we are passing data points, dimension,and other required values in LSH function. 
3. In LSH we create clusters according to the hash.
4. Hash for each point is created according to the formula (x1.h1 - h1.c1)/w
5. The data is again sorted and arranged according to clusters created.
6. Then in search_LSH function we are passing query_pt which is to be found as testing.
8. After that, search_LSH function is trying to find the value.
9. If search_LSH does not found point, we are trying to find nearest point in cluster




********** STEPS TO TEST THE PROGRAM *************

1. Open the program in c IDE
2. According to the requirement, change the values "dim", "ndata", "m" in the main function.
3. The definition of "dim", "ndata", "m" is mentioned by the professor.
4. The output is displaying hashs, cluster, start index of cluster and nearest distance from query. Also total clusters.
