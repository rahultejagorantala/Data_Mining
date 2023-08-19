#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//This function compares to arrays to check if they have all the values same.
int compareArray(double *a,double *b, int size)	{

    for(int i=0;i<size;i++){
        if(a[i]!=b[i]){
             return 1;
        }
    }
    return -1;
}

//This routine creates clusters according to the hash. hash for each point is created according to the formula (x1.h1 - h1.c1)/w.
//The data is again sorted and arranged according to clusters created.
int LSH(int dim, int ndata, double *data,
        int m, double W, double h[][m], double *b,
        int *cluster_start, int *cluster_size, double cluster_hashval[][m], double hyper[]){
    // initializing the variables.
    int hashpoint = 0;
    int index;
    int matchcounter;
    double point1[m];
    int max_clus_val = 10000;
    double cluster_hashval_dummy[max_clus_val][m];
    double centroid[dim];
    double point[dim];
    point[0] = -1;
    int nclusters = 0;
    double *data_dummy = malloc(ndata*dim*sizeof(double));
    for(int k=0; k<max_clus_val; k++) {
        for (int l = 0; l < m; l++) {
            cluster_hashval_dummy[k][l] = -99;
        }
    }
    for(int i=0;i < dim; i++){
        centroid[i] = 0;
    }

    for (int i = 0; i < ndata*dim; i++){
        centroid[i%dim]+=data[i];
    }
    //calculating the centroid.
    for (int i = 0; i< dim; i++){
        centroid[i] = centroid[i]/ndata;
    }
    // finding the xi.hi for each data point
    for (int i = 0; i <= ndata*dim; i++){
        index = 0;
        if (point[0] != -1 && i%dim == 0 ){

            if (i< dim+1){
                hashpoint=0;
            }
            else{
                hashpoint+=1;
            }
            for (int j = 0; j < m*dim; j++){
                if (j%dim == 0 && j>0){
                    index +=1;
                }
                h[hashpoint][index]+= (point[j%dim]*hyper[j]);
            }
        }
        point[i%dim] = data[i];
    }
    //calculating the bias.
    point[0] = -1;
    index= 0 ;
    for(int i = 0;i<= m*dim; i++){
        if (i%dim == 0 && point[0] != -1){
            if (i > dim+1){
                index+=1;
            }
            for (int j=0; j<dim; j++){
                b[index]+= (point[j]*centroid[j]);
            }
        }
        point[i%dim] = hyper[i];
    }
    //finding the hash for each point and creating a unique cluster_hashval array which has all the cluster hash values.
    for(int i=0; i< ndata; i++){
        for (int j=0; j < m; j++){
            h[i][j] = floor((h[i][j] - b[j])/W);
            point1[j] = h[i][j];
        }
        for(int k=0; k<1000; k++){
            if (cluster_hashval_dummy[k][0] == -99){
                for(int l=0; l<m;l++){
                    cluster_hashval_dummy[k][l] = point1[l];
                }
                nclusters = k;
                break;
            }
            matchcounter = 0;
            for(int l=0; l<m; l++){
                if (point1[l] == cluster_hashval_dummy[k][l]){
                    matchcounter+=1;
                }
            }
            if (matchcounter == m){
                break;
            }
        }

    }
    for(int k=0; k<=nclusters; k++) {
        for (int l = 0; l < m; l++) {
            cluster_hashval[k][l] = cluster_hashval_dummy[k][l];
        }
    }
    //sorting the data according to the clusters found and assigning back to the same pointer.
    int ind = -1;
    int counter;
    int *cluster_start1 = malloc(nclusters*sizeof(double));
    int *cluster_size1  = malloc(nclusters*sizeof(double));
    for(int k =0;k<=nclusters; k++){
        counter = 0;
        for(int i=0; i< ndata; i++){

            if (compareArray(h[i],cluster_hashval[k], m) == -1){
                for(int j =i*dim; j < (i*dim)+ dim;j++){
                    counter++;
                    ind++;
                    if (counter ==1){
                        cluster_start1[k] = ind;
                    }
                    data_dummy[ind] = data[j];
                }
            }
        }
        cluster_size1[k] = counter;
    }
    data = data_dummy;
    for (int i = 0; i<= nclusters; i++){
        cluster_start[i] = cluster_start1[i];
        cluster_size[i] = cluster_size1[i];
    }

    return nclusters;
}


//This routine finds the hash for the given point and uses the start indexes and sizes to determine the matching hash and nearest point in the cluster.
int search_LSH(int dim, int ndata, double *data,
               int m, double W, double h[][m], double *b,
               int nclusters, int *cluster_start, int *cluster_size, double cluster_hashval[][m],
               double *query_pt, double *result_pt, double hyper[]) {
    //Initializing the variables.
    double hash[m];
    int index = 0;
    double mindist = INT_MAX;
    double dist = 0;
    double point[dim];
    //Finding hash for the given point.
    for (int j=0; j < m*dim; j++){
        if (j%dim == 0 && j>0){
            index +=1;
        }
        hash[index]+= query_pt[j%dim]*hyper[j];
    }
    printf("HASH for the given point is\n");
    for (int j = 0; j<m; j++){
        hash[j] = floor((hash[j] - b[j])/W);
        printf("%f\n", hash[j]);
    }
//    for(int k=0; k<=nclusters; k++) {
//        for (int l = 0; l < m; l++) {
//            //printf(" inside LSH Search clusterhashval %lf %d %d\n", cluster_hashval[k][l], k,l);
//            // print this to see the unique cluster hash values.
//        }
//    }
    //Finding the matching cluster and nearest point in the cluster to the given point.
    for(int k =0;k<=nclusters; k++) {

        if (compareArray(hash,cluster_hashval[k], m) == -1){
            printf("the given query point belongs to cluster %d and start index is %d \n", k+1, cluster_start[k]);
            for (int i = cluster_start[k]; i <= cluster_start[k] + cluster_size[k]; i++){

                if(i%dim == 0 && i > cluster_start[k]){
                    //printf(" compare %f\n", data[i]);
                    if (sqrt(dist) < mindist){
                        mindist = sqrt(dist);
                        for(int j = 0; j < dim; j++){
                            result_pt[j] = point[j];
                        }
                    }
                    dist = 0;
                }
                if(i<cluster_start[k] + cluster_size[k]){
                    dist += (query_pt[i%dim] - data[i])* (query_pt[i%dim] - data[i]);
                    point[i%dim] = data[i];
                }
            }
            break;
        }
    }
    if (mindist == INT_MAX){
        printf("Given query point generates a hash that does not belong to any cluster.");
        return 0;
    }
    printf("nearest distance from the query point is %f\n", mindist);
    printf("and the point is ");

    for(int i = 0; i < dim; i++){
        printf("%lf   "   , result_pt[i]);
    }

    return 0;
}
int main(void) {
    //Initializing the variables.
    int dim = 16;
    int ndata = 100000;
    int m = 5;
    double *data = malloc(ndata*dim*sizeof(double));
    double *hyper = malloc(m*dim*sizeof(double));
    double* h = malloc((ndata* m) * sizeof(double));
    double b[m];
    double W = 1;
    int nclusters = 10000;  //initialized the number of clusters to 10000 as we don't know the number beforehand. setting this to a maximum possible value for 1 million points(as per my trials).
    // If w is not 1 then nclusters should be increased to 100000
    int *cluster_start = malloc(nclusters*sizeof(int));
    int *cluster_size = malloc(nclusters*sizeof(int));
    double cluster_hashval[nclusters][m];

    double query_pt[dim];
    double result_pt[dim];
    // creating dataset for points and hyperplanes
    for (int i = 0; i < dim; i++){
        query_pt[i] = ((double)rand() / (double)RAND_MAX) ;
    }

    for (int i = 0; i < ndata*dim; i++){
        data[i] = (double)rand() / (double)RAND_MAX;
    }
    for (int i = 0; i < m*dim; i++){
        hyper[i] = (double)rand() / (double)RAND_MAX;
    }
    //Calling LSH to create the Clusters and new sorted dataset.
    nclusters = LSH(dim, ndata, data,
                    m, W, h, b, cluster_start, cluster_size, cluster_hashval, hyper);
    // Calling the search LSH to find the nearest point to the given point.
    search_LSH(dim, ndata, data, m, W, h, b, nclusters, cluster_start, cluster_size,
               cluster_hashval, query_pt, result_pt, hyper);
    printf("\n number of clusters: %d", nclusters );
    return 0;
}
