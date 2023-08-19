#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// This routine finds the initial centroids before kmeans.
int initial_centers(int dim, int i0, int im, double *data, int k, double *cluster_centroid) {
    int j=2;
    int m = -1;
    double point[dim];
    double new_point[dim];
    double cent_dist[k];
    double dist = 1000000;
    double tmp=0;

    while(j<k){  // looping to find the centroids
        float max = -99999;
        //traversing through each point of the data.
        for (int i = 0; i < dim*(im+1); i++){
            point[i%dim] = data[i];
            if (i%dim == dim-1){
                m++;
                for(int l = 0; l < j*dim; l++){
                    tmp += (point[l%dim]-cluster_centroid[l]) * (point[l%dim]-cluster_centroid[l]);
                    if (l%dim == dim-1){
                        if(tmp < dist){
                            dist = tmp;
                        }
                        tmp = 0;
                    }
                }

                if (dist > max){
                    max = dist;
                    for(int z=0;z<dim;z++){
                        new_point[z] = point[z];
                    }
                }
                dist = 100000;
            }
        }
        // assigning all the centroids to cluster_centroid pointer
        for(int z=0;z<dim;z++){
            cluster_centroid[(j*dim)+z]=new_point[z];
        }
        j++;

    }
    return 0;
}


//This function returns the distance between given two points.
double find_distance(double x[], double y[], int dim){
    double dist = 0;
    for(int i=0; i < dim;i++){
        dist += (x[i]-y[i])*(x[i]-y[i]);
    }
    return sqrt(dist);
}

//This function returns the nearest cluster and point from the given query point.
int search_kmeans(int dim, int ndata, double *data, int k, int *cluster_start, int *cluster_size, double *cluster_radius, double *cluster_centroid, double *query_pt, double *finalPt){
    double least_dist_from_pt = (double)RAND_MAX;
    double min_dis_from_prev = 0;
    int findingPt = 0;
    int minCluster = -1;
    double centroid[dim];
    double point[dim];
    double *cluster_dis = malloc(k*sizeof(double));

    for(int q=0; q < k; q++){
        centroid[q%dim] = cluster_centroid[q];
        point[q%dim] = data[q];
       cluster_dis[q] = find_distance(query_pt, point, dim);
    }
    int pause = 0;
    while(pause == 0 ){
        double min = (double)RAND_MAX;
        int index = -100;
        for(int z=0;z < k; z++){
            if(cluster_dis[z] < min){
                min = cluster_dis[z];
                index =z;
                minCluster = z;
            }
        }
        if(index == -100)
            break;
        if(index != -100){
                if( cluster_radius[index] < cluster_dis[index] - least_dist_from_pt){
                   cluster_dis[index] = (double)RAND_MAX;
                    break;
                }

        }
        for(int i=cluster_start[index];i<cluster_size[index]+cluster_start[index];i++){
            double distance = find_distance(&data[i*dim], query_pt, dim);
            if(distance<least_dist_from_pt){
                least_dist_from_pt = distance;
            }

            if(distance == least_dist_from_pt){
                least_dist_from_pt = distance;
                min_dis_from_prev = distance;
                for(int j=0; j < dim; j++){
                    finalPt[j] = data[i*dim+j];
                }
            }
        }
        cluster_dis[index] = (double)RAND_MAX;
        findingPt = findingPt + cluster_size[index];
    }
    printf("minimum distance: %f is found in cluster %d\n", least_dist_from_pt, minCluster);
    return findingPt;
}

//This function forms the clusters and its centroids.
double kmeans(int dim, int i0, int im, double *data, int k,
              int *cluster_assign, int *cluster_start, int *cluster_size,
              double *cluster_radius, double *cluster_centroid){

    short stop_iteration = 0;
    int iteration_count = 0;
    double err=0;
    double point[dim];
    double centroid[dim];
    double cluster_centroid_dummy[im][dim];
    short cluster_counter= -1;
    int point_counter = -1;
    while(stop_iteration == 0){  //looping untill we get same clusters and centroids.
        err = 0.0;
        for(int c =0; c < k; c++){
            cluster_size[c] = 0;
        }
        int count_cluster_change = 0;
        for(int p=i0; p < im*dim;p++){  //looping through all the data points
            double minDistance = 9999999;
            short isSameCluster = -1;
            point[p%dim] = data[p];
            if (p%dim == dim-1){
                point_counter++;
                for(int c=0; c < k*dim;c++){  //looping through all the cluster centroids
                    centroid[c%dim] = cluster_centroid[c];
                    if (c%dim == dim-1){
                        cluster_counter+=1;
                        double distance = find_distance(centroid, point, dim);
                        if(distance<minDistance){
                            minDistance = distance;
                        }
                        if(distance == minDistance){
                            minDistance = distance;
                            isSameCluster = cluster_counter;
                        }
                    }
                }
                cluster_counter = -1;
                if(cluster_assign[point_counter] == -1 || cluster_assign[point_counter] != isSameCluster){
                    count_cluster_change++;
                }
                cluster_assign[point_counter] = isSameCluster;   //assigning the cluster values or each data point.

                err = err + minDistance/(im-i0);
            }

        }
        point_counter = -1;
        // Initialized a 2-d array for cluster centroids.
        for(int c=0; c < k;c++){
            for(int d=0; d < dim; d++){
                cluster_centroid_dummy[c][d] = 0;
            }
        }
        //Calculating the new cluster centroids from the clusters formed from the above step.
        for(int p=i0; p < im;p++){
            cluster_size[cluster_assign[p]]++;
            for(int d=0; d < dim; d++){
                cluster_centroid_dummy[cluster_assign[p]][d] += data[p*dim+d];
            }
        }
        for(int c=0; c < k;c++){
            if(c == 0){
                cluster_start[c] = 0;
            } else {
                cluster_start[c] =  cluster_start[c-1] + cluster_size[c-1];
            }
            for(int d=0; d < dim; d++){
                if(cluster_size[c] > 0){
                    cluster_centroid_dummy[c][d] = cluster_centroid_dummy[c][d]/cluster_size[c];

                }
            }
        }
        //Reassigning the cluster centroid values back to cluster_centroid pointer.
        for (int i = 0; i < k; i++){
            for (int j = 0; j < dim; j++){
                cluster_centroid[(i*dim) + j] = cluster_centroid_dummy[i][j];
            }
        }
        if(count_cluster_change == 0){
            stop_iteration = 1;
        }
        iteration_count++;
    }
    //Finding the radius of a cluster by taking the distance between farthest point and centroid.
    int curr_loc = 0;
    for(int i=0; i < k; i++){
        double far_point_distance = 0;
        int far_point=0;
        for(int j=i0; j < im; j++){
            if(cluster_assign[j]==i){
                for(int k=0; k < dim; k++){
                    double temp = data[j*dim+k];
                    data[j*dim+k] = data[curr_loc*dim+k];
                    data[curr_loc*dim+k] = temp;
                }
                cluster_assign[curr_loc] = k;
                double distance = find_distance(&data[j*dim], cluster_centroid_dummy[i], dim);
                if(far_point_distance<distance){
                    far_point_distance = distance;
                }

                if(distance == far_point_distance){
                    far_point_distance = distance;
                    far_point = curr_loc;
                }
                curr_loc++;
            }
        }
        cluster_radius[i]=far_point_distance;
    }
    return err;
}


// Main has the function calls to kmeans and kmeans search functions.
int main(void) {
    int dim = 16;
    int ndata = 10000;
    int total_centroid = 100;
    double *data = malloc(ndata*dim*sizeof(double));
    double centroids[total_centroid*dim];
    double centroid_one[dim];
    int *cluster_size = malloc(total_centroid*sizeof(int));
    int *cluster_assign = malloc(ndata*sizeof(int));
    int *cluster_start = malloc(total_centroid*sizeof(int));
    double *cluster_radius = malloc(total_centroid*sizeof(double));
    double centroid_two[dim];
    double point[dim];
    //Creating Random Data
    for (int i = 0; i < ndata*dim; i++){
        data[i] = (double)rand() / (double)RAND_MAX;
    }

    // first centroid
    for (int i = 0; i < dim; i++){
        centroid_one[i] = data[i];
    }

    // finding 2nd centroid
    double dis=-100;
    for (int i = 0; i < ndata*dim; i++){
        double tmp=0;
        for(int j=0;j<dim;j++){
            tmp += (data[i]-centroid_one[j]) * (data[i]-centroid_one[j]);
            point[j] = data[i];
            i++;
        }
        i--;
        tmp = sqrt(tmp);
        if(tmp>dis){
            dis=tmp;
            for(int k=0;k<dim;k++){
                centroid_two[k] = point[k];
            }
        }
    }

    for (int i = 0; i < dim; i++){
        centroids[i] = centroid_one[i];
    }

    for (int i = dim; i < 2*dim; i++){
        centroids[i] = centroid_two[i%dim];
    }
  // Initializing cluster assign pointer
    for(int i=0; i < ndata; i++){
        cluster_assign[i] = -1;
    }
  // Calling intital_centers to find the initial centroids of the clusters.
    initial_centers(dim,0,ndata-1,data,total_centroid,centroids);

    double err = kmeans(dim, 0, ndata, data, total_centroid, cluster_assign, cluster_start, cluster_size, cluster_radius, centroids);
    printf("clusters built \n");
    for(int i = 0; i < total_centroid; i ++){
        printf("cluster %d start %d size %d , radius %f\n", i, cluster_start[i], cluster_size[i], cluster_radius[i]);
    }
    printf("Error from kmeans is %f\n", err);

    printf("Query point is : ");
    double *queryPoints = malloc(dim*sizeof(double));
    int i=0;
    while(i<dim){
      queryPoints[i] = (double)rand() / (double)RAND_MAX;
        printf("%f ", queryPoints[i]);
      i++;
    }

    printf("\n");
    // pointer for nearest point
    double *finalPt = malloc(dim*sizeof(double));
    int points_explored = search_kmeans(dim, ndata, data, total_centroid, cluster_start, cluster_size, cluster_radius, centroids, queryPoints,finalPt);
    // print nearest point

    printf("Nearest point from the given query point is: ");
    for(int i=0; i< dim; i++){
        printf("%f ", finalPt[i]);
    }
    printf("\nand points checked are %d", points_explored);
    printf("\n");
    return 0;
}

