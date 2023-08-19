#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


// This function takes a cluster as input and divided it into two clusters based on the chosen dimension which has the highest variance among all the domensions
int bipartition(int dim, int i0, int im, double *data, int chosen_dim, double mean, int *cluster_start, int *cluster_size,double **cluster_bdry, double **cluster_centroid,short *cluster_assign){
    int l = 0;
    int first_cls = 0;
    int second_cls = im-i0-dim;
    int a = 0;
    double *cls_data = malloc((im-i0)*sizeof(double));
    cluster_size[0]=0;
    cluster_size[1]=0;
  // Initializing varaibles
    for(int i=0; i < dim; i++){
        cluster_bdry[0][i*2]=999999;
        cluster_bdry[0][i*2+1]=0;
        cluster_bdry[1][i*2]=999999;
        cluster_bdry[1][i*2+1]=0;
        cluster_centroid[0][i]=0;
        cluster_centroid[1][i]=0;
    }

    // Going through all the datapoints between start and end of a cluster and dividing all the points based on dimension.
    for(int k=i0; k < im; k = k+dim){
      a = (data[k+chosen_dim] > mean) ? 1:0;
      for(int q=k; q < k + dim; q++){
        if(data[q] > cluster_bdry[a][2*(q%dim)+1]) {cluster_bdry[a][2*(q%dim)+1] = data[q];}
        if(data[q] < cluster_bdry[a][2*(q%dim)])   {cluster_bdry[a][2*(q%dim)]   = data[q];}
        cluster_centroid[a][q%dim] += data[q];
        cluster_size[a]++;
        if (a == 1){
          cls_data[second_cls]=data[q];
          second_cls++;
        }
        else{
          cls_data[first_cls]=data[q];
          first_cls++;
        }
      }
      if (a==1){
        second_cls = second_cls-2*dim;
      }
      l++;
    }

    cluster_start[0]=i0;
    cluster_start[1]=i0+cluster_size[0];
  // assigning the value of centroid for each dimension.
    for(int i=0; i < dim; i++){
        cluster_centroid[0][i] = cluster_centroid[0][i]/cluster_size[0]*dim;
        cluster_centroid[1][i] = cluster_centroid[1][i]/cluster_size[1]*dim;
    }
    for(int i =i0; i<im; i++){
        data[i] = cls_data[i-i0];
    }
    return 0;
}
// This Function creates a KD Tree from the given dataset.
int kdtree(int dim, int ndata, double *data, int kk,int *cluster_start, int *cluster_size,double **cluster_bdry, double **cluster_centroid,short *cluster_assign){
    // initializing the cluster to 1.
    int jj = 1;
  //Looping through all the clusters
    while (jj < kk)
    {
       for (int j = 0; j < jj; j++)
       {
            // Initializing the variables
            int next_ind = (kk/(jj))*j;
            int cluster1_pointer = next_ind;
            int cluster_start1 = cluster_start[cluster1_pointer];
            int cluster_end = cluster_start[cluster1_pointer]+cluster_size[cluster1_pointer]*dim;
            int cluster2_pointer = next_ind + kk/(2*jj);
            int cluster_ind = 0;
            double mean[dim];
            double variance[dim];
            int max_dim = 0;
            double max_var = 0;
            for(int i=0; i < dim; i++){
                mean[i] = 0;
            }
            int cluster_start_bi[2];
            int cluster_size_bi[2] ;
            double *cluster_bdry_bi[2];
            for(int i =0; i < 2; i++)
                cluster_bdry_bi[i]=(double*)malloc(2*dim*sizeof(double));
            double *cluster_centroid_bi[2];
            for(int i =0; i < 2; i++)
                cluster_centroid_bi[i]=(double*)malloc(dim*sizeof(double));
            short *cluster_assign_bi =  malloc(ndata*sizeof(short));
            for(int k=cluster_start1; k < cluster_end; k++){
              if (cluster_ind == dim){
                cluster_ind = 0;
              }
              mean[cluster_ind] += data[k];
              cluster_ind++;
            }
         //Calculating the mean of the cluster
            for(int i=0; i < dim; i++){
                mean[i] /= cluster_size[cluster1_pointer];
            }
         //calculating the variance of each dimension to find the highest variance
            cluster_ind = 0;
            for(int k=cluster_start1; k < cluster_end; k++){
              if (cluster_ind == dim){
                cluster_ind = 0;
              }
              variance[cluster_ind] += (data[k]-mean[cluster_ind])*(data[k]-mean[cluster_ind]);
              cluster_ind++;
            }

            for(int i=0; i < dim; i++){
                if(variance[i] > max_var){
                    max_var = variance[i];
                    max_dim = i;
                }
            }
            printf("for cluster %d max variance is %f and dim is : %d\n",jj, max_var, max_dim);
            // calling bypartition method to divide the cluster into two partitions
            bipartition(dim, cluster_start[cluster1_pointer], cluster_start[cluster1_pointer] + cluster_size[cluster1_pointer]*dim , data, max_dim,mean[max_dim], cluster_start_bi, cluster_size_bi, cluster_bdry_bi, cluster_centroid_bi, cluster_assign_bi);
            for(int i = 0; i<2 ;i++){
              int cls;
              cls = i == 0 ? cluster1_pointer: cluster2_pointer;
              cluster_start[cls] = cluster_start_bi[i];
              cluster_size[cls] = cluster_size_bi[i]/dim;
            }
         // reassigning the output values from bipartition procedure to variables.
            for(int i=0; i < dim; i++){
                cluster_bdry[cluster1_pointer][i*2]=cluster_bdry_bi[0][i*2];
                cluster_bdry[cluster1_pointer][i*2+1]=cluster_bdry_bi[0][i*2+1];
                cluster_bdry[cluster2_pointer][i*2]=cluster_bdry_bi[1][i*2];
                cluster_bdry[cluster2_pointer][i*2+1]=cluster_bdry_bi[1][i*2+1];
                cluster_centroid[cluster1_pointer][i]=cluster_centroid_bi[0][i];
                cluster_centroid[cluster2_pointer][i]=cluster_centroid_bi[1][i];
            }
       }
       jj = 2*jj;
    }

    return 0;
}
// This Routine searches the given data point in the kdtree, if the value is not found it finds the nearest value in the cluster.
int search_kdtree(int dim, int ndata, double *data, int kk,int *cluster_start, int *cluster_size, double **cluster_bdry,double *query_pt, double *result_pt){
    int found = -1;
    double nearest = (double)RAND_MAX;
    for(int i=0; i < kk; i++){
        for(int j=0; j< dim; j++){   // If the value is not found
            if(cluster_bdry[i][2*j]>query_pt[j] || cluster_bdry[i][2*j+1]<query_pt[j]){
               if((fabs(cluster_bdry[i][2*j]-query_pt[j]) < nearest)){
                    found = i;
                    nearest = fabs(cluster_bdry[i][2*j]-query_pt[j]);
                } else if(fabs(cluster_bdry[i][2*j+1]-query_pt[j]) < nearest){
                    found = i;
                    nearest = fabs(cluster_bdry[i][2*j+1]-query_pt[j]);
                }

            }
          else{   //If the value is found
            found = i; }

        }
    }
    double *min_dis_point =  malloc(dim*sizeof(double));
// This finds the distance between the closet point in the cluster
 double min_dis = (double)RAND_MAX;
  for(int i=cluster_start[found]; i < cluster_start[found] + cluster_size[found]*dim;i = i + dim){
        double distance = 0;
        for(int j=0; j < dim; j++){
            distance += (data[i+j]-query_pt[j])*(data[i+j]-query_pt[j]);
        }
        distance = sqrt(distance);
        if(distance < min_dis){
            min_dis = distance;
            for(int j=0; j < dim; j++){
                min_dis_point[j] = data[i+j];
            }
        }
    }
    printf("nearest distance within the cluster is %f\n", min_dis);
    return found;
}

int main()
{
    int dim = 5;
    int ndata = 10000;
    int kk = 32;
    double *data = malloc(ndata*dim*sizeof(double));
    double *cluster_centroid[kk];
    int cluster_start[kk];
    int cluster_size[kk];
    double *cluster_bdry[kk];
    short *cluster_assign = malloc(ndata*sizeof(short));
    for(int i =0; i < kk; i++) {
        cluster_bdry[i] = (double *) malloc(2 * dim * sizeof(double));
        cluster_centroid[i] = (double *) malloc(dim * sizeof(double));
    }
    // Generating a random dataset using rand() function.
    for (int i = 0; i < ndata*dim; i++){
        data[i] = (double)rand() / (double)RAND_MAX;
    }

  printf("\n");
    //testing the KDtree funcionality
    cluster_size[0]=kk;
    cluster_start[0]=0;

    kdtree(dim, ndata, data, kk, cluster_start, cluster_size, cluster_bdry, cluster_centroid, cluster_assign);
    for(int i = 0; i < kk; i++){
       printf("cluster %d starts at %d and has a size %d\n", i, cluster_start[i], cluster_size[i]);
    }
    // demo of searching operation
    double *query_pt = malloc(dim*sizeof(double));
    int i=0;
    while(i<dim){
        query_pt[i] = 0.8;
        i++;
    }
    double *result_pt = malloc(dim*sizeof(double));
    int cluster = search_kdtree(dim, ndata, data, kk, cluster_start, cluster_size, cluster_bdry, query_pt,result_pt);
    printf("nearest point to given point in the cluster %d\n", cluster);
    return 0;
}

