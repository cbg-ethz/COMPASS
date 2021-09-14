/* Copyright@Louxin Zhang, National University of Singapore, 2020 
 *   
  This porgram computes the pairwise CASet measures between a 
     set of 1-labeled trees with the same root, which is the first
     root in the imput trees. 
     and their pairwise distances are output
     to <dist_output_file>. In particular, if the two trees are given,
     the CASet distance between them will be computed.

   Tree file format:
     Input rooted trees are listed one by one. each tree is given
     by a head line followed by the edges.
     The headline give the id of a tree and the number of tree edges
    For example, the header line "#tree 10: 28" represents the 10-th tree of
    28 tree edges.

   a tree edge line has two parts. For example, "a b_c_d" represents an edge from
   a node that is labeled with a to a node that is labeled with
   a subset {b, c, d} of three labels b, c, d;
   the root is assumed to the the first node.

 For example, the following five four lines specify the tree
  with 5 nodes which are labled with 0, 1, {2, 5}, 3 and 4  and
 whose root is 0:
 #tree 0:  4
 0 1
 1 3
 3 2_5
 3 4

 LIMITATION: (1) the program uses array as data strucutre to
 store trees; (2) the max. size of an input tree is set to be 40
 nodes; (3) If a node is labeled with a set of mutations, they are seprared
 with '_'. (4) the program is only applied to lableld trees with the
 with integer labels from 0 to 40. (5) The number of input trees is set to 20000.



 *   Complie command: gcc -lm Pairwise_Triplets_v2.c -o PW_Triplets
 *   The run command: ./PW_Triplets  <inputtrees_file> <dist_output_file> no_nodes
 *
 */
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NODES_NO 30
#define MAXSIZE  31 
#define MAX_NO_TREES 20000


 
typedef struct {
    unsigned int first;
    unsigned int second;
} edge;




void Decode(char *str, char a, int *num, char *names[]){
   int len, len1;
   int i, j,k;
   char temp[20];


   /*
    printf("%s  %d\n", str, *num);
   */

   k=0;
   len=strlen(str);
   len1=0;
   for (i=0; i<len; i++) {
     if (str[i]!='_') { temp[len1]=str[i]; len1 +=1;}
     if (str[i]=='_' || i==len-1) {
            temp[len1]='\0';
            names[*num]=(char *)malloc(len1); strcpy(names[*num], temp);
            *num +=1;
            len1=0;
     }

   } /* for */
}

 
edge *prufer_tree(const unsigned int *code, size_t len)
{
    unsigned int i, j, e = 0;
    const unsigned int n = len + 2; /* Nodes */
    edge *tree = malloc((n - 1) * sizeof(edge));
    unsigned int *degrees = calloc(n + 1, sizeof(unsigned int)); /* 1-based */
    if (tree == NULL || degrees == NULL) {
        free(tree);
        free(degrees);
        return NULL;
    }
    /* Set the degrees from the occurrences in the code */
    for (i = 0; i < len; i++) {
        degrees[ code[i]]++;
    }
    /* Add 1 to them all */
    for (i = 1; i <= n; i++) {
        degrees[i]++;
    }
    /* Add edges to nodes in the code */
    for (i = 0; i < len; i++) {
        /* Find the lowest-numbered node with degree 1 */
        for (j = 1; degrees[j] != 1; j++);
        /* Add the edge */
        tree[e].first = j;
        tree[e].second = code[i];
        degrees[j]--;
        degrees[ code[i]]--;
        e++;
    }
    /* Find the last 2 degree-1 nodes */
    for (i = 1; degrees[i] != 1; i++);
    for (j = i + 1; degrees[j] != 1; j++);
    /* Add the last edge */
    tree[e].first = i;
    tree[e].second = j;
    free(degrees);
 
    return tree;
}




int Check_Name(char *node_strings[], int no_nodes, char *str1){
  int i;

  if (str1==NULL) return -1;
  for (i=0; i<no_nodes; i++) {
	if (strcmp(str1, node_strings[i])==0) return i;
  } 
  return -1;
}



int min_of_two(int a, int b){
  if (b==0) return a; 
  else if (a==0)  return b;
  
  if (a< b) return a; else return b;
}



void Adjacent_matrix(int no1, int start1[], int end1[], int adjacent[][MAXSIZE]){
 int i, j, k;

 for (i=0; i<no1; i++) {
   for (j=0; j<no1; j++) {adjacent[i][j]=0;}
 } 

 for (i=0; i<no1; i++) {
     for (k=0; k<no1-1; k++) {	   
       if (start1[k]==i) { adjacent[i][end1[k]]=1; } 
       if (end1[k]==i) { adjacent[i][start1[k]]=1; }  
     }
 }


} /* Adjacent */


int Check_Map_Image( int j,  int no1,  int map[]) {
   int i; 
   for (i=0; i<no1; i++) {
	 if (map[i]==j) return i;
   }
   return -1;
}


void Distance_Comput(int no1, int start1[], int end1[], int distance[][MAXSIZE]){

 int i, j, k;
 /*
 int distance[MAXSIZE][MAXSIZE];
 */
 int adjacent[MAXSIZE][MAXSIZE];
 int dist;
 int total;

 Adjacent_matrix(no1, start1,  end1,  adjacent);


 total=0;
 for (i=0; i<no1; i++) {
   for (j=0; j<no1; j++) { 
	   distance[i][j]=adjacent[i][j];
	  if (distance[i][j]!=0) total=1+total; 
   }
 }


 while (total < (no1*no1 - no1)){
/* printf("total[%d]\n", total); */
  for (i=0; i<no1; i++) {
     for (j=0; j<no1; j++) { 
       if (i!=j && distance[i][j]==0) {
	 for (k=0; k<no1; k++) {
	   if (distance[i][k]!=0 && adjacent[k][j]!=0) 
	       distance[i][j]=min_of_two(distance[i][j], 
			      distance[i][k]+adjacent[k][j]);
	 }
	 if (distance[i][j]!=0) { total=1+total; }
       } 
   } /* j */
  } /* for */
 } /* while */
 /* obtain distance matrix */
}


void   Edge_Partitions(int no1,  int start1[], int end1[],  int distance[][MAXSIZE], int partition1[][MAXSIZE]){
 int i, j, k;
 /*
 int distance[MAXSIZE][MAXSIZE];
 int adjacent[MAXSIZE][MAXSIZE];
 */
 int dist;
 int total;



 for (i=0; i<no1-1; i++){
    /* printf("Edge[%d]:(%2d,%2d)\n", i, start1[i], end1[i]); */ 
       partition1[i][0]=no1;
      /*  printf("%2d:  ", partition1[i][0]); */
   for (j=0; j<no1; j++) {
      if (distance[start1[i]][j]>distance[end1[i]][j])
	    partition1[i][j+1]=0;
      else partition1[i][j+1]=1;
      /* printf("%2d ", partition1[i][j+1]);  */
   }
    /* printf("\n"); */ 
 }

} /* edge_partition */


int Comput_size(int edge, int v1,  int p1[][MAXSIZE], int no1){
   int i;
   int s1, s0;
   
   s1=0; 
   for (i=1; i<=no1; i++){
     if (p1[edge][i]==v1) s1=1+s1;
   }  
   return s1;
}

int Is_Same(int ind1, int ind2, int v1, int v2, int p1[][MAXSIZE], int p2[][MAXSIZE], int no1, int map[]) {
  int i, j, k;

  for (i=1; i<=no1; i++) {
     if (p1[ind1][i]==v1 && p2[ind2][map[i]]!=v2) return 0;
  }
  return 1;
}	
 

void Comput_map(int map[], int no1, char *tree1_names[], int no2,  char *tree2_names[]){
  int i, j;

  
  for (i=0; i<no1; i++){
    map[i+1]=i+1;
    /*
    map[i+1]=-1;
    for (j=0; j<no2; j++) {
      if (strcmp(tree1_names[i], tree2_names[j])==0) 
      {  map[i+1]=j+1; break;}
    }
    */
  }
}






void Ancester_Dist(double *a_dist, int no_mut1, int no_mut2, 
		int part1[][MAXSIZE], int part2[][MAXSIZE], int map[]){
    int i, j, k; 
    double dist, no_common;
    int common, as[2];
    int no_ans1[MAXSIZE], no_ans2[MAXSIZE];

    dist=0.0;
    no_common=0;
    for (i=0; i<no_mut1; i++) {
	no_ans1[i]=0;
	if (map[i]!=-1) no_common +=1;
	for (j=0; j<no_mut1; j++) { if (part1[i][j]==1) no_ans1[i] +=1;  }
    }
    for (i=0; i<no_mut2; i++) {
	no_ans2[i]=0;
	for (j=0; j<no_mut2; j++) { if (part2[i][j]==1) no_ans2[i] +=1;  }
    }

    for (i=0; i<no_mut1; i++) {
            common=0; 
	if (map[i]!=-1) {
            for (k=0; k<no_mut1; k++) {
             if (part1[i][k]==1) {
		 if (map[k]!=-1 && part2[map[i]][map[k]]==1) common +=1;
             } 
            }
	  if (no_ans1[i]+no_ans2[map[i]]!=0){
           dist=dist+(double)(no_ans1[i]+no_ans2[map[i]]-2*common)/(no_ans1[i]+no_ans2[map[i]]-common);
	  } 
       }
    }


    *a_dist =(dist+no_mut1+no_mut2-2*no_common)/(double)(no_mut1+no_mut2- no_common);
} /* anscestor */


void  Triplets_Dist1(int *cat, double *top_dist, double *low_dist, int no_mut1, 
		int no_mut2, int part1[][MAXSIZE], int part2[][MAXSIZE], int map[]){
     int i, j, k;
    double t_dist, l_dist;
	    int no_common;
    int ij, map_ij, both_ij;
    int inotj,  map_inotj,  both_inotj; 
    int jnoti,  map_jnoti,  both_jnoti; 
    int no_ans1[MAXSIZE], no_ans2[MAXSIZE];
    int both_ances1[MAXSIZE][MAXSIZE][MAXSIZE];
    int both_ances2[MAXSIZE][MAXSIZE][MAXSIZE];
    int diff_ances1[MAXSIZE][MAXSIZE][MAXSIZE];
    int diff_ances2[MAXSIZE][MAXSIZE][MAXSIZE];


    t_dist=0.0; l_dist=0.0;



    for (i=0; i<no_mut1; i++) {
	 for (k=0; k<no_mut1; k++)  both_ances1[i][i][k]=0;
         for (j=i+1; j<no_mut1; j++) {
	      for (k=0; k<no_mut1; k++) {
		    if (part1[i][k]==1 && part1[j][k]==1){
	                 both_ances1[i][j][k]=1; both_ances1[j][i][k]=1; 
		    } else {
	                 both_ances1[i][j][k]=0; both_ances1[j][i][k]=0; 
	            }
	      } 
         }
    } /* for */

    for (i=0; i<no_mut2; i++) {
	 for (k=0; k<no_mut1; k++)  both_ances2[i][i][k]=0;
         for (j=i+1; j<no_mut2; j++) {
	      for (k=0; k<no_mut2; k++) {
		    if (part2[i][k]==1 && part2[j][k]==1){
	                 both_ances2[i][j][k]=1; both_ances2[j][i][k]=1; }     
		    else { both_ances2[i][j][k]=0; both_ances2[j][i][k]=0; } 
	      } 
         }
    }

    for (i=0; i<no_mut1; i++) {
	  for (j=i+1; j<no_mut1; j++) {
             ij=0;  map_ij=0;  both_ij=0;
	        for (k=0; k<no_mut1; k++){
	          if (both_ances1[i][j][k]==1) { 
		     ij +=1;
                     if (map[k]!=-1 && both_ances2[map[i]][map[j]][map[k]]==1) both_ij +=1;
		  }
		}
		for (k=0; k<no_mut2; k++) {
	          if (both_ances2[map[i]][map[j]][k]==1) map_ij +=1;
	        }
		if (ij+map_ij!=0) 
		  t_dist = t_dist + (double)(ij+map_ij-2*both_ij)/(double)(ij+map_ij-both_ij);

         } /* j loop */
    } /* done for first tree side*/


    
    *top_dist= 2*t_dist/((no_mut1)*(no_mut1-1));
    

    *cat=floor(*top_dist/0.025);

 } /* triple dist1 */



	

/* Parts1 is a 2-dim array  parts1[i][k] 
 * indicates whether k is an acnestor or not
 * 1 if it k is an ancestor of only i. 
 * 0 if it is neither.
 * */
void Triplets(int no1, int start1[], int end1[], int distances1[][MAXSIZE],
int parts1[][MAXSIZE]){
  int i, j, k;

  for (i=0; i<no1; i++) {
      for (k=0; k<no1; k++) parts1[i][k]=0;
      /* parts1[i][i]=1; */
      for (k=0; k<no1; k++) { 
	 if (i!=k && distances1[i][k] +distances1[0][k]==distances1[0][i])
              parts1[i][k]=1; 
      } /* k lopp */ 
  } /* i loop */

} /* triplelets */

void Triplets_rooting(int no1, int start1[], int end1[], int distances1[][MAXSIZE], int parts1[][MAXSIZE], int r){
  int i, j, k;

  for (i=0; i<no1; i++) {
      for (k=0; k<no1; k++) parts1[i][k]=0;
      /* parts1[i][i]=1; */
      for (k=0; k<no1; k++) { 
	 if (i!=k && distances1[i][k] +distances1[r][k]==distances1[r][i])
              parts1[i][k]=1; 
      } /* k lopp */ 
  } /* i loop */

} /* triplelets */

void  Remove_Space(char h_line[], int len){
      int i;

      for (i=0; i<len-1; i++) h_line[i]=h_line[i+1];
      h_line[len-1]='\0';
}

void Print_Tree(edge *tree, int size){
   int i;

   for (i=0; i<size-1; i++) {
     printf("%d %d\n", tree[i].first, tree[i].second);
   }
   printf("\n");
}


void Add_One(int count[][20], int sub){
     int i, pos, done;

     done=0; pos=0;
     while (done ==0) {
	     if (count[sub][pos]<9) { 
		     count[sub][pos]=count[sub][pos]+1; done=1;
		 } else  { count[sub][pos]=0; pos=pos+1; }
     }
}	


void Print_Out_Count(FILE *file, int count[], int size){
	int i, k;

	 i=size-1;
	 while (count[i]==0 && i>=0) i=i-1;
	 if (i==-1) fprintf(file, "0");
	 for (k=i; k>=0; k--) { fprintf(file, "%d", count[k]); }
	 fprintf(file, "\n");
}


short Check_names(char *label_names[], short  node_count,  short *node_ind, 
      char *str){
    int i;

    for (i=0; i<node_count; i++) {
        if (strcmp(label_names[i], str)==0) {
           *node_ind=i;
           return 1;
        }
    } 
    *node_ind=node_count;
   return 0;
}


void Adjacent_matrix_Jan(short  no1, short graph[][2],
  short adjacent[][MAXSIZE]){
 short i, j, k;

 for (i=0; i<no1; i++) {
   for (j=0; j<no1; j++) {adjacent[i][j]=0;}
 }

 for (i=1; i<=no1-1; i++) {
   adjacent[graph[i][0]][graph[i][1]]=1;
   adjacent[graph[i][1]][graph[i][0]]=1;
 }


} /* Adjacent */

void Distance_Comput_Jan(short graph[][2],  short distance[][MAXSIZE]){

 short i, j, k;
 short adjacent[MAXSIZE][MAXSIZE];
 short  dist;
 short  total;
 short no1;


 no1=graph[0][0]+1;
 Adjacent_matrix_Jan(no1, graph,  adjacent);


 total=0;
 for (i=0; i<no1; i++) {
   for (j=0; j<no1; j++) {
          distance[i][j]=adjacent[i][j];
   }
 }
 total=2*(no1-1);

 while (total < (no1*no1 - no1)){
/* printf("total[%d]\n", total); */
  for (i=0; i<no1; i++) {
     for (j=0; j<no1; j++) {
       if (i!=j && distance[i][j]==0) {
         for (k=0; k<no1; k++) {
           if (distance[i][k]!=0 && adjacent[k][j]!=0)
               distance[i][j]=min_of_two(distance[i][j],
                              distance[i][k]+adjacent[k][j]);
         }
         if (distance[i][j]!=0) { total=1+total; }
       }
   } /* j */
  } /* for */
 } /* while */
 /* obtain distance matrix */
}

short  Orientation_Jan(char *r,  char *node_labels[],  short graph[][2], short distance[][MAXSIZE]){
   short i,  k;
   short no_edges;
   short root;


   root=0;
   no_edges=graph[0][0];
   for (i=0; i<=no_edges; i++) {
     if (strcmp(r, node_labels[i])==0) { root=i; break;}
   }

   for (i=1; i<=no_edges; i++) {
     if (distance[root][graph[i][0]]> distance[root][graph[i][1]]) {
            k=graph[i][0]; graph[i][0]=graph[i][1]; graph[i][1]=k;
     }
   }
   return  root;
}




void Print_Trees(short graph[][2], char *label_names[]){
  short  j;

  for (j=1; j<=graph[0][0]; j++)
      printf("%s %s\n", label_names[graph[j][0]], label_names[graph[j][1]]);
        printf("\n");
}


void Initialization(short no_nodes, short no_mut, char partition[][MAXSIZE]){
    short i, j;

    for (i=0; i<no_nodes; i++){
       for (j=0; j<no_mut; j++) { partition[i][j]='a'; }
       partition[i][no_mut]='\0';
    }
}



void   Comput_Mutation(short no_nodes, char *node_labels[], short no_muts, 
   char node_muts[][MAXSIZE], short mut_to_node[]){
   short i, j, k;
   int num;
   char *mutations[MAXSIZE];
  


   for (j=0; j<no_muts; j++) {
       mut_to_node[j]=-1;
   }

   for (j=0; j<no_nodes; j++) {
     for (i=0; i<no_muts; i++) { node_muts[j][i]='a'; }
     num=0;
     Decode(node_labels[j], '_', &num, mutations);
     for (i=0; i<num; i++) { 
       node_muts[j][atoi(mutations[i])]='1';  
       mut_to_node[atoi(mutations[i])]=j;
     }
     node_muts[j][no_muts]='\0';
     /*
     printf("%s\n", node_muts[j]);
      */
   }


}

/* for a pair of mutations, find their common ancestor  */
void Comput_Mut_Ancestors(short root,  short no_nodes, short no_muts,
   short dist[][MAXSIZE], char node_muts[][MAXSIZE], short mut_to_node_map[], 
   char *common_ancestors[][MAXSIZE]) {

   short i, j, k, m;
   char non_muts[MAXSIZE];
   short node1, node2;

   for (i=0; i<no_muts; i++) non_muts[i]='a'; non_muts[no_muts]='\0';

   for (i=0; i<no_muts; i++) {
      node1=mut_to_node_map[i];
      /*
      printf("mut(%d): node1(%d) \n", i, node1);
      */
      if (node1==-1) {
         for (j=i+1; j<no_muts; j++) {
           common_ancestors[i][j]=(char *)malloc(no_muts+1);
           strcpy(common_ancestors[i][j], non_muts);
         } 
      } else {
       for (j=i+1; j<no_muts; j++) {
        node2=mut_to_node_map[j]; 
        /*
         printf("mut-j(%d): node2(%d) \n", j, node2);
         */
        common_ancestors[i][j]=(char *)malloc(no_muts+1);
        strcpy(common_ancestors[i][j], non_muts);
        if (node2!=-1) {
         for (k=0; k<no_nodes; k++) {
           if (dist[root][k]+dist[k][node1]==dist[root][node1] &&
               dist[root][k]+dist[k][node2]==dist[root][node2]) {
             if (k!=node1 && k!=node2)
              for (m=0; m<no_muts; m++) {
                if (node_muts[k][m]=='1') common_ancestors[i][j][m]=node_muts[k][m];
              }
           } 
         } /* for */
       } /* node 2 */
     } /* j */
    } /* node 1*/
   } /* i */
   
} /* commut mut ancestor */


short CASet(float *t_dist, short no_muts, char *c_ances1[][MAXSIZE],
  char *c_ances2[][MAXSIZE]) {
   short i, j, m;
   short diff, cup;
   float score;

   

   score=0;
   for (i=0; i<no_muts; i++) {
     for (j=i+1; j<no_muts; j++) {
         diff=0;  cup=0;
         /* printf("%s\n%s\n", c_ances1[i][j], c_ances2[i][j]); */
         for (m=0; m<no_muts; m++) {
          if (c_ances1[i][j][m]=='1' && c_ances2[i][j][m]=='a') {diff++; cup++; } 
          if (c_ances1[i][j][m]=='a' && c_ances2[i][j][m]=='1') {diff++; cup++; } 
          if (c_ances1[i][j][m]=='1' && c_ances2[i][j][m]=='1') {  cup++; } 
         }
         if (cup!=0) score=score+diff/(float)cup; 
         /* printf("%d -- %d\n", diff, cup); */
     }
   }
   if (no_muts!=0 && no_muts!=1) { 
       score=2*score; score=score/no_muts; score=score/(no_muts-1);
   }
   *t_dist=score;
   return 0;
}


void Compute_LCA(short r, short  no_nodes,  short  dist[][MAXSIZE], short LCA[][MAXSIZE]){
     short a, b, k; 
     short least_dist, ancestor;


   for (a=0; a<no_nodes; a++){
     
      /* printf("inside loop: %d\n", a); */
     LCA[a][a]=a;
     for (b=a+1; b<no_nodes; b++) {
      /* printf("inside b loop: %d\n", b); */
     least_dist=dist[r][a]; ancestor=r;
       for (k=0; k<no_nodes; k++) {
	   if (dist[r][k]+dist[k][a]==dist[r][a] && 
			   dist[r][k]+dist[k][b] ==dist[r][b]){ 
           
	    if (dist[k][a]<least_dist) { least_dist=dist[k][a]; ancestor=k; } 
           }
       } /* k */
        LCA[a][b]=ancestor;
        LCA[b][a]=ancestor;
     } /* b */
   }/* a*/

} /* end of compute lca */


/* 19 nicses: */
short  Compute_triplet_type(short LCA[][MAXSIZE], short a, short b, short  c){
  
   if (a!=b && b!=c && c!=a) {
   if (LCA[a][b]==LCA[a][c] && LCA[a][b]==LCA[b][c]) {
       if (LCA[b][c]==a) return 17;  else if (LCA[a][c]==b) return 18;  
       else if (LCA[b][a]==c) return 19;  else return 7;
   } else if (LCA[a][b]==LCA[a][c]) {

       if (LCA[a][b]==a) { 
	      if (LCA[b][c]==b) return 1; else if (LCA[b][c]==c) return 2; 
       } else {
	      if (LCA[b][c]==b) return 8; else if (LCA[b][c]==c) return 9;
	      else return 14;
       }
   } else if (LCA[a][b]==LCA[b][c]){
	 if  (LCA[a][b]==b){ 
	    if (LCA[a][c]==a) return 3; else if (LCA[a][c]==c) return 4; 
         } else {
	   if (LCA[a][c]==a) return 10; else if (LCA[a][c]==c) return 11;
	   else return 15;
	 }
   } else {
        if (LCA[a][c]==c)  {
	    if (LCA[a][b]==a) return 5; else if (LCA[a][b]==b) return 6; 
	} else {
 	   if (LCA[a][b]==a) return 12; else if (LCA[a][b]==b) return 13;
	   else return 16;
        }
   }
   }

   if (a==b && b==c) return 29;

   if (a==b && b!=c) { 
       if (LCA[a][c]==a) return 22; else if (LCA[a][c]==c) return 25;
       else return 28;
   }

   if (c==b && b!=a) { 
       if (LCA[a][c]==c) return 20;
       else if (LCA[a][c]==a) return 23; else return 26;
   }

   if (c==a && b!=a) { 
       if (LCA[a][b]==a) return 21;
       else if (LCA[a][b]==b) return 24; else return 27;
   }
} /* compute_triplets_type */



double  Compute_Triplet_Dist(short no_mut1,  short no_mut2, short no_muts_ij,  short m_to_v1[], 
    short m_to_v2[],  short LCA1[][MAXSIZE], short LCA2[][MAXSIZE]){

    short  a, b, c;
    short v;
    int t1, t2;
    short  type1;
    short  type2;
    short x, y, z;
    

    v=0; t1=0; t2=0;
   for (a=0; a<no_mut1; a++) {
      for (b=a+1; b<no_mut1; b++) {
         for (c=b+1; c<no_mut1; c++) {
          if ( m_to_v1[a]!=-1 && m_to_v1[b]!=-1 && m_to_v1[c]!=-1){  
             x= m_to_v1[a]; y= m_to_v1[b]; z= m_to_v1[c];
             type1=Compute_triplet_type(LCA1, x, y, z);
             t1++;
          } else type1=-1;
          if ( m_to_v2[a]!=-1 && m_to_v2[b]!=-1 && m_to_v2[c]!=-1){  
              x= m_to_v2[a]; y= m_to_v2[b]; z= m_to_v2[c];
              type2 = Compute_triplet_type(LCA2, x, y, z);
             t2++;
          } else  type2=-1; 
          
          if (type1==type2 && (type1!=-1 || type2!=-1)) v=v+1;
       }
     }
   }

   /* printf("t1, t2 %d %d\n", t1, t2); */
   if (t1 > t2) return ((double)v)/t1; else return ((double)v)/t2;

}


short Comput_mut_in(short no_muts,  short map1[], short map2[]){
   short i, j, k;

   j=0;
   for (i=0; i<no_muts; i++) {
     if ( map1[i]!=-1 ||  map2[i]!=-1) j++;
   }
   return j;
}



/*
 *   Compile command: gcc -lm Pairwise_Triplets_v2.c -o PW_Triplet
 *   Run command:     ./PW_Triplet <input_file> <output_file> <no_muts>
 *   where input_file contains a number of trees given in the requested format
 *   and the pairwside CASet distance is writen in the <output_file>,
 *   no_muts is an integre, which means, each integer is less then it will be a
 *   mutation label.
 *        */


void main(int argc, char *argv[]){
FILE *In;
char file_name[30];
int i, j, k;

int  no1, no_edges1, no2, no_edges2; 
int u1, u2;
int  start1[MAXSIZE], end1[MAXSIZE];
int  start2[MAXSIZE], end2[MAXSIZE];
char *tree1_names[MAXSIZE], *tree2_names[MAXSIZE];
int  partit1[MAXSIZE][MAXSIZE];
int  partit2[MAXSIZE][MAXSIZE];
int  distances1[MAXSIZE][MAXSIZE], distances2[MAXSIZE][MAXSIZE]; 
int  map[MAXSIZE];
int order;
int node1, node2, node3, node4, node5;
int code[MAXSIZE];

int choice;
FILE *t_OUT;
edge *tree;
edge *all_trees[MAX_NO_TREES];
short graphs[MAX_NO_TREES][MAXSIZE][2];
char *node_labels[MAX_NO_TREES][MAXSIZE];
short distance[MAX_NO_TREES][MAXSIZE][MAXSIZE];
unsigned int count;
int lp1, lp2;
float t_dist, low_dist, a_dist;
int s, t;
char line[200];
char h_line[200];
short interval;
char str1[10], str2[10];
int  no_edges;
int head, tail;
short edge_count;
short node_count;
short node_ind;

short no_mutations[MAX_NO_TREES];
char  *mutations[MAX_NO_TREES][MAXSIZE];
short root;
short no_node_muts[MAX_NO_TREES][MAXSIZE];
char  node_muts[MAX_NO_TREES][MAXSIZE][MAXSIZE];
short mut_to_node_map[MAX_NO_TREES][MAXSIZE];
char  *common_ancestors[MAX_NO_TREES][MAXSIZE][MAXSIZE];
short LCA[MAX_NO_TREES][MAXSIZE][MAXSIZE];
short total_no_muts;
short no_muts_ij;


   if (argc!=4) { 
    printf("Command: pairwise_CASet <inut_file> <output_file> <no_muts>\n");
	 exit(100);
   }
   /* tree processing */
   no_edges=-1;
   In=fopen(argv[1], "r");
   if (In ==NULL) { printf("Tree_file_name is not readable\n"); exit(100);}

   count=-1;
   while (fgets(line, 200, In)>0) {
     line[strlen(line)-1]='\0';
     strcpy(h_line, line);
     if (line[0]=='#') {
         count++;
         edge_count=1;
         node_count=0;
         sscanf(h_line, "%s %s  %d", str1, str2, &no_edges);  
         graphs[count][0][0]=no_edges;
     } else {
         sscanf(h_line, "%s %s", str1, str2);
         if (Check_names(node_labels[count], node_count,  &node_ind, str1)==0){
            node_labels[count][node_ind]=(char *)malloc(strlen(str1)+1);
            strcpy(node_labels[count][node_ind], str1);
            node_count++;
         } 
         graphs[count][edge_count][0]=node_ind; 
         if (Check_names(node_labels[count], node_count,  &node_ind, str2)==0){
            node_labels[count][node_ind]=(char *)malloc(strlen(str2)+1);
            strcpy(node_labels[count][node_ind], str2);
            /* in this case, node_inde=node_count; */
            node_count++;
         } 
         graphs[count][edge_count][1]=node_ind;
         edge_count++;
   }
 }
 count++;


 total_no_muts=(short)atoi(argv[3])+1;

 for (i=0; i<count; i++) {
  
  /*
  printf("----------- %d\n", i);
  */

  Distance_Comput_Jan(graphs[i],  distance[i]);
  root=Orientation_Jan("00", node_labels[i],  graphs[i], distance[i]);
  /*
  Print_Trees(graphs[i], node_labels[i]);
   */
  
  /*
  k=0;
  for (j=0; j<=graphs[i][0][0]; j++) { 
        Decode(node_labels[i][j], '_', &k, mutations[i]);
  }
  no_mutations[i]=k;
  */

 /*
  printf("mutaitons \n");
  for (j=0; j<no_mutations[i]; j++) {
     printf("%s ", mutations[i][j]);
  }
   printf("\n");
 */
  

 
  Comput_Mutation(graphs[i][0][0]+1, node_labels[i], total_no_muts,
 node_muts[i], mut_to_node_map[i]); 

  /*
  for (j=0; j<=graphs[i][0][0]; j++) { 
        printf("id:%2d label(%s): %s\n", j, node_labels[i][j], node_muts[i][j]);
  }
  printf("\n");
  */


  /*
  Comput_Mut_Ancestors(root, graphs[i][0][0]+1, MAXSIZE-1, distance[i], 
   node_muts[i], mut_to_node_map[i],  common_ancestors[i]); 
   */

  Compute_LCA(root, graphs[i][0][0]+1, distance[i], LCA[i]);


  /*
  printf("---------------hhhh----\n");
  for (j=0; j<MAXSIZE-1; j++) {
     for (k=j+1; k<MAXSIZE-1; k++) {
         printf("mut(%d, %d): %s\n", j, k, common_ancestors[i][j][k]);
     }
  }
  */
 }



   t_OUT=fopen(argv[2], "w");
   for (i=0; i<count; i++) {
      for (j=i+1; j<count; j++) {     
        /*
        printf("%d %d: ", i, j);
        */
        no_muts_ij=Comput_mut_in(total_no_muts,
              mut_to_node_map[i], mut_to_node_map[j]);

       t_dist=Compute_Triplet_Dist(total_no_muts,total_no_muts, no_muts_ij, 
          mut_to_node_map[i], mut_to_node_map[j],  LCA[i], LCA[j]);
       fprintf(t_OUT, "%f\n", 1-t_dist);
       /*
       fprintf(t_OUT, "%f\n", t_dist);
       */
   } /* j loop */
    /*
    if (i==0) break;
    */
      
   } /* i loop */
   fclose(t_OUT);
}
