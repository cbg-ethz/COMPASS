/* @Copyright  Louxin Zhang, National University of Singapore, 2019 
 * This porgrm computes the pairwise distances between trees
 * given thorgh the list of edges after header line in the input file:
 *
 *
 *  Tree file format:
 *    Input rooted trees are listed one by one. each tree is given 
 *    by a head line followed by the edges.
 *    The headline give the id of a tree and the number of tree edges 
 *    For example, the header line "#tree 10: 28" represents the 10-th tree of 28 tree
 *    edges.
 *
 *    a tree edge line has two parts. For example, "a b_c_d" represents an edge from
 *    a node that is labeled with a to a node that is labeled with 
 *   a subset {b, c, d} of three labels b, c, d; 
 *   the root is assumed to the the first node.
 *
 *  For example, the following five four lines specify the tree 
      with 5 nodes which are labled with 0, 1, {2, 5}, 3 and 4  and 
      whose root is 0:
  *  #tree 0:  4 
  *  0 1
  *  1 3
  *  3 2_5
  *  3 4
  *
  *
  *  LIMITATION: (1) the program uses array as data strucutre to 
  *     store trees; (2) the max. size of an input tree is set to be 40
  *    nodes; (3) If a node is labeled with a set of mutations, they are seprared
  *    with '_'. (4) the program is only applied to lableld trees with the 
  *    with integer labels from 0 to 40. (5) The number of input trees is set to 20000.
 *
  *  Compile command: gcc  Pairwise_Bourque_DISTS_v3.c -o PW_BD 
  *  Run commpnd:  PW_BD <input_tree_file> 
  *
 */
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>

#define MAXSIZE  40 

/* deinfe the sise of bipartitie graphs used to compute 1 or 2-BD */
#define INFINITY 199
#define EDGE_NO 1000
#define NODE_NO 70 

#define MAX_NO_TREES 20000

#define MAX_MUTS 40


typedef struct {
    unsigned int first;
    unsigned int second;
} edge;


void swap(char str1[], char str2[]) {

    char temp[MAXSIZE];

    strcpy(temp, str1);
    strcpy(str1, str2);
    strcpy(str2, temp);

}


void sort(char str[][MAXSIZE], short len) {

    short  i = 1;
    short  j;
    
   i=1;
    while (i <= len) {

        j = i+1;
        while (j <= len) {
            if  (strcmp(str[i], str[j]) > 0) {
                swap(str[i], str[j]);
            } 
            j++;
        }

        i++;
    }

} /* sorting */





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
    for (i = 0; i < len; i++) { degrees[ code[i]]++; }
    /* Add 1 to them all */
    for (i = 1; i <= n; i++) { degrees[i]++; }
   /* Add edges to nodes in the code */
    for (i = 0; i < len; i++) {
        for (j = 1; degrees[j] != 1; j++);
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


short  Check_Name(char *node_strings[], short  no_nodes, char *str1){
  short i;

  if (str1==NULL) return -1;
  for (i=0; i<no_nodes; i++) {
	if (strcmp(str1, node_strings[i])==0) return i;
  } 
  return -1;
}




short  min_of_two(short  a, short b){
  if (b==0) return a; 
  else if (a==0)  return b;
  
  if (a< b) return a; else return b;
}



void Adjacent_matrix(short  no1, short start1[], short end1[], short adjacent[][MAXSIZE]){
 short i, j, k;

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


void Distance_Comput(short no1, short start1[], short end1[], short distance[][MAXSIZE]){

 short i, j, k;
 short adjacent[MAXSIZE][MAXSIZE];
 short  dist;
 short  total;

 Adjacent_matrix(no1, start1,  end1,  adjacent);


 total=0;
 for (i=0; i<no1; i++) {
   for (j=0; j<no1; j++) { distance[i][j]=adjacent[i][j];
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


void Adjacent_matrix_Jan(short  no1, short graph[][2], 
  short adjacent[][MAXSIZE]){
 short i, j, k;

 for (i=0; i<no1; i++) {
   for (j=0; j<no1; j++) {adjacent[i][j]=0;}
 }

/*
 for (i=0; i<no1; i++) {
     for (k=0; k<no1-1; k++) {
       if (graph[k][0]==i) { adjacent[i][end1[k]]=1; }
       if (end1[k]==i) { adjacent[i][start1[k]]=1; }
     }
 }
*/

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



short Check_Map_Image( short j,  short no1,  short map[]) {
   short i; 
   for (i=1; i<=no1; i++) { if (map[i]==j) return i; }
   return 0;
}




void Comput_map1(short map[], short no1, char *tree1_names[], short no2,  
    char *tree2_names[], short no_common, char *common_names[]){
  short i, j;

  for (i=0; i<no1; i++){
    map[i+1]=-1;
    if (Check_Name(common_names, no_common, tree1_names[i])>=0) {
       for (j=0; j<no2; j++) {
          if (strcmp(tree1_names[i], tree2_names[j])==0) {  map[i+1]=j+1; break;}
       }
    }
  }
} /* end comput_map */




short Get_Index(char *mut_name, short  no_mut, char *mut[]){
    short  k;

    for (k=0; k<no_mut; k++) {
         if (strcmp(mut_name, mut[k])==0) return k;
    }
    return -1;
}

void Decode(char *str, char a, short  *num, char *names[]){
   short len, len1;
   short i, j,k;
   char temp[20];


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


void   Edge_Partitions_Rooted(short  no1,  short start1[],short end1[],  
     short  distance[][MAXSIZE], short  partition1[][MAXSIZE], short r){
 short i, j, k;
 short  dist;
 short head, tail;



 for (i=0; i<no1-1; i++){
       partition1[i][0]=no1;
   for (j=0; j<no1; j++) {
      if (distance[r][start1[i]] < distance[r][end1[i]]){
              head=start1[i]; tail=end1[i]; /* start is more close root. */
      } else { head=end1[i]; tail=start1[i];}

      if (distance[head][j]<distance[tail][j])
            partition1[i][j+1]=0;
      else partition1[i][j+1]=1;  /* all the nodes below an edge */
   }
 }

} /* edge_partition_rooted */


/*
void   Edge_Partitions_Rooted(short no1,  short start1[], short end1[],
      char *tree1_names[],  short distance[][MAXSIZE],  short no_mut,
      char *mut_names[], short partition1[][MAXSIZE], short r){
*/
void  Edge_Partitions_Rooted_Jan(short graph[][2], char *label_names[], 
  short distance[][MAXSIZE],  char partition1[][MAXSIZE],   short  r){
   short i, j, k;
   short  dist;
   short  total;
   short  head, tail;
   short no_mut1;
   char *node_mut[MAXSIZE];
   short ind;
   short no_edges, no_nodes;


  no_edges=graph[0][0];
  no_nodes=no_edges+1;
  
  for (i=1; i<=no_edges; i++){
    /*
    partition1[i][0]=no_mut;
     */
     if (distance[r][graph[i][0]] < distance[r][graph[i][1]]){
          head=graph[i][0]; tail=graph[i][1]; /* start is more close root.*/
     } else { head=graph[i][1]; tail=graph[i][0];}


     for (j=0; j<no_nodes; j++) {
            no_mut1=0;
            Decode(label_names[j], '_', &no_mut1, node_mut);
        if (distance[head][j]<distance[tail][j]){
            /* partition1[i][j+1]=0; */
            for (k=0; k<no_mut1; k++) {
                ind=atoi(node_mut[k]); partition1[i][ind]='0';
            }
        } else {
          /* partition1[i][j+1]=1; */  /* all the nodes below an edge */
            for (k=0; k<no_mut1; k++) {
                ind=atoi(node_mut[k]); partition1[i][ind]='1';
            }
        }
     } /* j loop */

 } /* i loop */

} /* edge_partition_rooted */


short BD_Rooted_Jan(short no_nodes, short no_edges1, short no_edges2, 
   char p1[][MAXSIZE], char p2[][MAXSIZE], short no_mut1, short no_mut2){
short i, j, k;
short  common,  common_mut[MAXSIZE];
char p1_common[MAXSIZE][MAXSIZE], p2_common[MAXSIZE][MAXSIZE];
short n1, n2, Equal, E0, E1;
char *zeros, *ones;


   common=0;
   for (i=0; i<no_nodes; i++) {
     if (p1[1][i]!='a' && p2[1][i] !='a') { common_mut[common]=i; common +=1; }
   }

 if (common >0){
  zeros=(char *)malloc(common+1); ones=(char *)malloc(common+1);
  for (i=0; i<common; i++) { zeros[i]='0'; ones[i]='1'; }
  zeros[common]='\0'; ones[common]='\0';
 }

   for (i=1; i<=no_edges1; i++) {
     for (j=0; j<common; j++) p1_common[i][j]=p1[i][common_mut[j]];
     p1_common[i][common]='\0';
   }

   for (i=1; i<=no_edges2; i++) {
     for (j=0; j<common; j++) p2_common[i][j]=p2[i][common_mut[j]];
     p2_common[i][common]='\0';
   }

   sort(p1_common, no_edges1);
   /*
   for (i=1; i<=no_edges1; i++) { printf("%s\n",  p1_common[i]); }
   */
   sort(p2_common, no_edges2);
   /*
   for (i=1; i<=no_edges2; i++) { printf("%s\n",  p2_common[i]); }
   */
   
   i=1; j=1; k=0;
   if (common>0) {
   while (i<=no_edges1 || j<= no_edges2 ) {
     Equal=strcmp(p1_common[i], p2_common[j]); 
     E0=strcmp(p1_common[i], zeros); E1=strcmp(p1_common[i], ones);
     if (Equal==0) {
       n1=1; while (strcmp(p1_common[i+n1], p1_common[i])==0) n1++; i=i+n1;
       n2=1; while (strcmp(p2_common[j+n2], p2_common[j])==0) n2++; j=j+n2;
       if (E0!=0 && E1!=0) k=k+ min_of_two(n1, n2);
     } else if (Equal <0) {
       n1=1; while (strcmp(p1_common[i+n1], p1_common[i])==0) n1++; i=i+n1;
     } else {
       n2=1; while (strcmp(p2_common[j+n2], p2_common[j])==0) n2++; j=j+n2;
     }
   }
   }

  if (common==no_mut1 && common==no_mut2) return no_edges1+no_edges2- 2*k; 
  else return no_edges1+no_edges2-k;

} /* BD_Rooted */
 

short H_BD_Rooted_Jan(char p1[][MAXSIZE],  short mut_nos1,
      char *mut_pos1[], short edge_nos1, short edges1[], char p2[][MAXSIZE], 
        short mut_nos2, char *mut_pos2[], short edge_nos2, short edges2[]){
 
  short i, j, k;
  short  common,  common_mut[MAXSIZE];
  char p1_common[MAXSIZE][MAXSIZE], p2_common[MAXSIZE][MAXSIZE];
  char *zeros, *ones;
  short n1, n2, Equal, E0,E1;
  short x[MAXSIZE];

  for (i=0; i<MAXSIZE; i++)  x[i]=0; 
  for (i=0; i<mut_nos1; i++) { k=(short)atoi(mut_pos1[i]);  x[k]++; }
  for (i=0; i<mut_nos2; i++) x[atoi(mut_pos2[i])]++;
 
  common=0; 
  for (i=0; i<MAXSIZE; i++) { if (x[i]==2) {common_mut[common]=i; common++;}}
  /*
  printf("                         common %d\n", common);
  */

  if (common==0) k=0; 
  else {

  zeros=(char *)malloc(common+1); ones=(char *)malloc(common+1);
  for (i=0; i<common; i++) { zeros[i]='0'; ones[i]='1'; }
  zeros[common]='\0'; ones[common]='\0';

  for (i=1; i<=edge_nos1; i++) {
     for (j=0; j<common; j++) p1_common[i][j]=p1[edges1[i-1]][common_mut[j]];
     p1_common[i][common]='\0';
     /*
     printf("-------        %d, %s\n",  i, p1[edges1[i-1]]);
     printf("-------       %s\n",  p1_common[i]);
     */
  }
  
  for (i=1; i<=edge_nos2; i++) {
     for (j=0; j<common; j++) p2_common[i][j]=p2[edges2[i-1]][common_mut[j]];
     p2_common[i][common]='\0';
     /*
     printf("------        %d, %s\n",  i, p2[edges2[i-1]]);
     printf("------        %s\n",  p2_common[i]);
     */
  }

  sort(p1_common, edge_nos1); sort(p2_common, edge_nos2);
  
  i=1; j=1; k=0;
   while (i<=edge_nos1 && j<= edge_nos2 ) {
     Equal=strcmp(p1_common[i], p2_common[j]);
     E0=strcmp(p1_common[i], zeros); E1=strcmp(p1_common[i], ones);
      if (Equal==0) {
       n1=1; 
       
       while (i+n1<= edge_nos1 && strcmp(p1_common[i+n1], p1_common[i])==0) n1++;
        i=i+n1;
        n2=1; 
       while (j+n2<= edge_nos2 && strcmp(p2_common[j+n2], p2_common[j])==0) n2++;
        j=j+n2;
       if (E0!=0 && E1!=0) k=k+ min_of_two(n1, n2);
     } else if (Equal <0) {
       n1=1; 
       while (i+n1<= edge_nos1 && strcmp(p1_common[i+n1], p1_common[i])==0) 
       n1++; i=i+n1;
     } else {
       n2=1; 
       while (j+n2<= edge_nos2 && strcmp(p2_common[j+n2], p2_common[j])==0) n2++;
       j=j+n2;
     }
   }
  } /* common */

  /*
  printf( "edges_no1 %d, no2 %d, k value: %d\n", edge_nos1, edge_nos2,  k);
   */

  if (common==mut_nos1 && common==mut_nos2) return edge_nos1+edge_nos2- 2*k;
  /*
  else if (common==0) { 
      return edge_nos1+edge_nos2+1;
  } */
  else return edge_nos1+edge_nos2-k;
} /* H_BD */



short Min_Matching(short BG[][4], short no, short edge_no);

short  Compute_High_BD(short no_nodes1, short no_nodes2, char part1[][MAXSIZE],
        short mut_nos1[], char *mut_pos1[][MAXSIZE], 
        short edge_nos1[],  short edges1[][MAXSIZE], 
        char part2[][MAXSIZE], short mut_nos2[], char *mut_pos2[][MAXSIZE], 
         short edge_nos2[], short edges2[][MAXSIZE], char *label_names1[],
         char *label_names2[]){

     short lp1, lp2, ind;
     short  BG[EDGE_NO][4];

       if (no_nodes1<=no_nodes2) {
         ind=0;
         for (lp1=0; lp1<no_nodes2; lp1++) {
          for (lp2=0; lp2<no_nodes2; lp2++) {
            BG[ind][0]=lp1; BG[ind][1]=no_nodes1+lp2;
             if (lp1<no_nodes1){
              BG[ind][2]=
                 H_BD_Rooted_Jan(part1,  mut_nos1[lp1], mut_pos1[lp1], 
                edge_nos1[lp1], edges1[lp1], 
                 part2, mut_nos2[lp2], mut_pos2[lp2], edge_nos2[lp2], edges2[lp2]);
    /* printf("%s   %s: %d\n", label_names1[lp1], label_names2[lp2], BG[ind][2]);
 */
           } else { BG[ind][2]=edge_nos2[lp2]; }
            ind +=1;
         }
       }
        return Min_Matching(BG, no_nodes2, no_nodes2*no_nodes2);
      } else if (no_nodes1 > no_nodes2) {
         ind=0;
         for (lp1=0; lp1<no_nodes1; lp1++) {
             for (lp2=0; lp2<no_nodes1; lp2++) {
                   BG[ind][0]=lp1; BG[ind][1]=no_nodes1+lp2;
                  if (lp2 <no_nodes2)
                   BG[ind][2]= H_BD_Rooted_Jan(part1,  mut_nos1[lp1], mut_pos1[lp1],
                   edge_nos1[lp1], edges1[lp1],
                 part2, mut_nos2[lp2], mut_pos2[lp2], edge_nos2[lp2], edges2[lp2]);
                  else BG[ind][2]=edge_nos1[lp1];
                  ind +=1;
             }
        }
         return Min_Matching(BG, no_nodes1, no_nodes1*no_nodes1);
   }
} /* COmpute high_-BD */


void Compute_Subgraphs(short gph[][2], char *label_names[], short dist[][MAXSIZE],
 short mut_nos[], char *mut_positions[][MAXSIZE], short edge_nos[], short edges[][MAXSIZE], short order, short r){

short num, num1;
short i, j;
short no_local_muts;
short array[MAXSIZE];

   num1=gph[0][0];
   num=num1+1;

  for (i=0; i<num; i++) {
    no_local_muts=0; /* no of mutations in the neighbore */
    for (j=0; j<num; j++) {
        array[j]=0;
        if  (dist[i][j]<= order
                && dist[i][j]+dist[r][i]==dist[r][j]) {
           array[j]=1;
           Decode(label_names[j], '_', &no_local_muts, mut_positions[i]);
        }
    }
    mut_nos[i]=no_local_muts;
    
    edge_nos[i]=0;
    for (j=1; j<=num1; j++) {
       if (array[gph[j][0]]==1 && array[gph[j][1]]==1) {
           edges[i][edge_nos[i]]=j;  edge_nos[i]=edge_nos[i]+1;
       } 
    }
  }
} 

/* part1[i] hodl the parttion in the ith subgragh */
/* part1[i][j] holds the jth partition, j>0 in the ith subgraph */
/* part1[i][0] holds the mutations in the subgraph */
/* roote is the r-th node of the r-th subgraph */
void Sub_Partition_Rooted_Jan(char parts1[][MAXSIZE][MAXSIZE], 
   short  graph[][2], char *label_names[], short distances1[][MAXSIZE], 
   short order, short r,  short nbr_size1[]){
  short i,j, k, t, w;
  short num, ind, head, tail;
  short no_local_muts;
  char *local_mut[MAXSIZE];
  short no_nbs, nbrs[MAXSIZE], index, nbr_edges[MAXSIZE];
  short no1;
  short edge_ind;

 /*
 Sub_Partition_Rooted(parts1, no1, start1, end1, tree1_names, distances1, order, r1, no_mut1, mut1, nbr_size1);
 */

  no1=graph[0][0]+1;

  for (i=0; i<no1; i++) {
    no_nbs=0;  /* no of nodes in the order-neighbors of node i */
    no_local_muts=0; /* no of mutations in the neighbore */
    for (j=0; j<no1; j++) {
        if  (distances1[i][j]<= order
                && distances1[i][j]+distances1[r][i]==distances1[r][j]) {
           nbrs[no_nbs]=j; no_nbs =1+no_nbs;
           Decode(label_names[j], '_', &no_local_muts, local_mut);
        }
    }
    nbr_size1[i]=no_nbs;
  }

  for (i=0; i<no1; i++) {
    /* num=1; */
    for (j=0; j<no_local_muts; j++) {
           /*
           ind=Get_Index(local_mut[j], no_mut1, mut1);
           */
           index=atoi(local_mut[t]);
           parts1[i][0][index]='1';
           /* num=1+num; */
    }
    /* parts1[i][0][0]=no_local_muts; */
     /* parts1[i][0][0]: No. of muts in the neighbor*/


    edge_ind=0; 
    for (j=1; j<=no1-1; j++){  /* edge by edges */
       head=graph[j][0]; tail=graph[j][1];  /* head is more close to the root r */

       if (distances1[i][head]<=order &&  distances1[i][tail]<=order
           && distances1[i][head]+distances1[r][i]+1==distances1[r][tail] ) {
           /* egde is in the neighor and head is more close to i */
           nbr_edges[edge_ind]=j;
           edge_ind=edge_ind+1;
       }
    }



     for (j=0; j<edge_ind; j++) {
          head=graph[nbr_edges[j]][0]; tail=graph[nbr_edges[j]][1];
          for (k=0; k<no_nbs; k++) {
             no_local_muts=0;
             Decode(label_names[nbrs[k]], '_', &no_local_muts, local_mut);
             if (distances1[tail][nbrs[k]]+1== distances1[head][nbrs[k]]){
                    /* parts1[i][ind][k]=1; */
                   for (t=0; t<no_local_muts; t++) {
                     index=atoi(local_mut[t]);
                    /*
                    for (w=1; w<=num-1; w++) {
                      if (parts1[i][0][w]==index+1) { parts1[i][j+1][w]=1; }
                    }
                    */
                     parts1[i][j+1][index]='1';
                   }
              }  else {  /* parts1[i][ind][k]=0; */
                   for (t=0; t<no_local_muts; t++) {
                       /*
                        index=Get_Index(local_mut[t], no_mut1, mut1);
                       */
                        index=atoi(local_mut[t]);
                        parts1[i][j+1][index]='0';
                   }
             }
          }
    } /* j for */
  } /* i for */
} /* end of subpatition_rooted */



void  Remove_Space(char h_line[], short len){
      short i;

      for (i=0; i<len-1; i++) h_line[i]=h_line[i+1];
      h_line[len-1]='\0';
}

void dijkstra(short G[NODE_NO][NODE_NO], short n, short startnode,
     short distance[], short pred[]) {

        short cost[NODE_NO][NODE_NO];
        short visited[NODE_NO],count,mindistance,nextnode,i,j;

        for(i=0;i<n;i++) {
          for(j=0;j<n;j++) {
            if (j!=i)   cost[i][j]=G[i][j];
            else if (i==j) { cost[i][j]=INFINITY; }
          }
        }

        /*  initialize pred[],distance[] and visited[] */
        for(i=0;i<n;i++) {
          distance[i]=cost[startnode][i]; pred[i]=startnode; visited[i]=0;
        }
      distance[startnode]=0; visited[startnode]=1; count=1;

      while(count<n-1) {
          mindistance=INFINITY;

          /* nexnode gives the node at minimum distance */
          for(i=0;i<n;i++) {
            if((distance[i]< mindistance ) && !visited[i] ){
                mindistance=distance[i]; nextnode=i;
            }
          }


          /* check if a better path exists through nextnode */
          visited[nextnode]=1;
          for(i=0;i<n;i++) {
             if(!visited[i])
                if(mindistance+cost[nextnode][i]<distance[i]) {
                  distance[i]=mindistance+cost[nextnode][i];
                  pred[i]=nextnode;
             }
          }
              count++;
        } /* end of while */

} /* disji */


/* output graph size  no+2*/
/* covered_l[0] is the no of nodes covered by matching
 *  *  * followed by the indices of covered nodes */
void  Construct_ResidualGraph(short G[][NODE_NO], short BG[][4],
     short  edge_no,  short cost_l[], short cost_r[], short no){

     short  i,j,  covered[NODE_NO],   num;
     short s, t;

     s=2*no; t=s+1;
     for (i=0; i<=t; i++) {
        G[i][i]=0;
        for (j=i+1; j<=t; j++) { G[i][j]=INFINITY;G[j][i]=INFINITY; }
     }


     for (i=0;  i<s; i++) covered[i]=0; /* uncvered  */
     for (i=0;  i<edge_no; i++) {
       if (BG[i][3]==1) { covered[BG[i][0]]=1; covered[BG[i][1]]=1; }
     }

    for (i=0; i<no; i++) { if (covered[i]==0) G[s][i]=0; }
    for (i=no; i<s; i++) { if (covered[i]==0) G[i][t]=0; }

     for (i=0; i<edge_no; i++) {
        if (BG[i][3]==1) /* matching edge */
           G[BG[i][1]][BG[i][0]]=0;
        if (BG[i][3]==0) /* no-matching edges */
           G[BG[i][0]][BG[i][1]]=BG[i][2]+cost_l[BG[i][0]]-cost_r[BG[i][1]];
     }
} /* end contrstruction of residual graph */



/* modified matches from a shortest path found */
void  update_node_costs(short distance[], short cost_l[], short cost_r[], 
      short no){
    short i;

    for (i=0; i<no; i++){
       cost_l[i]=cost_l[i]+distance[i];
       cost_r[i+no]=cost_r[i+no]+distance[no+i];
    }
} /* update_node_costs */


void     update_BG( short pred[], short BG[][4], short edge_no, short no){
    short v1, v2; /* (v1, v2) is an edge */
    short i, s;

      s=2*no; v2=pred[s+1];
      do {
           v1=pred[v2];
        if (v2 >=no) {
           for (i=0; i<edge_no; i++) {
               if (BG[i][0]==v1 && BG[i][1]==v2) {BG[i][3]=1; break;}
           }
        } else {
           if (v1!=s){
             for (i=0; i<edge_no; i++) {
               if (BG[i][0]==v2 && BG[i][1]==v1) {BG[i][3]=0; break;}
             }
           }
        }
        v2=v1;
      } while(v2!=s);
} /* update BG */



/* BG[i] has left node BG[i][0], right node BG[i][1], and cost
  BG[i][2]; BG[i][3]=0 if it is not a matching edge and
  1 otherwise
*/
short Min_Matching(short BG[][4], short no, short edge_no){
    short i, j, k;
    short cost_l[NODE_NO], cost_r[NODE_NO];
    short G[NODE_NO][NODE_NO]; 
    short size, weight;
    short distance[NODE_NO], pred[NODE_NO];

   for (i=0; i<no; i++) { cost_l[i]=0; cost_r[no+i]=0;}
   for (i=0; i<edge_no; i++) { BG[i][3]=0; }
   size=0;
   while (size < no) {
     Construct_ResidualGraph(G, BG, edge_no, cost_l, cost_r, no);
     dijkstra(G,2*no+2, 2*no,  distance,  pred);
     update_node_costs(distance, cost_l, cost_r, no);
     update_BG(pred, BG, edge_no, no);
     size +=1;
   }

   weight=0;
   for (i=0; i<edge_no; i++){
      if (BG[i][3]==1)  { weight=weight+BG[i][2]; }
   }
   return weight;
}  /*mini_matching  */


short  Compute_root(short no_edges, short start[], short end[]){
    short i, j;
    short *in;
    short *out;
   

    in=(short *)malloc((no_edges+1)*sizeof(short));
    out=(short *)malloc((no_edges+1)*sizeof(short));
 
    for (i=0; i<no_edges+1; i++) { in[i]=0; out[i]=0; }
    for (i=0; i<no_edges; i++) { in[end[i]] +=1; out[start[i]] +=1; }

    for (i=0; i<no_edges+1; i++) { 
      if (in[i]==0 && out[i]>0)  { j=i; break; } 
    }
     free(in); free(out);
    return j;

}

/*
void Orientation(short r, short no, short start[], short end[], short distance[][MAXSIZE]){
   short i, j, k;

   for (i=0; i<no-1; i++) {
     if (distance[r][start[i]]> distance[r][end[i]]) {
            k=start[i]; start[i]=end[i]; end[i]=k;
     }
   } 
}
*/

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


void Remove(short  graph[][2], short node){
    short i, j, k;
    short nb;
    short no_edges;
    short t;

    no_edges=graph[0][0];
    for  (k=1; k<=no_edges; k++) {
        if (graph[k][0]==node) { nb=graph[k][1]; t=k; break; }
        else if (graph[k][1]==node) { nb=graph[k][0]; t=k;  break;}
    }

    for (k=t+1; k<=no_edges; k++) {
       if (graph[k][0]==node) { 
              graph[k-1][0]=nb;  graph[k-1][1]=graph[k][1];  
       }  else if (graph[k][1]==node) { 
              graph[k-1][1]=nb;  graph[k-1][0]=graph[k][0];  
       } else {
               graph[k-1][0]=graph[k][0];  graph[k-1][1]=graph[k][1];  
       }

    }
    graph[0][0]= graph[0][0]-1;
}

void Perturbation(edge *all_trees[], short graph[][MAXSIZE][2], short size, 
   short no_mut[], char *mutations[], short no_trees){
     int i, j, k;
     int no_removed;

      /*
      printf("No_trees: %d, size: %d\n", no_trees, size);
       */
     srand(time(0));
     for (i=0; i<no_trees; i++) {
        for (k=0; k<size; k++) {
           graph[i][k+1][0]=(short )all_trees[i][k].first-1;
           graph[i][k+1][1]= (short )all_trees[i][k].second-1;
        }

        graph[i][0][0]=size;
        no_mut[i]=size+1;
        no_removed=rand()%200; 
        if (no_removed==0) {
           Remove(graph[i], size); Remove(graph[i], size-1);
           Remove(graph[i], size-2);
           mutations[i][size]='0'; mutations[i][size-1]='0';
           mutations[i][size-2]='0';
           no_mut[i]=size-3;
        } else if (no_removed==1 || no_removed==2) {
          /* removed 2 */
           /* printf("%d: remove two nodes\n", i); */
           Remove(graph[i], size); Remove(graph[i], size-1);
           mutations[i][size]='0'; mutations[i][size-1]='0';
           no_mut[i]=size-2;
        } else if ( no_removed >2 && no_removed <5) {
           /* printf("%d: remove one nodes\n", i); */
           Remove(graph[i], size);
           mutations[i][size]='0';  
           no_mut[i]=size-1;
        }   
     }
}


char Integer_to_str(short a, short i) {
     if (i==1) { return '0'+(a/10); } else { return '0'+(a%10); }
}

short Finding_Label_Table(short graph[][2], short  hash[], short map[]){
      short size;
      short start;
      short k;
       
       start=0;
       for (k=0; k<MAXSIZE; k++) { hash[k]=-1; map[k]=-1; }
       size=graph[0][0];

       for (k=1; k<=size; k++) {
           if (hash[graph[k][0]]==-1) { 
               hash[graph[k][0]]=start;  map[start]=graph[k][0]; start++;
           }
           if (hash[graph[k][1]]==-1) {
               hash[graph[k][1]]=start; map[start]=graph[k][1]; start++;
           }
       }
       return start;
} /* find labels */



void Comput_Labels(short size, short map[], short id,  char *label_names[]){
       short k;

       for (k=0; k<size; k++) {
           if (map[k]!=id) label_names[k]=(char *)malloc(3);
           else label_names[k]=(char *)malloc(9);
           label_names[k][0]=Integer_to_str(map[k], 1);
           label_names[k][1]=Integer_to_str(map[k], 2);
           label_names[k][2]='\0';
       }
} /* compute label */

void Merge_Labels(short graph[][2], short merge_size, char *label_names[])
{
   int i, j, k;
   short merges[3];
   short size;
   short no_labels;
   short hash[MAXSIZE];   /* give a label the index */
   short map[MAXSIZE];
   short t;


/*
  for (j=1; j<=graph[0][0]; j++)
       printf("%d %d\n", graph[j][0], graph[j][1]);
  printf("\n\n size: %d\n", graph[0][0]);
 */

   /*  srand(time(0)); */
    size= graph[0][0];
   merges[0]=(short)abs(rand()%(size+1));
   while (merges[0]==0) merges[0]= (short) abs(rand()%(size+1));
   merges[1]= (short)abs(rand()%(size+1));
   while (merges[1]==0 || merges[1]==merges[0]) 
        merges[1]= (short) abs(rand()%(size+1));
   merges[2]= (short) abs(rand()%(size+1));
   while (merges[2]==0 || merges[2]==merges[0] || merges[2]==merges[1]) 
     merges[2]= (short) abs(rand()%(size+1));

   /*
   printf("-------%d  %d  %d------\n", merges[0],  merges[1], merges[2]); 
   */


   if (merge_size==3) {
       Remove(graph,  merges[1]);
       Remove(graph,  merges[2]);
       
       no_labels=Finding_Label_Table(graph, hash, map);
       Comput_Labels(no_labels, map, merges[0],label_names);
       t=hash[merges[0]];
     
       label_names[t][2]='_';
       label_names[t][3]=Integer_to_str(merges[1], 1);
       label_names[t][4]=Integer_to_str(merges[1], 2);
       label_names[t][5]='_';
       label_names[t][6]=Integer_to_str(merges[2], 1);
       label_names[t][7]=Integer_to_str(merges[2], 2);
       label_names[t][8]='\0';

   } else if (merge_size==2) {
       Remove(graph,  merges[1]);

       no_labels=Finding_Label_Table(graph, hash, map);
       Comput_Labels(no_labels, map, merges[0],label_names);
       t=hash[merges[0]];

       label_names[t][2]='_';
       label_names[t][3]=Integer_to_str(merges[1], 1);
       label_names[t][4]=Integer_to_str(merges[1], 2);
       label_names[t][5]='\0';
   } else {
       no_labels=Finding_Label_Table(graph, hash, map);
       Comput_Labels(no_labels, map, merges[0],label_names);
   }

   size=graph[0][0];

   for (i=1; i<=size; i++) {
       graph[i][0]=hash[graph[i][0]]; graph[i][1]=hash[graph[i][1]];
       /*
       printf("%s,  %s\n", label_names[graph[i][0]], label_names[graph[i][1]]);
       */
   }
   /*
   for (i=0; i<no_labels; i++) {
       printf("i(%d),  %s\n", i, label_names[i]);
   }
    */
   
   
}

void Multiple_Labels(short graph[][MAXSIZE][2], char *label_names[][MAXSIZE],
 int no_trees){
     int i, j, k;
     int no_removed;
     int size;

        srand(time(0));
     for (i=0; i<no_trees; i++) {
         /*printf("%d-th tree\n", i); */
        no_removed=rand()%150;
        if (no_removed==2) {
           /* printf("three_labels are merged\n"); */
           Merge_Labels(graph[i], 3, label_names[i]);
        } else if (no_removed==1 || no_removed==5) {
           /* printf("two_labels are merged\n"); */
           Merge_Labels(graph[i], 2,  label_names[i]);
        } else {
           Merge_Labels(graph[i], 0,  label_names[i]);
        }

     }
}

void Initialiation_Parts(short no_mut,  short graph[][2], char parts1[][MAXSIZE][MAXSIZE]){
 short i, j, k;
 short no_edges, no_nodes; 

 no_edges= graph[0][0]; 
 no_nodes= no_edges+1;

 for (i=0; i<no_nodes; i++) {
   for (k=0; k<no_nodes; k++) parts1[i][0][k]='a'; parts1[i][j][no_nodes]='\0';
   for (j=1; j<=no_edges; j++) {
    for (k=0; k<no_nodes; k++) parts1[i][j][k]='a'; parts1[i][j][no_nodes]='\0';
   }
 }
}

void Initialization(short no_edges, short no_mut, char partition[][MAXSIZE]){
    short i, j;

    for (i=0; i<=no_edges; i++){
       for (j=0; j<no_mut; j++) { partition[i][j]='a'; }
       partition[i][no_mut]='\0';
    }
}

void Print_Trees(short graph[][2], char *label_names[]){
  short  j;

  for (j=1; j<=graph[0][0]; j++) 
      printf("%s %s\n", label_names[graph[j][0]], label_names[graph[j][1]]); 
      /*
      printf("\n");
       */
}

void Printf_Partition(short graph[][2], char *names[], short  no_edges, 
  char partit1[][MAXSIZE], short t){
    short i;

   for (i=1; i<=no_edges; i++) printf("%s  (%s %s) %d\n",  
       partit1[i], names[graph[i][0]], names[graph[i][1]], t);
   printf("\n");
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


/*  compile command: gcc Pairwise_Bourque_DISTS_v3.c -o PW_BD
 *  Run command:   PW_BD  <input_trees_file>
 */
void main(int argc, char *argv[]){
FILE *In, *BD1_OUT, *BD2_OUT, *BD_OUT;
char file_name[100];

short i, j, k, t;
short    u1, u2, order;
edge *all_trees[MAX_NO_TREES];
short graphs[20001][MAXSIZE][2]; 
     /* graph[i][0][0] denote the number of nodes of ith trees */
int  code[40], node1, lp1, lp2,  BD_value, no_trees, count;
short  r1, r2;
char line[200], h_line[200];
char *label_names[MAX_NO_TREES][MAXSIZE];
char *mutations[MAX_NO_TREES][MAX_MUTS];
char all_muts[MAX_MUTS];
short no_mut[MAX_NO_TREES];
short distance[MAX_NO_TREES][MAXSIZE][MAXSIZE];
short no_nodes1, no_nodes2;
char o_partition[MAX_NO_TREES][MAXSIZE][MAXSIZE];
char *node_sub_edge_parti[MAX_NO_TREES][MAXSIZE][MAXSIZE];
short  gph1_mut_nos[MAX_NO_TREES][MAXSIZE];
char   *gph1_mut_positions[MAX_NO_TREES][MAXSIZE][MAXSIZE];
short   gph1_edges[MAX_NO_TREES][MAXSIZE][MAXSIZE];
short  gph1_edge_nos[MAX_NO_TREES][MAXSIZE];

short  gph2_mut_nos[MAX_NO_TREES][MAXSIZE];
char   *gph2_mut_positions[MAX_NO_TREES][MAXSIZE][MAXSIZE];
short   gph2_edges[MAX_NO_TREES][MAXSIZE][MAXSIZE];
short  gph2_edge_nos[MAX_NO_TREES][MAXSIZE];
short root;
int  no_edges;
short edge_count, node_count, node_ind;
char str1[100], str2[100];
/*
char *node_labels[MAX_NO_TREES][MAXSIZE];
*/
short no_muts;
char outfile_names[3][200];





   if (argc!=2) {
    printf("Command: PW_BD  <inut_file>\n");
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
         graphs[count][0][0]=(short)no_edges;
     } else {
         sscanf(h_line, "%s %s", str1, str2);
         if (Check_names(label_names[count], node_count,  &node_ind, str1)==0){
            label_names[count][node_ind]=(char *)malloc(strlen(str1)+1);
            strcpy(label_names[count][node_ind], str1);
            node_count++;
         }
         graphs[count][edge_count][0]=node_ind;
         if (Check_names(label_names[count], node_count,  &node_ind, str2)==0){
            label_names[count][node_ind]=(char *)malloc(strlen(str2)+1);
            strcpy(label_names[count][node_ind], str2);
            /* in this case, node_inde=node_count; */
            node_count++;
         }
         graphs[count][edge_count][1]=node_ind;
         edge_count++;
   }
 }
 count++;



 no_trees=count;


   for (i=0; i<MAX_MUTS-1; i++) all_muts[i]='a'; all_muts[MAX_MUTS-1]='\0'; 
   for (i=0; i<no_trees; i++) {
     k=0;
     for (j=0; j<=graphs[i][0][0]; j++) {
        Decode(label_names[i][j], '_', &k, mutations[i]);
     }
     no_mut[i]=k;
     for (j=0; j<k; j++) {
       all_muts[atoi(mutations[i][j])]='1';
     }
   }


   no_muts=0;
   for (i=0; i<MAX_MUTS; i++){
      if (all_muts[i]=='1') no_muts=i+1;

   }




   for (i=0; i<no_trees; i++) { 
      /*
      printf("#tree %d:  %d\n", i, graphs[i][0][0]);
      Print_Trees(graphs[i], label_names[i]); 
      */
 
      Distance_Comput_Jan(graphs[i],  distance[i]);
      root=Orientation_Jan("00", label_names[i],  graphs[i], distance[i]);
      /*
      Print_Trees(graphs[i], label_names[i]); 
      */
      Initialization(graphs[i][0][0], no_muts, o_partition[i]);
        /* the biparttion of edge i>=1 is in a string partit1[i] */
      Edge_Partitions_Rooted_Jan(graphs[i], label_names[i], distance[i], 
         o_partition[i], root);

      Compute_Subgraphs(graphs[i], label_names[i], distance[i], gph1_mut_nos[i], 
       gph1_mut_positions[i],  gph1_edge_nos[i], gph1_edges[i], 1, root); 

      Compute_Subgraphs(graphs[i], label_names[i], distance[i], gph2_mut_nos[i], 
       gph2_mut_positions[i],  gph2_edge_nos[i], gph2_edges[i], 2, root); 
     
 }

   
   for (i=0; i<3; i++) {
      strcpy(outfile_names[i], argv[1]);
      outfile_names[i][strlen(argv[1])]='_';
      outfile_names[i][strlen(argv[1])+1]='o';
      outfile_names[i][strlen(argv[1])+2]='u';
      outfile_names[i][strlen(argv[1])+3]='t';
      outfile_names[i][strlen(argv[1])+4]='_';
      outfile_names[i][strlen(argv[1])+5]='B';
      outfile_names[i][strlen(argv[1])+6]='D';
   }
      outfile_names[0][strlen(argv[1])+7]='\0';
      outfile_names[1][strlen(argv[1])+7]='1';
      outfile_names[2][strlen(argv[1])+7]='2';
      outfile_names[1][strlen(argv[1])+8]='\0';
      outfile_names[2][strlen(argv[1])+8]='\0';

   printf("------------here\n");
   
   BD_OUT=fopen(outfile_names[0], "w");
   BD1_OUT=fopen(outfile_names[1], "w");
   BD2_OUT=fopen(outfile_names[2], "w");
   printf("---------+%d---here\n", no_trees);

   for (i=0; i<no_trees -1; i++) {
     for (j=i+1; j<no_trees; j++) {
        /*
         printf("            in j(%d) loop \n", j);
         Print_Trees(graph[j], label_names[j]);
         */
       
        fprintf(BD_OUT, "%d\n",   
          BD_Rooted_Jan(no_muts, graphs[i][0][0], graphs[j][0][0],
              o_partition[i], o_partition[j], no_mut[i], no_mut[j]));


       no_nodes1=graphs[i][0][0]+1;
       no_nodes2=graphs[j][0][0]+1;
       BD_value=Compute_High_BD(no_nodes1, no_nodes2, o_partition[i],  
        gph1_mut_nos[i], gph1_mut_positions[i], 
        gph1_edge_nos[i], gph1_edges[i], o_partition[j], 
        gph1_mut_nos[j], gph1_mut_positions[j], 
        gph1_edge_nos[j], gph1_edges[j], label_names[i], label_names[j]);
       fprintf(BD1_OUT, "%d\n", BD_value);

       BD_value=Compute_High_BD(no_nodes1, no_nodes2, o_partition[i],  
        gph2_mut_nos[i], gph2_mut_positions[i], 
        gph2_edge_nos[i], gph2_edges[i], o_partition[j], 
        gph2_mut_nos[j], gph2_mut_positions[j], 
        gph2_edge_nos[j], gph2_edges[j], label_names[i], label_names[j]);
       fprintf(BD2_OUT, "%d\n", BD_value);

     } /* j loop */
     /*
     if (i==0) break;
      */
   } /* i loop */

   fclose(BD_OUT);
   fclose(BD1_OUT);
   fclose(BD2_OUT);

}
